"""Catalog designer entry point.

    python -m igi_catalog.designer --design design.yaml --paths paths.yaml --dataset IGI-SYN-SEQ-01 --out output/
"""
import argparse, csv, json, os, random, sys, time
from types import SimpleNamespace

import yaml

from .annotation import load_or_build, representative_transcripts
from .genome import Genome
from .expression import Expression
from .germline import Germline
from .context import ContextAnnotator
from .clones import CloneModel
from .hla_loss import HlaLoss
from .binding import NetMHCpan
from .snv_indel import SnvIndelDesigner
from .fusions import FusionDesigner
from .structural import StructuralDesigner
from .expressed_classes import ExpressedClassDesigner


def build_env(paths, design, ds_name):
    dcfg = design["datasets"][ds_name]
    t0 = time.time()
    tx = load_or_build(paths["gtf"], paths.get("annotation_cache"))
    rep = representative_transcripts(tx)
    genome = Genome(paths["reference_fasta"])
    expr = Expression(paths["expression_tsv"], tx)
    germline = Germline(paths["germline_vcf"][dcfg["baseline"]], sex=dcfg["sex"])
    ctx = ContextAnnotator(genome, germline, paths["exome_bed"], paths["rmsk_bed"], paths["segdups_bed"])
    clones = CloneModel(dcfg, paths["arms_bed"])
    work = os.path.join(paths["workdir"], ds_name)
    os.makedirs(work, exist_ok=True)
    netmhc = NetMHCpan(paths["netmhcpan_cmd"], os.path.join(work, "netmhcpan"), os.path.join(work, "netmhcpan_cache.json"), threads=int(paths.get("netmhcpan_threads", 6)))
    hla = sorted(set(dcfg["hla"]))
    # owner decision D4: every event is scored twice, once over all alleles and once over the alleles the
    # tumour still has after HLA loss, and carries a flag when its best allele is one that was lost
    hla_loss = HlaLoss(dcfg, clones)
    hla_retained, hla_lost = hla_loss.retained(), hla_loss.lost()
    print(f"[env] {ds_name}: {len(tx)} transcripts, {len(rep)} representative CDS, loaded in {time.time() - t0:.0f}s", flush=True)
    if hla_lost:
        print(f"[env] {ds_name}: HLA lost in clone {hla_loss.region['clone']} "
              f"({hla_loss.region['label']}): {', '.join(hla_lost)}; retained: {', '.join(hla_retained)}",
              flush=True)
    return SimpleNamespace(tx=tx, rep=rep, genome=genome, expr=expr, germline=germline, ctx=ctx,
                           clones=clones, netmhc=netmhc, hla=hla, hla_loss=hla_loss,
                           hla_retained=hla_retained, hla_lost=hla_lost, work=work)


def write_table(events, path):
    cols = []
    for e in events:
        for k in e:
            if k not in cols:
                cols.append(k)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for e in events:
            w.writerow({k: ("" if v is None else v) for k, v in e.items()})
    return path


def write_cards(cards, path, header):
    with open(path, "w") as fh:
        fh.write(header + "\n\n")
        fh.write("\n".join(cards))


def write_outputs(events, cards, out_dir, ds_name, design, summary_extra=None, fusions=None, fusion_cards=None,
                  extra_tables=None):
    os.makedirs(out_dir, exist_ok=True)
    path = write_table(events, os.path.join(out_dir, f"{ds_name}.snv_indel.tsv"))
    write_cards(cards, os.path.join(out_dir, f"{ds_name}.snv_indel.diffcards.txt"),
                f"# {ds_name} designed SNV/indel events: haplotype-aware diff cards (design v{design['version']}, seed {design['seed']})")
    if fusions:
        write_table(fusions, os.path.join(out_dir, f"{ds_name}.fusions.tsv"))
        write_cards(fusion_cards or [], os.path.join(out_dir, f"{ds_name}.fusions.diffcards.txt"),
                    f"# {ds_name} designed gene fusions: junction cards (design v{design['version']}, seed {design['seed']})")
    for name, rows in (extra_tables or {}).items():
        if rows:
            write_table(rows, os.path.join(out_dir, f"{ds_name}.{name}.tsv"))
    # summary counts
    from collections import Counter
    summ = {
        "dataset": ds_name, "n_events": len(events),
        "by_class": dict(Counter(e["class"] for e in events)),
        "by_subclass": dict(Counter(e["subclass"] for e in events)),
        "by_consequence": dict(Counter(e["consequence"] for e in events)),
        "by_clonality_tier": dict(Counter(e["clonality_tier"] for e in events)),
        "by_expression_tier": dict(Counter(e["expression_tier"] for e in events)),
        "by_binding_tier": dict(Counter(e["binding_tier"] for e in events)),
        "by_context": dict(Counter(e["context"] for e in events)),
        "chr1to6_fraction": round(sum(1 for e in events if e["chr1to6"]) / max(1, len(events)), 3),
        "flagposts": sum(1 for e in events if e["flagpost"]),
        "by_binding_tier_retained": dict(Counter(e.get("binding_tier_retained", "na") for e in events)),
        "best_allele_is_lost": sum(1 for e in events if e.get("best_allele_is_lost")),
        "tier_changed_by_hla_loss": sum(1 for e in events
                                       if e.get("binding_tier") != e.get("binding_tier_retained")),
        "grid_cells_filled": len({(e["clonality_tier"], e["expression_tier"], e["binding_tier"]) for e in events if e["subclass"] == "grid_missense"}),
        "grid_cells_total": len(design["tiers"]["clonality"]) * len(design["tiers"]["expression"]) * len(design["tiers"]["binding"]),
    }
    for name, rows in (extra_tables or {}).items():
        if rows:
            summ[name] = {"n": len(rows), "by_subclass": dict(Counter(r.get("subclass", "") for r in rows))}
    if fusions:
        summ["fusions"] = {
            "n": len(fusions),
            "by_subclass": dict(Counter(f["subclass"] for f in fusions)),
            "by_mechanism": dict(Counter(f["mechanism"] for f in fusions)),
            "in_frame": sum(1 for f in fusions if f["in_frame"]),
            "wes_visible": sum(1 for f in fusions if f["wes_visible"]),
            "by_binding_tier": dict(Counter(f["binding_tier"] for f in fusions)),
        }
    if summary_extra:
        summ.update(summary_extra)
    with open(os.path.join(out_dir, f"{ds_name}.snv_indel.summary.json"), "w") as fh:
        json.dump(summ, fh, indent=2)
    return path, summ


def main(argv=None):
    ap = argparse.ArgumentParser(description="IGI-SYN-SEQ catalog designer")
    ap.add_argument("--design", required=True)
    ap.add_argument("--paths", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--scale", type=float, default=1.0, help="multiply all counts (smoke tests)")
    a = ap.parse_args(argv)
    design = yaml.safe_load(open(a.design))
    paths = yaml.safe_load(open(a.paths))
    if a.scale != 1.0:
        c = design["counts"]
        for k, v in list(c.items()):
            if isinstance(v, int) and k != "flagposts":
                c[k] = max(1, int(v * a.scale))
            elif isinstance(v, dict):
                c[k] = {kk: max(1, int(vv * a.scale)) if isinstance(vv, int) else vv for kk, vv in v.items()}
    design["pool_scale"] = a.scale
    rng = random.Random(f"{design['seed']}:{a.dataset}")
    env = build_env(paths, design, a.dataset)
    dcfg = design["datasets"][a.dataset]
    des = SnvIndelDesigner(env, a.dataset, dcfg, design, rng)
    t0 = time.time()
    log = lambda m: print(m, flush=True)
    events = des.design(log=log)
    clones = list(dcfg["clones"])
    fus = FusionDesigner(env, a.dataset, dcfg, design, rng, des.next_id)
    fusions = fus.design(design["counts"]["fusions"], clones, log=log)
    sv = StructuralDesigner(env, a.dataset, dcfg, design, rng, des.next_id)
    svs = sv.design(design["counts"]["svs"], clones, log=log)
    exc = ExpressedClassDesigner(env, a.dataset, dcfg, design, rng, des.next_id, paths)
    exc.design_ctas(design["counts"]["ctas"], clones, log=log)
    exc.design_ervs(design["counts"]["ervs"], clones, log=log)
    exc.design_splice(design["counts"]["splice"], clones, des, log=log)
    path, summ = write_outputs(events, des.cards, a.out, a.dataset, design,
                               {"runtime_s": round(time.time() - t0), "scale": a.scale},
                               fusions=fusions, fusion_cards=fus.cards,
                               extra_tables={"svs": svs, "viruses": sv.viral, "expressed": exc.events,
                                             "hla_loh": env.hla_loss.rows()})
    print(json.dumps(summ, indent=2))
    print(f"[done] {len(events)} events -> {path}")


if __name__ == "__main__":
    main()
