"""Catalog designer entry point.

    python -m igi_catalog.designer --design design.yaml --paths paths.yaml --dataset IGI-SYN-SEQ-01 --out output/
"""
import argparse, csv, json, os, random, sys, time
from collections import Counter, defaultdict
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
        # `events` is the SNV/indel table only, so this counts that table's share. The design target is
        # across every class, and promote_flagposts spreads them, so the total is reported separately or
        # the figure reads as a shortfall against the target when there is none.
        "flagposts": sum(1 for e in events if e["flagpost"]),
        "by_binding_tier_retained": dict(Counter(e.get("binding_tier_retained", "na") for e in events)),
        "best_allele_is_lost": sum(1 for e in events if e.get("best_allele_is_lost")),
        "tier_changed_by_hla_loss": sum(1 for e in events
                                       if e.get("binding_tier") != e.get("binding_tier_retained")),
        "grid_cells_filled": len({(e["clonality_tier"], e["expression_tier"], e["binding_tier"]) for e in events if e["subclass"] == "grid_missense"}),
        "grid_cells_total": len(design["tiers"]["clonality"]) * len(design["tiers"]["expression"]) * len(design["tiers"]["binding"]),
    }
    all_rows = list(events) + list(fusions or [])
    for rows in (extra_tables or {}).values():
        all_rows += list(rows or [])
    summ["flagposts_all_classes"] = {
        "total": sum(1 for e in all_rows if e.get("flagpost")),
        "target": design["counts"].get("flagposts"),
        "by_class": dict(Counter(e.get("class", "") for e in all_rows if e.get("flagpost"))),
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


def promote_flagposts(events, target, log=print):
    """Bring the flagpost count up to the design target across every class.

    A flagpost is a positive control: an event a caller has no excuse for missing. Only the twelve
    hotspot substitutions and four fusions were ever marked, sixteen against a target of forty, so the
    design's figure described an intent nothing implemented. The shortfall is filled from the events
    that are already unambiguous rather than by designing new ones: truncal, strongly binding, well
    expressed and in clean sequence context. Ranking is deterministic, so the same run marks the same
    events.
    """
    have = [e for e in events if e.get("flagpost")]
    if len(have) >= target:
        return have
    rank = {"strong": 0, "weak": 1, "non": 2, "na": 3}
    expr = {"T1000": 0, "T100": 1, "T10": 2, "T1": 3, "T0": 4}

    def tier_of(e):
        # Each class records its level under a different name, and a single lookup on expression_tier
        # silently excluded every CTA, ERV and splice event from promotion. ERVs carry no tier at all,
        # only a target TPM, so theirs is derived the same way every other tier is.
        t = e.get("expression_tier") or e.get("target_expression_tier")
        if t:
            return t
        tpm = e.get("target_tumor_tpm") or e.get("target_tumor_junction_tpm") or e.get("gene_tpm")
        try:
            return Expression.tier(float(tpm)) if tpm not in (None, "") else ""
        except (TypeError, ValueError):
            return ""

    def score(e):
        return (rank.get(e.get("binding_tier"), 3),
                expr.get(tier_of(e), 4),
                0 if e.get("clone") == "T" else 1,
                0 if e.get("context", "clean") == "clean" else 1,
                -float(e.get("expected_vaf_tumor") or e.get("expected_vaf_dna") or 0),
                str(e.get("event_id")))

    pool = [e for e in events
            if not e.get("flagpost")
            and e.get("binding_tier") in ("strong", "weak")
            and tier_of(e) in ("T10", "T100", "T1000")]
    pool.sort(key=score)
    # A floor per class before ranking takes over. Ranked purely on binding and expression the list
    # fills with substitutions, and a caller that handles SNVs but not ERVs would still score perfectly
    # on the positive controls. Each class that has a qualifying event gets at least PER_CLASS_MIN.
    PER_CLASS_MIN = 2
    added = 0
    by_class = defaultdict(list)
    for e in pool:
        by_class[e.get("class", "")].append(e)
    already = Counter(e.get("class", "") for e in have)
    for cls in sorted(by_class):
        for e in by_class[cls][:max(0, PER_CLASS_MIN - already.get(cls, 0))]:
            if len(have) + added >= target:
                break
            e["flagpost"] = True
            added += 1
    for e in pool:
        if len(have) + added >= target:
            break
        if e.get("flagpost"):
            continue
        e["flagpost"] = True
        added += 1
    total = len(have) + added
    log(f"  flagposts: {len(have)} designed as such, {added} promoted -> {total} of {target}"
        + ("" if total >= target else "; not enough qualifying events to reach the target"))
    return [e for e in events if e.get("flagpost")]


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
    promote_flagposts(events + fusions + exc.events, design["counts"].get("flagposts", 0), log=log)
    path, summ = write_outputs(events, des.cards, a.out, a.dataset, design,
                               {"runtime_s": round(time.time() - t0), "scale": a.scale},
                               fusions=fusions, fusion_cards=fus.cards,
                               extra_tables={"svs": svs, "viruses": sv.viral, "expressed": exc.events,
                                             "hla_loh": env.hla_loss.rows()})
    print(json.dumps(summ, indent=2))
    print(f"[done] {len(events)} events -> {path}")


if __name__ == "__main__":
    main()
