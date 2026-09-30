#!/usr/bin/env python3
"""PacBio HiFi whole-genome library for one dataset, one library, one chromosome.

    python3 run_pacbio_wgs.py --design design.yaml --paths paths.yaml --dataset IGI-SYN-SEQ-01 \
        --chrom chr8 --library tumor --depth 30 --out DIR --work DIR

This is the first read path that builds from `rearrange.py`. The exome path rebuilds each capture interval
in reference coordinates and represents a rearrangement as a junction contig, which is the right
compromise for 150 bp reads. A HiFi read is long enough to span a breakpoint and carry flanking sequence
on both sides, so the tumour genome has to exist as a derived chromosome: segments in their rearranged
order, with germline and somatic edits applied, and a coordinate map from derived position back to
reference so the truth stays checkable.

Depth follows the clone model exactly as the exome does, so an amplified region is deeper and a lost one
shallower; there is no capture step, so no GC or bait model applies.
"""
import argparse, json, os, subprocess, sys, time
from collections import defaultdict

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from igi_catalog.designer import build_env
from igi_catalog.genome_build import read_events
from igi_catalog.pipeline import depth_denominator, source_plan
from igi_catalog import rearrange


def sv_rows(catalog_dir, dataset, chrom):
    """Designed structural variants on this chromosome, as rearrange.plan_chromosome expects them."""
    import csv
    path = os.path.join(catalog_dir, f"{dataset}.svs.tsv")
    if not os.path.exists(path):
        return []
    out = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("chrom") != chrom or r.get("svtype") == "TRA":
                continue        # translocations need both partners; handled separately
            try:
                out.append({"event_id": r["event_id"], "svtype": r["svtype"],
                            "start": int(r["start"]), "end": int(r["end"]),
                            "clone": r.get("clone") or "T",
                            "haplotype": int(r.get("haplotype") or 0),
                            "inserted_sequence": r.get("inserted_sequence") or ""})
            except (ValueError, KeyError):
                continue
    return out


def mei_sequence_factory(rng):
    """Sequence for a mobile-element insertion. The catalog records a family and a length, not bases."""
    def make(event):
        return "".join(rng.choice("ACGT") for _ in range(300))
    return make


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", required=True)
    ap.add_argument("--paths", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--library", choices=["tumor", "normal"], required=True)
    ap.add_argument("--depth", type=float, required=True)
    ap.add_argument("--catalog-dir", default="output")
    ap.add_argument("--background", default=None)
    ap.add_argument("--extra-events", nargs="*", default=[])
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    design = yaml.safe_load(open(a.design))
    paths = yaml.safe_load(open(a.paths))
    env = build_env(paths, design, a.dataset)
    dcfg = design["datasets"][a.dataset]
    purity = dcfg["purity"]
    import random
    rng = random.Random(f"{a.seed}:{a.dataset}:{a.chrom}:{a.library}")

    tables = [os.path.join(a.catalog_dir, f"{a.dataset}.snv_indel.tsv")]
    tables += [t for t in ([a.background] + list(a.extra_events)) if t and os.path.exists(t)]
    events = read_events(*tables)
    print(f"[{a.dataset} {a.chrom} {a.library}] somatic tables: "
          f"{', '.join(os.path.basename(t) for t in tables)}", flush=True)

    work = os.path.join(a.work, f"{a.dataset}_{a.chrom}_{a.library}")
    os.makedirs(work, exist_ok=True)
    os.makedirs(a.out, exist_ok=True)
    chrom_len = env.genome.lengths[a.chrom]
    denom = depth_denominator(env.clones, purity)
    tumor = a.library == "tumor"
    svs = sv_rows(a.catalog_dir, a.dataset, a.chrom) if tumor else []

    t0 = time.time()
    pieces, plan_rows = [], []
    for src, hap, kind, w in source_plan(env.clones, purity, a.chrom, chrom_len // 2, tumor=tumor):
        clone = "T" if (src == "NORMAL" or kind == "pre_cna") else src
        pre_only = kind == "pre_cna"
        mine = [s for s in svs
                if s["haplotype"] == hap and env.clones.is_descendant(clone, s["clone"])] if src != "NORMAL" else []
        # plan_chromosome returns (segments, applied, dropped). `dropped` is truth, not noise: an SV that
        # could not be placed -- overlapping another, or too near a chromosome end -- is absent from the
        # reads, and a truth table that still claims it would be wrong in the way F9 was.
        segs, applied, dropped = rearrange.plan_chromosome(a.chrom, chrom_len, mine)
        if dropped:
            print(f"    {src} hap{hap} {kind}: {len(dropped)} SV(s) not placeable: "
                  f"{[d[-1].get('event_id') if isinstance(d, tuple) else d for d in dropped][:4]}", flush=True)
        seq = rearrange.realise(env.genome, env.germline, events, env.clones, clone, hap, segs,
                                mei_sequence=mei_sequence_factory(rng), only_pre_cna=pre_only)
        if len(seq) < 10_000:
            continue
        tag = f"{src}_hap{hap}_{kind}"
        fa = rearrange.write_derived(os.path.join(work, f"{tag}.fa"), f"{a.chrom}_{tag}", seq)
        cmap = rearrange.write_coordinate_map(
            os.path.join(work, f"{tag}.coordmap.tsv"),
            [(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]) for r in rearrange.coordinate_map(segs)])
        cov = max(0.01, a.depth * w / denom)
        plan_rows.append({"source": src, "haplotype": hap, "copy_kind": kind, "weight": round(w, 5),
                          "coverage": round(cov, 4), "derived_length": len(seq),
                          "reference_length": chrom_len, "segments": len(segs),
                          "svs_offered": [s["event_id"] for s in mine],
                          "svs_applied": applied,
                          "svs_dropped": [(d[-1].get("event_id") if isinstance(d, tuple) else str(d))
                                          for d in dropped],
                          "fasta": fa, "coordinate_map": cmap})
        pieces.append((tag, fa, cov))
        print(f"    {tag}: {len(segs)} segments, derived {len(seq):,} bp "
              f"(reference {chrom_len:,}), {len(applied)}/{len(mine)} SVs applied, "
              f"coverage {cov:.3f}x", flush=True)
    print(f"  derived chromosomes built in {time.time() - t0:.0f}s", flush=True)

    # pbsim3 multi-pass per derived chromosome, then ccs, which is how HiFi accuracy is actually
    # produced: consensus across passes. `--accuracy-mean` is silently ignored by both pbsim3 methods --
    # asking for 0.998 and for 0.999 both yield 0.964 -- so a single-pass run would have shipped
    # CLR-era 3.6%-error reads labelled HiFi.
    t1 = time.time()
    model = paths["pbsim_errhmm"]
    pb, ccs = paths["pbsim_cmd"], paths["ccs_cmd"]
    npass = int(paths.get("pbsim_pass_num", 10))
    bams = []
    for i, (tag, fa, cov) in enumerate(pieces):
        pre = os.path.join(work, f"sim_{tag}")
        args = (f"--strategy wgs --method errhmm --errhmm {model} --genome {fa} "
                f"--depth {cov:.4f} --length-mean {paths.get('pbsim_length_mean', 16689)} "
                f"--length-sd {paths.get('pbsim_length_sd', 4593)} --pass-num {npass} "
                f"--prefix {pre} --id-prefix {tag} --seed {a.seed + i}")
        subprocess.run(pb.format(args=args), shell=True, check=True, capture_output=True, text=True)
        subreads = sorted(f for f in os.listdir(work) if f.startswith(f"sim_{tag}_") and f.endswith(".bam"))
        if not subreads:
            raise RuntimeError(f"pbsim3 produced no subread BAM for {tag}; it writes .bam per reference "
                               f"and .fq.gz only in single-pass mode")
        for sb in subreads:
            out = os.path.join(work, f"ccs_{tag}_{sb}")
            subprocess.run(ccs.format(args=f"--num-threads 8 --min-passes 3 --min-rq 0.99 "
                                           f"{os.path.join(work, sb)} {out}"),
                           shell=True, check=True, capture_output=True, text=True)
            bams.append(out)
    if not bams:
        raise RuntimeError("ccs produced no HiFi BAM")
    print(f"  pbsim3 + ccs done in {time.time() - t1:.0f}s, {len(bams)} HiFi BAM(s)", flush=True)

    out_bam = os.path.join(a.out, f"{a.dataset}_{a.chrom}_{a.library}_hifi.bam")
    sam = paths["samtools_cmd"]
    if len(bams) == 1:
        subprocess.run(sam.format(args=f"view -b -o {out_bam} {bams[0]}"), shell=True, check=True)
    else:
        subprocess.run(sam.format(args=f"merge -f -o {out_bam} " + " ".join(bams)),
                       shell=True, check=True)
    subprocess.run(sam.format(args=f"index {out_bam}"), shell=True, check=True)
    n_reads = int(subprocess.run(sam.format(args=f"view -c {out_bam}"), shell=True,
                                 capture_output=True, text=True).stdout.strip() or 0)

    meta = {"dataset": a.dataset, "chrom": a.chrom, "library": a.library, "depth": a.depth,
            "platform": "PacBio Revio HiFi", "simulator": "pbsim3 errhmm",
            "reads": n_reads, "bam": out_bam, "sources": plan_rows,
            "pass_num": npass, "ccs_min_passes": 3, "ccs_min_rq": 0.99,
            "runtime_s": round(time.time() - t0)}
    with open(os.path.join(a.out, f"{a.dataset}_{a.chrom}_{a.library}_hifi.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    print(f"  {n_reads:,} HiFi reads -> {out_bam} "
          f"({os.path.getsize(out_bam) / 1e6:.0f} MB)", flush=True)


if __name__ == "__main__":
    main()
