#!/usr/bin/env python3
"""Illumina short-read whole-genome library for one dataset, one library, one chromosome.

    python3 run_wgs.py --design design.yaml --paths paths.yaml --dataset IGI-SYN-SEQ-01 \
        --chrom chr8 --library tumor --depth 30 --out DIR --work DIR

This is the exome path without the capture step. `WesBuilder` already rebuilds each interval on every
clone and haplotype that retains it, applies germline and somatic edits in reference coordinates, weights
depth by absolute copy number, replaces an interval with a junction contig where a rearrangement sits, and
names reads uniquely across chromosomes. None of that is capture-specific, so WGS reuses it and differs in
three ways:

  * intervals are fixed windows tiling the whole chromosome rather than bait regions, so there is no
    on-target or off-target distinction and no flanking bands;
  * the GC model is flat. The measured curve is *capture efficiency*, fitted from how much sequence the
    baits pulled down at each GC, and applying it to WGS would impose a bias that has no cause here. Real
    WGS does carry a milder PCR-driven GC bias, which this does not model; that is recorded as a known
    limitation rather than approximated with the wrong curve;
  * depth is uniform across the genome, so the copy-number signal is the only thing modulating it.

Windows are 5 Mb so a chromosome is never held in memory twice; the builder fetches and edits per window.
"""
import argparse, json, os, subprocess, sys, time

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from igi_catalog.designer import build_env
from igi_catalog.genome_build import read_events
from igi_catalog.junctions import fusion_junctions, sv_junctions, viral_junctions
from igi_catalog.pipeline import GcBias, WesBuilder, quality_model, write_record_map
from igi_catalog.readnames import NAME_SPACE_STRIDE, shuffle_and_rename

WINDOW = 5_000_000
ORDINAL = {f"chr{i}": i for i in range(1, 23)}
ORDINAL.update({"chrX": 23, "chrY": 24, "chrM": 25})


def whole_chromosome(chrom, length, window=WINDOW):
    """Fixed windows tiling a chromosome, 1-based inclusive."""
    out, pos = [], 1
    while pos <= length:
        out.append((pos, min(length, pos + window - 1)))
        pos += window
    return {chrom: out}


def rows(path):
    import csv
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", required=True)
    ap.add_argument("--paths", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--library", choices=["tumor", "normal"], required=True)
    ap.add_argument("--depth", type=float, required=True)
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--background", default=None)
    ap.add_argument("--extra-events", nargs="*", default=[])
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--max-intervals", type=int, default=None)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    design = yaml.safe_load(open(a.design))
    paths = yaml.safe_load(open(a.paths))
    env = build_env(paths, design, a.dataset)
    purity = design["datasets"][a.dataset]["purity"]
    tables = [t for t in ([a.catalog, a.background] + list(a.extra_events)) if t and os.path.exists(t)]
    events = read_events(*tables)
    print(f"  somatic tables applied: {', '.join(os.path.basename(t) for t in tables)}", flush=True)

    length = env.genome.lengths[a.chrom]
    windows = whole_chromosome(a.chrom, length)
    work = os.path.join(a.work, f"{a.dataset}_{a.chrom}_{a.library}_wgs")
    os.makedirs(work, exist_ok=True)
    os.makedirs(a.out, exist_ok=True)

    t0 = time.time()
    # a flat GC model: the measured curve is capture efficiency and does not apply without baits
    wb = WesBuilder(env, events, purity, work, gc_bias=GcBias(None))
    wb.record_prefix = f"{ORDINAL.get(a.chrom, 0):02d}"
    if a.library == "tumor":
        cat_dir = os.path.dirname(a.catalog)
        jn = (sv_junctions(env.genome, rows(os.path.join(cat_dir, f"{a.dataset}.svs.tsv")))
              + fusion_junctions(env.genome, rows(os.path.join(cat_dir, f"{a.dataset}.fusions.tsv")))
              + viral_junctions(env.genome, paths.get("viral_fasta"),
                                rows(os.path.join(cat_dir, f"{a.dataset}.viruses.tsv"))))
        jn = [j for j in jn if j.get("chrom") == a.chrom or j.get("end_chrom") == a.chrom]
        wb.index_junctions(jn)
        print(f"  {len(jn)} junction contigs on {a.chrom}", flush=True)

    sources, n_iv = wb.write_source_fastas(windows, tumor=(a.library == "tumor"),
                                           max_intervals=a.max_intervals)
    print(f"[{a.dataset} {a.chrom} {a.library} WGS] {n_iv} windows of {WINDOW // 10**6} Mb, "
          f"{len(sources)} sources, {time.time() - t0:.0f}s", flush=True)

    qual = quality_model(paths)
    print(f"  quality model: {qual}", flush=True)
    pieces = wb.simulate(sources, paths["art_cmd"], a.depth,
                         f"{a.dataset}_{a.chrom}_{a.library}_wgs", seed=a.seed, qual=qual)
    print(f"  ART done in {time.time() - t0:.0f}s", flush=True)

    r1 = os.path.join(a.out, f"{a.dataset}_{a.chrom}_{a.library}_wgs_R1.fastq.gz")
    r2 = os.path.join(a.out, f"{a.dataset}_{a.chrom}_{a.library}_wgs_R2.fastq.gz")
    rmap = os.path.join(a.out, f"{a.dataset}_{a.chrom}_{a.library}_wgs_readmap.tsv.gz")
    recmap = os.path.join(a.out, f"{a.dataset}_{a.chrom}_{a.library}_wgs_records.tsv.gz")
    cat1, cat2 = os.path.join(work, "all_1.fq"), os.path.join(work, "all_2.fq")
    for idx, out in ((0, cat1), (1, cat2)):
        with open(out, "w") as fh:
            subprocess.run(["cat"] + [p[idx] for p in pieces if os.path.exists(p[idx])],
                           stdout=fh, check=True)
    # the WGS name space must not collide with the exome's: the exome uses ordinal x stride, so WGS is
    # offset by a further half-stride, which keeps both inside the flowcell's coordinate range
    offset = ORDINAL.get(a.chrom, 0) * NAME_SPACE_STRIDE + NAME_SPACE_STRIDE // 2
    n_reads = shuffle_and_rename(cat1, cat2, r1, r2, rmap, work, seed=a.seed, index_offset=offset)
    write_record_map(wb.record_map, recmap)
    for f in (cat1, cat2):
        if os.path.exists(f):
            os.remove(f)
    # The seed is recorded so a checker cannot be told a different one; see the note in run_wes.py.
    meta = {"dataset": a.dataset, "chrom": a.chrom, "library": a.library, "assay": "WGS",
            "seed": a.seed,
            "depth": a.depth, "window_bp": WINDOW, "windows": n_iv,
            "gc_bias": "flat (the measured curve is capture efficiency)",
            "quality_model": qual, "pairs_written": n_reads,
            "r1": r1, "r2": r2, "read_map": rmap, "record_map": recmap,
            "runtime_s": round(time.time() - t0)}
    with open(os.path.join(a.out, f"{a.dataset}_{a.chrom}_{a.library}_wgs.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    for p in pieces:
        for f in p[:2]:
            if os.path.exists(f):
                os.remove(f)
    print(f"  {n_reads:,} pairs -> {r1} ({os.path.getsize(r1) / 1e6:.0f} MB)", flush=True)


if __name__ == "__main__":
    main()
