#!/usr/bin/env python3
"""Build one exome library (tumour or normal) for one dataset and one chromosome.

    python3 run_wes.py --design design.yaml --paths paths.yaml --dataset IGI-SYN-SEQ-01 \
        --chrom chr1 --library tumor --depth 150 --out /path/out
"""
import argparse, json, os, time
import yaml

from igi_catalog.designer import build_env
from igi_catalog.genome_build import read_events
from igi_catalog.pipeline import WesBuilder, merged_capture, concat_fastqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", required=True)
    ap.add_argument("--paths", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--library", choices=["tumor", "normal"], required=True)
    ap.add_argument("--depth", type=float, required=True)
    ap.add_argument("--catalog", required=True, help="the dataset's snv_indel.tsv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--max-intervals", type=int, default=None)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    design = yaml.safe_load(open(a.design))
    paths = yaml.safe_load(open(a.paths))
    env = build_env(paths, design, a.dataset)
    purity = design["datasets"][a.dataset]["purity"]
    events = read_events(a.catalog)
    capture = merged_capture(paths["exome_bed"], {a.chrom})
    work = os.path.join(a.work, f"{a.dataset}_{a.chrom}_{a.library}")
    os.makedirs(work, exist_ok=True)
    os.makedirs(a.out, exist_ok=True)

    t0 = time.time()
    wb = WesBuilder(env, events, purity, work)
    sources, n_iv = wb.write_source_fastas(capture, tumor=(a.library == "tumor"),
                                           max_intervals=a.max_intervals)
    print(f"[{a.dataset} {a.chrom} {a.library}] {n_iv} intervals, {len(sources)} sources, "
          f"{time.time() - t0:.0f}s", flush=True)
    for (src, hap), (path, w) in sorted(sources.items()):
        print(f"    {src} hap{hap}: weight {w:.4f}  {os.path.getsize(path) / 1e6:.1f} MB", flush=True)

    t1 = time.time()
    pieces = wb.simulate(sources, paths["art_cmd"], a.depth,
                         f"{a.dataset}_{a.chrom}_{a.library}", seed=a.seed)
    print(f"  ART done in {time.time() - t1:.0f}s", flush=True)

    r1 = os.path.join(a.out, f"{a.dataset}_{a.chrom}_{a.library}_R1.fastq.gz")
    r2 = os.path.join(a.out, f"{a.dataset}_{a.chrom}_{a.library}_R2.fastq.gz")
    concat_fastqs(pieces, r1, r2)
    meta = {"dataset": a.dataset, "chrom": a.chrom, "library": a.library, "depth": a.depth,
            "intervals": n_iv, "sources": {f"{s}_hap{h}": w for (s, h), (_p, w) in sources.items()},
            "coverage_per_source": {f"{s}_h{h}": c for _r1, _r2, s, h, c in pieces},
            "r1": r1, "r2": r2, "runtime_s": round(time.time() - t0)}
    with open(os.path.join(a.out, f"{a.dataset}_{a.chrom}_{a.library}.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    for p in pieces:
        for f in p[:2]:
            if os.path.exists(f):
                os.remove(f)
    print(f"  wrote {r1} ({os.path.getsize(r1) / 1e6:.0f} MB)", flush=True)


if __name__ == "__main__":
    main()
