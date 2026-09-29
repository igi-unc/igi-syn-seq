#!/usr/bin/env python3
"""Build one exome library (tumour or normal) for one dataset and one chromosome.

    python3 run_wes.py --design design.yaml --paths paths.yaml --dataset IGI-SYN-SEQ-01 \
        --chrom chr1 --library tumor --depth 150 --out /path/out
"""
import argparse, json, os, subprocess, time
import yaml

from igi_catalog.designer import build_env
from igi_catalog.genome_build import read_events
from igi_catalog.pipeline import WesBuilder, merged_capture, write_record_map
from igi_catalog.readnames import shuffle_and_rename


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
    for key, (path, w) in sorted(sources.items()):
        print(f"    {key[0]} hap{key[1]} {key[2]} cn-profile {key[3]}: weight {w:.4f}  "
              f"{os.path.getsize(path) / 1e6:.1f} MB", flush=True)

    t1 = time.time()
    pieces = wb.simulate(sources, paths["art_cmd"], a.depth,
                         f"{a.dataset}_{a.chrom}_{a.library}", seed=a.seed)
    print(f"  ART done in {time.time() - t1:.0f}s", flush=True)

    # concatenate, then shuffle and give Illumina-style names in one disk-based pass, writing the map
    # from each read name back to its source record (owner decision D2)
    r1 = os.path.join(a.out, f"{a.dataset}_{a.chrom}_{a.library}_R1.fastq.gz")
    r2 = os.path.join(a.out, f"{a.dataset}_{a.chrom}_{a.library}_R2.fastq.gz")
    rmap = os.path.join(a.out, f"{a.dataset}_{a.chrom}_{a.library}_readmap.tsv.gz")
    recmap = os.path.join(a.out, f"{a.dataset}_{a.chrom}_{a.library}_records.tsv.gz")
    cat1 = os.path.join(work, "all_1.fq")
    cat2 = os.path.join(work, "all_2.fq")
    for idx, out in ((0, cat1), (1, cat2)):
        with open(out, "w") as fh:
            subprocess.run(["cat"] + [p[idx] for p in pieces if os.path.exists(p[idx])],
                           stdout=fh, check=True)
    n_reads = shuffle_and_rename(cat1, cat2, r1, r2, rmap, work, seed=a.seed)
    write_record_map(wb.record_map, recmap)
    for f in (cat1, cat2):
        if os.path.exists(f):
            os.remove(f)
    print(f"  {n_reads:,} pairs, Illumina names, maps in {os.path.basename(rmap)} and "
          f"{os.path.basename(recmap)}", flush=True)

    meta = {"dataset": a.dataset, "chrom": a.chrom, "library": a.library, "depth": a.depth,
            "intervals": n_iv, "sources": {f"{k[0]}_hap{k[1]}_{k[2]}_cn{k[3]}": w for k, (_p, w) in sources.items()},
            "coverage_per_source": [{"source": s, "hap": h, "coverage": c} for _r1, _r2, s, h, c in pieces],
            "r1": r1, "r2": r2, "read_map": rmap, "record_map": recmap,
            "pairs_written": n_reads, "runtime_s": round(time.time() - t0)}
    with open(os.path.join(a.out, f"{a.dataset}_{a.chrom}_{a.library}.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    for p in pieces:
        for f in p[:2]:
            if os.path.exists(f):
                os.remove(f)
    print(f"  wrote {r1} ({os.path.getsize(r1) / 1e6:.0f} MB)", flush=True)


if __name__ == "__main__":
    main()
