#!/usr/bin/env python3
"""Fit capture efficiency as a function of interval GC, from a real exome alignment.

Input is `samtools bedcov <capture.bed> <wes.bam>`: chrom, start, end, summed base depth over the
interval. Output is the JSON that `gc_bias_curve` in the site paths file points at.

The factors are normalised so their captured-base-weighted mean is exactly 1, which is what makes the
curve redistribute depth across GC without changing a library's total yield. Bins with no measurement are
left out, and `pipeline.GcBias.factor()` treats a missing bin as 0 rather than 1: the unmeasured bins are
the extreme GC tails where real coverage is near zero, so falling back to 1.0 there would invent depth the
capture does not produce.

    python3 fit_gc_bias.py --bedcov bedcov.txt --reference GRCh38.fa \
        --source "IPISRC044 BostonGene WES normal, 12M read pairs, bwa mem to GRCh38 noalt" \
        --out ipisrc044_wes_gc_bias.json
"""
import argparse
import json
from collections import defaultdict

import pysam


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bedcov", required=True, help="samtools bedcov output: chrom start end sum_depth")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--source", default="", help="provenance string recorded in the output")
    ap.add_argument("--bin-width", type=int, default=5, help="GC bin width in percent")
    ap.add_argument("--min-bases", type=int, default=0,
                    help="optionally drop a bin holding fewer captured bases than this. Off by default, "
                         "so the fit reproduces the shipped curve exactly. Note that a dropped bin "
                         "becomes factor 0.0 in the builder, not 1.0, so dropping is only safe in the "
                         "extreme GC tails where real coverage is already near zero.")
    a = ap.parse_args()

    fa = pysam.FastaFile(a.reference)
    w = a.bin_width
    depth_sum = defaultdict(float)   # summed base-depth per bin
    base_sum = defaultdict(int)      # captured bases per bin
    n = 0
    total_depth = 0.0
    total_bases = 0

    for line in open(a.bedcov):
        f = line.split()
        if len(f) < 4:
            continue
        chrom, start, end, dsum = f[0], int(f[1]), int(f[2]), float(f[3])
        length = end - start
        if length <= 0:
            continue
        seq = fa.fetch(chrom, start, end).upper()
        acgt = seq.count("A") + seq.count("C") + seq.count("G") + seq.count("T")
        if acgt == 0:
            continue  # an all-N interval has no GC and no meaningful depth
        gc = (seq.count("G") + seq.count("C")) / acgt
        b = min(100 // w - 1, int(gc * 100) // w)
        depth_sum[b] += dsum
        base_sum[b] += length
        total_depth += dsum
        total_bases += length
        n += 1

    mean_depth = total_depth / total_bases

    # Per-bin mean depth relative to the library mean. This is the quantity the builder needs: it says
    # how much of baseline coverage a capture interval at this GC actually receives.
    rel = {}
    for b in sorted(base_sum):
        if base_sum[b] < a.min_bases:
            continue
        rel[b] = (depth_sum[b] / base_sum[b]) / mean_depth

    kept_bases = sum(base_sum[b] for b in rel)
    share = {b: base_sum[b] / kept_bases for b in rel}

    # Normalise to a captured-base-weighted mean of 1. Without this the curve would scale total yield as
    # a side effect of changing its shape, and a depth argument would no longer mean what it says.
    wmean = sum(rel[b] * share[b] for b in rel)
    factors = {b: rel[b] / wmean for b in rel}

    def label(b):
        return f"{b * w}-{(b + 1) * w}"

    out = {
        "source": a.source,
        "mean_depth_measured": round(mean_depth, 3),
        "n_intervals": n,
        "bin_width_pct": w,
        "weighted_mean_before_normalisation": round(wmean, 6),
        "factor_by_gc_bin": {label(b): round(factors[b], 4) for b in sorted(factors)},
        "captured_base_share_by_gc_bin": {label(b): round(share[b], 5) for b in sorted(share)},
    }
    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=2)
        fh.write("\n")

    check = sum(factors[b] * share[b] for b in factors)
    print(f"{n:,} intervals, mean depth {mean_depth:.3f}x, {len(factors)} bins kept")
    print(f"weighted mean after normalisation {check:.6f} (must be 1.000000)")
    for b in sorted(factors):
        print(f"  {label(b):>7}  factor {factors[b]:.4f}  share {share[b]:.5f}")
    assert abs(check - 1.0) < 1e-6, "normalisation failed"


if __name__ == "__main__":
    main()
