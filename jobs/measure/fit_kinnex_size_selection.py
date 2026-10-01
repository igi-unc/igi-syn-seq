#!/usr/bin/env python3
"""Fit Kinnex library size selection: how a cDNA's length changes its chance of reaching a read.

Kinnex gives one read per molecule, and TPM is already length-normalised, so the molecule sampling weight
is abundance alone -- NOT abundance x length, which is the short-read weight. Measured on the release RNA
manifest the two differ by 2.3x in mean cDNA length (1,788 bp against 4,197 bp), so using the wrong one
would badly over-represent long transcripts.

Even with the right weight the simulated molecules come out shorter than real segments: p10 498 bp against
1,355, p50 1,173 against 1,974. Real library prep loses short cDNA -- bead cleanups, and short fragments
ligating into arrays less efficiently. This script fits that loss as an efficiency per length bin, the same
shape of correction as the WES GC-bias curve: the ratio of the real length density to the simulated
molecule density, normalised so the abundance-weighted mean efficiency is 1 and total yield is unchanged.

One honest caveat, recorded in the output: the real segments are HG002 and the simulated molecules come from
a TCGA-BRCA basal baseline, so the ratio conflates size selection with genuine transcriptome differences. At
the short end selection dominates -- a 500 bp cDNA is not three times rarer in one transcriptome than
another -- but the curve should not be read as pure prep efficiency.

    python3 fit_kinnex_size_selection.py --real-lengths seglen.txt \
        --manifest <dataset>_full_rna_transcripts.tsv --out kinnex_size_selection.json
"""
import argparse
import csv
import json
from collections import defaultdict

# Bin edges in bp. Fine where most of the mass and most of the disagreement sits, coarse in the long tail.
EDGES = [0, 250, 500, 750, 1000, 1250, 1500, 2000, 2500, 3000, 4000, 5000, 7000, 10000, 15000, 10**9]


def binof(n):
    for i in range(len(EDGES) - 1):
        if EDGES[i] <= n < EDGES[i + 1]:
            return i
    return len(EDGES) - 2


def label(i):
    hi = EDGES[i + 1]
    return f"{EDGES[i]}-{'inf' if hi >= 10**9 else hi}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--real-lengths", required=True, help="one segment length per line, from a real BAM")
    ap.add_argument("--manifest", required=True, help="RNA transcript manifest: length + abundance columns")
    ap.add_argument("--out", required=True)
    ap.add_argument("--source", default="")
    ap.add_argument("--max-efficiency", type=float, default=6.0,
                    help="cap, so a sparsely populated bin cannot dominate the library")
    a = ap.parse_args()

    real = defaultdict(float)
    n_real = 0
    for line in open(a.real_lengths):
        line = line.strip()
        if not line:
            continue
        real[binof(int(line))] += 1
        n_real += 1

    sim = defaultdict(float)
    n_sim = 0.0
    for r in csv.DictReader(open(a.manifest), delimiter="\t"):
        try:
            L, w = int(r["length"]), float(r["abundance"])
        except (KeyError, ValueError, TypeError):
            continue
        if L <= 0 or w <= 0:
            continue
        sim[binof(L)] += w          # abundance alone: one molecule, one read
        n_sim += w

    eff, rows = {}, []
    for i in range(len(EDGES) - 1):
        pr = real.get(i, 0.0) / n_real
        ps = sim.get(i, 0.0) / n_sim
        if ps <= 0:
            continue                # nothing simulated here, so no ratio to take
        e = min(a.max_efficiency, pr / ps)
        eff[i] = e
        rows.append((i, pr, ps, e))

    # normalise so the abundance-weighted mean efficiency is 1: the curve redistributes which molecules are
    # sequenced without changing how many reads the library yields
    wsum = sum(eff[i] * (sim[i] / n_sim) for i in eff)
    for i in eff:
        eff[i] /= wsum

    out = {
        "source": a.source,
        "note": ("Efficiency per cDNA length bin, applied as a multiplier on abundance when sampling "
                 "molecules. Normalised to an abundance-weighted mean of 1. Conflates library size "
                 "selection with transcriptome differences between the real sample and the simulated one; "
                 "see fit_kinnex_size_selection.py."),
        "real_segments": n_real,
        "bin_edges_bp": EDGES[:-1] + ["inf"],
        "efficiency_by_length_bin": {label(i): round(eff[i], 4) for i in sorted(eff)},
        "real_density_by_length_bin": {label(i): round(real.get(i, 0.0) / n_real, 5) for i in sorted(eff)},
        "simulated_density_by_length_bin": {label(i): round(sim[i] / n_sim, 5) for i in sorted(eff)},
    }
    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=2)
        fh.write("\n")

    print(f"{n_real:,} real segments, {n_sim:,.0f} units of simulated abundance")
    print(f"{'bin':>12} {'real':>8} {'sim':>8} {'efficiency':>11}")
    for i, pr, ps, _e in rows:
        print(f"{label(i):>12} {pr:8.4f} {ps:8.4f} {eff[i]:11.3f}")
    chk = sum(eff[i] * (sim[i] / n_sim) for i in eff)
    print(f"\nabundance-weighted mean efficiency after normalisation: {chk:.6f} (must be 1.000000)")
    assert abs(chk - 1.0) < 1e-6
    # How hard the curve is working. This is the number to look at, because the obvious "check" -- that
    # the corrected density matches the real one -- is an identity, not a test: the efficiency IS that
    # ratio, so it reproduces the target by construction. Only a built library can validate the curve.
    strong = [(label(i), eff[i]) for i in sorted(eff) if eff[i] < 0.1 or eff[i] > 2.5]
    print("\nbins where the curve moves density by more than 2.5x:")
    for lab, e in strong:
        print(f"{lab:>12} {e:8.3f}")
    print("\nThese are large corrections and they are not pure library chemistry: the real segments come "
          "from HG002\nand the simulated molecules from a TCGA-BRCA basal baseline, so the ratio also "
          "carries the difference\nbetween two transcriptomes, and the low-abundance short tail that a "
          "0.01 TPM floor keeps but a real\nlibrary barely samples. The curve corrects the observable; it "
          "does not explain it.")


if __name__ == "__main__":
    main()
