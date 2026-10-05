#!/usr/bin/env python3
"""Measure the real np / ec distributions and write catalog/resources/pacbio_np_model.json.

    fit_pacbio_np.py <hifi_wgs.bam> <kinnex_segmented.bam> [n_reads]

np used to be a single constant in every synthetic read while rq varied per read, so a tool filtering on
pass count saw no variation across a whole library. It is now drawn from the empirical CDF this writes.

np is NOT derived from each read's own rq, although that would be more self-consistent: Badread's per-read
accuracy spread is narrower than Revio's, so deriving it gives a far too tight np distribution. The
marginal is taken from real data and rq stays computed from the read's own qualities.
"""
import collections
import json
import os
import statistics as st
import subprocess
import sys


def tags(bam, n_reads):
    npv, ratio = [], []
    p = subprocess.Popen(["samtools", "view", bam], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for i, line in enumerate(p.stdout):
        if i >= n_reads:
            break
        n = e = None
        for f in line.rstrip("\n").split("\t")[11:]:
            if f.startswith("np:i:"):
                n = int(f[5:])
            elif f.startswith("ec:f:"):
                e = float(f[5:])
        if n:
            npv.append(n)
            if e:
                ratio.append(e / n)
    p.stdout.close()
    p.terminate()
    return npv, ratio


def fit(bam, n_reads, label):
    npv, ratio = tags(bam, n_reads)
    c = collections.Counter(npv)
    tot = sum(c.values())
    cum, pmf = 0.0, []
    for k in sorted(c):
        cum += c[k] / tot
        pmf.append([k, round(cum, 8)])
    return {"source": os.path.basename(bam), "n_reads": tot,
            "np_mean": round(st.mean(npv), 3), "np_median": st.median(npv), "np_cdf": pmf,
            "ec_over_np": round(st.median(ratio), 5) if ratio else None}


def main():
    hifi, kinnex = sys.argv[1], sys.argv[2]
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 300_000
    out = {
        "_provenance": "np / ec measured directly from the real BAMs named in each entry; "
                       "see docs/data-sources.md",
        "_why": ("np was a single constant (8) in every synthetic read, so a tool filtering on pass count "
                 "saw no variation while rq varied per read. np is now drawn from this empirical CDF and "
                 "ec is np times the measured ec/np ratio, which is how the two relate in real data."),
        "hifi_wgs": fit(hifi, n, "hifi_wgs"),
        "kinnex": fit(kinnex, n, "kinnex"),
    }
    p = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "..", "..", "catalog", "resources", "pacbio_np_model.json")
    with open(os.path.abspath(p), "w") as fh:
        json.dump(out, fh, indent=1)
    for k in ("hifi_wgs", "kinnex"):
        d = out[k]
        print(f"  {k}: n={d['n_reads']:,} np mean {d['np_mean']} median {d['np_median']} "
              f"ec/np {d['ec_over_np']}")


if __name__ == "__main__":
    main()
