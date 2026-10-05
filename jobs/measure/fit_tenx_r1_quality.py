#!/usr/bin/env python3
"""Measure real 10x read-1 per-cycle qualities and write catalog/resources/tenx_r1_quality.json.

    fit_tenx_r1_quality.py <gex_R1.fastq.gz> <tcr_R1.fastq.gz> [n_reads]

R1 qualities were drawn as a flat two-level string -- the top bin 98 % of the time, one lower bin
otherwise -- which gave 663 distinct quality strings in 20,000 reads. Real NovaSeq X R1 is four-level
binned (`#*9I`, so Q2/Q9/Q24/Q40) with a per-cycle shape: the first cycle is the worst at about 83 % top
bin and later cycles run above 93 %.

Note on the target. R2 of the same library has 19,461 distinct quality strings in 20,000 reads, which is
NOT what R1 should match: R1 is only 26 four-level-binned cycles, and the real library has 2,819. That is
the number the model is checked against.
"""
import collections
import gzip
import json
import os
import sys


def fit(path, n_max):
    per = collections.defaultdict(collections.Counter)
    n = 0
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 3:
                for c, ch in enumerate(line.rstrip("\n")):
                    per[c][ch] += 1
                n += 1
                if n >= n_max:
                    break
    cycles = []
    for c in range(max(per) + 1):
        tot = sum(per[c].values())
        cycles.append(sorted(((ch, round(k / tot, 6)) for ch, k in per[c].items()), key=lambda t: -t[1]))
    return {"n_reads": n, "read_length": len(cycles), "cycles": cycles,
            "alphabet": sorted({ch for c in cycles for ch, _ in c}),
            "source": os.path.basename(path)}


def main():
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 500_000
    out = {
        "_provenance": "per-cycle R1 quality distributions measured from the real IPISRC044 10x libraries "
                       "named below",
        "_why": ("R1 qualities were drawn as a flat two-level string, which gave 663 distinct quality "
                 "strings in 20,000 reads where the real library has 2,819. Real NovaSeq X R1 is "
                 "four-level binned with a per-cycle shape, so it is measured and drawn per cycle."),
        "gex": fit(sys.argv[1], n),
        "tcr": fit(sys.argv[2], n),
    }
    p = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "..", "..", "catalog", "resources", "tenx_r1_quality.json")
    with open(os.path.abspath(p), "w") as fh:
        json.dump(out, fh, indent=1)
    for k in ("gex", "tcr"):
        d = out[k]
        print(f"  {k}: n={d['n_reads']:,} len={d['read_length']} alphabet={''.join(d['alphabet'])}")


if __name__ == "__main__":
    main()
