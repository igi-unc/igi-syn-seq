#!/usr/bin/env python3
"""Measure read length and accuracy from a real long-read BAM or FASTQ.

This produces the `pacbio_*` and `ont_*` length and identity constants in the site paths file, and is
also how a simulated output is checked against its target.

Two traps are handled here deliberately:

1. **Accuracy must be averaged as error probability, not as Phred.** Phred is a logarithm, so the mean of
   Phred scores is not the Phred of the mean error rate -- Jensen's inequality. Averaging Phred scores
   and converting once gave 98.92 % where the correct answer is 92.15 %. This script sums 10**(-Q/10)
   over every base and converts at the end.
2. **Badread's realized identity is offset from its requested identity, and the offset is not constant:**
   about +0.43 points for PacBio and -0.22 for ONT cDNA. So the measured value here is a *target*, and the
   value written into the paths file is the target corrected by the offset for that model. Re-measure the
   simulated output with this same script to confirm where it actually landed.

PacBio HiFi BAMs carry a per-read `rq` tag, which is the instrument's own accuracy estimate and is the
right figure to compare against for HiFi. Base qualities are reported alongside it; for ONT there is no
`rq`, so base qualities are all there is.

    python3 measure_longread.py --bam reads.bam --label "HG002 Revio HiFi" --n 200000
    python3 measure_longread.py --fastq sim_R1.fastq.gz --label "simulated chr21"
"""
import argparse
import gzip
import math
import shlex
import subprocess
import sys


def stats(lengths, err_sum, n_bases, rq_sum, rq_n, label):
    if not lengths:
        print(f"{label}: no reads read")
        return
    lengths.sort()
    n = len(lengths)
    mean = sum(lengths) / n

    def pct(p):
        return lengths[min(n - 1, int(p * n))]

    sd = (sum((x - mean) ** 2 for x in lengths) / n) ** 0.5
    print(f"=== {label}")
    print(f"  reads              {n:,}")
    print(f"  length mean        {mean:.0f}")
    print(f"  length sd          {sd:.0f}")
    print(f"  length min/med/max {lengths[0]} / {pct(0.5)} / {lengths[-1]}")
    print(f"  length p10/p90     {pct(0.10)} / {pct(0.90)}")
    if n_bases:
        # mean error PROBABILITY, then convert once -- see the module docstring
        e = err_sum / n_bases
        q = -10 * math.log10(e) if e > 0 else float("inf")
        print(f"  basecall accuracy  {1 - e:.5f}  (Q{q:.1f}, from {n_bases:,} bases)")
    if rq_n:
        print(f"  read quality rq    {rq_sum / rq_n:.5f}  (instrument estimate, n={rq_n:,})")


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--bam")
    g.add_argument("--fastq")
    ap.add_argument("--label", default="reads")
    ap.add_argument("--n", type=int, default=200000, help="reads to sample")
    ap.add_argument("--samtools", default="samtools",
                    help="samtools command; may be a multi-word wrapper such as "
                         "'singularity exec -B /mnt img.sif samtools'")
    a = ap.parse_args()

    lengths, err_sum, n_bases, rq_sum, rq_n = [], 0.0, 0, 0.0, 0

    if a.bam:
        p = subprocess.Popen(shlex.split(a.samtools) + ["view", a.bam], stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL, text=True, bufsize=1 << 20)
        for line in p.stdout:
            if len(lengths) >= a.n:
                break
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            lengths.append(len(f[9]))
            for ch in f[10]:
                err_sum += 10 ** (-(ord(ch) - 33) / 10.0)
            n_bases += len(f[10])
            for t in f[11:]:
                if t.startswith("rq:f:"):
                    rq_sum += float(t[5:])
                    rq_n += 1
                    break
        p.stdout.close()
        p.terminate()
    else:
        op = gzip.open if a.fastq.endswith(".gz") else open
        with op(a.fastq, "rt") as fh:
            while len(lengths) < a.n:
                if not fh.readline():
                    break
                seq = fh.readline().rstrip("\n")
                fh.readline()
                qual = fh.readline().rstrip("\n")
                if not qual:
                    break
                lengths.append(len(seq))
                for ch in qual:
                    err_sum += 10 ** (-(ord(ch) - 33) / 10.0)
                n_bases += len(qual)

    stats(lengths, err_sum, n_bases, rq_sum, rq_n, a.label)


if __name__ == "__main__":
    main()
