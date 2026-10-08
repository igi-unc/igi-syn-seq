#!/usr/bin/env python3
"""Measure a simulated long-read library's TRUE identity from its alignment, and compare it with what
the reads claim about themselves.

    alignment_identity.py --bam aligned.bam --out out.json --label name [--max-reads 20000]

Why this exists. Every long-read identity setting in this project was calibrated with
jobs/measure/measure_longread.py, which derives accuracy from the base-quality string (and from the `rq`
tag for real HiFi BAMs). For REAL reads that is the right measurement: the quality string is the
instrument's own estimate and is reasonably calibrated. For SIMULATED reads it is not a measurement at
all. Badread generates sequence from an error model and qualities from a SEPARATE qscore model, so the
quality string is the simulator's claim about the read, not a property of it. The two can disagree, and
for the PacBio model they disagree by a factor of five: reads whose qualities claim Q27.8 carry Q20.4 of
actual error (F24, CHANNEL entry 23).

Reads simulated from an unmodified reference make the alignment ground truth -- there are no germline or
designed variants to confuse with error -- so NM over the aligned length IS the error rate.

Reported per library:

    true_error_per_bp      NM summed over primary alignments / aligned length (M+I+D)
    claimed_error_per_bp   mean error PROBABILITY from the quality strings, converted once at the end
    ratio                  true / claimed; 1.0 means the reads' qualities tell the truth
    indel_fraction         share of the true error that is inserted or deleted bases

Accuracy is averaged as an error PROBABILITY and converted to Phred once. Phred is a logarithm, so a mean
of Phred scores is not the Phred of the mean error rate.
"""
import argparse
import json
import math

import pysam

CONSUMES_QUERY = {0, 1, 4, 7, 8}     # M I S = X
CONSUMES_REF = {0, 2, 3, 7, 8}       # M D N = X


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--label", default="library")
    ap.add_argument("--max-reads", type=int, default=20000)
    a = ap.parse_args()

    nm = aligned = ins = dele = matches = 0
    q_err = 0.0
    q_bases = 0
    n = n_skipped = clipped = 0
    bam = pysam.AlignmentFile(a.bam, "rb")
    for rec in bam.fetch(until_eof=True):
        if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
            n_skipped += 1
            continue
        if not rec.has_tag("NM"):
            raise SystemExit("alignment has no NM tag; minimap2 writes one, so this BAM is not what "
                             "this measurement expects")
        n += 1
        nm += rec.get_tag("NM")
        for op, length in rec.cigartuples:
            if op in (0, 7, 8):
                aligned += length
                matches += length
            elif op == 1:
                aligned += length
                ins += length
            elif op == 2:
                aligned += length
                dele += length
            elif op in (4, 5):
                clipped += length
        q = rec.query_qualities
        if q is not None:
            # error probability per base, summed; converted once below
            q_err += sum(10.0 ** (-v / 10.0) for v in q)
            q_bases += len(q)
        if n >= a.max_reads:
            break
    bam.close()

    if not n:
        raise SystemExit("no primary alignments")
    true_err = nm / aligned
    claimed_err = (q_err / q_bases) if q_bases else None
    out = {
        "label": a.label,
        "reads": n,
        "secondary_or_unmapped_skipped": n_skipped,
        "aligned_bases": aligned,
        "soft_hard_clipped_bases": clipped,
        "clipped_fraction": round(clipped / (aligned + clipped), 5) if aligned + clipped else None,
        "nm": nm,
        "true_error_per_bp": round(true_err, 8),
        "true_identity": round(1 - true_err, 6),
        "true_q": round(-10 * math.log10(true_err), 2) if true_err > 0 else None,
        "indel_bases": ins + dele,
        "inserted": ins,
        "deleted": dele,
        # NM counts mismatches plus inserted plus deleted bases, so the indel share of the error is
        # (I+D)/NM. 91 % for the delivered HiFi library, which is what makes the overstatement matter
        # for a homopolymer tier.
        "indel_fraction_of_error": round((ins + dele) / nm, 4) if nm else None,
    }
    if claimed_err is not None:
        out.update(
            claimed_error_per_bp=round(claimed_err, 8),
            claimed_identity=round(1 - claimed_err, 6),
            claimed_q=round(-10 * math.log10(claimed_err), 2) if claimed_err > 0 else None,
            ratio_true_over_claimed=round(true_err / claimed_err, 3) if claimed_err > 0 else None,
        )
    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=2)
        fh.write("\n")
    print(json.dumps(out, indent=2), flush=True)


if __name__ == "__main__":
    main()
