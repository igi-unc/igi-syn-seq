#!/usr/bin/env python3
"""Self-consistency checks on a generated catalog.

    python3 validate_catalog.py --table output/IGI-SYN-SEQ-02.snv_indel.tsv --reference /path/GRCh38.fa
    python3 validate_catalog.py --table output/IGI-SYN-SEQ-02.snv_indel.tsv   # skips reference checks

Exits non-zero if any check fails, so it can be used as a regression gate.
"""
import argparse
import csv
import statistics
import sys

TIER_ORDER = ["T_amp", "T_LOH", "T_het", "A", "B", "A1"]


def check(rows, fa):
    fails = []

    def fail(msg):
        fails.append(msg)

    # 1. REF alleles agree with the reference genome
    if fa is not None:
        bad = [r for r in rows
               if fa.fetch(r["chrom"], int(r["pos"]) - 1, int(r["pos"]) - 1 + len(r["ref"])).upper() != r["ref"]]
        if bad:
            fail(f"{len(bad)} rows whose REF allele disagrees with the reference, e.g. "
                 f"{bad[0]['chrom']}:{bad[0]['pos']} {bad[0]['ref']}")

    # 2. the mutant haplotype window carries the ALT allele where the card marks it
    checked = off_bad = 0
    for r in rows:
        w = r.get("mut_hap_window") or ""
        if not w:
            continue
        checked += 1
        if w[20:20 + len(r["alt"])] != r["alt"]:
            off_bad += 1
    if off_bad:
        fail(f"{off_bad}/{checked} mutant windows do not carry ALT at the marked offset")

    # 3. window lengths track the indel size
    bad_len = [r for r in rows
               if r.get("mut_hap_window") and r.get("wt_hap_window")
               and len(r["mut_hap_window"]) != len(r["wt_hap_window"]) + len(r["alt"]) - len(r["ref"])]
    if bad_len:
        fail(f"{len(bad_len)} rows whose mutant window length is inconsistent with the indel")

    # 4. expected VAF is a probability
    bad_vaf = [r for r in rows if not 0.0 <= float(r["expected_vaf_tumor"]) <= 1.0]
    if bad_vaf:
        fail(f"{len(bad_vaf)} rows with an expected VAF outside [0, 1]")

    # 5. VAF decreases down the clonality ladder
    med = {}
    for r in rows:
        med.setdefault(r["clonality_tier"], []).append(float(r["expected_vaf_tumor"]))
    present = [t for t in TIER_ORDER if t in med]
    medians = [statistics.median(med[t]) for t in present]
    if medians != sorted(medians, reverse=True):
        fail("median expected VAF is not monotonically decreasing across "
             + " > ".join(present) + f" (got {[round(m, 3) for m in medians]})")

    # 6. a binding tier implies a peptide, and no tier implies none
    for r in rows:
        if r["binding_tier"] in ("strong", "weak", "non") and not r["best_peptide"]:
            fail(f"{r['event_id']}: binding tier {r['binding_tier']} without a peptide")
            break
        if r["binding_tier"] == "na" and r["best_peptide"]:
            fail(f"{r['event_id']}: tier 'na' but a peptide is recorded")
            break

    # 7. the recorded tier matches the recorded rank
    for r in rows:
        if not r["best_rank_el"]:
            continue
        rank = float(r["best_rank_el"])
        want = "strong" if rank <= 0.5 else "weak" if rank <= 2.0 else "non"
        if r["binding_tier"] != want:
            fail(f"{r['event_id']}: rank {rank} labelled {r['binding_tier']}, expected {want}")
            break

    # 8. a wild-type counterpart must differ from the mutant peptide and be the same length
    for r in rows:
        wt, mut = r.get("wt_peptide"), r.get("best_peptide")
        if wt and mut and (wt == mut or len(wt) != len(mut)):
            fail(f"{r['event_id']}: wild-type counterpart {wt} is not a same-length variant of {mut}")
            break

    return fails, present, medians


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--reference", help="GRCh38 FASTA; reference-allele checks are skipped without it")
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.table), delimiter="\t"))
    fa = None
    if a.reference:
        import pysam
        fa = pysam.FastaFile(a.reference)
    fails, tiers, medians = check(rows, fa)
    print(f"{a.table}: {len(rows)} rows")
    for t, m in zip(tiers, medians):
        print(f"  median expected VAF {t:6s} {m:.3f}")
    if fails:
        print("\nFAILED:")
        for f in fails:
            print(f"  - {f}")
        sys.exit(1)
    print("\nall checks passed")


if __name__ == "__main__":
    main()
