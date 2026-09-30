#!/usr/bin/env python3
"""Left-align the indel coordinates in existing truth tables (F14).

    python3 normalise_truth.py --paths paths.yaml --tables output/*.snv_indel.tsv output/*.background.tsv

A caller's VCF is normalised and a truth table built from the designer's own coordinates is not. Inside a
repeat the same deletion has several equally valid representations, so an unnormalised row does not just
fail to match: it is scored as a false negative for the event that is really there and a false positive
for the call that found it. About 27 % of the designed and background indels needed shifting.

The shift is a change of representation, not of sequence. The builder applies whichever form the table
carries and produces identical bases either way, so a table normalised after the fact still describes the
reads that were generated from it. Every row is checked against the reference before and after.
"""
import argparse, csv, os, shutil, sys

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from igi_catalog.genome import Genome, left_align


def normalise(path, genome, dry_run=False):
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols, rows = reader.fieldnames, list(reader)
    n_indel = n_shift = n_bad = 0
    for r in rows:
        ref, alt = r.get("ref", ""), r.get("alt", "")
        if not ref or not alt or len(ref) == len(alt) == 1:
            continue
        n_indel += 1
        pos, chrom = int(r["pos"]), r["chrom"]
        if genome.seq(chrom, pos - 1, pos - 1 + len(ref)) != ref:
            n_bad += 1
            continue
        np_, nr, na = left_align(genome, chrom, pos, ref, alt)
        if (np_, nr, na) == (pos, ref, alt):
            continue
        # the shifted form must describe the same edit: same reference bases, same resulting sequence
        assert genome.seq(chrom, np_ - 1, np_ - 1 + len(nr)) == nr, f"{r['event_id']}: shifted ref mismatch"
        before = genome.seq(chrom, np_ - 6, pos - 1 + len(ref) + 5)
        off_old, off_new = pos - (np_ - 5), np_ - (np_ - 5)
        seq_old = before[:off_old] + alt + before[off_old + len(ref):]
        seq_new = before[:off_new] + na + before[off_new + len(nr):]
        assert seq_old == seq_new, f"{r['event_id']}: shift changes the sequence"
        r["pos"], r["ref"], r["alt"] = str(np_), nr, na
        n_shift += 1
    if n_shift and not dry_run:
        shutil.copy2(path, path + ".prenorm")
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
    return n_indel, n_shift, n_bad


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--paths", required=True)
    ap.add_argument("--tables", nargs="+", required=True)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)
    genome = Genome(yaml.safe_load(open(a.paths))["reference_fasta"])
    total = 0
    for path in a.tables:
        if not os.path.exists(path):
            continue
        n_indel, n_shift, n_bad = normalise(path, genome, a.dry_run)
        total += n_shift
        flag = " DRY RUN" if a.dry_run else ""
        print(f"{os.path.basename(path)}: {n_indel} indels, {n_shift} left-aligned"
              + (f", {n_bad} reference mismatch" if n_bad else "") + flag)
    print(f"[done] {total} row(s) normalised")


if __name__ == "__main__":
    main()
