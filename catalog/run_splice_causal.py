#!/usr/bin/env python3
"""Turn the causal splice-site variants of tumour-specific splice events into a somatic edit table.

    python3 run_splice_causal.py --design design.yaml --paths paths.yaml \
        --dataset IGI-SYN-SEQ-01 --expressed output/IGI-SYN-SEQ-01.expressed.tsv --out output/

A tumour-specific splice event is defined by a somatic base change at a canonical donor or acceptor
dinucleotide. The designer records that change as a string in the expressed-class table, which no read
path consumes, so the junction it is supposed to explain had no genomic cause in the sequence. This
writes the same changes in the column layout `genome_build.read_events` accepts, and checks each one
against the reference before doing so.
"""
import argparse, csv, json, os, sys
from collections import Counter

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from igi_catalog.genome import Genome

CANONICAL = {"GT", "GC", "AT"}          # donor dinucleotides seen in human introns
CANONICAL_ACC = {"AG", "AC"}


def parse_variant(s):
    """'chr1:12345G>A' -> (chrom, pos, ref, alt)."""
    loc, change = s.split(":", 1)
    i = 0
    while i < len(change) and change[i].isdigit():
        i += 1
    pos = int(change[:i])
    ref, alt = change[i:].split(">", 1)
    return loc, pos, ref, alt


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", required=True)
    ap.add_argument("--paths", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--expressed", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    paths = yaml.safe_load(open(a.paths))
    genome = Genome(paths["reference_fasta"])
    rows, skipped, dinuc = [], [], Counter()
    for r in csv.DictReader(open(a.expressed), delimiter="\t"):
        if r.get("class") != "splice" or not r.get("causal_variant"):
            continue
        chrom, pos, ref, alt = parse_variant(r["causal_variant"])
        obs = genome.seq(chrom, pos - 1, pos)
        if obs != ref:
            skipped.append(f"{r['event_id']}: reference has {obs} at {chrom}:{pos}, table says {ref}")
            continue
        if alt == ref:
            skipped.append(f"{r['event_id']}: alt equals ref at {chrom}:{pos}")
            continue
        pair = genome.seq(chrom, pos - 1, pos + 1)
        dinuc[pair] += 1
        # the clone decides the timing exactly as it does for designed SNVs: truncal changes predate the
        # copy-number events, subclonal ones follow them
        clone = r["clone"]
        rows.append({
            "event_id": r["event_id"], "dataset": a.dataset, "class": "snv", "subclass": "splice_causal",
            "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
            "haplotype": int(r["haplotype"]), "clone": clone, "ccf": r["ccf"],
            "timing": "pre_cna" if clone == "T" else "post_cna",
            "clonality_tier": r["clonality_tier"],
            "gene": r["gene"], "transcript": r["transcript"], "mechanism": r["mechanism"],
            "exon_index": r["exon_index"], "consequence": "splice_site",
            "reference_dinucleotide": pair,
            "canonical_site": pair in CANONICAL or pair in CANONICAL_ACC,
            "chr1to6": r["chr1to6"],
        })
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, f"{a.dataset}.splice_causal.tsv")
    if rows:
        cols = list(rows[0])
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
    summ = {"dataset": a.dataset, "n": len(rows), "skipped": skipped,
            "reference_dinucleotides": dict(dinuc),
            "non_canonical": sum(1 for r in rows if not r["canonical_site"]),
            "by_clone": dict(Counter(r["clone"] for r in rows)),
            "by_mechanism": dict(Counter(r["mechanism"] for r in rows))}
    with open(os.path.join(a.out, f"{a.dataset}.splice_causal.summary.json"), "w") as fh:
        json.dump(summ, fh, indent=2)
    print(json.dumps(summ, indent=2))
    if skipped:
        raise SystemExit(f"{len(skipped)} causal variant(s) did not match the reference; see above")
    print(f"[done] {len(rows)} causal splice-site variants -> {path}")


if __name__ == "__main__":
    main()
