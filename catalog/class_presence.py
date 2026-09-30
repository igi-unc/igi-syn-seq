#!/usr/bin/env python3
"""Per-class presence report: what the catalog designed against what reached the reads.

    python3 class_presence.py --dataset IGI-SYN-SEQ-01 --catalog-dir output \
        --records merged/..._tumor_records.tsv.gz --junctions merged/..._tumor_junctions.tsv \
        --rna-manifest release/rna/..._full_rna_transcripts.tsv --out report.json

Owner decision D7 condition 2. "The class is present" is not evidence the path works: a class present at
one record out of thirty says the path is broken, not working, and an aggregate pass hides that. So each
class is reported as designed / reachable / observed, where reachable means the assay can physically see
it at all -- a deep intronic variant is not reachable in an exome, and an off-capture viral integration is
not reachable in either, per D8.
"""
import argparse, csv, gzip, json, os
from collections import Counter, defaultdict


def rows(path):
    if not path or not os.path.exists(path):
        return []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--catalog-dir", default="output")
    ap.add_argument("--records", help="merged record map from the exome build")
    ap.add_argument("--junctions", help="merged junctions table from the exome build")
    ap.add_argument("--rna-manifest", help="RNA transcript manifest with observed_pairs")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    cd, ds = a.catalog_dir, a.dataset

    cat = {
        "snv": [r for r in rows(f"{cd}/{ds}.snv_indel.tsv") if r.get("class") == "snv"],
        "indel": [r for r in rows(f"{cd}/{ds}.snv_indel.tsv") if r.get("class") == "indel"],
        "fusion": rows(f"{cd}/{ds}.fusions.tsv"),
        "sv": rows(f"{cd}/{ds}.svs.tsv"),
        "virus": rows(f"{cd}/{ds}.viruses.tsv"),
        "cta": [r for r in rows(f"{cd}/{ds}.expressed.tsv") if r.get("class") == "cta"],
        "erv": [r for r in rows(f"{cd}/{ds}.expressed.tsv") if r.get("class") == "erv"],
        "splice": [r for r in rows(f"{cd}/{ds}.expressed.tsv") if r.get("class") == "splice"],
        "background": rows(f"{cd}/{ds}.background.tsv"),
        "splice_causal": rows(f"{cd}/{ds}.splice_causal.tsv"),
    }

    exome = {}
    placed = {r["event_id"].split(":")[0] for r in rows(a.junctions)}
    for cls, rws in cat.items():
        designed = len(rws)
        if cls in ("snv", "indel", "background", "splice_causal"):
            # small variants are reachable in the exome when on target; the acceptance per-chromosome
            # checks measure whether they are observed, so this reports reachability only
            reach = sum(1 for r in rws if str(r.get("ctx_on_target", "True")) == "True")
            exome[cls] = {"designed": designed, "reachable": reach, "observed": None,
                          "note": "observed measured per chromosome by acceptance.py"}
        elif cls in ("fusion", "sv", "virus"):
            ids = {r["event_id"] for r in rws}
            exome[cls] = {"designed": designed, "reachable": None,
                          "observed": len(ids & placed),
                          "note": "reaches the exome only as a junction contig at a captured breakpoint"}
        else:
            exome[cls] = {"designed": designed, "reachable": 0, "observed": 0,
                          "note": "expression class; not visible in DNA except via its own variants"}
    if cat["virus"]:
        exome["virus"]["note"] += "; D8: the integration site is off-capture by design, so 0 is correct"

    rna = {}
    man = rows(a.rna_manifest)
    if man:
        bysrc = defaultdict(lambda: [0, 0])
        for r in man:
            got = float(r.get("observed_pairs") or 0) > 0
            bysrc[r["source"]][0] += 1
            bysrc[r["source"]][1] += 1 if got else 0
        want = {"reference": None, "fusion": len(cat["fusion"]), "erv": len(cat["erv"]),
                "splice_isoform": len(cat["splice"]), "cta": len(cat["cta"]),
                "virus": len(cat["virus"])}
        for src, designed in want.items():
            n, obs = bysrc.get(src, [0, 0])
            rna[src] = {"designed_events": designed, "records": n, "records_with_reads": obs,
                        "all_records_produced_reads": (n > 0 and obs == n)}

    out = {"dataset": ds, "exome": exome, "rna": rna,
           "verdict": {
               "exome_junction_classes_observed": sum(1 for c in ("fusion", "sv", "virus")
                                                      if (exome[c]["observed"] or 0) > 0),
               "rna_classes_with_every_record_read": sorted(
                   k for k, v in rna.items() if v["all_records_produced_reads"]),
               "rna_classes_missing_entirely": sorted(k for k, v in rna.items() if v["records"] == 0),
           }}
    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
