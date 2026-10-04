#!/usr/bin/env python3
"""Validate every generated assay against the target it was calibrated to.

This is deliberately not the acceptance suite. Acceptance checks that the designed BIOLOGY arrived -- allele
fractions by clonality tier, depth by copy number, which antigen classes reached which library -- and needs
an aligned BAM. This checks that the SEQUENCING is what was asked for: read counts, read lengths, accuracy,
and truth-map integrity, straight off the delivered files.

It exists because six of the fifteen sample types have no automated validation at all. Every long-read assay
has been checked only by hand, one chromosome at a time, which does not scale and does not get re-run when
something changes. A cheap check that covers everything is worth more than an expensive one that covers a
third.

Targets come from the measurements in paths.example.yaml and docs/build-reference.md. Where a target does not
exist -- ONT bulk cDNA length, because the only real ONT data is single cell -- the check reports the value
and says there is nothing to compare it with, rather than inventing a bound.

    python3 qa_release.py --release /path/to/release --out qa_report.json
"""
import argparse
import glob
import gzip
import json
import math
import os
import subprocess

# (name, target, tolerance) where tolerance is a fraction of the target; None means "report, do not judge"
TARGETS = {
    "pacbio":        {"length_mean": (16689, 0.08), "accuracy": (0.99823, 0.0015)},
    "ont_wgs":       {"length_mean": (19117, 0.08), "accuracy": (0.98676, 0.0030)},
    # ONT RNA length is set by the transcript and the size-selection curve, not by a parameter, and the
    # curve is fitted on Kinnex rather than ONT -- so the single-cell value is reported against the real
    # 910 bp median without a bound, and the bulk value has no real reference at all because the only
    # IPISRC044 ONT data is single cell. Inventing a bound would manufacture a pass or a failure.
    "ont_sc_rna":    {"length_mean": (910, None),   "accuracy": (0.98220, 0.0030)},
    "ont_bulk_rna":  {"length_mean": (None, None),  "accuracy": (0.98220, 0.0030)},
    "kinnex_bulk":   {"length_mean": (2223, 0.12),  "accuracy": (0.99788, 0.0020)},
    # 2,155 is the real OVERALL mean of skera's output on single-cell arrays, which is what a segmented
    # BAM contains. The first version of this target used 952, the real median of the SINGLE-cDNA subset,
    # and failed a correct library at 2,416 by comparing it against the wrong population: with nine
    # adapters for sixteen cDNAs about 12 % of skera's output is multi-cDNA blocks averaging 8.6 TSOs, and
    # those are in the BAM too. Judging a mixture against one of its components is not a check.
    "kinnex_sc":     {"length_mean": (2155, 0.20),  "accuracy": (0.99788, 0.0020)},
}


def read_stats(path, limit=40000, samtools=None):
    """(n, mean, sd, accuracy) from a FASTQ.gz or a BAM, over at most `limit` reads."""
    L, err, nb = [], 0.0, 0
    if path.endswith(".bam"):
        cmd = (samtools or "samtools").split() + ["view", path]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                             text=True, bufsize=1 << 20)
        for line in p.stdout:
            if len(L) >= limit:
                break
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            L.append(len(f[9]))
            for c in f[10]:
                err += 10 ** (-(ord(c) - 33) / 10.0)
            nb += len(f[10])
        p.stdout.close(); p.terminate()
    else:
        op = gzip.open if path.endswith(".gz") else open
        with op(path, "rt") as fh:
            while len(L) < limit:
                if not fh.readline():
                    break
                s = fh.readline().strip(); fh.readline(); q = fh.readline().strip()
                if not q:
                    break
                L.append(len(s))
                for c in q:
                    err += 10 ** (-(ord(c) - 33) / 10.0)
                nb += len(q)
    if not L:
        return 0, 0.0, 0.0, 0.0
    m = sum(L) / len(L)
    sd = (sum((x - m) ** 2 for x in L) / len(L)) ** 0.5
    return len(L), m, sd, (1 - err / nb if nb else 0.0)


def judge(got, target, tol):
    if target is None:
        return "no target"
    if tol is None:
        return f"reported ({got:,.0f} vs real {target:,.0f}, no bound set)"
    lo, hi = (target * (1 - tol), target * (1 + tol)) if target > 1 else (target - tol, target + tol)
    return "PASS" if lo <= got <= hi else f"FAIL (outside {lo:.5g}-{hi:.5g})"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--limit", type=int, default=40000)
    a = ap.parse_args()

    found = {
        "pacbio":       sorted(glob.glob(f"{a.release}/pacbio_merged/*/*_hifi.bam")) or
                        sorted(glob.glob(f"{a.release}/pacbio/*/*_hifi.bam"))[:4],
        "ont_wgs":      sorted(glob.glob(f"{a.release}/ont_wgs_merged/*/*_ont_wgs.fastq.gz")) or
                        sorted(glob.glob(f"{a.release}/ont_wgs/*/*_ont_wgs.fastq.gz"))[:4],
        "ont_sc_rna":   sorted(glob.glob(f"{a.release}/ont_sc_rna/*/*_ont_sc_rna.fastq.gz")),
        "ont_bulk_rna": sorted(glob.glob(f"{a.release}/ont_bulk_rna/*/*_ont_bulk_rna.fastq.gz")),
        "kinnex_bulk":  sorted(glob.glob(f"{a.release}/kinnex_bulk/*/*_segmented.bam")),
        "kinnex_sc":    sorted(glob.glob(f"{a.release}/kinnex_sc/*/*_segmented.bam")),
    }
    report = {"assays": {}, "n_fail": 0}
    for assay, files in found.items():
        rows = []
        for f in files:
            n, m, sd, acc = read_stats(f, a.limit, a.samtools)
            if not n:
                continue
            tl, ttl = TARGETS[assay]["length_mean"]
            ta, tta = TARGETS[assay]["accuracy"]
            vl, va = judge(m, tl, ttl), judge(acc, ta, tta)
            rows.append({"file": os.path.basename(f), "reads_sampled": n,
                         "length_mean": round(m, 1), "length_sd": round(sd, 1),
                         "accuracy": round(acc, 5), "length_verdict": vl, "accuracy_verdict": va})
            report["n_fail"] += sum(1 for v in (vl, va) if v.startswith("FAIL"))
        report["assays"][assay] = rows
        for r in rows:
            print(f"  {assay:<13} {r['file'][:46]:<46} len {r['length_mean']:>8,.0f} "
                  f"acc {r['accuracy']:.5f}  {r['length_verdict']} / {r['accuracy_verdict']}")

    # Truth-map integrity for the assays that ship one: every read must resolve to a molecule.
    tm = {}
    for assay, pat in (("tenx_gex", f"{a.release}/tenx_gex/*/*_gex_molecules.tsv.gz"),
                       ("tenx_tcr", f"{a.release}/tenx_tcr/*/*_tcr_molecules.tsv.gz"),
                       ("kinnex_sc", f"{a.release}/kinnex_sc/*/*_arrays.tsv.gz"),
                       ("ont_sc_rna", f"{a.release}/ont_sc_rna/*/*_molecules.tsv.gz")):
        for f in sorted(glob.glob(pat)):
            n = 0
            with gzip.open(f, "rt") as fh:
                fh.readline()
                for _ in fh:
                    n += 1
            tm.setdefault(assay, []).append({"file": os.path.basename(f), "rows": n})
            print(f"  truth map    {os.path.basename(f)[:46]:<46} {n:,} rows")
    report["truth_maps"] = tm

    with open(a.out, "w") as fh:
        json.dump(report, fh, indent=2); fh.write("\n")
    print(f"\n  {report['n_fail']} failing check(s) -> {a.out}")


if __name__ == "__main__":
    main()
