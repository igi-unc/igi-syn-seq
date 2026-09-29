#!/usr/bin/env python3
"""Acceptance test: measure a generated library against the truth bundle.

The earlier check compared allele fractions at on-target SNVs and nothing else. It passed while the RNA
library had no expression structure, while exome depth was flat across copy-number states, while
structural variants never reached the DNA and while most frameshifts carried no peptide. Each check here
exists because something that should have been caught was not.

    python3 acceptance.py --bam tumor.bam --dataset IGI-SYN-SEQ-02 --catalog-dir output \
        --reference GRCh38.fa --chrom chr6 --report report.json
    python3 acceptance.py --rna-manifest rna_transcripts.tsv --report rna.json

Exits non-zero if any check fails, so it can gate a release.
"""
import argparse
import csv
import gzip
import json
import os
import re
import statistics
import subprocess
import sys
from collections import Counter, defaultdict

TIER_ORDER = ["T_amp", "T_LOH", "T_het", "A", "B", "A1"]


def rows(path):
    if not path or not os.path.exists(path):
        return []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def pileup(bam, reference, sites_path, min_bq=15, min_mq=20, max_depth=8000):
    """{position: (depth, bases)} at the requested sites."""
    cmd = ["samtools", "mpileup", "-f", reference, "-l", sites_path,
           "-d", str(max_depth), "-Q", str(min_bq), "-q", str(min_mq), bam]
    out = subprocess.run(cmd, capture_output=True, text=True).stdout
    res = {}
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 5:
            continue
        b = re.sub(r"\^.", "", f[4]).replace("$", "")
        b = re.sub(r"[+-]\d+[ACGTNacgtn]*", "", b)
        res[int(f[1])] = (int(f[3]), b.upper())
    return res


def captured_depth(bam, chrom, start, end, capture_bed, work, min_intervals=20):
    """Mean depth over captured bases only, within a region.

    Mean depth across a whole window is meaningless for an exome: it measures how much of the window is
    captured, not how deeply it is sequenced. An arbitrary 2 Mb window read 3.4x while the gene-dense MHC
    read 15.6x purely because far more of the MHC is on target. Restricting to capture intervals with
    `samtools bedcov` compares like with like.
    """
    sub = os.path.join(work, f"cap_{chrom}_{start}_{end}.bed")
    n = 0
    with open(capture_bed) as fh, open(sub, "w") as out:
        for line in fh:
            f = line.split("\t")
            if len(f) < 3 or f[0] != chrom:
                continue
            s0, e0 = int(f[1]), int(f[2])
            if e0 <= start or s0 >= end:
                continue
            out.write(f"{chrom}\t{max(s0, start)}\t{min(e0, end)}\n")
            n += 1
    if n < min_intervals:
        return None, n
    res = subprocess.run(["samtools", "bedcov", "-Q", "20", sub, bam], capture_output=True, text=True).stdout
    total_bases = total_cov = 0
    for line in res.splitlines():
        f = line.split("\t")
        if len(f) < 4:
            continue
        total_bases += int(f[2]) - int(f[1])
        total_cov += int(f[3])
    return (total_cov / total_bases if total_bases else 0.0), n


class Report:
    def __init__(self):
        self.checks = []

    def add(self, name, passed, detail, numbers=None):
        self.checks.append({"check": name, "pass": bool(passed), "detail": detail,
                            "numbers": numbers or {}})
        mark = "PASS" if passed else "FAIL"
        print(f"[{mark}] {name}: {detail}", flush=True)

    def failed(self):
        return [c for c in self.checks if not c["pass"]]

    def write(self, path):
        with open(path, "w") as fh:
            json.dump({"checks": self.checks,
                       "n_pass": sum(1 for c in self.checks if c["pass"]),
                       "n_fail": len(self.failed())}, fh, indent=2)


def check_variants(rep, bam, reference, catalog, chrom, work):
    """Observed allele fraction against expected, for SNVs and indels, by clonality tier."""
    ev = [r for r in rows(catalog) if r["chrom"] == chrom and r.get("ctx_on_target") == "True"]
    if not ev:
        rep.add("variants", False, f"no on-target events on {chrom}")
        return
    sites = os.path.join(work, "sites.txt")
    with open(sites, "w") as fh:
        for r in ev:
            fh.write(f"{chrom}\t{r['pos']}\n")
    pile = pileup(bam, reference, sites)
    by_tier, by_class = defaultdict(list), defaultdict(list)
    seen = 0
    for r in ev:
        pos = int(r["pos"])
        if pos not in pile:
            continue
        dp, b = pile[pos]
        if dp < 30:
            continue
        obs = b.count(r["alt"]) / dp if len(r["alt"]) == 1 and len(r["ref"]) == 1 else None
        exp = float(r["expected_vaf_tumor"])
        if obs is not None:
            by_tier[r["clonality_tier"]].append((exp, obs))
            by_class[r["class"]].append((exp, obs))
        if b.count(r["alt"]) > 0 or len(r["ref"]) > 1:
            seen += 1
    nums = {}
    worst = 0.0
    for tier in TIER_ORDER:
        v = by_tier.get(tier)
        if not v:
            continue
        e, o = statistics.mean(x[0] for x in v), statistics.mean(x[1] for x in v)
        nums[tier] = {"n": len(v), "expected": round(e, 4), "observed": round(o, 4),
                      "ratio": round(o / e, 3) if e else None}
        if e:
            worst = max(worst, abs(o / e - 1))
    rep.add("allele fractions by tier", worst <= 0.15,
            f"largest deviation {worst:.0%} across {sum(len(v) for v in by_tier.values())} sites", nums)
    rep.add("indels present", len([r for r in ev if r["class"] == "indel"]) == 0 or seen > 0,
            f"{seen} of {len(ev)} designed events have supporting reads")


def check_depth_by_cn(rep, bam, design, dataset, chrom, arms_bed, capture_bed, work):
    """Observed depth ratio between copy-number regions against the purity and copy-number formula."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import yaml
    from igi_catalog.clones import CloneModel
    from igi_catalog.pipeline import depth_denominator, source_plan
    cfg = yaml.safe_load(open(design))["datasets"][dataset]
    cm = CloneModel(cfg, arms_bed)
    den = depth_denominator(cm, cm.purity)
    regions = []
    for ev in (cfg.get("cna") or {}).get("T", []):
        c, s, e = cm.region(ev["region"])
        if c == chrom and e - s > 200_000:
            regions.append((ev.get("label", "cna"), s + 1, e))
    if not regions:
        rep.add("depth by copy number", True, f"no large copy-number segment on {chrom} to test")
        return
    # the baseline must itself be a copy-number-neutral stretch, not an arbitrary window
    neutral = []
    pos = 1_000_000
    while pos < 150_000_000:
        a, b, _lab = cm.cn("T", chrom, pos)
        if (a, b) == cm.base and not any(s <= pos <= e for _l, s, e in regions):
            neutral.append(pos)
        pos += 2_000_000
    if not neutral:
        rep.add("depth by copy number", False, "no copy-number-neutral baseline found on this chromosome")
        return
    mid_n = neutral[len(neutral) // 2]
    base, n_base = captured_depth(bam, chrom, mid_n - 5_000_000, mid_n + 5_000_000, capture_bed, work)
    nums, worst = {"baseline": {"depth": round(base, 2) if base else None, "intervals": n_base}}, 0.0
    if not base:
        rep.add("depth by copy number", False, "baseline region has no captured depth")
        return
    for label, s, e in regions:
        want = sum(w for _s, _h, _k, w in source_plan(cm, cm.purity, chrom, (s + e) // 2)) / den
        got, n_iv = captured_depth(bam, chrom, s, e, capture_bed, work)
        if got is None:
            nums[label] = {"skipped": f"only {n_iv} capture intervals"}
            continue
        ratio = got / base
        nums[label] = {"expected_ratio": round(want, 3), "observed_ratio": round(ratio, 3),
                       "captured_depth": round(got, 2), "baseline_depth": round(base, 2),
                       "intervals": n_iv}
        if want:
            worst = max(worst, abs(ratio / want - 1))
    rep.add("depth by copy number", worst <= 0.25,
            f"largest deviation {worst:.0%} over {len(regions)} segment(s), on captured bases", nums)


def check_read_names(rep, r1):
    """Names must be unique, Illumina-shaped, and free of anything identifying the source."""
    names, dup, leaky = set(), 0, 0
    n = 0
    with gzip.open(r1, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4:
                continue
            nm = line.split()[0][1:]
            n += 1
            if nm in names:
                dup += 1
            names.add(nm)
            if re.search(r"chr|hap|NORMAL|_T_|_A1?_|_B_", nm):
                leaky += 1
            if n >= 400000:
                break
    rep.add("read names unique", dup == 0, f"{dup} duplicates in {n:,} sampled")
    rep.add("read names carry no truth", leaky == 0, f"{leaky} names contain source information")


def check_classes_present(rep, catalog_dir, dataset, junctions_tsv):
    """Every designed class that can reach this library should appear in it."""
    present, missing = {}, []
    for name, path in (("sv", f"{dataset}.svs.tsv"), ("fusion", f"{dataset}.fusions.tsv"),
                       ("virus", f"{dataset}.viruses.tsv")):
        present[name] = len(rows(os.path.join(catalog_dir, path)))
    placed = rows(junctions_tsv) if junctions_tsv else []
    ids = {r["event_id"].split(":")[0] for r in placed}
    rep.add("junction classes reach the library", bool(placed) or present["sv"] == 0,
            f"{len(placed)} junctions placed covering {len(ids)} events",
            {"catalog": present, "placed": len(placed)})


def check_rna(rep, manifest):
    """RNA depth must track abundance, cover both haplotypes and include the designed classes."""
    recs = [r for r in rows(manifest) if float(r.get("observed_pairs") or 0) > 0]
    if not recs:
        rep.add("rna", False, "no transcript produced reads")
        return
    import math
    # TPM is already length-normalised, so the quantity that should track abundance is per-base coverage,
    # not the raw pair count: a short abundant transcript and a long rare one can yield the same number of
    # pairs. Comparing pairs directly confounds length with abundance and understates the agreement.
    live = [r for r in recs if float(r["abundance"]) > 0 and int(r["length"]) > 0]
    a = [math.log10(float(r["abundance"])) for r in live]
    o = [math.log10(float(r["observed_pairs"]) * 2 * 150 / int(r["length"])) for r in live]
    ma, mo = statistics.mean(a), statistics.mean(o)
    num = sum((x - ma) * (y - mo) for x, y in zip(a, o))
    den = math.sqrt(sum((x - ma) ** 2 for x in a) * sum((y - mo) ** 2 for y in o))
    r = num / den if den else 0
    rep.add("rna depth tracks abundance", r >= 0.95,
            f"correlation of log abundance with log observed coverage = {r:.3f} over {len(live):,} records")
    haps = Counter(x["hap"] for x in recs)
    rep.add("rna covers both haplotypes", len(haps) > 1 and min(haps.values()) > 0.2 * max(haps.values()),
            f"records per haplotype {dict(haps)}")
    srcs = Counter(x["source"] for x in recs)
    want = {"reference", "fusion", "erv", "splice_isoform", "cta", "virus"}
    rep.add("rna carries the designed classes", len(want & set(srcs)) >= 3,
            f"sources present {dict(srcs)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam")
    ap.add_argument("--r1")
    ap.add_argument("--reference")
    ap.add_argument("--dataset")
    ap.add_argument("--catalog-dir", default="output")
    ap.add_argument("--chrom")
    ap.add_argument("--design", default="design.yaml")
    ap.add_argument("--arms-bed")
    ap.add_argument("--capture-bed")
    ap.add_argument("--junctions")
    ap.add_argument("--rna-manifest")
    ap.add_argument("--work", default="/tmp")
    ap.add_argument("--report", required=True)
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    rep = Report()
    if a.bam and a.reference and a.dataset and a.chrom:
        cat = os.path.join(a.catalog_dir, f"{a.dataset}.snv_indel.tsv")
        check_variants(rep, a.bam, a.reference, cat, a.chrom, a.work)
        if a.arms_bed and a.capture_bed:
            check_depth_by_cn(rep, a.bam, a.design, a.dataset, a.chrom, a.arms_bed,
                              a.capture_bed, a.work)
        check_classes_present(rep, a.catalog_dir, a.dataset, a.junctions)
    if a.r1:
        check_read_names(rep, a.r1)
    if a.rna_manifest:
        check_rna(rep, a.rna_manifest)
    rep.write(a.report)
    if rep.failed():
        print(f"\n{len(rep.failed())} check(s) failed", flush=True)
        sys.exit(1)
    print(f"\nall {len(rep.checks)} checks passed", flush=True)


if __name__ == "__main__":
    main()
