#!/usr/bin/env python3
"""Designed allele fraction against the allele fraction the generator actually wrote.

    python3 truth_af.py --dataset IGI-SYN-SEQ-02 --catalog-dir output --design design.yaml \
        --records merged/..._tumor_records.tsv.gz --readmap merged/..._tumor_readmap.tsv.gz \
        --report truth_af.json [--chrom chr5,chr6]

Every other allele-fraction measurement in this release reads a BAM, which answers a different question
than the one it appears to answer. Two failed checks on IGI-SYN-SEQ-02's WGS arm are what forced the
distinction: `allele fractions by tier` compared the design against a samtools pileup of a bwa alignment
and reported the tumour het tier 13 to 15 % low, which reads as a statement about the release and is
really a statement about the whole path. The release was correct. Measured here, the delivered reads sat
within 0.2 % of the design in aggregate and no site was off by more than 3 %; the deficit appeared during
alignment and pileup, where roughly an eighth of the alt-carrying reads at those sites are not recovered.

This check takes the aligner out of the path. The truth map is the release's own account of what it put
in each library: records.tsv.gz lists every simulated interval copy as (chrom, start, end, population,
haplotype, copy_kind) and readmap.tsv.gz maps every read to the copy it came from. For a site, the copies
that carry its alt allele are those of the populations descended from the variant's clone, on the
variant's haplotype, and -- when the variant followed the copy-number change -- only the copy that
postdates it. Reads from those copies over reads from every copy of the interval is the allele fraction
in the library, exactly, with no alignment, no mapping quality and no base quality in the way.

What each of the two checks is good for:

  * this one fails when the GENERATOR and the design disagree, which is a release defect.
  * the pileup check fails when a CONSUMER cannot recover what the generator wrote, which is worth
    knowing and is not a defect. Its bound was never drawn with that distinction in mind, so it is the
    one that will go red on a correct release.

Each copy of an interval shares its coordinates with the others, so read counts are directly comparable
and any site covered by more than one interval is reported rather than silently averaged.
"""
import argparse, csv, gzip, json, math, os, shlex, shutil, subprocess, sys
from bisect import bisect_left
from collections import defaultdict

import yaml

TIERS = ["T_LOH", "T_amp", "T_het", "A", "B", "A1"]


def rows(path):
    if not path or not os.path.exists(path):
        return []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def stream(path):
    """Rows one at a time. The record map for one exome arm is 8.3 million rows; slurping it is GBs."""
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            yield r


def descendants(clones):
    """Every population that carries a clone's variants, itself included.

    A variant in clone A is in every cell of A and of A's children, and in no cell of T or B. Reading the
    tree matters: hard-coding "the tumour populations" is right only for truncal variants, which is three
    of the six tiers.
    """
    kids = defaultdict(list)
    for name, c in clones.items():
        if c.get("parent"):
            kids[c["parent"]].append(name)
    out = {}
    for name in clones:
        seen, stack = set(), [name]
        while stack:
            cur = stack.pop()
            if cur in seen:
                continue
            seen.add(cur)
            stack += kids.get(cur, [])
        out[name] = seen
    return out


def carries(rec, site, desc):
    """Does this interval copy carry the site's alt allele?"""
    if rec["source"] not in desc.get(site["clone"], set()):
        return False
    if rec["haplotype"] != site["haplotype"]:
        return False
    # "pre_cna" is the copy that already existed when the copy-number change duplicated it, so it carries
    # only what predated the change; "all" carries everything. A post-CNA variant is therefore on one
    # copy however many the locus has, which is what the catalog's multiplicity column says.
    if rec["copy_kind"] == "all":
        return True
    return rec["copy_kind"] == "pre_cna" and site.get("timing") == "pre_cna"


def record_counts(readmap, wanted=None, prefer_awk=True):
    """Reads per record; every record when `wanted` is None.

    The exome readmap is 68 million lines and the WGS one is about a billion: a Python loop over the
    latter takes longer than the alignment the check exists to avoid, and the work per line is one split
    and one dict bump. awk does it an order of magnitude faster, so it gets first refusal and Python is
    the fallback when awk or zcat is missing. Both paths are asserted to agree in the tests.
    """
    gz = readmap.endswith(".gz")
    if prefer_awk and shutil.which("awk") and (shutil.which("zcat") if gz else True):
        reader = "zcat" if gz else "cat"
        cmd = (f"{reader} -- {shlex.quote(readmap)} | "
               "awk -F'\\t' 'NR>1{c[$2]++} END{for (k in c) print k\"\\t\"c[k]}'")
        p = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
        if p.returncode == 0 and p.stdout:
            out = defaultdict(int)
            for line in p.stdout.splitlines():
                rec, n = line.rsplit("\t", 1)
                if wanted is None or rec in wanted:
                    out[rec] = int(n)
            return out
    counts = defaultdict(int)
    op = gzip.open if gz else open
    with op(readmap, "rt") as fh:
        next(fh, None)
        for line in fh:
            rec = line.rstrip("\n").split("\t")[1]
            if wanted is None or rec in wanted:
                counts[rec] += 1
    return counts


def cached_record_counts(readmap, wanted, cache):
    """The same counts, computed once per library rather than once per invocation.

    14_accept runs the acceptance per chromosome, so a check that streams the readmap itself streams it
    ten times over for one library -- twenty minutes each on the WGS arm, for a table that does not
    depend on which chromosome is being asked about. The cache is keyed only by being newer than the
    readmap it came from, the same currency rule the release uses everywhere else (section 20.2): a
    cache older than its input is a stale answer wearing a current one's clothes, so it is rebuilt.
    """
    if cache and os.path.exists(cache) and os.path.getmtime(cache) >= os.path.getmtime(readmap):
        out = defaultdict(int)
        with open(cache) as fh:
            for line in fh:
                rec, n = line.rstrip("\n").rsplit("\t", 1)
                if rec in wanted:
                    out[rec] = int(n)
        return out
    every = record_counts(readmap, None)
    if cache:
        tmp = f"{cache}.{os.getpid()}"
        with open(tmp, "w") as fh:
            for rec, n in sorted(every.items()):
                fh.write(f"{rec}\t{n}\n")
        os.replace(tmp, cache)
    return defaultdict(int, {k: v for k, v in every.items() if k in wanted})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--catalog-dir", default="output")
    ap.add_argument("--design", default="design.yaml")
    ap.add_argument("--records", required=True, help="merged record map for this library")
    ap.add_argument("--readmap", required=True, help="merged read -> record map for this library")
    ap.add_argument("--chrom", default="", help="comma separated; default every chromosome in the catalog")
    ap.add_argument("--counts", default="",
                    help="cache of reads per record, built if missing or older than the readmap and "
                         "reused by every later invocation against the same library")
    ap.add_argument("--report", required=True)
    ap.add_argument("--bound", type=float, default=0.05,
                    help="largest tolerated relative deviation of a tier's pooled fraction (default 0.05)")
    a = ap.parse_args()

    clones = yaml.safe_load(open(a.design))["datasets"][a.dataset]["clones"]
    desc = descendants(clones)
    chroms = set(a.chrom.split(",")) if a.chrom else None

    sites = []
    for r in rows(os.path.join(a.catalog_dir, f"{a.dataset}.snv_indel.tsv")):
        if chroms and r["chrom"] not in chroms:
            continue
        if r["clonality_tier"] not in TIERS or not r.get("expected_vaf_tumor"):
            continue
        sites.append(r)
    if not sites:
        print(f"no designed variants to measure for {a.dataset}"
              + (f" on {a.chrom}" if a.chrom else ""), flush=True)
        sys.exit(2)

    # Index the sites by chromosome and position so each record is matched by a binary search over the
    # sites it could contain. The obvious nested loop is 8.3 million records times 873 sites for one
    # exome arm -- 7.2 billion comparisons, measured at 66 minutes -- and the acceptance would pay it
    # per dataset. Streaming the records matters as much: as a list of dicts, that file is gigabytes.
    by_chrom = defaultdict(list)
    for v in sites:
        by_chrom[v["chrom"]].append((int(v["pos"]), v))
    for c in by_chrom:
        by_chrom[c].sort(key=lambda t: t[0])
    at = {c: [p for p, _ in lst] for c, lst in by_chrom.items()}

    want, spans = defaultdict(list), defaultdict(set)
    for rec in stream(a.records):
        lst = by_chrom.get(rec["chrom"])
        if not lst or (chroms and rec["chrom"] not in chroms):
            continue
        s, e = int(rec["start"]), int(rec["end"])
        lo = bisect_left(at[rec["chrom"]], s)
        hi = bisect_left(at[rec["chrom"]], e)
        if lo == hi:
            continue
        for _, v in lst[lo:hi]:
            key = f"{v['chrom']}:{v['pos']}"
            want[rec["record"]].append((key, carries(rec, v, desc)))
            spans[key].add((s, e))

    counts = cached_record_counts(a.readmap, want, a.counts)

    per = defaultdict(lambda: [0, 0])
    for rec, lst in want.items():
        c = counts.get(rec, 0)
        for key, is_alt in lst:
            per[key][1] += c
            if is_alt:
                per[key][0] += c

    by = {f"{v['chrom']}:{v['pos']}": v for v in sites}
    tiers, out_sites = defaultdict(lambda: [0.0, 0.0, 0, 0]), []
    for key, (alt, tot) in per.items():
        if tot == 0:
            continue
        v = by[key]
        exp, obs = float(v["expected_vaf_tumor"]), alt / tot
        t = tiers[v["clonality_tier"]]
        t[0] += exp * tot
        t[1] += alt
        t[2] += tot
        t[3] += 1
        out_sites.append({"site": key, "tier": v["clonality_tier"], "haplotype": v["haplotype"],
                          "timing": v.get("timing"), "multiplicity": v.get("multiplicity"),
                          "designed": round(exp, 4), "truth": round(obs, 4),
                          "ratio": round(obs / exp, 4) if exp else None,
                          "alt_reads": alt, "reads": tot})

    nums, failures = {}, []
    for tier in TIERS:
        if tier not in tiers:
            continue
        e_reads, alt, tot, n = tiers[tier]
        exp, obs = e_reads / tot, alt / tot
        # The counts are reads, so the sampling error here is the generator's own draw, not an aligner's.
        # It is tiny at release depth -- which is the point of measuring this way -- so the bound is a
        # flat relative tolerance with the sampling error reported beside it rather than folded in.
        sigma = math.sqrt(max(exp * (1 - exp) / tot, 1e-12)) / exp
        dev = abs(obs - exp) / exp if exp else 0.0
        nums[tier] = {"sites": n, "reads": tot, "designed": round(exp, 4), "truth": round(obs, 4),
                      "ratio": round(obs / exp, 4) if exp else None,
                      "relative_sigma": round(sigma, 4), "deviation": round(dev, 4),
                      "bound": a.bound}
        if dev > a.bound:
            failures.append(f"{tier} {dev:.1%} > {a.bound:.0%} ({dev / sigma:.0f} sigma, {n} sites)")

    multi = sorted(k for k, v in spans.items() if len(v) > 1)
    report = {"dataset": a.dataset, "chrom": a.chrom or "all", "n_sites": len(out_sites),
              "tiers": nums, "failures": failures,
              "sites_covered_by_more_than_one_interval": multi,
              "sites": sorted(out_sites, key=lambda s: (s["site"].split(":")[0],
                                                        int(s["site"].split(":")[1])))}
    with open(a.report, "w") as fh:
        json.dump(report, fh, indent=2)

    for tier in TIERS:
        if tier in nums:
            d = nums[tier]
            print(f"  {tier:>6}: {d['sites']:>4} sites, {d['reads']:>9} reads, designed {d['designed']:.4f}, "
                  f"written {d['truth']:.4f}, ratio {d['ratio']:.3f}", flush=True)
    if multi:
        print(f"\n{len(multi)} site(s) lie under more than one interval; their fractions pool copies of "
              f"different intervals and are not exact: {', '.join(multi[:5])}", flush=True)
    if not out_sites:
        print("\nno site had a single read in the map: nothing was measured", flush=True)
        sys.exit(2)
    if failures:
        print(f"\nthe generator and the design disagree: {'; '.join(failures)}", flush=True)
        sys.exit(1)
    print(f"\nall {len(nums)} tier(s) within {a.bound:.0%} of design across {len(out_sites)} sites",
          flush=True)


if __name__ == "__main__":
    main()
