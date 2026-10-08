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
import csv
import os
import re
import subprocess

# (name, target, tolerance) where tolerance is a fraction of the target; None means "report, do not judge"
TARGETS = {
    "pacbio":        {"length_mean": (16689, 0.08), "accuracy": (0.99823, 0.0015)},
    "ont_wgs":       {"length_mean": (19117, 0.08), "accuracy": (0.98676, 0.0030)},
    # Single-cell ONT RNA length is reported against the real 910 bp IPISRC044 median without a bound,
    # because it is set by the transcript and the size-selection curve rather than by a parameter and the
    # curve is fitted on Kinnex. Bulk DOES have a real reference now: 2,091 bp, the mean over 400,000
    # reads of the HG002 Kinnex FLNC BAM, which is a genuine full-length cDNA library. It used to be
    # generated with the single-cell model, at 910 bp, in a sample the design calls "cDNA, full length" --
    # so this target is the one that would have caught that.
    "ont_sc_rna":    {"length_mean": (910, None),   "accuracy": (0.98220, 0.0030)},
    "ont_bulk_rna":  {"length_mean": (2091, 0.15),  "accuracy": (0.98220, 0.0030)},
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



ILLUMINA_NAME = re.compile(r"^[A-Z]+:\d+:[A-Z0-9]+:\d+:\d+:\d+:\d+$")
PACBIO_NAME = re.compile(r"^m\d+_\d+_\d+_s\d+/\d+/ccs(/\d+_\d+)?$")
HEX8 = re.compile(r"^[0-9a-f]{8}(/.*)?$")


def _fq_head(path, n_lines):
    out = subprocess.run(f"pigz -dc {path} 2>/dev/null | head -{n_lines}", shell=True,
                         capture_output=True, text=True).stdout
    return out.splitlines()


def ont_deliverables(release):
    """The ONT FASTQs that are actually SHIPPED, which is not everything matching ont_*.

    `ont_wgs/` holds the per-chromosome intermediates the merge consumes -- 48 files a dataset -- and
    `ont_wgs_merged/` holds the four libraries. Globbing `ont_*/*/*.fastq.gz` picks up the intermediates
    and reports 192 failures for files nobody receives, which is the same mistake as checking the
    builder's output directory instead of the report directory: a check that fires on scratch is worse
    than no check, because it buries the real result.
    """
    wgs = sorted(glob.glob(f"{release}/ont_wgs_merged/*/*_ont_wgs.fastq.gz"))
    if not wgs:
        wgs = sorted(glob.glob(f"{release}/ont_wgs/*/*_ont_wgs.fastq.gz"))
    return (wgs
            + sorted(glob.glob(f"{release}/ont_sc_rna/*/*_ont_sc_rna.fastq.gz"))
            + sorted(glob.glob(f"{release}/ont_bulk_rna/*/*_ont_bulk_rna.fastq.gz")))


def check_ont_format(release, results):
    """No ONT read may carry the simulator's description, and the map that replaced it must be there.

    Badread writes the source reference, strand and coordinates into a genomic read's header and the
    source molecule id into an RNA read's. A benchmark whose reads carry their own answer is not a
    benchmark, so the description is stripped into `*_read_map.tsv.gz` and this is the check that it was.
    """
    for f in ont_deliverables(release):
        heads = [l for i, l in enumerate(_fq_head(f, 40000)) if i % 4 == 0]
        if not heads:
            continue
        leaky = [h for h in heads if " " in h]
        mp = f.replace(".fastq.gz", "_read_map.tsv.gz")
        results.append({
            "check": "ont_header_has_no_description", "file": os.path.basename(f),
            "value": f"{len(leaky)}/{len(heads)} headers carry a description",
            "verdict": "PASS" if not leaky else f"FAIL ({leaky[0][:70]})"})
        results.append({
            "check": "ont_read_map_present", "file": os.path.basename(mp),
            "value": f"{os.path.getsize(mp) / 1e6:.0f} MB" if os.path.exists(mp) else "absent",
            "verdict": "PASS" if os.path.exists(mp) and os.path.getsize(mp) > 1000 else "FAIL"})


def check_pacbio_format(release, results, samtools="samtools"):
    """Every PacBio BAM must be one the PacBio tools will read.

    The delivered HiFi BAMs came from `samtools import`: eleven fields, no tags, UUID read names, no read
    group, no index. Each clause here is a thing some tool refuses -- the @RG ID must be eight hex digits
    or pbindex aborts with "ERROR: stoul", the movie must be in PU or skera names every segment with a
    leading slash, and np must vary or a pass-count filter has nothing to filter on.
    """
    # pacbio_merged, not pacbio: the latter holds 48 per-chromosome intermediates a dataset. Same reason
    # as ont_deliverables.
    bams = sorted(glob.glob(f"{release}/pacbio_merged/*/*_hifi.bam"))
    if not bams:
        bams = sorted(glob.glob(f"{release}/pacbio/*/*_hifi.bam"))[:4]
    bams += sorted(glob.glob(f"{release}/kinnex_*/*/*_hifi_reads.bam"))
    bams += sorted(glob.glob(f"{release}/kinnex_*/*/*_segmented.bam"))
    for f in [b for b in bams if "non_passing" not in b and "rebuild" not in b]:
        base = os.path.basename(f)
        hdr = subprocess.run(f"{samtools} view -H {f}", shell=True, capture_output=True,
                             text=True).stdout
        rg = next((l for l in hdr.splitlines() if l.startswith("@RG")), "")
        fields = dict(x.split(":", 1) for x in rg.split("\t")[1:] if ":" in x)
        results.append({"check": "rg_id_is_8_hex", "file": base, "value": fields.get("ID", "(none)"),
                        "verdict": "PASS" if HEX8.match(fields.get("ID", "")) else "FAIL"})
        results.append({"check": "rg_pu_is_movie", "file": base, "value": fields.get("PU", "(none)"),
                        "verdict": "PASS" if fields.get("PU", "").startswith("m") else "FAIL"})
        results.append({"check": "rg_pm_is_instrument", "file": base, "value": fields.get("PM", "(none)"),
                        "verdict": "PASS" if fields.get("PM") else "FAIL"})

        recs = subprocess.run(f"{samtools} view {f} 2>/dev/null | head -20000",
                              shell=True, capture_output=True, text=True).stdout.splitlines()
        if not recs:
            continue
        names = [r.split("\t", 1)[0] for r in recs]
        bad = [n for n in names if not PACBIO_NAME.match(n)]
        results.append({"check": "pacbio_read_name", "file": base,
                        "value": f"{len(bad)}/{len(names)} malformed; e.g. {names[0]}",
                        "verdict": "PASS" if not bad else f"FAIL ({bad[0][:60]})"})
        tags = {}
        for r in recs:
            for fd in r.split("\t")[11:]:
                if fd[:5] in ("zm:i:", "np:i:", "ec:f:", "rq:f:", "qs:i:", "qe:i:"):
                    tags.setdefault(fd[:2], set()).add(fd[5:])
        missing = [t for t in ("zm", "np", "ec", "rq", "qs", "qe") if t not in tags]
        results.append({"check": "hifi_tags_present", "file": base,
                        "value": "missing " + ",".join(missing) if missing else "zm np ec rq qs qe",
                        "verdict": "PASS" if not missing else "FAIL"})
        n_np = len(tags.get("np", ()))
        results.append({"check": "np_varies", "file": base,
                        "value": f"{n_np} distinct np over {len(recs):,} reads",
                        "verdict": "PASS" if n_np >= 5 else f"FAIL ({n_np} distinct)"})
        pbi = f + ".pbi"
        results.append({"check": "pbi_present", "file": base,
                        "value": f"{os.path.getsize(pbi) / 1e6:.0f} MB" if os.path.exists(pbi)
                                 else "absent",
                        "verdict": "PASS" if os.path.exists(pbi) and os.path.getsize(pbi) > 128
                                   else "FAIL"})


def check_tenx_format(release, results):
    """10x reads must be named like Illumina reads, and R1 qualities must have real variety.

    R1 was a flat two-level draw giving 663 distinct quality strings in 20,000 reads. The bound here is
    the REAL library's 2,819, not R2's 19,461: R1 is 26 four-level-binned cycles and genuinely has low
    variety, so R2 was never the right comparison.
    """
    for f in sorted(glob.glob(f"{release}/tenx_*/*/*_R1_001.fastq.gz")):
        base = os.path.basename(f)
        lines = _fq_head(f, 80000)
        heads = [l[1:].split()[0] for i, l in enumerate(lines) if i % 4 == 0 and l]
        quals = [l for i, l in enumerate(lines) if i % 4 == 3]
        bad = [h for h in heads if not ILLUMINA_NAME.match(h)]
        results.append({"check": "tenx_read_name", "file": base,
                        "value": f"{len(bad)}/{len(heads)} malformed; e.g. {heads[0] if heads else ''}",
                        "verdict": "PASS" if heads and not bad else f"FAIL ({bad[0][:40] if bad else 'no reads'})"})
        d = len(set(quals))
        results.append({"check": "r1_quality_variety", "file": base,
                        "value": f"{d} distinct over {len(quals):,} reads (real library: 2,819)",
                        "verdict": "PASS" if 1200 <= d <= 6000 else f"FAIL ({d})"})


def check_tcr_truth(release, results):
    """CDR3s must be anchored, free of internal cysteines, and tied to the V and J the table names."""
    for f in sorted(glob.glob(f"{release}/tenx_tcr/*/*_tcr_clonotypes.tsv")):
        base = os.path.basename(f)
        with open(f) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        if not rows:
            continue
        n = len(rows)
        anch = sum(1 for r in rows for k in ("trb", "tra")
                   if r.get(k, "").startswith("C") and r.get(k, "").endswith("F"))
        icys = sum(1 for r in rows for k in ("trb", "tra") if "C" in r.get(k, "")[1:-1])
        vj = sum(1 for r in rows if r.get("trb_v") and r.get("trb_j")
                 and r.get("tra_v") and r.get("tra_j"))
        fp = sum(1 for r in rows if r.get("flagpost_antigen"))
        results.append({"check": "cdr3_anchored_C_to_F", "file": base,
                        "value": f"{anch}/{2 * n} chains", "verdict": "PASS" if anch == 2 * n else "FAIL"})
        results.append({"check": "cdr3_internal_cysteine", "file": base,
                        "value": f"{icys}/{2 * n} chains ({100 * icys / (2 * n):.1f}%)",
                        "verdict": "PASS" if icys <= 0.03 * 2 * n else f"FAIL ({icys})"})
        results.append({"check": "clonotype_names_v_and_j", "file": base,
                        "value": f"{vj}/{n} clonotypes", "verdict": "PASS" if vj == n else "FAIL"})
        results.append({"check": "flagpost_antigen_linked", "file": base, "value": f"{fp} clonotypes",
                        "verdict": "PASS" if fp >= 5 else f"FAIL ({fp})"})


def check_alphabet(release, results, samtools="samtools", limit=400000):
    """Every delivered base must be A, C, G, T or N.

    This is the check that was missing. A VCF ALT column carries symbolic alleles as well as sequence --
    `*` for an allele removed by a spanning deletion, `<DEL>`, breakend notation -- and the generator
    appended whatever allele the genotype selected straight into the sequence. 2.65 % of the HG002 Q100
    records carry a `*`, so literal `*` characters reached the delivered IGI-SYN-SEQ-01 reads of every
    arm built from the genome, and the same records silently deleted 960,223 reference bases across
    74,017 sites. The aligners tolerated the character and Cell Ranger refused it, which is how it
    surfaced: downstream, in someone else's run, rather than here.

    A length check cannot see this and an accuracy check cannot either. One grep over the reads can.
    """
    for f in (sorted(glob.glob(f"{release}/*/*/*_R1.fastq.gz"))
              + sorted(glob.glob(f"{release}/*/*/*_R1_001.fastq.gz"))
              + ont_deliverables(release)):
        seqs = [l for i, l in enumerate(_fq_head(f, limit)) if i % 4 == 1]
        if not seqs:
            continue
        bad = sorted({c for l in seqs for c in l.upper()} - set("ACGTN"))
        results.append({"check": "bases are ACGTN", "file": os.path.basename(f),
                        "value": f"{len(seqs):,} reads, offending characters: {bad or 'none'}",
                        "verdict": "PASS" if not bad else f"FAIL ({''.join(bad)})"})
    for f in (sorted(glob.glob(f"{release}/kinnex_*/*/*_segmented.bam"))
              + sorted(glob.glob(f"{release}/pacbio_merged/*/*_hifi.bam"))):
        out = subprocess.run(f"{samtools} view {f} 2>/dev/null | head -{limit // 4} | cut -f10",
                             shell=True, capture_output=True, text=True).stdout
        seqs = out.splitlines()
        if not seqs:
            continue
        bad = sorted({c for l in seqs for c in l.upper()} - set("ACGTN"))
        results.append({"check": "bases are ACGTN", "file": os.path.basename(f),
                        "value": f"{len(seqs):,} records, offending characters: {bad or 'none'}",
                        "verdict": "PASS" if not bad else f"FAIL ({''.join(bad)})"})


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

    # Format conformance: not "is the sequencing what was asked for" but "will the tools read it at all".
    # Eight defects found by review on 2026-10-05 were all of this kind, and none of them would have been
    # caught by a length or accuracy check.
    fmt = []
    check_ont_format(a.release, fmt)
    check_pacbio_format(a.release, fmt, a.samtools)
    check_tenx_format(a.release, fmt)
    check_tcr_truth(a.release, fmt)
    check_alphabet(a.release, fmt, a.samtools)
    report["format"] = fmt
    report["n_fail"] += sum(1 for r in fmt if r["verdict"].startswith("FAIL"))
    print()
    for r in fmt:
        print(f"  {r['check']:<30} {r['file'][:40]:<40} {str(r['value'])[:46]:<46} {r['verdict']}")

    with open(a.out, "w") as fh:
        json.dump(report, fh, indent=2); fh.write("\n")
    print(f"\n  {report['n_fail']} failing check(s) -> {a.out}")


if __name__ == "__main__":
    main()
