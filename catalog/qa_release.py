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
import datetime
import glob
import gzip
import json
import math
import csv
import os
import re
import shlex
import subprocess

# (name, target, tolerance) where tolerance is a fraction of the target; None means "report, do not judge"
TARGETS = {
    # Accuracy targets for the three PacBio arms are the SIMULATOR'S measured true accuracy, 0.99695,
    # not real HiFi's 0.99823. read_stats derives accuracy from the quality strings, and those are now
    # calibrated to the alignment (pacbio_qual_error_scale), so this check asks "do the reads describe
    # themselves correctly" -- which is answerable -- rather than "did Badread reach real HiFi", which it
    # cannot: its pacbio2021 model floors at 3.05e-3/bp against real 1.77e-3. The gap is recorded in
    # build-reference.md and is a simulator limit, not a build defect, so it must not sit in a check that
    # fails on every release.
    "pacbio":        {"length_mean": (16689, 0.08), "accuracy": (0.99695, 0.0015)},
    "ont_wgs":       {"length_mean": (19117, 0.08), "accuracy": (0.98676, 0.0030)},
    # Single-cell ONT RNA length is reported against the real 910 bp IPISRC044 median without a bound,
    # because it is set by the transcript and the size-selection curve rather than by a parameter and the
    # curve is fitted on Kinnex. Bulk DOES have a real reference now: 2,091 bp, the mean over 400,000
    # reads of the HG002 Kinnex FLNC BAM, which is a genuine full-length cDNA library. It used to be
    # generated with the single-cell model, at 910 bp, in a sample the design calls "cDNA, full length" --
    # so this target is the one that would have caught that.
    "ont_sc_rna":    {"length_mean": (910, None),   "accuracy": (0.98220, 0.0030)},
    "ont_bulk_rna":  {"length_mean": (2091, 0.15),  "accuracy": (0.98220, 0.0030)},
    "kinnex_bulk":   {"length_mean": (2223, 0.12),  "accuracy": (0.99695, 0.0020)},
    # 2,155 is the real OVERALL mean of skera's output on single-cell arrays, which is what a segmented
    # BAM contains. The first version of this target used 952, the real median of the SINGLE-cDNA subset,
    # and failed a correct library at 2,416 by comparing it against the wrong population: with nine
    # adapters for sixteen cDNAs about 12 % of skera's output is multi-cDNA blocks averaging 8.6 TSOs, and
    # those are in the BAM too. Judging a mixture against one of its components is not a check.
    "kinnex_sc":     {"length_mean": (2155, 0.20),  "accuracy": (0.99695, 0.0020)},
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


_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(x):
    """Local, so qa_release stays importable without the catalog package on sys.path."""
    return x.translate(_COMP)[::-1]


IUPAC_SET = set("RYSWKMBDHV")


# One awk pass per deliverable: count the records, and tally every character that is not ACGTN. `gsub`
# strips the legal alphabet so the per-character loop only ever runs on an offence, which makes a clean
# library cost one gsub per read and nothing else. Memory is bounded by the number of DISTINCT offending
# characters, not by the library, so this streams a 67 GB BAM in constant space.
# The regex pre-test is not a micro-optimisation: `gsub` rebuilds the string on every read it touches,
# and skipping it on the ~99.9999 % of reads that are clean took one exome library from 105 s to 11 s --
# 9.5x, same counts. Over a 1.3 TB release that is the difference between a check that runs on every
# build and one that gets skipped because it is too slow, which is how a check stops being a check.
#
# `ph` (placeholder) is set for BAM input only. In SAM, a SEQ field of exactly `*` means "no sequence
# stored" -- it is not a base at all, and counting it as a symbolic allele was a false positive in this
# scan's first version: ds-02's chr1to6 Kinnex BAM reported `* x7` with zero `*` embedded in any sequence.
# A record with no sequence is still wrong in a deliverable, so it is counted and reported, but under its
# own name, because it needs a different fix from a `*` written into the middle of a read.
_TALLY = (r'{ n++; if (ph && $0 == "*") { e++; next }'
          r'  if ($0 ~ /[^ACGTNacgtn]/) { s = $0; gsub(/[ACGTNacgtn]/, "", s);'
          r'  for (i = 1; i <= length(s); i++) c[substr(s, i, 1)]++ } }'
          r'END { printf "RECORDS %d\nEMPTY %d\n", n, e; for (k in c) printf "CHAR %s %d\n", k, c[k] }')


def _scan_alphabet(path, samtools="samtools", complete=True, limit=400000):
    """Count non-ACGTN characters in a deliverable's sequence.

    Returns (n_records, {char: count}, n_empty), where n_empty counts records carrying no sequence at
    all -- SAM's `*` placeholder, which is not a base and must not be reported as one.

    COMPLETE BY DEFAULT, and that is the whole point of this function. The first version of this check
    read 400,000 lines -- 100,000 reads -- off the front of each file and declared the alphabet clean if
    it saw nothing. The exome arm carries IUPAC codes at roughly ONE READ PER MILLION, because only 94
    reference positions are ambiguous and an exome covers few of them, so a 100,000-read sample misses
    them with probability ~0.9 and reports PASS on a library that aborts razers3 on load. A sample can
    witness contamination; it cannot establish absence at a rate below its own resolution, and absence is
    exactly what every consumer of these files depends on.

    `limit` applies only when `complete` is False, which exists for a fast pre-flight and never for a
    verdict -- the value string says which was done, so no reader has to guess.
    """
    if path.endswith(".bam"):
        src = f"{samtools} view {shlex.quote(path)} | cut -f10"
        ph = 1
    else:
        src = f"pigz -dc {shlex.quote(path)} | awk 'NR % 4 == 2'"
        ph = 0
    if complete:
        # No `head` in the pipeline, so pipefail is safe and a decompression failure cannot masquerade
        # as a clean library -- the difference between "nothing offended" and "nothing was read".
        cmd = f"set -o pipefail; {src} | awk -v ph={ph} '{_TALLY}'"
    else:
        cmd = f"{src} | head -{limit // 4} | awk -v ph={ph} '{_TALLY}'"
    r = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
    if r.returncode != 0:
        return None, {}, 0
    n, counts, empty = 0, {}, 0
    for line in r.stdout.splitlines():
        p = line.split()
        if p[:1] == ["RECORDS"]:
            n = int(p[1])
        elif p[:1] == ["EMPTY"]:
            empty = int(p[1])
        elif p[:1] == ["CHAR"]:
            counts[p[1].upper()] = counts.get(p[1].upper(), 0) + int(p[2])
    return n, counts, empty


def _alphabet_verdict(bad, path=""):
    """Verdict for a set of offending characters, naming the cause so the fix is unambiguous.

    EVERY non-ACGTN character FAILS, wherever it appears. This was not always so, and the history is the
    argument for it: the codes were first reported as WARN outside the 10x arms, on the reasoning that
    minimap2, pbmm2 and bwa all tolerate ambiguity and only Cell Ranger refuses a run outright. That
    reasoning was an ENUMERATION OF THE CONSUMERS I HAPPENED TO KNOW, and it was wrong within a day:
    OptiType's razers3 aborts on load with `seqan::ParseError: Unexpected character 'Y' found`, because
    SeqAn's reader accepts only A, C, G, T and N. It killed HLA typing on the normal exome through all
    seven retries while every aligner in the same run stayed quiet.

    The lesson generalises past razers3. A delivered base either is a base or it is not, and a reference
    dataset cannot know what will read it -- callers, typers, assemblers, counters, things not written
    yet. Ambiguity codes are normalised to N at source now (genome.normalize_bases), so a clean alphabet
    costs nothing to guarantee and the only arms that fail here are the ones built before that, which
    need rebuilding regardless. Tolerating a defect because the tools we tested happen to survive it is
    how this reached someone else's pipeline in the first place.
    """
    if not bad:
        return "PASS"
    tag = "".join(sorted(bad))
    chars = set(bad)
    if chars <= IUPAC_SET:
        return (f"FAIL ({tag}; IUPAC codes -- razers3/SeqAn aborts on load and Cell Ranger refuses the "
                f"run, rebuild this arm under genome.normalize_bases)")
    if chars & IUPAC_SET:
        return f"FAIL ({tag}; IUPAC codes AND a symbolic allele -- see check_alphabet)"
    return f"FAIL ({tag}; symbolic allele written into the sequence -- generator defect)"


def delivered_sequence_files(release):
    """Every shipped file whose sequence a consumer reads -- BOTH MATES, not just R1.

    The earlier globs ended at `*_R1.fastq.gz` and `*_R1_001.fastq.gz`, so half of every paired library
    went unexamined, including the 10x R2 that carries the cDNA and is the only mate Cell Ranger aligns.
    """
    pats = ("*_R1.fastq.gz", "*_R2.fastq.gz", "*_R1_001.fastq.gz", "*_R2_001.fastq.gz")
    out = []
    for p in pats:
        out += glob.glob(f"{release}/*/*/{p}")
    out += ont_deliverables(release)
    out += glob.glob(f"{release}/kinnex_*/*/*_segmented.bam")
    out += glob.glob(f"{release}/pacbio_merged/*/*_hifi.bam")
    return sorted(set(out))


def check_alphabet(release, results, samtools="samtools", complete=True, limit=400000):
    """Every delivered base must be A, C, G, T or N. Counted completely, over every shipped file.

    This is the check that was missing, and it has now caught two different defects that every other
    check in this release passed over:

      *, <DEL>, breakend notation  -- a symbolic VCF allele written into the sequence. 2.65 % of the
                                      HG002 Q100 records carry a `*`, so literal `*` characters reached
                                      the delivered reads of every ds-01 arm built from the genome, and
                                      the same records silently deleted 960,223 reference bases across
                                      74,017 sites. A generator defect; EditSet.add now rejects these at
                                      the point every edit passes through.
      R Y S W K M B D H V          -- IUPAC ambiguity codes, which GRCh38 itself contains at 94
                                      positions (36 on chr10) and the viral reference at 497. Not
                                      invented here, but no instrument emits them and strict consumers
                                      abort on them, so genome.normalize_bases turns them into N at
                                      every point sequence enters.

    Both surfaced downstream, in someone else's run, rather than here -- the first in Cell Ranger, the
    second in razers3. A length check cannot see either and an accuracy check cannot either. One
    complete pass over the bases can, which is why this one does not sample.
    """
    for f in delivered_sequence_files(release):
        n, counts, empty = _scan_alphabet(f, samtools, complete, limit)
        if n is None:
            results.append({"check": "bases are ACGTN", "file": os.path.basename(f),
                            "value": "UNREADABLE -- decompression or samtools failed",
                            "verdict": "FAIL (could not read the file; this is not a clean alphabet)"})
            continue
        if not n:
            continue
        unit = "records" if f.endswith(".bam") else "reads"
        scope = "complete" if complete else f"first {limit // 4:,} sampled"
        detail = ", ".join(f"{c}x{counts[c]:,}" for c in sorted(counts)) or "none"
        results.append({"check": "bases are ACGTN", "file": os.path.basename(f),
                        "value": f"{n:,} {unit} ({scope}), offending: {detail}",
                        "verdict": _alphabet_verdict(sorted(counts), f)})
        if empty:
            results.append({"check": "every record has a sequence", "file": os.path.basename(f),
                            "value": f"{empty:,} of {n:,} {unit} carry no sequence",
                            "verdict": f"FAIL ({empty:,} records with SEQ='*' -- a delivered read with "
                                       f"no bases; not an alphabet problem, see check_alphabet)"})



# Which directory a deliverable is assembled FROM. The merged libraries and the chr1to6 subsets are built
# out of the per-chromosome arms, so each of these is a "newer than its inputs" relation that must hold.
ASSEMBLED_FROM = {
    "merged":          "wes",
    "wgs_merged":      "wgs",
    "ont_wgs_merged":  "ont_wgs",
    "pacbio_merged":   "pacbio",
    "wes_chr1to6":     "wes",
    "wgs_chr1to6":     "wgs",
    "ont_wgs_chr1to6": "ont_wgs",
    "pacbio_chr1to6":  "pacbio",
}


def _ts(epoch):
    return datetime.datetime.fromtimestamp(epoch).strftime("%Y-%m-%d %H:%M")


def check_tcr_r2_direction(release, results, limit=400000, min_reads=200):
    """R2 must be ANTISENSE, not merely consistent -- and the TCR truth map can prove it.

    check_tenx_r2_strand establishes that reads from one molecule agree with each other, which catches an
    unstranded library but cannot tell antisense from sense: a library flipped wholesale the wrong way
    would pass it. Direction matters, because Cell Ranger decides chemistry on it (SC5P-R2 antisense,
    SC3Pv2 sense) and would then mis-call the assay rather than refuse it.

    The TCR arm ships `cdr3_nt` in its truth map, so the direction is checkable against delivered data
    with no reference and no aligner: a read covering the CDR3 must contain the REVERSE COMPLEMENT of that
    sequence, never the sequence itself. Measured after the orientation fix:

        ds-01 TCR   sense 0   antisense 61,976
        ds-02 TCR   sense 0   antisense 59,727

    The GEX arm cannot be checked this way -- its truth map carries a record id rather than sequence -- so
    for GEX the direction rests on tenx.orient_antisense plus the consistency check. That asymmetry is
    worth knowing rather than papering over: TCR verifies direction on delivered bases, GEX does not.
    """
    for mp in sorted(glob.glob(f"{release}/tenx_tcr/*/*_tcr_molecules.tsv.gz")):
        ds_dir = os.path.dirname(mp)
        tag = "chr1to6" if "chr1to6" in os.path.basename(mp) else "full"
        r2s = [f for f in glob.glob(f"{ds_dir}/*_R2_001.fastq.gz")
               if (tag == "chr1to6") == ("chr1to6" in os.path.basename(f))]
        if not r2s:
            continue
        sense = anti = 0
        with gzip.open(r2s[0], "rt") as fq, gzip.open(mp, "rt") as fm:
            fm.readline()
            for _ in range(limit):
                h = fq.readline()
                if not h:
                    break
                seq = fq.readline().strip(); fq.readline(); fq.readline()
                row = fm.readline()
                if not row:
                    break
                f = row.rstrip("\n").split("\t")
                if len(f) < 6 or len(f[5]) < 20:
                    continue
                probe = f[5][:20]
                if probe in seq:
                    sense += 1
                elif revcomp(probe) in seq:
                    anti += 1
        n = sense + anti
        if n < min_reads:
            verdict = f"SKIP (only {n} reads covering a CDR3)"
        elif sense == 0:
            verdict = "PASS"
        else:
            verdict = (f"FAIL ({sense}/{n} reads carry the CDR3 in SENSE orientation; 5' chemistry is "
                       f"antisense, see check_tcr_r2_direction)")
        results.append({"check": "R2 is antisense to the transcript", "file": os.path.basename(r2s[0]),
                        "value": f"{anti} antisense / {sense} sense of {n} CDR3-covering reads",
                        "verdict": verdict})


def check_currency(release, results):
    """Every assembled deliverable must be NEWER than the per-chromosome reads it was assembled from.

    This check exists because of the question "so all of the ds-01 files are now valid?", which I could
    not answer from the build logs. Nine arms had been regenerated overnight and every one reported
    COMPLETED with no failures, so the arms were current -- but no merge had run afterwards, and the
    deliverables are the MERGED libraries, which the staged inputs symlink to. Every file a consumer could
    open still predated the rebuild, the exome merge by nine days, and it still carried ~9,000 `*`
    characters a day after that defect was fixed and the reads were regenerated.

    The individual stages do compare mtimes and would have refused to ship a stale output, so nothing
    would have been delivered wrongly. But a per-stage guard only fires when that stage runs, and the
    failure mode here was a stage NOT RUNNING. Nothing looked at the release as a whole and asked whether
    it was current, so "all arms completed" read as "the release is ready" -- to me, in my own reporting.

    A fix that has not reached the file the consumer opens has not reached the consumer, and that is a
    property of the release tree rather than of any one job, so this is the layer that has to check it.

    **A PASS here is not a statement that the deliverable is correct.** It says only that it was assembled
    from the reads currently on disk. `merged/IGI-SYN-SEQ-02` passes precisely because its exome inputs are
    also from Sep 30 -- and that library is the one that aborted razers3. Currency and validity are
    different questions, and check_alphabet answers the second. Together "current and clean" means
    something; either one alone does not.
    """
    for out_dir, in_dir in sorted(ASSEMBLED_FROM.items()):
        for ds_path in sorted(glob.glob(f"{release}/{out_dir}/*")):
            if not os.path.isdir(ds_path):
                continue
            ds = os.path.basename(ds_path)
            outs = [f for f in glob.glob(f"{ds_path}/*") if os.path.isfile(f)]
            ins = [f for f in glob.glob(f"{release}/{in_dir}/{ds}/*") if os.path.isfile(f)]
            if not outs or not ins:
                continue
            newest_in = max(os.path.getmtime(f) for f in ins)
            oldest_out = min(os.path.getmtime(f) for f in outs)
            lag_h = (newest_in - oldest_out) / 3600.0
            if lag_h > 0:
                verdict = (f"FAIL (stale by {lag_h:,.1f} h -- re-run the stage that builds {out_dir}/; "
                           f"its inputs in {in_dir}/ are newer)")
            else:
                verdict = "PASS"
            results.append({"check": "deliverable newer than inputs", "file": f"{out_dir}/{ds}",
                            "value": f"output {_ts(oldest_out)}, newest input {_ts(newest_in)}",
                            "verdict": verdict})


def _seed_orientation(a, b, seed=25):
    """Same / opposite / unknown orientation of two reads known to come from the same window.

    Seeds rather than whole reads because ART's HS25 profile puts errors in most reads, and several
    offsets because two 90 bp reads drawn from a 450 bp window need not overlap at all.
    """
    for off in (32, 20, 44, 8, 56):
        s = a[off:off + seed]
        if len(s) < seed:
            continue
        if s in b:
            return "same"
        if revcomp(s) in b:
            return "opposite"
    return "unknown"


def check_tenx_r2_strand(release, results, limit=200000, min_pairs=50):
    """R2 reads from one molecule must share an orientation. 10x 5' chemistry is antisense throughout.

    The third defect of this release, and the second reported by a downstream consumer rather than caught
    here. Both 10x builders simulated R2 with `art_illumina ... -na`; ART draws either strand with equal
    probability and `-na` suppresses the ALN that records which, so the libraries were UNSTRANDED. Cell
    Ranger infers chemistry from R2 strand bias and refused the ds-02 run outright:

        Unable to distinguish between [SC5P-R2, SC3Pv2] chemistries based on the R2 read mapping
        Total Reads = 100000  Mapped reads = 92806  Sense reads = 45140  Antisense reads = 44751

    A coin flip matches no chemistry. Nothing here measured orientation, so nothing here could have seen
    it -- the same gap as the alphabet before check_alphabet existed.

    This check needs NO reference and no aligner, which is what lets it live in the release QA: reads
    sharing a truth-map record were drawn from the same 5' window, so they must agree in orientation.
    Under a random strand draw an informative pair disagrees half the time; in a correct library, never.
    So it detects the defect without knowing which direction is right; the direction itself is asserted
    separately, in tenx.orient_antisense, against Cell Ranger's own chemistry definitions.

    Measured on the delivered libraries when it was written -- all six at a coin flip, which also settled
    that TCR was affected and not merely suspected:

        ds-01 GEX full      2006 same / 1996 opposite
        ds-02 GEX full      1997 same / 2003 opposite
        ds-01 TCR           2184 same / 2233 opposite
        ds-02 TCR           2290 same / 2264 opposite
    """
    for arm, pat, mpat in (("tenx_gex", "*_R2_001.fastq.gz", "*_gex_molecules.tsv.gz"),
                           ("tenx_tcr", "*_R2_001.fastq.gz", "*_tcr_molecules.tsv.gz")):
        for r2 in sorted(glob.glob(f"{release}/{arm}/*/{pat}")):
            ds_dir = os.path.dirname(r2)
            tag = "chr1to6" if "chr1to6" in os.path.basename(r2) else "full"
            maps = [m for m in glob.glob(f"{ds_dir}/{mpat}") if (tag == "chr1to6") == ("chr1to6" in m)]
            if not maps:
                continue
            seqs = [l for i, l in enumerate(_fq_head(r2, limit * 4)) if i % 4 == 1]
            if not seqs:
                continue
            groups = {}
            with gzip.open(maps[0], "rt") as fh:
                fh.readline()
                for i, line in enumerate(fh):
                    if i >= len(seqs):
                        break
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 4:
                        continue
                    groups.setdefault(tuple(f[3:]), []).append(seqs[i])
            same = opp = 0
            for reads in groups.values():
                if len(reads) < 2:
                    continue
                for other in reads[1:]:
                    v = _seed_orientation(reads[0], other)
                    if v == "same":
                        same += 1
                    elif v == "opposite":
                        opp += 1
                if same + opp >= 4000:
                    break
            n = same + opp
            if n < min_pairs:
                verdict = f"SKIP (only {n} informative pairs in {len(seqs):,} reads)"
            elif opp == 0:
                verdict = "PASS"
            else:
                verdict = (f"FAIL ({opp}/{n} pairs from one molecule disagree -- R2 is unstranded; "
                           f"Cell Ranger cannot infer chemistry, see check_tenx_r2_strand)")
            results.append({"check": "R2 orientation is consistent", "file": os.path.basename(r2),
                            "value": f"{same} same / {opp} opposite of {n} informative pairs",
                            "verdict": verdict})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--limit", type=int, default=40000)
    # The alphabet check reads every base of every deliverable, which is minutes per library and the
    # only way to establish absence at ~1e-6 per read. This makes it sample instead, for a fast
    # pre-flight only -- a sampled PASS is not evidence of a clean alphabet and the report says so.
    ap.add_argument("--alphabet-sample", action="store_true",
                    help="sample the alphabet check instead of counting completely (pre-flight only)")
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
    check_alphabet(a.release, fmt, a.samtools, complete=not a.alphabet_sample)
    check_tenx_r2_strand(a.release, fmt)
    check_tcr_r2_direction(a.release, fmt)
    check_currency(a.release, fmt)
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
