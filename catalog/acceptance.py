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
import bisect
import csv
import gzip
import json
import math
import os
import re
import statistics
import yaml
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


def pileup(bam, reference, sites_path, min_bq=15, min_mq=20, max_depth=8000, region=None):
    """{position: (depth, bases)} at the requested sites.

    `region` matters at release scale. With `-l` alone, mpileup streams the whole BAM and filters, which on
    a 6.25 GB whole-genome alignment takes over a quarter of an hour per call; on the slice BAMs it was
    free, so the cost only appeared once the release existed. Every site in the list is on one chromosome
    by construction, so passing `-r` lets mpileup seek with the index instead.
    """
    cmd = ["samtools", "mpileup", "-f", reference, "-l", sites_path,
           "-d", str(max_depth), "-Q", str(min_bq), "-q", str(min_mq)]
    if region:
        cmd += ["-r", region]
    cmd += [bam]
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


def callable_runs(genome, chrom, start, end, min_run=1000):
    """Stretches of non-N reference within a region, 0-based half-open.

    This is the WGS analogue of a capture interval. chr21 is 18.6 % N, and a window straddling the
    acrocentric arm or the centromere reads at a fraction of its true depth simply because much of it
    cannot be sequenced: chr21:0-6,000,000 is only 11.5 % callable, so measured over the whole window it
    reads at 11.5 % of true depth, an apparent 8.7x deletion. Comparing such a window against a clean one
    (chr21:30-36 Mb is 100 % callable) would read as homozygous loss that is not there -- the same mistake
    as measuring exome depth over a whole window instead of over capture intervals, which once made an
    arbitrary 2 Mb window look 4.6x shallower than the MHC.

    Runs shorter than `min_run` are dropped: they are the ragged edges of assembly gaps, where mappability
    is poor for reasons that have nothing to do with copy number.
    """
    seq = genome.seq(chrom, start, end)
    runs, i, n = [], 0, len(seq)
    while i < n:
        if seq[i] == "N":
            i += 1
            continue
        j = i
        while j < n and seq[j] != "N":
            j += 1
        if j - i >= min_run:
            runs.append((start + i, start + j))
        i = j
    return runs


def captured_depth(bam, chrom, start, end, capture_bed, work, min_bases=5000, genome=None):
    """Mean depth over the region's measurable bases: captured bases, or non-N bases for WGS.

    Mean depth across a whole window is meaningless for an exome: it measures how much of the window is
    captured, not how deeply it is sequenced. An arbitrary 2 Mb window read 3.4x while the gene-dense MHC
    read 15.6x purely because far more of the MHC is on target. Restricting to capture intervals with
    `samtools bedcov` compares like with like.

    With `capture_bed=None` the assay has no capture, so the sub-intervals come from the reference's non-N
    runs instead (see `callable_runs`). `genome` is then required. Everything downstream is identical, so
    the copy-number check reads the same whichever assay it is given.

    The floor is on measurable bases, not on interval count. A focal amplicon is small by definition: the
    MYC amplicon holds 11 capture intervals but 15 kb of captured sequence, which at amplicon depth is
    several million read bases and plenty to estimate a mean from. An interval-count floor of 20 skipped
    it, so the amplified arm of the copy-number model went unmeasured while the check reported a pass.
    """
    sub = os.path.join(work, f"cap_{chrom}_{start}_{end}.bed")
    n = bases = 0
    if capture_bed is None:
        if genome is None:
            raise ValueError("captured_depth needs either a capture bed or a genome to find non-N runs")
        with open(sub, "w") as out:
            for lo, hi in callable_runs(genome, chrom, start, end):
                out.write(f"{chrom}\t{lo}\t{hi}\n")
                n += 1
                bases += hi - lo
    else:
        with open(capture_bed) as fh, open(sub, "w") as out:
            for line in fh:
                f = line.split("\t")
                if len(f) < 3 or f[0] != chrom:
                    continue
                s0, e0 = int(f[1]), int(f[2])
                if e0 <= start or s0 >= end:
                    continue
                lo, hi = max(s0, start), min(e0, end)
                out.write(f"{chrom}\t{lo}\t{hi}\n")
                n += 1
                bases += hi - lo
    if bases < min_bases:
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
    pile = pileup(bam, reference, sites, region=chrom)
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
            by_tier[r["clonality_tier"]].append((exp, obs, dp))
            by_class[r["class"]].append((exp, obs, dp))
        if b.count(r["alt"]) > 0 or len(r["ref"]) > 1:
            seen += 1
    nums = {}
    thin, failures, worst, worst_tier = [], [], 0.0, None
    # A flat tolerance is the wrong test. The precision of an allele-fraction estimate depends on the
    # expected number of alt reads, n x depth x VAF, not on the number of sites: the A1 tier sits at a VAF
    # of 0.015 to 0.065, so a dozen sites yield only a few dozen alt reads and one sigma is already 13 to
    # 20 %, while T_het at 0.35 has a sigma of 2 to 4 % and a 15 % bound would miss a real bias six times
    # over. Each tier is therefore compared against its own sampling error, with a floor so that a
    # systematic offset still fails where sigma is tiny.
    MIN_ALT_READS = 20
    SIGMA_MULT = 3.0
    FLOOR = 0.10
    for tier in TIER_ORDER:
        v = by_tier.get(tier)
        if not v:
            continue
        e = statistics.mean(x[0] for x in v)
        o = statistics.mean(x[1] for x in v)
        total_depth = sum(x[2] for x in v)
        alt_reads = e * total_depth
        # relative standard error of the mean observed fraction
        sigma = math.sqrt((1 - e) / (e * total_depth)) if e > 0 and total_depth else float("inf")
        dev = abs(o / e - 1) if e else 0.0
        bound = max(SIGMA_MULT * sigma, FLOOR)
        counted = alt_reads >= MIN_ALT_READS
        nums[tier] = {"n": len(v), "expected": round(e, 4), "observed": round(o, 4),
                      "ratio": round(o / e, 3) if e else None,
                      "expected_alt_reads": round(alt_reads),
                      "relative_sigma": round(sigma, 4),
                      "bound": round(bound, 4), "deviation": round(dev, 4),
                      "counted_towards_verdict": counted}
        if not counted:
            thin.append(f"{tier} alt~{alt_reads:.0f}")
            continue
        if dev > worst:
            worst, worst_tier = dev, tier
        if dev > bound:
            failures.append(f"{tier} {dev:.0%} > {bound:.0%} ({dev / sigma:.1f} sigma)")
    detail = (f"largest deviation {worst:.0%}"
              + (f" ({worst_tier}, {worst / max(nums[worst_tier]['relative_sigma'], 1e-9):.1f} sigma)"
                 if worst_tier else " (no tier with enough alt reads)")
              + f" across {sum(len(v) for v in by_tier.values())} sites")
    if failures:
        detail += "; beyond sampling error: " + ", ".join(failures)
    if thin:
        detail += f"; not counted: {', '.join(thin)}"
    rep.add("allele fractions by tier", not failures, detail, nums)
    rep.add("indels present", len([r for r in ev if r["class"] == "indel"]) == 0 or seen > 0,
            f"{seen} of {len(ev)} designed events have supporting reads")


def mean_gc_factor(genome, gc, chrom, start, end, capture_bed):
    """Captured-base-weighted mean GC factor over a region's capture intervals.

    Once capture efficiency varies with GC, the depth ratio between two regions is no longer copy number
    alone: it also carries the ratio of their GC compositions. Comparing a GC-rich amplicon against a
    GC-poor baseline without this would read as a copy-number error that is not there.

    The curve is CAPTURE efficiency, so without a capture there is no such correction to make: the WGS
    builders use a flat GcBias deliberately, and applying the exome curve here would invent a bias the
    assay does not have and then "correct" for it.
    """
    if capture_bed is None or not gc.enabled:
        return 1.0
    num = den = 0
    with open(capture_bed) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 3 or f[0] != chrom:
                continue
            s0, e0 = int(f[1]), int(f[2])
            if e0 <= start or s0 >= end:
                continue
            lo, hi = max(s0, start), min(e0, end)
            if hi - lo < 50:
                continue
            fac = gc.factor(gc.bin_of(genome.gc(chrom, lo, hi)))
            num += fac * (hi - lo)
            den += hi - lo
    return num / den if den else 1.0


def check_depth_by_cn(rep, bam, design, dataset, chrom, arms_bed, capture_bed, work, gc_curve=None):
    """Observed depth ratio between copy-number regions against the purity and copy-number formula."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import yaml
    from igi_catalog.clones import CloneModel
    from igi_catalog.genome import Genome
    from igi_catalog.pipeline import GcBias, depth_denominator, source_plan
    cfg = yaml.safe_load(open(design))["datasets"][dataset]
    cm = CloneModel(cfg, arms_bed)
    paths = yaml.safe_load(open(gc_curve)) if gc_curve and gc_curve.endswith((".yml", ".yaml")) else None
    gc = GcBias(paths.get("gc_bias_curve") if paths else gc_curve)
    genome = Genome(paths["reference_fasta"]) if paths else None
    # Without a capture bed the assay is WGS and depth is measured over non-N reference instead; that
    # needs the genome, so fail loudly rather than silently measuring over whole windows.
    unit = "captured" if capture_bed else "callable"
    if capture_bed is None and genome is None:
        rep.add("depth by copy number", False,
                "WGS depth check needs a reference genome (pass --paths) to find non-N runs")
        return
    den = depth_denominator(cm, cm.purity)
    regions = []
    for ev in (cfg.get("cna") or {}).get("T", []):
        c, s, e = cm.region(ev["region"])
        if c == chrom and e - s > 200_000:
            regions.append((ev.get("label", "cna"), s + 1, e))
    if not regions:
        rep.add("depth by copy number", True, f"no large copy-number segment on {chrom} to test")
        return
    # The baseline must be a copy-number-neutral stretch that is actually captured, and it does not have
    # to be on this chromosome. Requiring the same chromosome left the design's headline events
    # unverifiable: chr13p is acrocentric with zero capture intervals, so dataset 01's 13q LOH had no
    # usable local baseline, and chr17 has both arms altered in both datasets so it has no neutral region
    # at all. The depth model is genome-wide and GC is now corrected explicitly, so a neutral region
    # elsewhere is a valid reference; the chromosome it came from is reported either way.
    def neutral_baseline(order):
        for ch in order:
            if ch not in (genome.lengths if genome else {ch: None}):
                continue
            limit = min(genome.lengths[ch] if genome else 150_000_000, 250_000_000)
            pos = 1_000_000
            cands = []
            while pos < limit:
                a, b, _lab = cm.cn("T", ch, pos)
                on_region = ch == chrom and any(s <= pos <= e for _l, s, e in regions)
                if (a, b) == cm.base and not on_region:
                    cands.append(pos)
                pos += 2_000_000
            # try the middle first, then work outwards: the centre of a chromosome is the least likely
            # to be centromeric or telomeric and so the most likely to be captured
            for p in sorted(cands, key=lambda x: abs(x - (cands[len(cands) // 2] if cands else 0))):
                lo, hi = max(1, p - 5_000_000), p + 5_000_000
                # The whole window must be neutral, not just its midpoint. Checking only the centre put
                # dataset 01's chr13 baseline at 10-20 Mb, whose midpoint sits on the acrocentric 13p but
                # whose upper half is inside the 13q LOH; since 13p carries almost no capture, nearly every
                # captured base in that window came from the lost region, so the check compared LOH against
                # LOH and read the RB1 loss as ratio 1.04 against an expected 0.65.
                step = 500_000
                probes = list(range(lo, hi + 1, step)) + [hi]
                if any(cm.cn("T", ch, q)[:2] != cm.base for q in probes):
                    continue
                d, n_iv = captured_depth(bam, ch, lo, hi, capture_bed, work, genome=genome)
                if d:
                    return ch, lo, hi, d, n_iv
        return None, None, None, None, 0

    others = [c for c in (f"chr{i}" for i in list(range(1, 23))) if c != chrom]
    b_chrom, b_lo, b_hi, base, n_base = neutral_baseline([chrom] + others)
    if not base:
        rep.add("depth by copy number", False,
                f"no copy-number-neutral, {unit} baseline found on any chromosome")
        return
    gc_base = mean_gc_factor(genome, gc, b_chrom, b_lo, b_hi, capture_bed) if genome else 1.0
    nums, worst = {"baseline": {"depth": round(base, 2), "intervals": n_base,
                               "chrom": b_chrom, "start": b_lo, "end": b_hi,
                               "same_chromosome": b_chrom == chrom}}, 0.0
    for label, s, e in regions:
        want = sum(w for _s, _h, _k, w in source_plan(cm, cm.purity, chrom, (s + e) // 2)) / den
        # fold in the region's GC composition relative to the baseline's
        gc_reg = mean_gc_factor(genome, gc, chrom, s, e, capture_bed) if genome else 1.0
        gc_ratio = (gc_reg / gc_base) if gc_base else 1.0
        want *= gc_ratio
        got, n_iv = captured_depth(bam, chrom, s, e, capture_bed, work, genome=genome)
        if got is None:
            nums[label] = {"skipped": f"only {n_iv} {unit} interval(s), under the measurable-base floor"}
            continue
        ratio = got / base
        nums[label] = {"expected_ratio": round(want, 3), "observed_ratio": round(ratio, 3),
                       "captured_depth": round(got, 2), "baseline_depth": round(base, 2),
                       "intervals": n_iv, "gc_factor_ratio": round(gc_ratio, 4)}
        if want:
            worst = max(worst, abs(ratio / want - 1))
    rep.add("depth by copy number", worst <= 0.25,
            f"largest deviation {worst:.0%} over {len(regions)} segment(s), on {unit} bases", nums)


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


def check_read_map(rep, read_map, record_map):
    """Every read must trace back to a source record, and the interval sets must all be represented.

    The read map is what owner decision D2 promised: a way back from a read to the clone, haplotype and
    interval it came from. A name that resolves to nothing is not a cosmetic defect, it is a hole in the
    truth, and it is exactly the failure that a parse bug in the naming step produced for a tenth of the
    library without anything noticing.
    """
    records, sets = {}, Counter()
    with gzip.open(record_map, "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        iset = head.index("interval_set") if "interval_set" in head else None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            records[f[0]] = f[iset] if iset is not None else "on_target"
    unresolved, n = 0, 0
    with gzip.open(read_map, "rt") as fh:
        fh.readline()
        for line in fh:
            _name, src = line.rstrip("\n").split("\t")
            n += 1
            if src in records:
                sets[records[src]] += 1
            else:
                unresolved += 1
    rep.add("every read traces to a source record", unresolved == 0,
            f"{unresolved:,} of {n:,} reads have no source record",
            {"reads": n, "unresolved": unresolved, "by_interval_set": dict(sets),
             "on_target_fraction": round(sets.get("on_target", 0) / max(1, n), 4)})
    if len(sets) > 1:
        frac = sets.get("on_target", 0) / max(1, n)
        rep.add("off-target reads present and in range", 0.35 <= frac <= 0.85,
                f"{100 * frac:.1f}% of reads are on target",
                {"by_interval_set": dict(sets)})


def check_classes_present(rep, catalog_dir, dataset, junctions_tsv, chrom=None, capture_bed=None,
                          genome_lengths=None, seed=1, clones=None):
    """Junctions placed against the number this chromosome could place.

    The catalog-wide count is the wrong yardstick. A rearrangement reaches an exome only if a breakpoint
    falls inside a capture interval, and whether any do on one chromosome is a draw: dataset 02 has a
    single chr8 SV breakpoint, at 1.45 Mb and off target, and no chr8 fusion at all, so nothing being
    placed is the right answer rather than a failure. Comparing against the whole catalog reported that
    correct behaviour as a defect.
    """
    present = {}
    tables = {}
    for name, path in (("sv", f"{dataset}.svs.tsv"), ("fusion", f"{dataset}.fusions.tsv"),
                       ("virus", f"{dataset}.viruses.tsv")):
        tables[name] = rows(os.path.join(catalog_dir, path))
        present[name] = len(tables[name])
    placed_all = rows(junctions_tsv) if junctions_tsv else []
    # Restrict to this chromosome. The junctions file of a merged library spans every chromosome, so
    # comparing its whole event count against a per-chromosome `expected` is not a comparison at all:
    # on dataset 01 chr13 it put 87 events from 24 chromosomes against 3 expected on chr13, passed on
    # `>=`, and in doing so hid a genuine seed mismatch between the builder and this check. If the file
    # has no chrom column it is a per-chromosome file already and needs no filtering.
    placed = [r for r in placed_all if not r.get("chrom") or not chrom or r["chrom"] == chrom]
    ids = {r["event_id"].split(":")[0] for r in placed}

    expected, detail_extra, evs = None, "", set()
    cols = {"sv": [("chrom", "start"), ("end_chrom", "end")],
            "fusion": [("chrom_5p", "breakpoint_5p"), ("chrom_3p", "breakpoint_3p")],
            "virus": [("chrom", "integration_pos")]}

    def retained_here(r, pos):
        """Whether the clone still holds the haplotype this breakpoint sits on.

        The builder only emits a junction on a copy the clone has: a breakpoint on a lost haplotype
        cannot appear, so counting it would fail the check for behaviour that is correct.
        """
        if clones is None:
            return True
        clone, hap = r.get("clone") or "T", str(r.get("haplotype") or "")
        retained = clones.retained_haplotypes(clone, chrom, pos) or []
        return not (hap.isdigit() and int(hap) not in retained)

    if chrom and capture_bed is None:
        # WGS: there is no capture, so every designed breakpoint on this chromosome is placeable. That is
        # a strictly stronger expectation than the exome's, and skipping the check here instead -- which
        # is what happens if `capture_bed` is simply absent -- would leave the assay with the weakest
        # junction check of the three rather than the strongest.
        evs = set()
        for name, rws in tables.items():
            for r in rws:
                for cc, pc in cols[name]:
                    if r.get(cc) != chrom or not str(r.get(pc) or "").strip():
                        continue
                    if retained_here(r, int(r[pc])):
                        evs.add(r["event_id"])
        expected = len(evs)
        missed = sorted(evs - {r["event_id"].split(":")[0] for r in placed})
        detail_extra = f" of {expected} designed on {chrom} (WGS: all are placeable)"
        if missed:
            detail_extra += f"; not placed: {', '.join(missed[:4])}"
    elif chrom and capture_bed:
        # the same padded, gap-merged intervals the builder places junctions into. Re-deriving them
        # from the raw BED without merging under-counted what is placeable: a breakpoint inside a
        # merged gap is placed by the builder but was not counted here, and the check only passed
        # because it compares with >=.
        # Every interval set the builder generates reads from, not just the bait set. Since off-target
        # bands were added a junction can be placed in one of them -- the one placed on chr8 sits in a
        # 5 kb distal window -- so counting only capture intervals made "placeable" smaller than what
        # was actually placed. The check passed anyway because it compares with >=, which is how the
        # inconsistency stayed invisible.
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import random as _random
        from igi_catalog.pipeline import merged_capture, off_target_bands
        cap_by_chrom = merged_capture(capture_bed, {chrom})
        spans = list(cap_by_chrom.get(chrom, []))
        if genome_lengths and chrom in genome_lengths:
            rng = _random.Random(f"{seed}:{dataset}:{chrom}:offtarget")
            for _name, _rel, band in off_target_bands(cap_by_chrom, chrom, genome_lengths[chrom], rng):
                spans += band
        spans.sort()

        def inside(pos):
            i = bisect.bisect_right(spans, (pos, float("inf")))
            return any(s <= pos <= e for s, e in spans[max(0, i - 2):i + 1])
        evs = set()
        for name, rws in tables.items():
            for r in rws:
                for cc, pc in cols[name]:
                    if r.get(cc) != chrom or not str(r.get(pc) or "").strip():
                        continue
                    pos = int(r[pc])
                    if not inside(pos) or not retained_here(r, pos):
                        continue
                    evs.add(r["event_id"])
        expected = len(evs)
        placed_ids = {r["event_id"].split(":")[0] for r in placed}
        missed = sorted(evs - placed_ids)
        detail_extra = f" of {expected} placeable on {chrom}"
        if missed:
            detail_extra += f"; not placed: {', '.join(missed[:4])}"

    # Equality, not `>=`. Placing more events than are placeable is as much a defect as placing fewer:
    # it means the two sides disagree about which windows the library was built from. `>=` is what let a
    # seed mismatch sit undetected.
    ok = len(ids) == expected if expected is not None else (bool(placed) or present["sv"] == 0)
    extra = sorted(ids - evs) if expected is not None else []
    if extra:
        detail_extra += f"; placed but not expected: {', '.join(extra[:4])}"
    rep.add("junction classes reach the library", ok,
            f"{len(placed)} junctions placed covering {len(ids)} events{detail_extra}",
            {"catalog": present, "placed": len(placed), "placeable_on_chrom": expected,
             "events_on_chrom": len(ids), "unexpected": len(extra)})


def check_rna(rep, manifest, expected_sources=None):
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
    # "at least three of six" counted `reference` and `virus` towards the total, and an episomal virus is
    # not tied to a chromosome, so a slice with no ERV, splice isoform or CTA whatsoever satisfied it. The
    # chr6 slice has exactly none of those three and passed. Expect the classes the catalog actually
    # places on this chromosome, and say which are missing.
    expect = set(expected_sources or ()) or {"reference"}
    missing = sorted(expect - set(srcs))
    rep.add("rna carries the designed classes", not missing,
            f"sources present {dict(srcs)}"
            + (f"; expected but absent: {', '.join(missing)}" if missing else ""),
            {"expected": sorted(expect), "present": dict(srcs)})


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
    ap.add_argument("--capture-bed",
                    help="capture intervals, for a WES library. Omit it with --assay wgs; depth is then "
                         "measured over non-N reference instead")
    ap.add_argument("--assay", default="wes", choices=["wes", "wgs"],
                    help="wgs drops the capture bed from the depth and placeability checks. The GC curve "
                         "is capture efficiency and is not applied; every designed site is placeable")
    ap.add_argument("--junctions")
    ap.add_argument("--rna-manifest")
    ap.add_argument("--read-map")
    ap.add_argument("--record-map")
    ap.add_argument("--work", default="/tmp")
    ap.add_argument("--seed", type=int, default=1,
                    help="the seed the library was built with; the off-target windows depend on it. "
                         "Prefer --build-meta, which reads it from the builder's own output")
    ap.add_argument("--build-meta",
                    help="the builder's per-chromosome JSON (<dataset>_<chrom>_<library>.json). The seed "
                         "is taken from it, so the checker cannot be told a seed the library was not "
                         "built with. Overrides --seed and warns if the two disagree")
    ap.add_argument("--paths", default=None,
                    help="site paths yaml; supplies the GC bias curve and reference for the depth check")
    ap.add_argument("--report", required=True)
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    if a.build_meta:
        _m = json.load(open(a.build_meta))
        if "seed" not in _m:
            raise SystemExit(f"{a.build_meta} records no seed; it predates seed recording, so pass "
                             f"--seed explicitly and make sure it matches how the library was built")
        for k in ("dataset", "chrom", "library"):
            if a.__dict__.get(k) and _m.get(k) and k != "library" and str(_m[k]) != str(a.__dict__[k]):
                raise SystemExit(f"--build-meta is for {k}={_m[k]} but --{k} says {a.__dict__[k]}")
        if a.seed != 1 and a.seed != _m["seed"]:
            print(f"[acceptance] --seed {a.seed} disagrees with the build's {_m['seed']}; "
                  f"using {_m['seed']}", file=sys.stderr)
        a.seed = _m["seed"]
    rep = Report()
    if a.bam and a.reference and a.dataset and a.chrom:
        cat = os.path.join(a.catalog_dir, f"{a.dataset}.snv_indel.tsv")
        check_variants(rep, a.bam, a.reference, cat, a.chrom, a.work)
        cap = None if a.assay == "wgs" else a.capture_bed
        if a.assay == "wgs" and a.capture_bed:
            print("[acceptance] --assay wgs: ignoring --capture-bed", file=sys.stderr)
        if a.arms_bed and (cap or a.assay == "wgs"):
            check_depth_by_cn(rep, a.bam, a.design, a.dataset, a.chrom, a.arms_bed,
                              cap, a.work, gc_curve=a.paths)
        lengths = None
        if a.paths:
            import pysam as _pysam
            _p = yaml.safe_load(open(a.paths))
            _fa = _pysam.FastaFile(_p['reference_fasta'])
            lengths = dict(zip(_fa.references, _fa.lengths))
        cmodel = None
        if a.arms_bed:
            sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
            from igi_catalog.clones import CloneModel
            cmodel = CloneModel(yaml.safe_load(open(a.design))["datasets"][a.dataset], a.arms_bed)
        check_classes_present(rep, a.catalog_dir, a.dataset, a.junctions, a.chrom,
                              cap, lengths, a.seed, cmodel)
    if a.r1:
        check_read_names(rep, a.r1)
    if a.read_map and a.record_map:
        check_read_map(rep, a.read_map, a.record_map)
    if a.rna_manifest:
        # which designed classes this slice could contain, from the catalog restricted to its chromosomes
        exp = {"reference"}
        if a.chrom:
            chroms = set(a.chrom.split(","))
            for cls, src in (("erv", "erv"), ("splice", "splice_isoform"), ("cta", "cta")):
                for r in rows(os.path.join(a.catalog_dir, f"{a.dataset}.expressed.tsv")):
                    if r.get("class") == cls and r.get("chrom") in chroms:
                        exp.add(src)
                        break
            for r in rows(os.path.join(a.catalog_dir, f"{a.dataset}.fusions.tsv")):
                if r.get("chrom_5p") in chroms:
                    exp.add("fusion")
                    break
        check_rna(rep, a.rna_manifest, exp)
    rep.write(a.report)
    if rep.failed():
        print(f"\n{len(rep.failed())} check(s) failed", flush=True)
        sys.exit(1)
    print(f"\nall {len(rep.checks)} checks passed", flush=True)


if __name__ == "__main__":
    main()
