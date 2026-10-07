#!/usr/bin/env python3
"""F24: is the simulator's indel error rate in homopolymers the same as real HiFi's?

    homopolymer_indel_rate.py --bam aligned.bam --ref GRCh38.fa --out out.json
                              [--exclude-vcf truth.vcf.gz] [--min-run 6] [--max-reads 20000]

Real HiFi accuracy comes from consensus over multiple passes of one molecule, and its residual errors
concentrate in homopolymers and tandem repeats. Badread applies an error model per read with no pass
structure, so matching marginal identity (0.99830 against a real 0.99823) does not establish that errors
land in the same PLACES. That matters because the design places 40 homopolymer indels per dataset whose
stated purpose is to be PacBio-hard: if the simulated reads are clean through homopolymer runs, that tier
is populated but hollow.

The measurement is the indel rate inside reference homopolymer runs of length >= `--min-run` over the rate
outside them, computed the same way for simulated and real reads, and compared as a ratio.

**Germline indels have to be excluded or they dominate.** A real HiFi indel error rate is of order 1e-5 per
base; human germline indels are of order 5e-4 per base, fifty times larger. Measuring raw indels would
therefore mostly measure biology, and biology is itself homopolymer-enriched, which would flatter the
simulator. `--exclude-vcf` drops any indel within `--pad` of a known truth indel. For the simulated arm
the cleaner route is to simulate from the plain reference, so every indel is by construction an error and
no exclusion is needed.
"""
import argparse
import collections
import gzip
import json

import pysam


def homopolymer_mask(seq, min_run):
    """Boolean list: is this base inside a homopolymer run of at least `min_run`?"""
    n = len(seq)
    mask = bytearray(n)
    i = 0
    s = seq.upper()
    while i < n:
        j = i + 1
        while j < n and s[j] == s[i]:
            j += 1
        if s[i] in "ACGT" and (j - i) >= min_run:
            for k in range(i, j):
                mask[k] = 1
        i = j
    return mask


def truth_indels(vcf_path, chroms):
    """{chrom: sorted [pos]} for every indel in a truth VCF, 1-based."""
    if not vcf_path:
        return {}
    out = collections.defaultdict(list)
    op = gzip.open if str(vcf_path).endswith(".gz") else open
    with op(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 5)
            if len(f) < 5 or f[0] not in chroms:
                continue
            ref, alts = f[3], f[4].split(",")
            if any(len(a) != len(ref) for a in alts):
                out[f[0]].append(int(f[1]))
    for c in out:
        out[c].sort()
    return out


def near(sorted_pos, pos, pad):
    import bisect
    i = bisect.bisect_left(sorted_pos, pos - pad)
    return i < len(sorted_pos) and sorted_pos[i] <= pos + pad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--exclude-vcf", default=None)
    ap.add_argument("--min-run", type=int, default=6)
    ap.add_argument("--pad", type=int, default=10)
    ap.add_argument("--max-reads", type=int, default=20000)
    ap.add_argument("--min-mapq", type=int, default=0,
                    help="skip alignments below this MAPQ. Repeats are homopolymer- and tandem-rich, so "
                         "mismapped reads land preferentially in the homopolymer bucket and inflate the "
                         "enrichment. Comparing a genome-wide real alignment against reads simulated from "
                         "one clean window without this filter is not a like-for-like comparison.")
    ap.add_argument("--region", default=None,
                    help="chrom:start-end, to hold the reference context identical across arms")
    ap.add_argument("--label", default="")
    a = ap.parse_args()

    fa = pysam.FastaFile(a.ref)
    bam = pysam.AlignmentFile(a.bam, "rb")
    chroms = {c for c in fa.references if "_" not in c and c not in ("chrEBV",)}
    excl = truth_indels(a.exclude_vcf, chroms)

    # denominators are aligned reference bases; numerators are indel EVENTS and indel BASES
    base = {"hp": 0, "non": 0}
    ev = {"hp": 0, "non": 0}
    bp = {"hp": 0, "non": 0}
    n_reads = n_excluded = n_lowmapq = 0
    reg = None
    if a.region:
        c, span = a.region.split(":")
        lo, hi = span.split("-")
        reg = (c, int(lo), int(hi))
    runlen_hist = collections.Counter()

    itr = (bam.fetch(reg[0], reg[1], reg[2]) if reg else bam.fetch(until_eof=True))
    for rec in itr:
        if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
            continue
        if rec.reference_name not in chroms:
            continue
        if rec.mapping_quality < a.min_mapq:
            n_lowmapq += 1
            continue
        if reg and not (rec.reference_name == reg[0]
                        and rec.reference_start >= reg[1] and rec.reference_end <= reg[2]):
            continue
        n_reads += 1
        if n_reads > a.max_reads:
            break
        rstart, rend = rec.reference_start, rec.reference_end
        if rend is None or rend - rstart < 1000:
            continue
        ref = fa.fetch(rec.reference_name, rstart, rend)
        mask = homopolymer_mask(ref, a.min_run)
        ex = excl.get(rec.reference_name, [])

        pos = rstart
        for op, ln in rec.cigartuples or []:
            if op in (0, 7, 8):                     # M/=/X consume both
                for k in range(ln):
                    (base.__setitem__("hp", base["hp"] + 1) if mask[pos - rstart + k]
                     else base.__setitem__("non", base["non"] + 1))
                pos += ln
            elif op == 2:                           # D consumes reference
                key = "hp" if mask[min(pos - rstart, len(mask) - 1)] else "non"
                if ex and near(ex, pos + 1, a.pad):
                    n_excluded += 1
                else:
                    ev[key] += 1
                    bp[key] += ln
                    if key == "hp":
                        runlen_hist[ln] += 1
                pos += ln
            elif op == 1:                           # I consumes query only
                key = "hp" if mask[min(max(pos - rstart - 1, 0), len(mask) - 1)] else "non"
                if ex and near(ex, pos + 1, a.pad):
                    n_excluded += 1
                else:
                    ev[key] += 1
                    bp[key] += ln
                    if key == "hp":
                        runlen_hist[ln] += 1
            elif op in (3, 4, 5, 6):                # N/S/H/P
                if op == 3:
                    pos += ln

    def rate(num, den):
        return (num / den) if den else None

    res = {
        "label": a.label, "bam": a.bam, "min_run": a.min_run,
        "reads_used": min(n_reads, a.max_reads),
        "excluded_vcf": a.exclude_vcf, "indels_excluded_as_truth": n_excluded,
        "min_mapq": a.min_mapq, "reads_below_mapq": n_lowmapq, "region": a.region,
        "aligned_ref_bases": base,
        "indel_events": ev, "indel_bases": bp,
        "event_rate_per_base": {k: rate(ev[k], base[k]) for k in ev},
        "base_rate_per_base": {k: rate(bp[k], base[k]) for k in bp},
        "hp_indel_run_lengths": dict(sorted(runlen_hist.items())[:12]),
    }
    r_hp = res["event_rate_per_base"]["hp"]
    r_non = res["event_rate_per_base"]["non"]
    res["enrichment_hp_over_non"] = (r_hp / r_non) if (r_hp and r_non) else None
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=2)
        fh.write("\n")
    print(f"  {a.label or a.bam}")
    print(f"    reads {res['reads_used']:,}  aligned bases hp {base['hp']:,} non {base['non']:,}")
    print(f"    indel events   hp {ev['hp']:,}  non {ev['non']:,}  (truth-excluded {n_excluded:,})")
    print(f"    indel rate/bp  hp {r_hp:.3e}  non {r_non:.3e}" if r_hp and r_non else "    rate unavailable")
    print(f"    ENRICHMENT hp/non = {res['enrichment_hp_over_non']:.2f}x"
          if res["enrichment_hp_over_non"] else "    enrichment unavailable")


if __name__ == "__main__":
    main()
