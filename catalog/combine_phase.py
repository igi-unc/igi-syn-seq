#!/usr/bin/env python3
"""Combine SHAPEIT5 statistical phase (panel sites) with WhatsHap read-backed blocks (all sites).

For every heterozygous PASS site of the sample:
  1. site present in SHAPEIT5 output      -> take SHAPEIT5 haplotypes (chromosome-wide phase set)
  2. else site inside a WhatsHap block that contains >=1 SHAPEIT5-phased site
                                          -> orient the WhatsHap block to agree with SHAPEIT5 and place the site
  3. else                                 -> random orientation (seeded), flagged PHASE_SOURCE=random
Homozygous sites are emitted phased trivially. Haploid regions (chrX non-PAR, chrY in a male) are emitted
as haploid. Output INFO/PHASE_SOURCE in {shapeit5, whatshap_block, random, hom, haploid}.
Usage: combine_phase.py --whatshap all.phased.vcf.gz --shapeit5 chrN.bcf [...] --out out.vcf --seed 1 --sex male
"""
import argparse, random, sys
import pysam

PAR_GRCH38 = {"chrX": [(10001, 2781479), (155701383, 156030895)]}

def in_par(chrom, pos):
    return any(a <= pos <= b for a, b in PAR_GRCH38.get(chrom, []))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--whatshap", required=True)
    ap.add_argument("--shapeit5", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--sex", choices=["male", "female"], default="male")
    a = ap.parse_args()
    rng = random.Random(a.seed)

    # SHAPEIT5 haplotypes keyed by (chrom, pos, ref, alt) -> (h0, h1)
    s5 = {}
    for f in a.shapeit5:
        for rec in pysam.VariantFile(f):
            gt = rec.samples[0]["GT"]
            if gt is None or len(gt) != 2 or None in gt:
                continue
            s5[(rec.chrom, rec.pos, rec.ref, rec.alts[0])] = (gt[0], gt[1])

    wh = pysam.VariantFile(a.whatshap)
    # add the tag to the reader's own header: records parsed afterwards reference this header, so INFO can be set on them
    if "PHASE_SOURCE" not in wh.header.info:
        wh.header.info.add("PHASE_SOURCE", 1, "String",
                           "Origin of the phase: shapeit5, whatshap_block, random, hom, haploid")
    out = pysam.VariantFile(a.out, "w", header=wh.header)
    counts = {}

    # pass 1: block orientation from WhatsHap PS blocks vs SHAPEIT5
    votes = {}  # (chrom, ps) -> [agree, disagree]
    for rec in wh:
        s = rec.samples[0]
        gt = s["GT"]
        if gt is None or None in gt or len(gt) != 2 or gt[0] == gt[1] or not s.phased:
            continue
        key = (rec.chrom, rec.pos, rec.ref, rec.alts[0])
        if key in s5 and rec.chrom in wh.header.contigs:
            ps = s.get("PS")
            if ps is None:
                continue
            v = votes.setdefault((rec.chrom, ps), [0, 0])
            v[0 if tuple(gt) == s5[key] else 1] += 1
    flip = {k: (v[1] > v[0]) for k, v in votes.items()}
    wh.reset()

    for rec in wh:
        s = rec.samples[0]
        gt = s["GT"]
        src = None
        if gt is None or None in gt:
            out.write(rec); continue
        haploid = a.sex == "male" and ((rec.chrom == "chrX" and not in_par(rec.chrom, rec.pos)) or rec.chrom == "chrY")
        if haploid:
            allele = max(set(gt), key=gt.count) if len(set(gt)) > 1 else gt[0]
            s["GT"] = (allele,); src = "haploid"
        elif len(gt) == 2 and gt[0] == gt[1]:
            s["GT"] = (gt[0], gt[1]); s.phased = True; src = "hom"
        else:
            key = (rec.chrom, rec.pos, rec.ref, rec.alts[0])
            if key in s5:
                s["GT"] = s5[key]; s.phased = True; src = "shapeit5"
            else:
                ps = s.get("PS") if s.phased else None
                if ps is not None and (rec.chrom, ps) in flip:
                    g = tuple(gt)
                    s["GT"] = (g[1], g[0]) if flip[(rec.chrom, ps)] else g
                    s.phased = True; src = "whatshap_block"
                else:
                    g = tuple(gt)
                    s["GT"] = g if rng.random() < 0.5 else (g[1], g[0])
                    s.phased = True; src = "random"
            s["PS"] = 1  # every phased het joins the single chromosome-wide phase set
        rec.info["PHASE_SOURCE"] = src
        counts[src] = counts.get(src, 0) + 1
        out.write(rec)
    out.close()
    for k in sorted(counts):
        print(f"{k}\t{counts[k]}", file=sys.stderr)

if __name__ == "__main__":
    main()
