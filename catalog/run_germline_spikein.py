#!/usr/bin/env python3
"""Add the dataset's pathogenic germline alleles to its baseline VCF.

    python3 run_germline_spikein.py --design design.yaml --paths paths.yaml \
        --dataset IGI-SYN-SEQ-01 --out output/ --vcf-out /path/to/IGI-SYN-SEQ-01.germline.vcf.gz

Writes a new phased VCF containing the baseline genotype plus the designed alleles, a truth table, and a
diff card. The alleles are resolved from coding coordinates and checked against the expected codon,
residue and consequence, so a change of annotation release fails loudly instead of moving the variant.
"""
import argparse, json, os, subprocess, sys, tempfile

import pysam
import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from igi_catalog.designer import build_env, write_table
from igi_catalog.germline import Germline
from igi_catalog.germline_spikein import GermlineSpikein
from igi_catalog import diffcards


# GRCh38 class I windows, hardcoded rather than read from the GTF: the annotation is not an input to this
# script, and these three spans are stable coordinates.
HLA_CLASS_I = {"HLA-A": ("chr6", 29942532, 29945870),
               "HLA-B": ("chr6", 31353875, 31357188),
               "HLA-C": ("chr6", 31268749, 31272130)}


def check_hla_phasing(vcf, dcfg):
    """Assert the VCF just written can actually express the HLA genotype the design declares.

    This exists because it could not. IGI-SYN-SEQ-02's germline VCF was built from a baseline whose MHC had
    never been recovered -- alt-aware alignment had parked the HLA reads on GRCh38's 525 HLA-* contigs -- so
    chr6:29.6-31.4 Mb held 214 records and 6 heterozygous calls against ds-01's 12,212 and 6,597, and
    HLA-A, -B and -C held NONE AT ALL. Both simulated haplotypes were byte-identical bare GRCh38 at all
    three loci, measured through haplotype_seqs: the delivered reads were not B*08:01 or B*27:05, they were
    the reference's own allele, in every ds-02 arm. HLA typing could not return the designed genotype and
    the designed HLA LOH had no allelic signal to be seen in.

    The recovery had been run, verified, merged and re-phased 13 h before that VCF was built. Only the
    `germline_vcf_baseline` key still named the superseded file. The recovery job asserts its own output;
    nothing asserted that the output was ever USED, which is the gap this closes -- at the one step that
    holds both the declared genotype and the file being written.

    Phase is checked too, and it is not a formality. `Germline.haplotype_seqs` assigns alleles by GT ORDER
    and never consults the phase flag, so unphased 0/1 calls would put every alternate allele on hap1 and
    leave hap0 bare reference: one chimera carrying both real alleles' variants, and -- since design.yaml
    loses hap1 to HLA LOH -- a tumour reading as homozygous reference. The record count would look right.
    """
    by_locus = {}
    for allele in dcfg.get("hla") or []:
        by_locus.setdefault(allele.split("*")[0], []).append(allele)
    problems = []
    for locus, (chrom, lo, hi) in sorted(HLA_CLASS_I.items()):
        alleles = by_locus.get(locus) or []
        if len(alleles) != 2:
            continue                      # not a declared class I pair; nothing to assert
        het_expected = alleles[0] != alleles[1]
        reg = f"{chrom}:{lo}-{hi}"
        r = subprocess.run(["bcftools", "view", "-H", "-r", reg, vcf], capture_output=True, text=True)
        if r.returncode != 0:
            raise SystemExit(f"[germline] could not read {reg} from {vcf}: {r.stderr.strip()[-200:]}")
        n_het = n_phased = 0
        for line in r.stdout.splitlines():
            cols = line.split("\t")
            if len(cols) < 10:
                continue
            gt = cols[9].split(":")[0]
            a1, a2 = (gt.replace("|", "/").split("/") + ["."])[:2]
            if a1 != a2 and "." not in (a1, a2):
                n_het += 1
                n_phased += ("|" in gt)
        tag = "heterozygous by design" if het_expected else "homozygous by design"
        print(f"[germline] {locus} {'/'.join(alleles)} ({tag}): {n_het} heterozygous, {n_phased} phased",
              flush=True)
        if het_expected and n_het == 0:
            problems.append(f"{locus} is {'/'.join(alleles)} in the design, but the VCF has no "
                            f"heterozygous call in {reg}; both haplotypes would be bare reference there, "
                            f"so the reads would carry neither allele")
        elif het_expected and n_phased < n_het:
            problems.append(f"{locus} has {n_het - n_phased} unphased heterozygous call(s) in {reg}; "
                            f"haplotype_seqs reads GT order, so those alleles would all land on hap1")
    if problems:
        raise SystemExit("[germline] FAIL: this VCF cannot express the declared HLA genotype:\n  - "
                         + "\n  - ".join(problems))


def write_spike_vcf(base_vcf, rows, path):
    """A VCF with the baseline's header and sample carrying only the designed records."""
    src = pysam.VariantFile(base_vcf)
    hdr = src.header.copy()
    hdr.add_line('##INFO=<ID=IGISYN,Number=1,Type=String,Description="IGI-SYN-SEQ designed germline allele">')
    out = pysam.VariantFile(path, "wz", header=hdr)
    sample = src.header.samples[0]
    for r in sorted(rows, key=lambda x: (x["chrom"], x["pos"])):
        rec = out.new_record(contig=r["chrom"], start=r["pos"] - 1,
                             alleles=(r["ref"], r["alt"]), id=r["event_id"])
        rec.info["IGISYN"] = r["event_id"]
        rec.samples[sample]["GT"] = (1, 0) if r["haplotype"] == 0 else (0, 1)
        rec.samples[sample].phased = True
        out.write(rec)
    out.close()
    src.close()
    pysam.tabix_index(path, preset="vcf", force=True)
    return path


def concat(base_vcf, spike_vcf, out_vcf):
    subprocess.run(["bcftools", "concat", "-a", "-Oz", "-o", out_vcf, base_vcf, spike_vcf], check=True)
    subprocess.run(["bcftools", "index", "-t", "-f", out_vcf], check=True)
    return out_vcf


def card(r):
    lines = [f"  germline  {r['gene']} {r['cds_change']}  {r['protein_change']}  ({r['consequence']})",
             f"  DNA       {r['chrom']}:{r['pos']} {r['ref']} > {r['alt']}  GT {r['genotype']} "
             f"(haplotype {r['haplotype']})",
             f"  tumour    somatic {r['somatic_loh_label']} removes the other copy, leaving this allele alone"]
    w, m, k = r["wt_protein"], r["mut_protein"], r["aa_index"]
    lo = max(0, k - 10)
    prot = [f"  protein   wild type {r['wt_protein_len']} aa -> mutant {r['mut_protein_len']} aa",
            f"  wt        ...{w[lo:k + 25]}...",
            f"  mut       ...{m[lo:k + 25]}...",
            f"            {' ' * (3 + k - lo)}^ first changed residue (aa {k + 1})"]
    return diffcards.render(r["event_id"], f"{r['gene']} {r['cds_change']} germline", lines, prot,
                            {"note": r["note"], "consequence": r["consequence"]})


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", required=True)
    ap.add_argument("--paths", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--vcf-out", required=True, help="the spiked germline VCF (bgzipped)")
    a = ap.parse_args(argv)

    design = yaml.safe_load(open(a.design))
    paths = yaml.safe_load(open(a.paths))
    dcfg = design["datasets"][a.dataset]
    specs = dcfg.get("germline_spikein") or []
    if not specs:
        print(f"[germline] {a.dataset}: nothing to add")
        return
    env = build_env(paths, design, a.dataset, verify_hla=False)
    # always build from the untouched baseline, never from a VCF this script already wrote: `germline_vcf`
    # points at the spiked file once it exists, so working from that key would either double-spike or
    # report the allele this script added as a pre-existing conflict
    base = paths.get("germline_vcf_baseline", paths["germline_vcf"])[dcfg["baseline"]]
    env.germline = Germline(base, sex=dcfg["sex"])
    gs = GermlineSpikein(env, a.dataset, dcfg)
    rows = []
    for spec in specs:
        r = gs.resolve(spec)
        print(f"[germline] {r['event_id']}: {r['gene']} {r['cds_change']} -> {r['chrom']}:{r['pos']} "
              f"{r['ref']}>{r['alt']} GT {r['genotype']} {r['protein_change']} "
              f"({r['wt_protein_len']} aa -> {r['mut_protein_len']} aa)", flush=True)
        rows.append(r)

    os.makedirs(os.path.dirname(os.path.abspath(a.vcf_out)), exist_ok=True)
    with tempfile.TemporaryDirectory(dir=env.work) as td:
        spike = write_spike_vcf(base, rows, os.path.join(td, "spikein.vcf.gz"))
        concat(base, spike, a.vcf_out)
    n_base = int(subprocess.run(["bcftools", "index", "-n", base], capture_output=True, text=True).stdout or 0)
    n_out = int(subprocess.run(["bcftools", "index", "-n", a.vcf_out], capture_output=True, text=True).stdout or 0)
    print(f"[germline] {base} ({n_base:,} records) + {len(rows)} -> {a.vcf_out} ({n_out:,} records)")
    if n_out != n_base + len(rows):
        raise SystemExit(f"record count mismatch: expected {n_base + len(rows)}, got {n_out}")
    check_hla_phasing(a.vcf_out, dcfg)

    os.makedirs(a.out, exist_ok=True)
    write_table(rows, os.path.join(a.out, f"{a.dataset}.germline_spikein.tsv"))
    with open(os.path.join(a.out, f"{a.dataset}.germline_spikein.diffcards.txt"), "w") as fh:
        fh.write(f"# {a.dataset} designed germline alleles (design v{design['version']})\n\n")
        fh.write("\n".join(card(r) for r in rows))
    print(f"[done] {len(rows)} germline allele(s) -> {a.out}")


if __name__ == "__main__":
    main()
