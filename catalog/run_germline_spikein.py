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
    env = build_env(paths, design, a.dataset)
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

    os.makedirs(a.out, exist_ok=True)
    write_table(rows, os.path.join(a.out, f"{a.dataset}.germline_spikein.tsv"))
    with open(os.path.join(a.out, f"{a.dataset}.germline_spikein.diffcards.txt"), "w") as fh:
        fh.write(f"# {a.dataset} designed germline alleles (design v{design['version']})\n\n")
        fh.write("\n".join(card(r) for r in rows))
    print(f"[done] {len(rows)} germline allele(s) -> {a.out}")


if __name__ == "__main__":
    main()
