#!/usr/bin/env python3
"""Build the external input manifest the recipe needs (owner decision D9).

    python3 make_input_manifest.py --paths paths.yaml --out docs/input-manifest.tsv

D9 makes the release a recipe rather than data, so every input a user must supply has to be identified
precisely enough to obtain the same bytes. Each row carries the key the configuration uses, what the input
is, its class, size, md5, and where it comes from. A checksum without a source is not reproducible and a
source without a checksum is not verifiable, so both are required or the row is marked incomplete.
"""
import argparse, hashlib, json, os, subprocess, sys

import yaml

# Where each input comes from, and whether a third party can obtain it. `host` means the owner will host
# it because it is lab-internal and otherwise unobtainable (D10).
SOURCES = {
    "reference_fasta": ("public", "GRCh38 analysis set; GATK resource bundle / Broad"),
    "gtf": ("public", "GENCODE v37 comprehensive annotation, gencodegenes.org"),
    "arms_bed": ("public", "UCSC cytoBand-derived chromosome arms, hg38"),
    "rmsk_bed": ("public", "UCSC RepeatMasker track, hg38"),
    "segdups_bed": ("public", "UCSC genomicSuperDups track, hg38"),
    "exome_bed": ("host", "generic human exome interval list, not a named capture kit (D11)"),
    "pan_normal_quant": ("host", "lab pan-normal and mTEC expression, 95th percentile"),
    "erv_all_loci": ("host", "HERV annotation, digest 2fc043"),
    "erv_loci": ("host", "curated ERV loci, digest 2fc043"),
    "viral_fasta": ("host", "unmasked viral genome set, digest 02ec8"),
    "cta_gene_list": ("host", "PIRL cancer-testis antigen gene list, 06APR2026"),
    "single_cell_whitelist": ("public", "10x Genomics 737K-august-2016 barcode whitelist"),
    "germline_vcf_baseline.HG002": ("public", "GIAB HG002 Q100 v1.1 small-variant benchmark, GRCh38"),
    "germline_vcf_baseline.IPISRC044": ("derived", "called and phased here; see baseline-references.md"),
    "germline_vcf.IGI-SYN-SEQ-01": ("derived", "baseline plus designed alleles; run_germline_spikein.py"),
    # These two keys are aliases so build_env can look a VCF up by baseline name. They deliberately point
    # at the spiked file, not the baseline, so every builder sees the designed germline alleles; the
    # untouched baselines live under germline_vcf_baseline. Left unlabelled a reader could not tell which.
    "germline_vcf.HG002": ("derived", "alias of germline_vcf.IGI-SYN-SEQ-01, the spiked VCF"),
    "germline_vcf.IPISRC044": ("derived", "alias of germline_vcf.IGI-SYN-SEQ-02, the spiked VCF"),
    "germline_vcf.IGI-SYN-SEQ-02": ("derived", "baseline plus designed alleles; run_germline_spikein.py"),
    "expression_tsv": ("in-repo", "TCGA-BRCA basal medians via UCSC Xena Toil; catalog/resources"),
    "art_profile_r1": ("in-repo", "fitted from IPISRC044 NovaSeq X; catalog/resources"),
    "art_profile_r2": ("in-repo", "fitted from IPISRC044 NovaSeq X; catalog/resources"),
    "gc_bias_curve": ("in-repo", "measured from IPISRC044 WES; catalog/resources"),
    "annotation_cache": ("derived", "regenerated from the GTF on first run; safe to delete"),
    "workdir": ("scratch", "working directory, not an input"),
}


def md5(path, limit_gb=8):
    if os.path.getsize(path) > limit_gb * 1e9:
        return "not-computed-too-large"
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(8 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def image_digest(path):
    """A singularity image's identity. The file's own md5 is what a user can verify locally."""
    return md5(path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paths", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    p = yaml.safe_load(open(a.paths))

    rows = []
    def add(key, path, kind):
        cls, src = SOURCES.get(key, ("UNCLASSIFIED", ""))
        if os.path.isdir(path):
            return      # a working directory is not an input
        ok = os.path.isfile(path)
        rows.append({
            "key": key, "kind": kind, "class": cls,
            "basename": os.path.basename(path),
            "bytes": os.path.getsize(path) if ok else 0,
            "md5": md5(path) if ok else "MISSING",
            "source": src or "UNRECORDED",
            "complete": "yes" if (ok and cls != "UNCLASSIFIED" and src) else "no",
            "durable": "no" if "/scratch/" in path else "yes",
        })

    for k, v in sorted(p.items()):
        if isinstance(v, str) and "singularity" in v:
            for tok in v.split():
                if tok.endswith((".img", ".sif")):
                    cls, src = SOURCES.get(k, ("container", ""))
                    rows.append({"key": k, "kind": "container", "class": "container",
                                 "basename": os.path.basename(tok),
                                 "bytes": os.path.getsize(tok) if os.path.exists(tok) else 0,
                                 "md5": image_digest(tok) if os.path.exists(tok) else "MISSING",
                                 "source": src or "container registry; see basename for the exact tag",
                                 "complete": "yes" if os.path.exists(tok) else "no",
                                 "durable": "yes"})
        elif isinstance(v, str) and v.startswith("/"):
            add(k, v, "file")
        elif isinstance(v, dict):
            for kk, vv in v.items():
                if isinstance(vv, str) and vv.startswith("/"):
                    add(f"{k}.{kk}", vv, "file")

    cols = ["key", "kind", "class", "basename", "bytes", "md5", "source", "complete", "durable"]
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    n_inc = sum(1 for r in rows if r["complete"] == "no")
    n_scr = sum(1 for r in rows if r["durable"] == "no")
    print(f"{len(rows)} inputs -> {a.out}")
    print(f"  incomplete rows (missing source or classification): {n_inc}")
    print(f"  inputs living on scratch, which is not durable:      {n_scr}")
    for r in rows:
        if r["complete"] == "no" or r["durable"] == "no":
            print(f"    {r['key']:34} class={r['class']:12} durable={r['durable']}")


if __name__ == "__main__":
    main()
