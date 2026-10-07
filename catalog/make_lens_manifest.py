#!/usr/bin/env python3
"""Emit the LENS sample manifests for the IGI-SYN-SEQ datasets.

    make_lens_manifest.py --design design.yaml --out-dir <dir> [--release full]

Column contract from uselens.io "Manifest specifications" (lens-v1.9.2):

    required   Dataset, Patient_Name, Run_Name, File_Prefix, Sequencing_Method, Normal
    optional   Alleles (required for mouse, optional for human), Group
    Run_Name   a two-letter prefix then `-`: first letter a (abnormal) or n (normal),
               second r (RNA) or d (DNA)
    Alleles    six comma-separated four-digit alleles, two each of HLA-A, -B, -C
    Normal     TRUE or FALSE, uppercase
    File_Prefix  RAFT gathers every FASTQ matching it and pairs mates on anchored
                 _1/_2 or _R1/_R2 suffixes

Two files are written, one with Alleles and one without: with it, LENS uses the designed
germline haplotypes and HLA LOH is something it has to detect; without it, LENS types the
HLA itself from the normal, which is the harder test and the one that exercises OptiType.

Group is a single value per dataset. The design has one synthetic patient and one timepoint,
and LENS v2.0.0-dev discriminates assays by Sequencing_Method, so the assay does not need to
be encoded here.
"""
import argparse
import os

import yaml

# (sample label, Sequencing_Method, normal?, RNA?, file-prefix template)
#
# The prefixes are the staged symlink names, which are deliberately not all the release names:
# WES is linked with a _wes_ infix because "<ds>_tumor" would otherwise also match
# "<ds>_tumor_wgs_R1" and "<ds>_tumor_ont_wgs", and one row would claim three assays.
SHORT_READ = [
    ("wes_tumor",   "WES",       False, False, "{d}_tumor_wes"),
    ("wes_normal",  "WES",       True,  False, "{d}_normal_wes"),
    ("rna_tumor",   "RNA-Seq",   False, True,  "{d}_full_rna"),
    ("wgs_tumor",   "WGS",       False, False, "{d}_tumor_wgs"),
    ("wgs_normal",  "WGS",       True,  False, "{d}_normal_wgs"),
    ("gex_tumor",   "scRNA-Seq", False, True,  "{d}-GEX_S1_L001"),
    ("tcr_tumor",   "scTCR-Seq", False, True,  "{d}-TCR_S1_L001"),
]

# Long-read arms. ONT is delivered as FASTQ and goes in the same way. PacBio HiFi and Kinnex
# are delivered as unaligned BAM, which the manifest specification has no column for, so they
# are held out of the manifest rather than guessed at; see the README this writes.
LONG_READ = [
    # (label, Sequencing_Method, normal?, RNA?, File_Prefix template, Platform)
    ("ont_wgs_tumor",   "WGS",       False, False, "{d}_tumor_ont_wgs",            "ONT"),
    ("ont_wgs_normal",  "WGS",       True,  False, "{d}_normal_ont_wgs",           "ONT"),
    ("ont_bulk_rna",    "RNA-Seq",   False, True,  "{d}_full_ont_bulk_rna",        "ONT"),
    ("ont_sc_rna",      "scRNA-Seq", False, True,  "{d}_full_ont_sc_rna",          "ONT"),
    # BAM-delivered. File_Prefix is the FULL .bam filename, matching the
    # IPISRC044_T1_sclrs_live.seg.bam row in the v2.0.0-dev manifest. preflight_resolve_inputs
    # prefix-matches and then keeps only names ending in .bam, so the .pbi beside each one is
    # excluded, and it fails loudly if a prefix matches more than one BAM.
    ("hifi_wgs_tumor",  "WGS",       False, False, "{d}_tumor_hifi.bam",           "PacBio"),
    ("hifi_wgs_normal", "WGS",       True,  False, "{d}_normal_hifi.bam",          "PacBio"),
    ("kinnex_bulk",     "RNA-Seq",   False, True,  "{d}_full_kinnex_bulk_segmented.bam", "PacBio"),
    ("kinnex_sc",       "scRNA-Seq", False, True,  "{d}_full_kinnex_sc_segmented.bam",   "PacBio"),
]

COLS = ["Patient_Name", "Dataset", "Run_Name", "File_Prefix",
        "Sequencing_Method", "Normal", "Group", "Alleles",
        "Platform", "Read_Type", "Sample_Type"]


def rows(dataset, cfg, include_long_read=True):
    out = []
    alleles = ",".join(cfg["hla"])
    entries = [t + ("Illumina",) for t in SHORT_READ]
    if include_long_read:
        entries += LONG_READ
    for label, method, is_normal, is_rna, tmpl, platform in entries:
        pre = f"{'n' if is_normal else 'a'}{'r' if is_rna else 'd'}"
        out.append({
            "Patient_Name": dataset,
            "Dataset": dataset,
            "Run_Name": f"{pre}-{label.replace('_', '-')}",
            "File_Prefix": tmpl.format(d=dataset),
            "Sequencing_Method": method,
            "Normal": "TRUE" if is_normal else "FALSE",
            "Group": "T0",
            "Alleles": alleles,
            "Platform": platform,
            "Read_Type": "short-read" if platform == "Illumina" else "long-read",
            "Sample_Type": "blood-normal" if is_normal else "tumor",
        })
    return out


def write(path, recs, with_alleles):
    cols = [c for c in COLS if keep_col(c, with_alleles)]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in recs:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    return path


def keep_col(c, with_alleles):
    return with_alleles or c != "Alleles"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", default="design.yaml")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--no-long-read", action="store_true")
    a = ap.parse_args()
    design = yaml.safe_load(open(a.design))
    for ds, cfg in design["datasets"].items():
        d = os.path.join(a.out_dir, ds)
        os.makedirs(d, exist_ok=True)
        recs = rows(ds, cfg, include_long_read=not a.no_long_read)
        p1 = write(os.path.join(d, f"lens.{ds}.with-alleles.manifest"), recs, True)
        p2 = write(os.path.join(d, f"lens.{ds}.no-alleles.manifest"), recs, False)
        print(f"  {ds}: {len(recs)} rows")
        print(f"    {os.path.basename(p1)}  (Alleles: {','.join(cfg['hla'])})")
        print(f"    {os.path.basename(p2)}  (LENS types HLA from the normal)")


if __name__ == "__main__":
    main()
