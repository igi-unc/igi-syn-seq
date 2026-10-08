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
# The two File_Prefix templates are the FULL and the CHR1TO6 staged names. They have to differ, and
# neither may be a prefix of the other: preflight_resolve_inputs walks every search directory
# recursively and matches on the basename alone, so inputs/fastqs/<ds>/full and .../chr1to6 are one
# namespace, not two. That is what made "<ds>-TCR_S1_L001" match four FASTQs when 10x TCR was linked
# into both.
SHORT_READ = [
    ("wes_tumor",   "WES",       False, False, "{d}_tumor_wes",   "{d}_chr1to6_tumor_wes"),
    ("wes_normal",  "WES",       True,  False, "{d}_normal_wes",  "{d}_chr1to6_normal_wes"),
    ("rna_tumor",   "RNA-Seq",   False, True,  "{d}_full_rna",    "{d}_chr1to6_rna"),
    ("wgs_tumor",   "WGS",       False, False, "{d}_tumor_wgs",   "{d}_chr1to6_tumor_wgs"),
    ("wgs_normal",  "WGS",       True,  False, "{d}_normal_wgs",  "{d}_chr1to6_normal_wgs"),
    ("gex_tumor",   "scRNA-Seq", False, True,  "{d}-GEX_S1_L001", "{d}-chr1to6-GEX_S1_L001"),
    # 10x TCR is the one row whose two releases name the SAME file: design section 12 keeps the TCR
    # library complete in chr1to6 because TRA/TRB/TRG lie outside chr1-6. It is linked only under
    # full/, so the prefix resolves to one pair from either manifest.
    ("tcr_tumor",   "scTCR-Seq", False, True,  "{d}-TCR_S1_L001", "{d}-TCR_S1_L001"),
]

# Long-read arms. ONT is delivered as FASTQ and goes in the same way. PacBio HiFi and Kinnex
# are delivered as unaligned BAM, which the manifest specification has no column for, so they
# are held out of the manifest rather than guessed at; see the README this writes.
LONG_READ = [
    # (label, Sequencing_Method, normal?, RNA?, full prefix, chr1to6 prefix, Platform)
    ("ont_wgs_tumor",   "WGS",       False, False, "{d}_tumor_ont_wgs",
     "{d}_chr1to6_tumor_ont_wgs",                                                  "ONT"),
    ("ont_wgs_normal",  "WGS",       True,  False, "{d}_normal_ont_wgs",
     "{d}_chr1to6_normal_ont_wgs",                                                 "ONT"),
    ("ont_bulk_rna",    "RNA-Seq",   False, True,  "{d}_full_ont_bulk_rna",
     "{d}_chr1to6_ont_bulk_rna",                                                   "ONT"),
    ("ont_sc_rna",      "scRNA-Seq", False, True,  "{d}_full_ont_sc_rna",
     "{d}_chr1to6_ont_sc_rna",                                                     "ONT"),
    # BAM-delivered. File_Prefix is the FULL .bam filename, matching the
    # IPISRC044_T1_sclrs_live.seg.bam row in the v2.0.0-dev manifest. preflight_resolve_inputs
    # prefix-matches and then keeps only names ending in .bam, so the .pbi beside each one is
    # excluded, and it fails loudly if a prefix matches more than one BAM.
    ("hifi_wgs_tumor",  "WGS",       False, False, "{d}_tumor_hifi.bam",
     "{d}_chr1to6_tumor_hifi.bam",                                                 "PacBio"),
    ("hifi_wgs_normal", "WGS",       True,  False, "{d}_normal_hifi.bam",
     "{d}_chr1to6_normal_hifi.bam",                                                "PacBio"),
    ("kinnex_bulk",     "RNA-Seq",   False, True,  "{d}_full_kinnex_bulk_segmented.bam",
     "{d}_chr1to6_kinnex_bulk_segmented.bam",                                      "PacBio"),
    ("kinnex_sc",       "scRNA-Seq", False, True,  "{d}_full_kinnex_sc_segmented.bam",
     "{d}_chr1to6_kinnex_sc_segmented.bam",                                        "PacBio"),
]
RELEASES = ("full", "chr1to6")

COLS = ["Patient_Name", "Dataset", "Run_Name", "File_Prefix",
        "Sequencing_Method", "Normal", "Group", "Alleles",
        "Platform", "Read_Type", "Sample_Type"]


def rows(dataset, cfg, include_long_read=True, release="full"):
    """One row per sample type. Patient_Name stays the dataset; Dataset carries the release,
    per design section 2: Dataset = IGI-SYN-SEQ-01 or IGI-SYN-SEQ-01-chr1to6, Patient_Name =
    IGI-SYN-SEQ-01, one synthetic patient per dataset."""
    if release not in RELEASES:
        raise ValueError(f"release must be one of {RELEASES}, not {release!r}")
    out = []
    alleles = ",".join(cfg["hla"])
    entries = [t + ("Illumina",) for t in SHORT_READ]
    if include_long_read:
        entries += LONG_READ
    for label, method, is_normal, is_rna, tmpl_full, tmpl_sub, platform in entries:
        tmpl = tmpl_full if release == "full" else tmpl_sub
        pre = f"{'n' if is_normal else 'a'}{'r' if is_rna else 'd'}"
        out.append({
            "Patient_Name": dataset,
            "Dataset": dataset if release == "full" else f"{dataset}-{release}",
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
    ap.add_argument("--release", default="both", choices=("full", "chr1to6", "both"))
    a = ap.parse_args()
    design = yaml.safe_load(open(a.design))
    releases = RELEASES if a.release == "both" else (a.release,)
    for ds, cfg in design["datasets"].items():
        d = os.path.join(a.out_dir, ds)
        os.makedirs(d, exist_ok=True)
        for release in releases:
            recs = rows(ds, cfg, include_long_read=not a.no_long_read, release=release)
            stem = ds if release == "full" else f"{ds}.{release}"
            p1 = write(os.path.join(d, f"lens.{stem}.with-alleles.manifest"), recs, True)
            p2 = write(os.path.join(d, f"lens.{stem}.no-alleles.manifest"), recs, False)
            print(f"  {ds} [{release}]: {len(recs)} rows")
            print(f"    {os.path.basename(p1)}  (Alleles: {','.join(cfg['hla'])})")
            print(f"    {os.path.basename(p2)}  (LENS types HLA from the normal)")
        # No prefix may be a prefix of another, within or across releases: one namespace, see above.
        allp = [r["File_Prefix"] for rel in RELEASES
                for r in rows(ds, cfg, include_long_read=not a.no_long_read, release=rel)]
        allp = sorted(set(allp))
        bad = [(x, y) for x in allp for y in allp if x != y and y.startswith(x)]
        if bad:
            raise SystemExit(f"{ds}: File_Prefix values shadow each other: {bad}")


if __name__ == "__main__":
    main()
