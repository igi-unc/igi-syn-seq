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


def versioned(tmpl, version):
    """Insert a build version into a File_Prefix template, keeping that family's separator.

    Every template starts with `{d}` and the character after it is the separator its family uses:
    `_` for the general names, `-` for the 10x ones, whose `<sample>_S1_L001_R1_001` shape Cell Ranger
    parses and must not be disturbed.

    The version goes immediately after the dataset and NOT at the end, because versioned and unversioned
    builds share one namespace -- `preflight_resolve_inputs` walks recursively and matches on basename,
    so a version directory separates nothing. Appending would leave the OLD prefix `<d>_tumor_wes` still
    matching the new `<d>_tumor_wes_<ver>_R1`, and one row would claim four FASTQs, which is the same
    shadowing the TCR note above describes. Inserting at the dataset breaks the relation at the version.
    """
    if not version:
        return tmpl
    if not tmpl.startswith("{d}") or len(tmpl) < 4:
        raise ValueError(f"template does not start with a dataset token: {tmpl!r}")
    sep = tmpl[3]
    if sep not in "_-.":
        raise ValueError(f"unexpected separator {sep!r} after the dataset in {tmpl!r}")
    return "{d}" + sep + version + tmpl[3:]


# The off-the-shelf demo manifests in inputs/metadata/ots use a shorter, differently ordered column
# set than the per-dataset ones, and they carry several patients in one file. Matching them matters: a
# manifest dropped into that directory is read alongside its neighbours.
OTS_COLS = ["Patient_Name", "Run_Name", "Dataset", "File_Prefix",
            "Sequencing_Method", "Normal", "Group"]


def rows(dataset, cfg, include_long_read=True, release="full", version="", only=None, group="T0"):
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
    if only:
        missing = set(only) - {e[0] for e in entries}
        if missing:
            raise SystemExit(f"--only names no such sample: {sorted(missing)}. "
                             f"Available: {sorted(e[0] for e in entries)}")
        entries = [e for e in entries if e[0] in only]
    for label, method, is_normal, is_rna, tmpl_full, tmpl_sub, platform in entries:
        tmpl = tmpl_full if release == "full" else tmpl_sub
        pre = f"{'n' if is_normal else 'a'}{'r' if is_rna else 'd'}"
        out.append({
            "Patient_Name": dataset,
            "Dataset": dataset if release == "full" else f"{dataset}-{release}",
            "Run_Name": f"{pre}-{label.replace('_', '-')}",
            "File_Prefix": versioned(tmpl, version).format(d=dataset),
            "Sequencing_Method": method,
            "Normal": "TRUE" if is_normal else "FALSE",
            "Group": group,
            "Alleles": alleles,
            "Platform": platform,
            "Read_Type": "short-read" if platform == "Illumina" else "long-read",
            "Sample_Type": "blood-normal" if is_normal else "tumor",
        })
    return out


def write(path, recs, with_alleles, cols=None):
    cols = cols or [c for c in COLS if keep_col(c, with_alleles)]
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
    ap.add_argument("--out-dir", default="")
    ap.add_argument("--no-long-read", action="store_true")
    ap.add_argument("--release", default="both", choices=("full", "chr1to6", "both"))
    ap.add_argument("--combined", default="",
                    help="write ONE manifest at this path containing every dataset, instead of one "
                         "file per dataset. Each dataset gets its own Group, because Group is the "
                         "pairing key and a shared value would let one patient's tumour pair with "
                         "another's normal.")
    ap.add_argument("--ots-columns", action="store_true",
                    help="use the shorter column set and order of inputs/metadata/ots")
    ap.add_argument("--only", default="",
                    help="comma-separated sample labels to include, e.g. "
                         "wes_tumor,wes_normal,rna_tumor. Default: every sample.")
    ap.add_argument("--subset-name", default="",
                    help="infix for the manifest filename when --only is used, e.g. wes-rna")
    ap.add_argument("--alleles", default="both", choices=("both", "with", "no"),
                    help="which allele variants to write")
    ap.add_argument("--version", default="",
                    help="build version to embed in File_Prefix and the manifest filename, "
                         "e.g. 20261026. Omit for the unversioned names.")
    a = ap.parse_args()
    if not a.out_dir and not a.combined:
        raise SystemExit("one of --out-dir or --combined is required")
    design = yaml.safe_load(open(a.design))
    releases = RELEASES if a.release == "both" else (a.release,)
    only = [x.strip() for x in a.only.split(",") if x.strip()] or None
    if only and not a.subset_name:
        raise SystemExit("--only needs --subset-name, or the file would overwrite the full manifest")
    if a.combined:
        cols = OTS_COLS if a.ots_columns else [c for c in COLS if keep_col(c, a.alleles == "with")]
        allrecs = []
        for i, (ds, cfg) in enumerate(sorted(design["datasets"].items()), start=1):
            for release in releases:
                allrecs += rows(ds, cfg, include_long_read=not a.no_long_read, release=release,
                                version=a.version, only=only, group=str(i))
        os.makedirs(os.path.dirname(os.path.abspath(a.combined)), exist_ok=True)
        write(a.combined, allrecs, a.alleles == "with", cols=cols)
        print(f"  combined: {len(allrecs)} rows -> {a.combined}")
        for g, ds in enumerate(sorted(design["datasets"]), start=1):
            print(f"    Group {g} = {ds}")
        return

    for ds, cfg in design["datasets"].items():
        d = os.path.join(a.out_dir, ds)
        os.makedirs(d, exist_ok=True)
        for release in releases:
            recs = rows(ds, cfg, include_long_read=not a.no_long_read, release=release,
                        version=a.version, only=only)
            stem = ds if release == "full" else f"{ds}.{release}"
            if a.version:
                stem = f"{stem}.{a.version}"
            if a.subset_name:
                stem = f"{stem}.{a.subset_name}"
            print(f"  {ds} [{release}]: {len(recs)} rows"
                  + (f"  ({', '.join(sorted(only))})" if only else ""))
            if a.alleles in ("both", "with"):
                p1 = write(os.path.join(d, f"lens.{stem}.with-alleles.manifest"), recs, True)
                print(f"    {os.path.basename(p1)}  (Alleles: {','.join(cfg['hla'])})")
            if a.alleles in ("both", "no"):
                p2 = write(os.path.join(d, f"lens.{stem}.no-alleles.manifest"), recs, False)
                print(f"    {os.path.basename(p2)}  (LENS types HLA from the normal)")
        # No prefix may be a prefix of another, within or across releases: one namespace, see above.
        # Checked against the UNVERSIONED prefixes too: both builds live in one basename namespace,
        # so a new prefix that shadows an old file is just as broken as one that shadows a new file.
        allp = [r["File_Prefix"] for rel in RELEASES for v in {a.version, ""}
                for r in rows(ds, cfg, include_long_read=not a.no_long_read, release=rel, version=v,
                              only=only)]
        allp = sorted(set(allp))
        bad = [(x, y) for x in allp for y in allp if x != y and y.startswith(x)]
        if bad:
            raise SystemExit(f"{ds}: File_Prefix values shadow each other: {bad}")


if __name__ == "__main__":
    main()
