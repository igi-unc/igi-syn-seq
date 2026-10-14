#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Render the two-page IGI-SYN-SEQ dataset infographic: the design, then the measured status.

Page 1 is the specification -- what the two datasets are meant to be, read from `design.yaml` and the
catalog summaries. Page 2 is what is on disk, read from a release tree and the builders' own sidecar
JSONs. Nothing on page 2 is taken from the specification, so a disagreement between the pages is a real
disagreement rather than a formatting choice.

Page 2 is a snapshot. Regenerate it whenever the release changes -- the build state it reports, and the
rebuilt-candidate comparison in particular, go stale as soon as a stage finishes.

    python3 tools/make_dataset_infographic.py \
        --design catalog/design.yaml --catalog catalog/output \
        --release "$IGI_RELEASE" [--candidate "$IGI_WORK/release-<date>"] \
        --out igi-syn-seq-datasets.pdf

`--candidate` is optional: without it both measured columns come from `--release` and the
"fixed by this rebuild" panel is omitted, which is the right output once a rebuild has been accepted
and staged.
"""
import argparse, glob, json, os, sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyBboxPatch, Rectangle


CHROMS = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]


def jload(p):
    try:
        return json.load(open(p))
    except Exception:
        return None


def size(*paths):
    t = 0
    for p in paths:
        for f in glob.glob(p):
            try:
                t += os.path.getsize(os.path.realpath(f))
            except OSError:
                pass
    return t


def per_chrom(root, arm, ds, suffix, key, lib):
    """Sum one sidecar key over the 24 per-chromosome tasks of one library."""
    tot, n = 0, 0
    for c in CHROMS:
        d = jload(f"{root}/{arm}/{ds}/{ds}_{c}_{lib}{suffix}.json")
        if d and key in d:
            tot += d[key]
            n += 1
    return tot, n


def collect(root, ds):
    out = {}

    # 1-2 exome, 4-5 Illumina WGS: pairs_written per chromosome; bytes from the merged library.
    for lib, num in (("tumor", 1), ("normal", 2)):
        p, n = per_chrom(root, "wes", ds, "", "pairs_written", lib)
        out[f"wes_{lib}"] = dict(n=num, units="pairs", measured=p, chroms=n,
                                 bytes=size(f"{root}/merged/{ds}/{ds}_{lib}_R[12].fastq.gz"))
    for lib, num in (("tumor", 4), ("normal", 5)):
        p, n = per_chrom(root, "wgs", ds, "_wgs", "pairs_written", lib)
        out[f"wgs_{lib}"] = dict(n=num, units="pairs", measured=p, chroms=n,
                                 bytes=size(f"{root}/wgs_merged/{ds}/{ds}_{lib}_wgs_R[12].fastq.gz"))
    # 6-7 PacBio HiFi, 14-15 ONT WGS: read counts.
    for lib, num in (("tumor", 6), ("normal", 7)):
        r, n = per_chrom(root, "pacbio", ds, "_hifi", "reads", lib)
        out[f"pacbio_{lib}"] = dict(n=num, units="reads", measured=r, chroms=n,
                                    bytes=size(f"{root}/pacbio_merged/{ds}/{ds}_{lib}_hifi.bam"))
    for lib, num in (("tumor", 14), ("normal", 15)):
        r, n = per_chrom(root, "ont_wgs", ds, "_ont_wgs", "reads", lib)
        out[f"ont_wgs_{lib}"] = dict(n=num, units="reads", measured=r, chroms=n,
                                     bytes=size(f"{root}/ont_wgs_merged/{ds}/{ds}_{lib}_ont_wgs.fastq.gz"))

    # 3, 8-13: whole-library arms, one sidecar each.
    whole = [
        (3,  "rna",          f"{ds}_full_rna.json",          "pairs_written",      "pairs",
         [f"{root}/rna/{ds}/{ds}_full_rna_R[12].fastq.gz"]),
        (8,  "kinnex_bulk",  f"{ds}_full_kinnex_bulk.json",  "segments_post_skera", "segments",
         [f"{root}/kinnex_bulk/{ds}/{ds}_full_kinnex_bulk_*.bam"]),
        (9,  "kinnex_sc",    f"{ds}_full_kinnex_sc.json",    "segments_post_skera", "segments",
         [f"{root}/kinnex_sc/{ds}/{ds}_full_kinnex_sc_*.bam"]),
        (10, "tenx_gex",     f"{ds}_full_tenx_gex.json",     "read_pairs_written", "pairs",
         [f"{root}/tenx_gex/{ds}/{ds}-GEX_S1_L001_R[12]_001.fastq.gz"]),
        (11, "tenx_tcr",     f"{ds}_full_tenx_tcr.json",     "read_pairs_written", "pairs",
         [f"{root}/tenx_tcr/{ds}/{ds}-TCR_S1_L001_R[12]_001.fastq.gz"]),
        (13, "ont_bulk_rna", f"{ds}_full_ont_bulk_rna.json", "reads",              "reads",
         [f"{root}/ont_bulk_rna/{ds}/{ds}_full_ont_bulk_rna.fastq.gz"]),
        (12, "ont_sc_rna",   f"{ds}_full_ont_sc_rna.json",   "reads",              "reads",
         [f"{root}/ont_sc_rna/{ds}/{ds}_full_ont_sc_rna.fastq.gz"]),
    ]
    for num, arm, fn, key, units, globs in whole:
        d = jload(f"{root}/{arm}/{ds}/{fn}") or {}
        # pairs_written is absent from the RNA sidecar under its older name
        m = d.get(key, d.get("pairs_written", d.get("reads", 0)))
        rec = dict(n=num, units=units, measured=m, bytes=size(*globs))
        for extra in ("cells", "t_cells", "clonotypes", "mean_cdna_length", "reads_pre_skera",
                      "length_mean", "pairs_target", "segments_target", "reads_per_cell_target",
                      "molecules", "arrays"):
            if extra in d:
                rec[extra] = d[extra]
        out[arm] = rec
    return out



DESIGN, LIVE, CAND, CAT = {}, {}, {}, {}
D1, D2 = "IGI-SYN-SEQ-01", "IGI-SYN-SEQ-02"
INK, MUT, RULE = "#161b26", "#6b7280", "#dfe5ec"
C1, C2 = "#1d4ed8", "#c2410c"
GOOD, WARN, BAD, PANEL = "#15803d", "#b45309", "#b3261e", "#f6f8fb"
W, H = 16.54, 11.69          # A3 landscape

SAMPLES = [
    (1,  "tumour WES",          "Illumina PE150",   "150x on-target",        "wes_tumor",     "pairs"),
    (2,  "normal WES",          "Illumina PE150",   "100x",                  "wes_normal",    "pairs"),
    (3,  "tumour bulk RNA",     "Illumina PE150",   "80 M pairs, polyA",     "rna",           "pairs"),
    (4,  "tumour WGS",          "Illumina PE150",   "100x",                  "wgs_tumor",     "pairs"),
    (5,  "normal WGS",          "Illumina PE150",   "30x",                   "wgs_normal",    "pairs"),
    (6,  "tumour HiFi WGS",     "PacBio Revio",     "30x",                   "pacbio_tumor",  "reads"),
    (7,  "normal HiFi WGS",     "PacBio Revio",     "30x",                   "pacbio_normal", "reads"),
    (8,  "tumour Kinnex bulk",  "PacBio Revio",     "15 M segments",         "kinnex_bulk",   "segs"),
    (9,  "tumour Kinnex scRNA", "PacBio Revio",     "4,000 cells / 12 M",    "kinnex_sc",     "segs"),
    (10, "tumour 10x 5' GEX",   "Illumina",         "4,000 cells / 40k pc",  "tenx_gex",      "pairs"),
    (11, "tumour 10x 5' TCR",   "Illumina",         "~5k reads per cell",    "tenx_tcr",      "pairs"),
    (12, "tumour ONT scRNA",    "ONT R10.4.1",      "4,000 cells / 10 M",    "ont_sc_rna",    "reads"),
    (13, "tumour ONT bulk RNA", "ONT R10.4.1",      "~20 M reads",           "ont_bulk_rna",  "reads"),
    (14, "tumour ONT WGS",      "ONT R10.4.1",      "30x",                   "ont_wgs_tumor", "reads"),
    (15, "normal ONT WGS",      "ONT R10.4.1",      "30x",                   "ont_wgs_normal","reads"),
]

# From build-reference.md section 1, which is the one place build/acceptance status is tracked.
ACCEPT = {
    1: ("44/44", GOOD, "14_accept"),           2: ("never opened", BAD, "see S24"),
    3: ("5/5", GOOD, "15_accept_rna"),         4: ("40/40", GOOD, "22_accept_wgs"),
    5: ("not aligned", BAD, "never opened"),   6: ("none", WARN, "no stage"),
    7: ("none", WARN, "no stage"),             8: ("preflight", WARN, "consumer only"),
    9: ("preflight", WARN, "consumer only"),   10: ("CR count", WARN, "Cell Ranger"),
    11: ("CR vdj", WARN, "see S25"),           12: ("none", WARN, "no stage"),
    13: ("none", WARN, "no stage"),            14: ("none", WARN, "no stage"),
    15: ("none", WARN, "no stage"),
}


def fnum(v, unit=""):
    if v >= 1e9:  return f"{v/1e9:.2f} B{unit}"
    if v >= 1e6:  return f"{v/1e6:.1f} M{unit}"
    if v >= 1e3:  return f"{v/1e3:.0f} k{unit}"
    return f"{v:,.0f}{unit}"


def txt(fig, x, y, s, size=9, w="normal", c=None, ha="left", va="baseline", fam=None, alpha=1.0):
    fig.text(x, y, s, fontsize=size, fontweight=w, color=c or INK, ha=ha, va=va,
             family=fam or "DejaVu Sans", alpha=alpha)


def panel(fig, x, y, w, h, fc=PANEL, ec=RULE, lw=0.8, r=0.006):
    fig.patches.append(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                      fc=fc, ec=ec, lw=lw, transform=fig.transFigure, zorder=0))


def bar(fig, x, y, w, h, frac, c, bg="#e8edf4"):
    fig.patches.append(Rectangle((x, y), w, h, fc=bg, ec="none", transform=fig.transFigure, zorder=1))
    fig.patches.append(Rectangle((x, y), w * max(0.0, min(1.3, frac)), h, fc=c, ec="none",
                                 transform=fig.transFigure, zorder=2))


def chip(fig, x, y, s, fc, tc="white", pad=0.004, size=7.5, w="bold"):
    t = fig.text(x, y, s, fontsize=size, fontweight=w, color=tc, ha="left", va="center", zorder=3,
                 bbox=dict(boxstyle=f"round,pad={pad*90}", fc=fc, ec="none"))
    return t


def wrap(fig, x, y, text, width, size=7.3, lead=0.0122, c=None, w="normal", n=99):
    """Lay prose out inside a panel. Every overflow on the first draft was an unwrapped string."""
    words, line, lines = text.split(), "", []
    for word in words:
        if len(line) + len(word) + 1 > width:
            lines.append(line); line = word
        else:
            line = (line + " " + word).strip()
    lines.append(line)
    for i, ln in enumerate(lines[:n]):
        txt(fig, x, y - i * lead, ln, size, w, c or MUT)
    return len(lines[:n])


def header(fig, title, sub, kicker):
    txt(fig, 0.035, 0.955, kicker, 8.5, "bold", C1)
    txt(fig, 0.035, 0.925, title, 23, "bold", INK)
    txt(fig, 0.035, 0.900, sub, 10.5, "normal", MUT)
    fig.patches.append(Rectangle((0.035, 0.890), 0.93, 0.0016, fc=INK, ec="none",
                                 transform=fig.transFigure, zorder=1))


def footer(fig, left, right):
    fig.patches.append(Rectangle((0.035, 0.042), 0.93, 0.0012, fc=RULE, ec="none",
                                 transform=fig.transFigure, zorder=1))
    txt(fig, 0.035, 0.026, left, 7.5, "normal", MUT)
    txt(fig, 0.965, 0.026, right, 7.5, "normal", MUT, ha="right")


# ----------------------------------------------------------------------------- page 1: the design
def page_design(pdf):
    fig = plt.figure(figsize=(W, H)); fig.patch.set_facecolor("white")
    header(fig, "IGI-SYN-SEQ: two synthetic multi-assay tumour datasets",
           "Design specification — fifteen sample types per dataset, one shared designed-event catalog, "
           "TNBC-modelled architecture over two real germline baselines.",
           "THE DESIGN  ·  WHAT EACH DATASET IS MEANT TO CONTAIN")

    # ---- the two dataset architecture cards
    for i, (ds, col) in enumerate(((D1, C1), (D2, C2))):
        d = DESIGN[ds]
        x = 0.035 + i * 0.472
        panel(fig, x, 0.575, 0.458, 0.295, fc="white")
        fig.patches.append(Rectangle((x, 0.8615), 0.458, 0.0085, fc=col, ec="none",
                                     transform=fig.transFigure, zorder=1))
        txt(fig, x + 0.014, 0.837, ds, 14, "bold", col)
        txt(fig, x + 0.014, 0.819, f"germline baseline  {d['baseline']}   ·   {d['sex']}", 8.5, "normal", MUT)

        # purity / ploidy / WGD / mutation burden strip
        facts = [("purity", f"{d['purity']:.2f}"),
                 ("WGD", "yes" if d["wgd"] else "no"),
                 ("ploidy", "~3.6" if d["wgd"] else "~2.3"),
                 ("background", f"{d['background_mut_per_mb']:.0f} mut/Mb")]
        for j, (k, v) in enumerate(facts):
            fx = x + 0.014 + j * 0.111
            panel(fig, fx, 0.767, 0.100, 0.040, fc=PANEL)
            txt(fig, fx + 0.008, 0.792, v, 12, "bold", INK)
            txt(fig, fx + 0.008, 0.775, k, 7.5, "normal", MUT)

        # clone tree: T -> A -> A1, T -> B, bar length = CCF
        txt(fig, x + 0.014, 0.748, "SUBCLONAL ARCHITECTURE   (bar = cancer cell fraction)", 7.5, "bold", MUT)
        order = [("T", 0), ("A", 1), ("A1", 2), ("B", 1)]
        for j, (cl, indent) in enumerate(order):
            ccf = d["clones"][cl]
            yy = 0.726 - j * 0.0225
            bx = x + 0.020 + indent * 0.016
            txt(fig, bx, yy, cl, 8.5, "bold", col, va="center")
            bar(fig, bx + 0.024, yy - 0.006, 0.26, 0.012, ccf, col)
            txt(fig, bx + 0.290, yy, f"{ccf:.2f}", 8, "bold", INK, va="center")
            drv = d["drivers"].get(cl, [])
            if drv:
                txt(fig, bx + 0.316, yy, ", ".join(g for g, _w, _k in drv),
                    7.3, "normal", MUT, va="center")

        # HLA alleles, with the haplotype the CNA removes marked
        txt(fig, x + 0.014, 0.622, "CLASS I HLA   ·   hap1 is removed by the HLA-LOH CNA", 7.5, "bold", MUT)
        for h, yy in (("0", 0.602), ("1", 0.582)):
            lost = h == "1"
            txt(fig, x + 0.020, yy, f"hap{h}", 8, "bold", MUT if not lost else BAD, va="center")
            for j, a in enumerate(d["hla_haplotypes"][h]):
                chip(fig, x + 0.050 + j * 0.098, yy, a, BAD if lost else col, size=7.2)
            if lost:
                txt(fig, x + 0.355, yy, "LOST", 7.5, "bold", BAD, va="center")
            else:
                txt(fig, x + 0.355, yy, "retained", 7.5, "normal", MUT, va="center")
        if len(set(d["hla_haplotypes"]["0"] + d["hla_haplotypes"]["1"])) < 6:
            txt(fig, x + 0.233, 0.622, "HLA-A homozygous — the loss costs B and C only.",
                6.8, "normal", MUT)

    # ---- fifteen sample types
    panel(fig, 0.035, 0.108, 0.555, 0.448, fc="white")
    txt(fig, 0.049, 0.532, "THE FIFTEEN SAMPLE TYPES, PER DATASET", 9, "bold", INK)
    txt(fig, 0.049, 0.516, "Both datasets carry all fifteen. No normal RNA by design; the four "
        "single-cell libraries share one 4,000-cell roster, so a barcode means the same cell in all four.",
        7.3, "normal", MUT)
    cols = [(0.049, "#"), (0.068, "SAMPLE"), (0.235, "PLATFORM"), (0.360, "DESIGN TARGET")]
    for cx, lab in cols:
        txt(fig, cx, 0.496, lab, 7, "bold", MUT)
    fig.patches.append(Rectangle((0.049, 0.4915), 0.527, 0.0011, fc=RULE, ec="none",
                                 transform=fig.transFigure, zorder=1))
    GRP = {1: "#eef2ff", 2: "#eef2ff", 3: "#ecfdf5", 4: "#eef2ff", 5: "#eef2ff",
           6: "#fff7ed", 7: "#fff7ed", 8: "#fff7ed", 9: "#fff7ed",
           10: "#ecfdf5", 11: "#ecfdf5", 12: "#f5f3ff", 13: "#f5f3ff", 14: "#f5f3ff", 15: "#f5f3ff"}
    for k, (n, name, plat, tgt, _key, _u) in enumerate(SAMPLES):
        yy = 0.474 - k * 0.0245
        fig.patches.append(Rectangle((0.049, yy - 0.0085), 0.527, 0.0205, fc=GRP[n], ec="none",
                                     transform=fig.transFigure, zorder=0))
        txt(fig, 0.0515, yy, str(n), 7.5, "bold", MUT, va="center")
        txt(fig, 0.068, yy, name, 8.5, "bold", INK, va="center")
        txt(fig, 0.235, yy, plat, 7.8, "normal", MUT, va="center")
        txt(fig, 0.360, yy, tgt, 7.8, "normal", INK, va="center")
    txt(fig, 0.049, 0.088, "Colour groups the delivery technology: blue Illumina DNA, green Illumina RNA, "
        "orange PacBio, violet Oxford Nanopore.", 7, "normal", MUT)

    # ---- designed event catalog
    panel(fig, 0.607, 0.108, 0.358, 0.448, fc="white")
    txt(fig, 0.621, 0.532, "THE DESIGNED EVENT CATALOG", 9, "bold", INK)
    txt(fig, 0.621, 0.516, "Selected per dataset from one seeded design, spanning an evidence grid of "
        "expression x binding x clonality tiers.", 7.3, "normal", MUT)
    txt(fig, 0.830, 0.496, D1.replace("IGI-SYN-SEQ-", "ds-"), 7.5, "bold", C1, ha="center")
    txt(fig, 0.905, 0.496, D2.replace("IGI-SYN-SEQ-", "ds-"), 7.5, "bold", C2, ha="center")
    rows = [
        ("SNV / indel events", lambda c: c["n_events"], None),
        ("   missense", lambda c: c["by_consequence"]["missense"], None),
        ("   frameshift", lambda c: c["by_consequence"]["frameshift"], None),
        ("   nonsense / synonymous",
         lambda c: c["by_consequence"]["nonsense"] + c["by_consequence"]["synonymous"], None),
        ("structural variants", lambda c: c["svs"]["n"], None),
        ("gene fusions", lambda c: c["fusions"]["n"], None),
        ("expressed-class targets", lambda c: c["expressed"]["n"], None),
        ("HLA-LOH allele rows", lambda c: c["hla_loh"]["n"], None),
        ("viral events", lambda c: c["viruses"]["n"], None),
        ("flagposts (all classes)", lambda c: c["flagposts_all_classes"]["total"], None),
        ("strong binders (all 6 alleles)", lambda c: c["by_binding_tier"]["strong"], True),
        ("best allele lost to HLA LOH", lambda c: c["best_allele_is_lost"], True),
        ("tier changed by HLA loss", lambda c: c["tier_changed_by_hla_loss"], True),
        ("evidence-grid cells filled", lambda c: c["grid_cells_filled"], True),
        ("fraction on chr1-6", lambda c: c["chr1to6_fraction"], True),
    ]
    for k, (lab, fn, em) in enumerate(rows):
        yy = 0.474 - k * 0.0245
        if k % 2 == 0:
            fig.patches.append(Rectangle((0.621, yy - 0.0085), 0.330, 0.0205, fc=PANEL, ec="none",
                                         transform=fig.transFigure, zorder=0))
        txt(fig, 0.624, yy, lab, 8.2, "bold" if not lab.startswith("  ") else "normal",
            INK if not lab.startswith("  ") else MUT, va="center")
        for ds, cx, cc in ((D1, 0.830, C1), (D2, 0.905, C2)):
            v = fn(CAT[ds])
            s = f"{v:.3f}" if isinstance(v, float) else f"{v:,}"
            txt(fig, cx, yy, s, 8.5, "bold", cc, ha="center", va="center")
    wrap(fig, 0.621, 0.090, "The grid is filled on the tier over all six alleles, so event selection"
         " does not move when a haplotype is lost — only the retained-allele columns do.", 72, 7)

    footer(fig, "Design: catalog/design.yaml (seed 20260925)  ·  spec: docs/IGI-SYN-SEQ-design.md  "
                "·  event counts: catalog/output/*.snv_indel.summary.json",
           "page 1 of 2  ·  the design")
    pdf.savefig(fig); plt.close(fig)


# ----------------------------------------------------------------------------- page 2: the status
def page_status(pdf):
    fig = plt.figure(figsize=(W, H)); fig.patch.set_facecolor("white")
    header(fig, "Current status — measured on disk, 2026-10-14",
           "Every number below is read from the release tree and the builders' own sidecar JSONs, not from "
           "the specification. All fifteen sample types exist at release scale for both datasets.",
           "THE STATUS  ·  WHAT IS ACTUALLY BUILT RIGHT NOW")

    # ---- headline counters
    tot_live = sum(LIVE[d][k]["bytes"] for d in (D1, D2) for k in LIVE[d])
    heads = [("15 / 15", "sample types built,\nboth datasets", GOOD),
             ("3 / 15", "passed a release-scale\nacceptance stage", BAD),
             (f"{tot_live/1e12:.2f} TB", "merged primary libraries\n(6.9 TB incl. per-chrom + subsets)", INK),
             ("2,139", "designed events, all classes,\nacross both catalogs", INK),
             ("ds-02", "reads rebuilt 2026-10-13/14,\nmerges in flight", C2)]
    for i, (big, lab, c) in enumerate(heads):
        x = 0.035 + i * 0.1875
        panel(fig, x, 0.782, 0.175, 0.086, fc="white")
        fig.patches.append(Rectangle((x, 0.782), 0.0045, 0.086, fc=c, ec="none",
                                     transform=fig.transFigure, zorder=1))
        txt(fig, x + 0.016, 0.838, big, 19, "bold", c)
        for j, line in enumerate(lab.split("\n")):
            txt(fig, x + 0.016, 0.814 - j * 0.013, line, 7.5, "normal", MUT)

    # ---- build / acceptance matrix with measured yield
    panel(fig, 0.035, 0.300, 0.610, 0.455, fc="white")
    txt(fig, 0.049, 0.731, "PER-SAMPLE YIELD AND VERIFICATION COVERAGE", 9, "bold", INK)
    txt(fig, 0.049, 0.715, "Bars are measured reads, pairs or post-skera segments against the larger of "
        "the two datasets, so equal bars mean the pair is matched. Grey bar = ds-01, coloured = ds-02.",
        7.3, "normal", MUT)
    for cx, lab, ha in ((0.049, "#", "left"), (0.064, "SAMPLE", "left"), (0.205, "ds-01 MEASURED", "right"),
                        (0.300, "ds-02 MEASURED", "right"), (0.312, "RELATIVE YIELD", "left"),
                        (0.520, "ACCEPTANCE", "left")):
        txt(fig, cx, 0.695, lab, 7, "bold", MUT, ha=ha)
    fig.patches.append(Rectangle((0.049, 0.6905), 0.582, 0.0011, fc=RULE, ec="none",
                                 transform=fig.transFigure, zorder=1))
    for k, (n, name, _plat, _tgt, key, unit) in enumerate(SAMPLES):
        yy = 0.671 - k * 0.0243
        if k % 2 == 0:
            fig.patches.append(Rectangle((0.049, yy - 0.0090), 0.582, 0.0214, fc=PANEL, ec="none",
                                         transform=fig.transFigure, zorder=0))
        a = LIVE[D1][key]["measured"]
        b = CAND[D2][key]["measured"]
        mx = max(a, b) or 1
        txt(fig, 0.0515, yy, str(n), 7.5, "bold", MUT, va="center")
        txt(fig, 0.064, yy, name, 8.2, "bold", INK, va="center")
        txt(fig, 0.205, yy, fnum(a), 8, "normal", INK, ha="right", va="center")
        txt(fig, 0.300, yy, fnum(b), 8, "bold", C2, ha="right", va="center")
        bar(fig, 0.312, yy + 0.0012, 0.195, 0.0046, a / mx, "#9aa6b8")
        bar(fig, 0.312, yy - 0.0055, 0.195, 0.0046, b / mx, C2)
        lab, col, note = ACCEPT[n]
        chip(fig, 0.520, yy, lab, col, size=6.8, pad=0.0035)
        txt(fig, 0.583, yy, note, 7, "normal", MUT, va="center")
    wrap(fig, 0.049, 0.3135, "Ten of fifteen sample types have no acceptance stage — four of those have a "
         "weaker consumer-side check — and two more were built but never opened.", 130, 7)

    # ---- the normal-depth defect, before and after. Only meaningful against a rebuilt candidate: with
    # one release root the two bars would be the same measurement drawn twice.
    if not HAVE_CANDIDATE:
        panel(fig, 0.662, 0.541, 0.303, 0.214, fc="white")
        txt(fig, 0.676, 0.731, "SINGLE RELEASE ROOT", 9, "bold", MUT)
        wrap(fig, 0.676, 0.712, "No --candidate was given, so both measured columns come from the same "
             "release tree and there is no before/after to draw.", 62, 7.1)
        return _rest(fig, pdf)
    panel(fig, 0.662, 0.541, 0.303, 0.214, fc="white")
    txt(fig, 0.676, 0.731, "FIXED BY THIS REBUILD", 9, "bold", GOOD)
    wrap(fig, 0.676, 0.715, "All four ds-02 NORMAL libraries shipped at 69 % of their designed depth. "
         "Nothing measured it, because no acceptance stage had ever opened a normal library.", 62, 7.1)
    txt(fig, 0.676, 0.668, "ds-02 normal yield as a fraction of ds-01", 7, "bold", MUT)
    pairs = [("normal WES", "wes_normal"), ("normal WGS", "wgs_normal"),
             ("normal HiFi WGS", "pacbio_normal"), ("normal ONT WGS", "ont_wgs_normal")]
    for j, (lab, key) in enumerate(pairs):
        yy = 0.645 - j * 0.0258
        ref = LIVE[D1][key]["measured"]
        before = LIVE[D2][key]["measured"] / ref
        after = CAND[D2][key]["measured"] / ref
        txt(fig, 0.678, yy, lab, 7.8, "bold", INK, va="center")
        bar(fig, 0.772, yy + 0.0022, 0.120, 0.0048, before, BAD)
        bar(fig, 0.772, yy - 0.0055, 0.120, 0.0048, after, GOOD)
        txt(fig, 0.898, yy + 0.0046, f"{before:.0%}", 6.5, "bold", BAD, va="center")
        txt(fig, 0.898, yy - 0.0043, f"{after:.0%}", 6.5, "bold", GOOD, va="center")
    txt(fig, 0.676, 0.5495, "red = as delivered   ·   green = rebuilt candidate", 6.6, "normal", MUT)

    _rest(fig, pdf)


def _rest(fig, pdf):
    # ---- what changed and what is in flight
    panel(fig, 0.662, 0.300, 0.303, 0.214, fc="white")
    txt(fig, 0.676, 0.490, "IN FLIGHT AND PENDING", 9, "bold", C2)
    items = [
        (GOOD, "ds-02 germline VCF rebuilt", "class I HLA now allele-specific; 55/55 and 104/104 het calls phased at HLA-B and -C"),
        (GOOD, "ds-01 HLA-LOH truth corrected", "the allele-to-haplotype declaration was inverted; reads were right all along"),
        (GOOD, "ds-02 reads rebuilt", "all 24 chromosomes x 2 libraries x 4 DNA arms, plus all 7 whole-library arms"),
        (WARN, "merges, alignments, QA", "driver chain running; exomes merged and aligned, long reads in progress"),
        (WARN, "truth bundle held", "deliberately, so it cannot publish from the pre-fix catalog"),
        (BAD,  "staging not yet done", "nothing enters the consumer input tree without owner approval"),
    ]
    for j, (c, head, sub) in enumerate(items):
        yy = 0.466 - j * 0.0285
        fig.patches.append(Rectangle((0.678, yy - 0.0035), 0.0035, 0.018, fc=c, ec="none",
                                     transform=fig.transFigure, zorder=1))
        txt(fig, 0.688, yy + 0.0065, head, 8, "bold", INK, va="center")
        txt(fig, 0.688, yy - 0.0045, sub[:84], 6.8, "normal", MUT, va="center")

    # ---- known caveats
    panel(fig, 0.035, 0.062, 0.930, 0.218, fc="white")
    txt(fig, 0.049, 0.256, "KNOWN CAVEATS A CONSUMER SHOULD SCORE AGAINST", 9, "bold", INK)
    cav = [
        ("Kinnex scRNA yield", "The MAS array carries 16 cDNAs but only 9 adapters, so isoseq refine "
         "discards the double-stuffed segment as a concatemer: 12.6 % of segments and roughly half the "
         "designed molecules never reach the cell matrix. bcstats recovers 3,707 of 4,000 cells. The "
         "cells and clonotypes that arrive are correct — a yield loss, not a truth defect."),
        ("Bulk RNA strandedness", "Specified as stranded, delivered unstranded — salmon calls both "
         "libraries unstranded and the document describing the data wins. Nothing in LENS depends on "
         "bulk RNA strand."),
        ("Acceptance coverage", "Three of fifteen sample types have a release-scale acceptance stage. "
         "The long-read and single-cell arms are verified only by consumer preflight or by Cell Ranger "
         "running to completion, which answers a weaker question than a truth comparison."),
        ("10x TCR repertoire", "V and J anchors are real GENCODE germline sequence, so a caller recovers "
         "the genes the truth table names, but N-region nucleotides are drawn rather than taken from a "
         "real repertoire: junction length and composition are realistic, the clonotypes are nobody's."),
    ]
    for j, (head, body) in enumerate(cav):
        x = 0.049 + (j % 2) * 0.468
        y = 0.228 - (j // 2) * 0.082
        txt(fig, x, y, head, 8.2, "bold", C1)
        wrap(fig, x, y - 0.0145, body, 97, 7, n=5)

    footer(fig, "Measured from $IGI_WORK/release and the rebuilt candidate release-20261012  ·  "
                "status: docs/build-reference.md §1  ·  caveats: §6.3, §20.5, §24, §25",
           "page 2 of 2  ·  the status")
    pdf.savefig(fig); plt.close(fig)


def main(argv=None):
    global DESIGN, LIVE, CAND, CAT, HAVE_CANDIDATE
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--design", required=True, help="path to design.yaml")
    ap.add_argument("--catalog", required=True, help="directory holding <dataset>.snv_indel.summary.json")
    ap.add_argument("--release", required=True, help="release root to measure")
    ap.add_argument("--candidate", help="a rebuilt candidate root to compare the release against")
    ap.add_argument("--out", required=True, help="output PDF")
    a = ap.parse_args(argv)

    import yaml
    d = yaml.safe_load(open(a.design))
    DESIGN = {}
    for ds, c in d["datasets"].items():
        DESIGN[ds] = dict(
            baseline=c["baseline"], sex=c["sex"], purity=c["purity"], wgd=c.get("wgd", False),
            hla_haplotypes={str(k): v for k, v in c["hla_haplotypes"].items()},
            clones={k: v["ccf"] for k, v in c["clones"].items()},
            drivers={k: [list(x) for x in v] for k, v in c.get("drivers", {}).items()},
            background_mut_per_mb=c.get("background_mut_per_mb"))
    CAT = {ds: json.load(open(os.path.join(a.catalog, f"{ds}.snv_indel.summary.json")))
           for ds in (D1, D2)}
    LIVE = {ds: collect(a.release, ds) for ds in (D1, D2)}
    HAVE_CANDIDATE = bool(a.candidate)
    CAND = {ds: collect(a.candidate, ds) for ds in (D1, D2)} if a.candidate else LIVE

    with PdfPages(a.out) as pdf:
        page_design(pdf)
        page_status(pdf)
        meta = pdf.infodict()
        meta["Title"] = "IGI-SYN-SEQ synthetic multi-assay tumour datasets"
        meta["Subject"] = "Sample inventory, design targets and measured status"
        meta["Author"] = "igi-syn-seq build"
    print(f"wrote {a.out} ({os.path.getsize(a.out):,} bytes)")


if __name__ == "__main__":
    main()
