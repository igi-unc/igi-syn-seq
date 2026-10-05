#!/usr/bin/env python3
"""Build the 10x Chromium 5' V(D)J (TCR) library for one dataset.

The fifteenth and last sample type. Structure is the same as 10x gene expression -- R1 is the 16 bp barcode
plus 10 bp UMI, R2 is cDNA -- so it reuses `tenx` for both. What differs is which molecules exist:

- Only T cells carry a receptor. The roster assigns clonotypes to its CD8, CD4 and Treg cells, so the
  library covers those and no others; a barcode absent here but present in gene expression is a non-T cell,
  which is the correct relationship between the two assays.
- The transcripts are assembled V + CDR3 + J + C from real GRCh38 segments, read straight from the GTF
  because the annotation cache drops TR gene types (one of 202 survives into env.tx). A caller aligns reads
  to germline segments to infer the junction, so invented framework would make every read unassignable.
- Cells of one clonotype share a nucleotide junction, because the reverse translation is deterministic.
  That is what clonal expansion means and what the assay exists to test.
- Some cells carry two alpha chains. The roster already models that (`tra2`), and both are emitted.

Inherited limitation, from `cells.CellRoster._cdr3`: the CDR3s are synthetic peptides rather than real
recombinants, so this library tests clonotype clustering and expansion, not germline V-J gene assignment.

    python3 run_tenx_tcr.py --design design.yaml --paths paths.yaml --dataset IGI-SYN-SEQ-01 \
        --cells 4000 --reads-per-cell 5000 --out /path/out --work /scratch
"""
import argparse
import collections
import csv
import gzip
import json
import os
import random
import subprocess
import time

import yaml

from igi_catalog.cells import CellRoster, load_whitelist
from igi_catalog.designer import build_env
from igi_catalog.molecules import umi
from igi_catalog.pipeline import quality_model
from igi_catalog.tenx import (R1_LEN, R2_LEN, cellranger_names, five_prime_window,
                              r1_quality, r1_sequence, reads_for_molecule)
from igi_catalog.readnames import NAME_SPACE_STRIDE, illumina_name
from igi_catalog.vdj import assemble, genes




def link_flagposts(roster, catalog_dir, dataset, log, n=5):
    """Attach flagpost neoantigens to the largest clonotypes, and report how many landed.

    By rank, because an expanded clonotype is the one a user will actually recover, so it is the one worth
    asserting an antigen for. The antigens come from the catalog's own flagpost column rather than from a
    list repeated here, so the two cannot drift apart.
    """
    path = os.path.join(catalog_dir, f"{dataset}.snv_indel.tsv")
    if not os.path.exists(path):
        log(f"  no catalog at {path}; flagpost_antigen left empty")
        return 0
    flag = []
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for row in rd:
            if str(row.get("flagpost", "")).lower() in ("1", "true", "yes") and row.get("best_peptide"):
                flag.append((row.get("gene", ""), row.get("best_peptide", ""),
                             row.get("best_allele", "")))
    if not flag:
        log("  catalog has no flagpost peptides; flagpost_antigen left empty")
        return 0
    top = sorted(roster.clonotypes, key=lambda c: -c["size"])[:n]
    for i, c in enumerate(top):
        gene, pep, allele = flag[i % len(flag)]
        c["flagpost_antigen"] = f"{gene}:{pep}:{allele}"
    log(f"  linked {len(top)} flagpost clonotypes to antigens, largest first "
        f"({', '.join(c['clonotype'] for c in top)})")
    return len(top)


def cdr3_nt_for(roster, clonotype, prefix, cell):
    """The nucleotide junction for one chain of one clonotype.

    The roster's clonotype table carries it for the beta and the primary alpha. A second alpha is per
    cell, so it is recombined again from the same keyed stream the roster used, which reproduces it
    exactly rather than guessing.
    """
    from igi_catalog import vdj
    if prefix in ("trb", "tra"):
        for row in roster.clonotypes:
            if row["clonotype"] == clonotype:
                return row[f"{prefix}_cdr3_nt"]
    r = vdj.recombine("TRA", random.Random(f"{roster.seed}:vdj:{clonotype}:TRA2:{cell['barcode']}"))
    return r["cdr3_nt"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", default="design.yaml")
    ap.add_argument("--paths", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--catalog-dir", default="output",
                    help="where the designer wrote <dataset>.snv_indel.tsv; the flagpost "
                         "antigen link is read from its flagpost and best_peptide columns")
    ap.add_argument("--cells", type=int, default=4000)
    ap.add_argument("--reads-per-cell", type=int, default=5000,
                    help="10x's own recommendation for V(D)J, and far fewer than gene expression needs "
                         "because only two transcripts per cell are being sequenced")
    ap.add_argument("--molecules-per-cell", type=int, default=120,
                    help="receptor transcripts captured per T cell; TCR mRNA is abundant, so a modest "
                         "molecule count carries a high duplicate rate at 5,000 reads")
    ap.add_argument("--molecules-per-chunk", type=int, default=400_000)
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--seed", type=int, default=9100,
                    help="MUST be 9100: the roster is shared with every other single-cell assay")
    ap.add_argument("--label", default=None)
    a = ap.parse_args()

    design = yaml.safe_load(open(a.design))
    paths = yaml.safe_load(open(a.paths))
    env = build_env(paths, design, a.dataset)
    work = os.path.join(a.work, f"{a.dataset}_tenx_tcr")
    os.makedirs(work, exist_ok=True)
    os.makedirs(a.out, exist_ok=True)
    rng = random.Random(a.seed)
    t0 = time.time()
    log = lambda m: print(m, flush=True)

    cell_seed = f"{a.seed}:{a.dataset}"
    # Ordinary Illumina names, the same scheme WES, WGS and bulk RNA use, each library in
    # its own slice of the name space. These two arms were writing "@<dataset>:<n>", which is
    # not a name any instrument produces and not a name Cell Ranger or a QC tool expects.
    name_base = NAME_SPACE_STRIDE * 21 + (0 if a.dataset.endswith("01") else NAME_SPACE_STRIDE // 2)
    roster = CellRoster(a.cells, load_whitelist(paths["single_cell_whitelist"]),
                        cell_seed, env.clones)
    tcells = [c for c in roster.cells if c.get("clonotype")]
    log(f"  roster {len(roster.cells):,} cells (seed {cell_seed}); "
        f"{len(tcells):,} carry a clonotype")

    log(f"  germline: {len(genes('TRB', 'v'))} TRBV, {len(genes('TRB', 'j'))} TRBJ, "
        f"{len(genes('TRA', 'v'))} TRAV, {len(genes('TRA', 'j'))} TRAJ functional genes")

    # One transcript per chain per clonotype, so expanded cells share a junction exactly. The roster
    # already recombined each clonotype from the real germline anchors and recorded which V and J it used,
    # so the transcript here is assembled from those same genes: the V gene a caller recovers from the
    # framework is the V gene the truth table names, which is the whole point of doing it this way.
    by_clono = {}
    for c in tcells:
        k = c["clonotype"]
        if k in by_clono:
            continue
        chains = []
        for chain, pre in (("TRB", "trb"), ("TRA", "tra"), ("TRA", "tra2")):
            if not c.get(pre):
                continue
            rec = {"v": c[f"{pre}_v"], "j": c[f"{pre}_j"],
                   "c": "TRBC2" if chain == "TRB" else "TRAC",
                   "cdr3_aa": c[pre], "cdr3_nt": cdr3_nt_for(roster, k, pre, c)}
            seq = assemble(rec)
            if seq:
                chains.append((chain, c[pre], seq, rec))
        by_clono[k] = chains
    n_ch = sum(len(v) for v in by_clono.values())
    log(f"  {len(by_clono):,} clonotypes -> {n_ch:,} chain transcripts "
        f"({sum(1 for v in by_clono.values() for ch,_,_,_ in v if ch=='TRA')} alpha, "
        f"{sum(1 for v in by_clono.values() for ch,_,_,_ in v if ch=='TRB')} beta)")

    mols = []
    for c in tcells:
        chains = by_clono.get(c["clonotype"]) or []
        if not chains:
            continue
        for _ in range(a.molecules_per_cell):
            chain, pep, seq, used = chains[rng.randrange(len(chains))]
            mols.append((c["barcode"], umi(rng), chain, pep, seq, used))
    mean_reads = (a.reads_per_cell * len(tcells)) / max(1, len(mols))
    log(f"  {len(mols):,} molecules over {len(tcells):,} T cells; "
        f"{mean_reads:.1f} reads per molecule")

    qual = quality_model(paths)
    r1_out, r2_out = cellranger_names(a.out, f"{a.dataset}-TCR")
    mpath = os.path.join(a.out, f"{a.dataset}_{a.label or 'full'}_tcr_molecules.tsv.gz")
    art = paths["art_cmd"]
    n_reads = n_chunks = 0

    with gzip.open(r1_out, "wt") as f1, gzip.open(r2_out, "wt") as f2, gzip.open(mpath, "wt") as fm:
        fm.write("read\tbarcode\tumi\tchain\tcdr3_aa\tcdr3_nt\tv\tj\tc\n")
        chunk, ci = [], 0

        def flush(chunk, ci):
            nonlocal n_reads, n_chunks
            if not chunk:
                return
            fa = os.path.join(work, f"c{ci:05d}.fa")
            with open(fa, "w") as fh:
                for i, m in enumerate(chunk):
                    w = five_prime_window(m[4])
                    if len(w) >= R2_LEN:
                        fh.write(f">m{i}\n{w}\n")
            pre = os.path.join(work, f"c{ci:05d}_")
            cmd = art.format(args=(f"{qual} -i {fa} -l {R2_LEN} -f {R2_LEN/450.0*1.02:.4f} "
                                   f"-ss HS25 -rs {a.seed + ci} -na -o {pre}"))
            subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
            fq = f"{pre}.fq" if os.path.exists(f"{pre}.fq") else f"{pre}1.fq"
            if os.path.exists(fq):
                with open(fq) as fh:
                    k = 0
                    while k < len(chunk):
                        h = fh.readline()
                        if not h:
                            break
                        s = fh.readline().strip(); fh.readline(); q = fh.readline().strip()
                        bc, u, chain, pep, _seq, used = chunk[k]
                        nm = illumina_name(name_base + n_reads + k)
                        f1.write(f"@{nm} 1:N:0:1\n{r1_sequence(bc,u)}\n+\n{r1_quality(rng, kind='tcr')}\n")
                        f2.write(f"@{nm} 2:N:0:1\n{s}\n+\n{q}\n")
                        fm.write(f"{nm}\t{bc}\t{u}\t{chain}\t{pep}\t{used['cdr3_nt']}\t"
                                 f"{used['v']}\t{used['j']}\t{used['c']}\n")
                        k += 1
                n_reads += k
            for p in (fa, fq):
                if os.path.exists(p):
                    os.remove(p)
            n_chunks += 1

        for m in mols:
            for _ in range(reads_for_molecule(rng, mean_reads)):
                chunk.append(m)
                if len(chunk) >= a.molecules_per_chunk:
                    flush(chunk, ci); ci += 1; chunk = []
        flush(chunk, ci)

    log(f"  {n_reads:,} read pairs over {n_chunks} ART chunk(s) -> {os.path.basename(r2_out)}")

    # The design asks for five flagpost clonotypes "mapped in the truth bundle to specific flagpost
    # neoantigens". The catalog owns which events are flagposts, so the link is made here, where it is
    # readable, and written into the clonotype table rather than left for a user to guess at.
    n_linked = link_flagposts(roster, a.catalog_dir, a.dataset, log)
    roster.write_clonotypes(os.path.join(a.out, f"{a.dataset}_{a.label or 'full'}_tcr_clonotypes.tsv"))

    sizes = collections.Counter(c["clonotype"] for c in tcells)
    meta = {"dataset": a.dataset, "seed": a.seed, "assay": "10x 5' V(D)J (TCR)",
            "cells": len(roster.cells), "t_cells": len(tcells),
            "clonotypes": len(by_clono), "chain_transcripts": n_ch,
            "reads_per_cell_target": a.reads_per_cell,
            "molecules_per_cell": a.molecules_per_cell,
            "mean_reads_per_molecule": round(mean_reads, 2),
            "read_pairs_written": n_reads,
            "largest_clonotype": max(sizes.values()) if sizes else 0,
            "trbv_genes": len(genes("TRB", "v")), "trbj_genes": len(genes("TRB", "j")),
            "trav_genes": len(genes("TRA", "v")), "traj_genes": len(genes("TRA", "j")),
            "flagpost_clonotypes_linked": n_linked,
            "caveat": ("junctions are recombined from the real GENCODE V and J germline anchors -- the "
                       "conserved Cys at IMGT 104 and the Phe of the J FGXG motif -- with a trimmed, "
                       "GC-rich N region, so the V and J a caller recovers are the ones the truth table "
                       "names. N-region nucleotides are drawn rather than taken from a real repertoire, "
                       "so junction length and composition are realistic but the clonotypes are not any "
                       "individual's"),
            "r1": r1_out, "r2": r2_out, "molecule_map": mpath,
            "runtime_s": round(time.time() - t0)}
    with open(os.path.join(a.out, f"{a.dataset}_{a.label or 'full'}_tenx_tcr.json"), "w") as fh:
        json.dump(meta, fh, indent=2); fh.write("\n")
    log(f"  done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
