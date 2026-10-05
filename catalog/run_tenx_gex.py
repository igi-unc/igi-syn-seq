#!/usr/bin/env python3
"""Build the 10x Chromium 5' v2 gene-expression library for one dataset.

    cell roster + molecule pool (shared with Kinnex sc, 10x TCR, ONT sc)
      -> per molecule: the 5'-proximal cDNA window, written once per PCR duplicate
      -> ART over those windows with the fitted NovaSeq X profile  ->  R2
      -> R1 synthesised as barcode + UMI, which we know rather than simulate
      -> Cell Ranger-named FASTQ pair + a molecule truth table

R1 is not simulated from a template: it is the 16 bp barcode and 10 bp UMI this builder assigned. Putting
sequencing error into a barcode would move cells between barcodes, which is a different experiment from the
one being set up.

    python3 run_tenx_gex.py --design design.yaml --paths paths.yaml --dataset IGI-SYN-SEQ-01 \
        --chroms chr1,...,chrY --cells 4000 --reads-per-cell 40000 --out /path/out --work /scratch
"""
import argparse
import gzip
import json
import os
import random
import subprocess
import time

import yaml

from igi_catalog.cells import CellRoster, load_whitelist
from igi_catalog.designer import build_env
from igi_catalog.genome_build import read_events
from igi_catalog.molecules import MoleculePool
from igi_catalog.pipeline import quality_model
from igi_catalog.rna_assembly import assemble
from igi_catalog.simulate import clone_weights
from igi_catalog.tenx import (R1_LEN, R2_LEN, cellranger_names, five_prime_window,
                              r1_quality, r1_sequence, reads_for_molecule)
from igi_catalog.readnames import NAME_SPACE_STRIDE, illumina_name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", default="design.yaml")
    ap.add_argument("--paths", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--chroms", required=True)
    ap.add_argument("--catalog-dir", default="output")
    ap.add_argument("--cells", type=int, default=4000)
    ap.add_argument("--reads-per-cell", type=int, default=40_000)
    ap.add_argument("--mean-molecules", type=int, default=8000)
    ap.add_argument("--min-tpm", type=float, default=0.01)
    ap.add_argument("--molecules-per-chunk", type=int, default=400_000,
                    help="molecules per ART invocation; each contributes a 450 bp record")
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--seed", type=int, default=9100,
                    help="MUST match every other single-cell assay or the roster and molecule pool "
                         "diverge and a barcode stops meaning one cell across assays")
    ap.add_argument("--label", default=None)
    ap.add_argument("--max-genes", type=int, default=None)
    a = ap.parse_args()

    chroms = set(a.chroms.split(","))
    design = yaml.safe_load(open(a.design))
    paths = yaml.safe_load(open(a.paths))
    env = build_env(paths, design, a.dataset)
    ds_cfg = design["datasets"][a.dataset]
    work = os.path.join(a.work, f"{a.dataset}_tenx_gex")
    os.makedirs(work, exist_ok=True)
    os.makedirs(a.out, exist_ok=True)
    rng = random.Random(a.seed)
    t0 = time.time()
    log = lambda m: print(m, flush=True)

    events = read_events(os.path.join(a.catalog_dir, f"{a.dataset}.snv_indel.tsv"))
    weights = clone_weights(env.clones, ds_cfg["purity"])
    rb = assemble(env, ds_cfg, paths, a.catalog_dir, a.dataset, chroms, weights, events,
                  min_tpm=a.min_tpm, max_genes=a.max_genes, log=log)
    rb.assign_record_ids()

    # The roster seed is scoped to the dataset. Seeding it from --seed alone gave both datasets the SAME
    # 4,000 barcodes and the same barcode-to-cell-type map, differing only in clone labels: two independent
    # tumours sharing one cell suspension, and a barcode collision for anyone who pools the datasets. The
    # scope keeps what matters -- every assay of ONE dataset shares a roster, because they all derive it
    # from the same string -- while making the two datasets independent.
    cell_seed = f"{a.seed}:{a.dataset}"
    # Ordinary Illumina names, the same scheme WES, WGS and bulk RNA use, each library in
    # its own slice of the name space. These two arms were writing "@<dataset>:<n>", which is
    # not a name any instrument produces and not a name Cell Ranger or a QC tool expects.
    name_base = NAME_SPACE_STRIDE * 20 + (0 if a.dataset.endswith("01") else NAME_SPACE_STRIDE // 2)
    wl = load_whitelist(paths["single_cell_whitelist"])
    roster = CellRoster(a.cells, wl, cell_seed, env.clones)
    log(f"  roster: {len(roster.cells):,} cells (seed {a.seed}; shared with every single-cell assay)")

    tx = [(r["rec"], r.get("gene", ""), len(r["sequence"]), r["abundance"],
           r.get("chrom", ""), r.get("pos", 0), r.get("hap", 0),
           "" if r.get("source") == "reference" else r.get("clone", ""))
          for r in rb.records if r["abundance"] > 0 and r["sequence"]]
    seq_by_rec = {r["rec"]: r["sequence"] for r in rb.records}
    pool = MoleculePool(roster, tx, env.clones, cell_seed, mean_molecules=a.mean_molecules)

    # Draw the pool, then decide how many times each molecule is sequenced. Reads per cell divided by
    # molecules per cell is the duplicate rate the UMIs exist to collapse; one read per molecule would
    # make the library duplicate-free and the UMIs decorative.
    mols = []
    for c in roster.cells:
        for idx, u in pool.draw(c):
            mols.append((c["barcode"], u, tx[idx][0]))
    mean_reads = a.reads_per_cell / max(1, len(mols) / len(roster.cells))
    log(f"  pool: {len(mols):,} molecules over {len(roster.cells):,} cells "
        f"({len(mols)/len(roster.cells):.0f} per cell); {mean_reads:.1f} reads per molecule")

    # One FASTA record per read: ART then emits one read per record at -f matching R2_LEN/window.
    # Writing the window once per duplicate is what makes duplicates share a barcode and UMI.
    qual = quality_model(paths)
    r1_out, r2_out = cellranger_names(a.out, f"{a.dataset}-GEX")
    rmap = os.path.join(a.out, f"{a.dataset}_{a.label or 'full'}_gex_molecules.tsv.gz")
    n_reads = n_chunks = 0
    art = paths["art_cmd"]

    with gzip.open(r1_out, "wt") as f1, gzip.open(r2_out, "wt") as f2, gzip.open(rmap, "wt") as fm:
        fm.write("read\tbarcode\tumi\trecord\n")
        chunk, chunk_i = [], 0
        def flush(chunk, chunk_i):
            nonlocal n_reads, n_chunks
            if not chunk:
                return
            fa = os.path.join(work, f"c{chunk_i:05d}.fa")
            with open(fa, "w") as fh:
                for i, (_bc, _u, rec) in enumerate(chunk):
                    w = five_prime_window(seq_by_rec[rec])
                    if len(w) < R2_LEN:
                        continue
                    fh.write(f">m{i}\n{w}\n")
            pre = os.path.join(work, f"c{chunk_i:05d}_")
            cov = R2_LEN / 450.0 * 1.02   # ~1 read per record
            cmd = art.format(args=(f"{qual} -i {fa} -l {R2_LEN} -f {cov:.4f} -ss HS25 "
                                   f"-rs {a.seed + chunk_i} -na -o {pre}"))
            subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
            fq = f"{pre}.fq" if os.path.exists(f"{pre}.fq") else f"{pre}1.fq"
            if os.path.exists(fq):
                with open(fq) as fh:
                    k = 0
                    while True:
                        h = fh.readline()
                        if not h:
                            break
                        s = fh.readline().strip(); fh.readline(); q = fh.readline().strip()
                        if k >= len(chunk):
                            break
                        bc, u, rec = chunk[k]
                        name = illumina_name(name_base + n_reads + k)
                        f1.write(f"@{name} 1:N:0:1\n{r1_sequence(bc,u)}\n+\n{r1_quality(rng, kind='gex')}\n")
                        f2.write(f"@{name} 2:N:0:1\n{s}\n+\n{q}\n")
                        fm.write(f"{name}\t{bc}\t{u}\t{rec}\n")
                        k += 1
                n_reads += k
            for p in (fa, fq):
                if os.path.exists(p):
                    os.remove(p)
            n_chunks += 1

        for bc, u, rec in mols:
            for _ in range(reads_for_molecule(rng, mean_reads)):
                chunk.append((bc, u, rec))
                if len(chunk) >= a.molecules_per_chunk:
                    flush(chunk, chunk_i); chunk_i += 1; chunk = []
        flush(chunk, chunk_i)

    log(f"  {n_reads:,} read pairs over {n_chunks} ART chunk(s) -> {os.path.basename(r2_out)}")
    roster.write(os.path.join(a.out, f"{a.dataset}_{a.label or 'full'}_cells.tsv"))
    roster.write_clonotypes(os.path.join(a.out, f"{a.dataset}_{a.label or 'full'}_clonotypes.tsv"))

    meta = {"dataset": a.dataset, "seed": a.seed, "assay": "10x 5' v2 gene expression",
            "chemistry": "10x Chromium 5' v2", "cells": len(roster.cells),
            "reads_per_cell_target": a.reads_per_cell, "mean_molecules_per_cell": a.mean_molecules,
            "molecules": len(mols), "mean_reads_per_molecule": round(mean_reads, 2),
            "read_pairs_written": n_reads, "r1_len": R1_LEN, "r2_len": R2_LEN,
            "five_prime_window": 450, "min_tpm": a.min_tpm,
            "r1": r1_out, "r2": r2_out, "molecule_map": rmap,
            "runtime_s": round(time.time() - t0)}
    with open(os.path.join(a.out, f"{a.dataset}_{a.label or 'full'}_tenx_gex.json"), "w") as fh:
        json.dump(meta, fh, indent=2); fh.write("\n")
    log(f"  done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
