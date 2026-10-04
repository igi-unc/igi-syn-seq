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
from igi_catalog.vdj import assemble, cdr3_nucleotide, load_segments


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", default="design.yaml")
    ap.add_argument("--paths", required=True)
    ap.add_argument("--dataset", required=True)
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
    roster = CellRoster(a.cells, load_whitelist(paths["single_cell_whitelist"]),
                        cell_seed, env.clones)
    tcells = [c for c in roster.cells if c.get("clonotype")]
    log(f"  roster {len(roster.cells):,} cells (seed {cell_seed}); "
        f"{len(tcells):,} carry a clonotype")

    segs = load_segments(paths["gtf"], env.genome)
    log(f"  {len(segs)} TR segments from the GTF (the annotation cache keeps only 1 of 202)")

    # One transcript per chain per clonotype, so expanded cells share a junction exactly.
    by_clono = {}
    for c in tcells:
        k = c["clonotype"]
        if k in by_clono:
            continue
        crng = random.Random(f"{cell_seed}:vdj:{k}")
        chains = []
        for chain, pep in (("TRB", c.get("trb")), ("TRA", c.get("tra")), ("TRA", c.get("tra2"))):
            if not pep:
                continue
            seq, used = assemble(segs, chain, cdr3_nucleotide(pep), crng)
            if seq:
                chains.append((chain, pep, seq, used))
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
        fm.write("read\tbarcode\tumi\tchain\tcdr3_aa\tv\tj\tc\n")
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
                        nm = f"{a.dataset}:{n_reads + k + 1}"
                        f1.write(f"@{nm} 1:N:0:1\n{r1_sequence(bc,u)}\n+\n{r1_quality(rng)}\n")
                        f2.write(f"@{nm} 2:N:0:1\n{s}\n+\n{q}\n")
                        fm.write(f"{nm}\t{bc}\t{u}\t{chain}\t{pep}\t"
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
            "segments_available": len(segs),
            "caveat": ("CDR3s are synthetic peptides, not real recombinants; this library tests clonotype "
                       "clustering and expansion, not germline V-J assignment"),
            "r1": r1_out, "r2": r2_out, "molecule_map": mpath,
            "runtime_s": round(time.time() - t0)}
    with open(os.path.join(a.out, f"{a.dataset}_{a.label or 'full'}_tenx_tcr.json"), "w") as fh:
        json.dump(meta, fh, indent=2); fh.write("\n")
    log(f"  done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
