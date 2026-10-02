#!/usr/bin/env python3
"""Build the Kinnex single-cell RNA library (PacBio MAS-seq, 10x 5' v2 cDNA) for one dataset.

Same array machinery as the bulk builder, with two differences that matter and one that turned out not to
exist:

1. **The molecules are 10x 5' v2 cDNAs**, not bare transcripts: a TruSeq R1 primer, the cell barcode, the
   UMI, the 5' TSO, the cDNA, polyA and the SMART primer's reverse complement. Every offset was recovered
   from the real HG002 single-cell arrays -- see `kinnex.tenx_segment` -- rather than taken from a kit
   document.
2. **The molecules come from the shared per-cell pool**, so a barcode and UMI mean the same molecule here
   as in 10x gene expression, 10x TCR and ONT single cell. That only holds because `MoleculePool.draw` is a
   pure function of the cell; it used to draw from one pool-wide stream, which made the result depend on how
   many cells had been drawn before, so two assays asking for the same cell got different molecules.
3. **The array carries sixteen cDNAs**, which a real array read states directly: a median of 16 10x TSOs
   (mean 15.1, p90 16) across 2,000 reads averaging 17,092 bp, against a single-cDNA median of 952 bp. The
   adapter layout is not resolved -- nine adapters are identified and no tenth is detectable -- so the array
   is built to reproduce skera's observed output on real data rather than from an asserted chemistry. That
   limitation is the open item for this assay; see `igi_catalog/kinnex.py` and `docs/build-reference.md`.

    python3 run_kinnex_sc.py --design design.yaml --paths paths.yaml --dataset IGI-SYN-SEQ-01 \
        --chroms chr1,...,chrY --cells 4000 --segments 12000000 --out /path/out --work /scratch
"""
import argparse
import json
import os
import random
import subprocess
import time

import yaml

from igi_catalog.cells import CellRoster, load_whitelist
from igi_catalog.designer import build_env
from igi_catalog.genome_build import read_events
from igi_catalog.kinnex import MasArrays, SizeSelection, tenx_segment
from igi_catalog.longread import chunk_fasta, combine_fastq, run_chunks
from igi_catalog.molecules import MoleculePool
from igi_catalog.pacbio_bam import source_of, write_hifi_bam
from igi_catalog.rna_assembly import assemble
from igi_catalog.simulate import clone_weights


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", default="design.yaml")
    ap.add_argument("--paths", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--chroms", required=True)
    ap.add_argument("--catalog-dir", default="output")
    ap.add_argument("--cells", type=int, default=4000)
    ap.add_argument("--segments", type=int, default=12_000_000,
                    help="target segmented reads. One molecule yields one segment")
    ap.add_argument("--mean-molecules", type=int, default=8000,
                    help="molecules captured per cell before sequencing depth is applied")
    ap.add_argument("--min-tpm", type=float, default=0.01)
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--jobs", type=int, default=None)
    ap.add_argument("--cdna-per-array", type=int, default=16,
                    help="cDNAs per array. 16 is measured directly: a real single-cell array read carries a "
                         "median of 16 10x TSOs (mean 15.1) over 17,092 bp, against a single-cDNA median of "
                         "952 bp. Only nine adapters are identified, so the first eight boundaries are "
                         "bracketed and the rest join directly; see igi_catalog/kinnex.py")
    ap.add_argument("--arrays-per-chunk", type=int, default=2000)
    ap.add_argument("--fill-passes", type=int, default=3)
    ap.add_argument("--label", default=None)
    ap.add_argument("--max-genes", type=int, default=None)
    a = ap.parse_args()

    chroms = set(a.chroms.split(","))
    design = yaml.safe_load(open(a.design))
    paths = yaml.safe_load(open(a.paths))
    env = build_env(paths, design, a.dataset)
    ds_cfg = design["datasets"][a.dataset]
    work = os.path.join(a.work, f"{a.dataset}_kinnex_sc")
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

    # The roster and the molecule pool are the shared objects: every single-cell assay must build them the
    # same way from the same seed, or a barcode stops meaning one cell across assays.
    wl = load_whitelist(paths["single_cell_whitelist"])
    roster = CellRoster(a.cells, wl, a.seed, env.clones)
    log(f"  roster: {len(roster.cells):,} cells, {roster.summary()}")

    tx = [(r["rec"], r.get("gene", ""), len(r["sequence"]), r["abundance"],
           r.get("chrom", ""), r.get("pos", 0), r.get("hap", 0),
           "" if r.get("source") == "reference" else r.get("clone", ""))
          for r in rb.records if r["abundance"] > 0 and r["sequence"]]
    seq_by_rec = {r["rec"]: r["sequence"] for r in rb.records}
    pool = MoleculePool(roster, tx, env.clones, a.seed, mean_molecules=a.mean_molecules)

    sz = SizeSelection(paths.get("kinnex_sc_size_selection"))
    mas = MasArrays(paths["kinnex_adapters_fasta"], paths["kinnex_mas_profile"], size_selection=sz)
    log(f"  MAS profile: {mas.n_seg} segments/array, {len(mas.adapters)} adapters, size selection "
        f"{'on' if sz.enabled else 'OFF (no sc curve configured)'}")

    # Draw every cell's molecules, build each as a 10x segment, then subsample to the target depth. The
    # pool is the biology; the target is the sequencing run, and they are different numbers.
    all_mols = []
    for c in roster.cells:
        for idx, u in pool.draw(c):
            all_mols.append((c["barcode"], u, tx[idx][0]))
    log(f"  pool: {len(all_mols):,} molecules over {len(roster.cells):,} cells "
        f"({len(all_mols)/max(1,len(roster.cells)):.0f} per cell)")
    if a.segments < len(all_mols):
        all_mols = rng.sample(all_mols, a.segments)
    elif a.segments > len(all_mols):
        log(f"  note: asked for {a.segments:,} segments but the pool holds {len(all_mols):,}; "
            f"sequencing every molecule once rather than inventing more")
    rng.shuffle(all_mols)

    # size selection acts on the finished molecule, which includes the 10x structure
    if sz.enabled:
        keep = [m for m in all_mols if rng.random() < min(1.0, sz.factor(len(seq_by_rec[m[2]])))]
        log(f"  size selection kept {len(keep):,} of {len(all_mols):,} molecules")
        all_mols = keep

    fa_all = os.path.join(work, "arrays.fa")
    members, alen, nseg = {}, [], 0
    mol_objs = [{"rec": r, "sequence": tenx_segment(bc, u, seq_by_rec[r], rng), "bc": bc, "umi": u}
                for bc, u, r in all_mols]
    with open(fa_all, "w") as fh:
        for name, seq, mols in mas.build_sc(mol_objs, rng, n_cdna=a.cdna_per_array):
            members[name] = [(m["bc"], m["umi"], m["rec"]) for m in mols]
            alen.append(len(seq))
            nseg += len(mols)
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    alen.sort()
    n_a = len(alen)
    log(f"  {n_a:,} arrays, {nseg:,} segments, mean {nseg/n_a:.2f} each; "
        f"length median {alen[n_a//2]:,} p90 {alen[9*n_a//10]:,} max {alen[-1]:,} bp")

    br = paths["badread_cmd"]
    n_par = a.jobs or min(8, int(os.environ.get("SLURM_CPUS_PER_TASK", 8)))
    timeout_s = int(paths.get("pacbio_chunk_timeout_s", 3600))
    want_len = int(alen[-1] * 1.5) + 20_000

    def simulate(fasta, tag, seed0):
        jobs = []
        for j, sub in enumerate(chunk_fasta(fasta, a.arrays_per_chunk, work, tag)):
            out = os.path.join(work, f"{tag}_{j:04d}.fq")
            args = (f"simulate --reference {sub} --quantity 1x "
                    f"--length {want_len},{want_len // 10} "
                    f"--identity {paths.get('pacbio_identity', '99.4,99.8,0.5')} "
                    f"--error_model {paths.get('pacbio_error_model', 'pacbio2021')} "
                    f"--qscore_model {paths.get('pacbio_qscore_model', 'pacbio2021')} "
                    f"--seed {seed0 + j} --start_adapter_seq '' --end_adapter_seq '' "
                    f"--junk_reads 0 --random_reads 0 --chimeras 0")
            jobs.append((f"{tag}/{j}", out, br.format(args=args)))
        return run_chunks(jobs, n_par, timeout_s, log=log)

    def covered_from(parts):
        got = set()
        for p in parts:
            with open(p) as fh:
                for line in fh:
                    if line.startswith("@"):
                        s = source_of(line.split(None, 1)[1] if " " in line else "")
                        if s:
                            got.add(s)
        return got

    parts = simulate(fa_all, "br0", a.seed + 1000)
    covered = covered_from(parts)
    for p_i in range(a.fill_passes):
        missing = [n for n in members if n not in covered]
        if not missing:
            break
        log(f"  fill pass {p_i + 1}: {len(missing):,} array(s) with no read "
            f"({len(missing)/len(members):.2%})")
        sub_fa = os.path.join(work, f"refill{p_i}.fa")
        keep, w = set(missing), False
        with open(fa_all) as src, open(sub_fa, "w") as dst:
            for line in src:
                if line.startswith(">"):
                    w = line[1:].split()[0] in keep
                if w:
                    dst.write(line)
        more = simulate(sub_fa, f"br{p_i + 1}", a.seed + 2000 + p_i * 500)
        parts += more
        covered |= covered_from(more)
    lost = [n for n in members if n not in covered]
    if lost:
        log(f"  {len(lost):,} array(s) never produced a read ({len(lost)/len(members):.2%})")

    label = a.label or "full"
    cat = os.path.join(work, "all.fq")
    combine_fastq(parts, cat, log=log)
    pre_bam = os.path.join(a.out, f"{a.dataset}_{label}_kinnex_sc_hifi_reads.bam")
    n_pre, by_src = write_hifi_bam([cat], pre_bam, one_per_source=True,
                                   np_passes=int(round(float(
                                       json.load(open(paths["kinnex_mas_profile"]))
                                       .get("passes_np_mean") or 8))),
                                   sample=a.dataset, library=f"{a.dataset}_kinnex_sc")
    log(f"  pre-skera: {n_pre:,} reads (one per array) -> {os.path.basename(pre_bam)}")

    seg_bam = os.path.join(a.out, f"{a.dataset}_{label}_kinnex_sc_segmented.bam")
    subprocess.run(paths["skera_cmd"].format(
        args=f"split {pre_bam} {paths['kinnex_adapters_fasta']} {seg_bam}"), shell=True, check=True)
    n_seg = int(subprocess.run(paths["samtools_cmd"].format(args=f"view -c {seg_bam}"), shell=True,
                               capture_output=True, text=True).stdout.strip() or 0)
    log(f"  post-skera: {n_seg:,} segments -> {os.path.basename(seg_bam)}")

    import gzip
    mpath = os.path.join(a.out, f"{a.dataset}_{label}_kinnex_sc_arrays.tsv.gz")
    with gzip.open(mpath, "wt") as fh:
        fh.write("array\tread_name\tsegment_index\tbarcode\tumi\trecord\n")
        for arr in sorted(members):
            rn = by_src.get(arr, "")
            for i, (bc, u, rec) in enumerate(members[arr]):
                fh.write(f"{arr}\t{rn}\t{i}\t{bc}\t{u}\t{rec}\n")
    roster.write(os.path.join(a.out, f"{a.dataset}_{label}_cells.tsv"))
    roster.write_clonotypes(os.path.join(a.out, f"{a.dataset}_{label}_clonotypes.tsv"))
    rb.write_manifest(os.path.join(a.out, f"{a.dataset}_{label}_kinnex_sc_transcripts.tsv"))

    meta = {"dataset": a.dataset, "seed": a.seed, "assay": "Kinnex single-cell RNA",
            "platform": "PacBio Revio", "array": f"MAS {mas.n_seg}-mer",
            "chemistry": "10x 5' v2", "cells": len(roster.cells),
            "mean_molecules_per_cell": a.mean_molecules, "min_tpm": a.min_tpm,
            "segments_target": a.segments, "molecules": nseg, "arrays": n_a,
            "arrays_with_a_read": len(covered), "arrays_lost": len(lost),
            "array_length": {"median": alen[n_a // 2], "p90": alen[9 * n_a // 10],
                             "max": alen[-1], "requested_read_length": want_len},
            "size_selection": "measured" if sz.enabled else "flat",
            "reads_pre_skera": n_pre, "segments_post_skera": n_seg,
            "pre_skera_bam": pre_bam, "segmented_bam": seg_bam, "array_map": mpath,
            "runtime_s": round(time.time() - t0)}
    with open(os.path.join(a.out, f"{a.dataset}_{label}_kinnex_sc.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
        fh.write("\n")
    log(f"  done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
