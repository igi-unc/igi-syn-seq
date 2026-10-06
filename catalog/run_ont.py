#!/usr/bin/env python3
"""Build an Oxford Nanopore library for one dataset: WGS, bulk RNA, or single-cell RNA.

One driver for three sample types, because they differ only in what the simulator is pointed at:

    --assay wgs      derived chromosomes, as run_pacbio_wgs.py builds them, with the ONT model
    --assay bulk_rna full-length cDNA molecules, as the Kinnex bulk builder samples them, no MAS array
    --assay sc_rna   10x 5' cDNA molecules from the shared per-cell pool, no MAS array

ONT has no array concatenation, so the long-read RNA assays are the Kinnex builders minus the array step:
one molecule is one read. That is the whole difference, and it is why these are one driver rather than three.

Model constants are measured, not guessed, and the two RNA assays and WGS use different ones:

    ont_length_mean / _sd        910 / 461    IPISRC044 R10.4.1 SUP scRNA, 2025-07-17 (sc only)
    ont_bulk_rna_length_mean/_sd 2091/ 843    HG002 Kinnex FLNC, real full-length cDNA
    ont_identity                 98.22,...    same, calibrated for Badread's -0.22 point cDNA offset
    ont_wgs_length_mean / _sd    19117/15530  ONT open data, HG002 PAW70337 R10.4.1 SUP
    ont_wgs_identity_target      0.98676      same

GIAB's own HG002 ONT is R9.4-era -- 2D reads from 2016 basecalled with Guppy v2/v3 -- so it cannot
calibrate an R10.4.1 model, which is why the genomic parameters come from ONT's open-data bucket.

    python3 run_ont.py --design design.yaml --paths paths.yaml --dataset IGI-SYN-SEQ-01 \
        --assay wgs --chrom chr1 --library tumor --depth 30 --out /path --work /scratch
"""
import argparse
import gzip
import json
import os
import shutil
import random
import subprocess
import time

import yaml

from igi_catalog.cells import CellRoster, load_whitelist
from igi_catalog.designer import build_env
from igi_catalog.genome_build import read_events
from igi_catalog.derived import build_derived
from igi_catalog.kinnex import SizeSelection, ont_cdna
from igi_catalog.longread import (_framed, chunk_fasta, combine_fastq, run_chunks,
                                  split_fasta_bp)
from igi_catalog.molecules import MoleculePool
from igi_catalog.rna_assembly import assemble
from igi_catalog.simulate import clone_weights


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", default="design.yaml")
    ap.add_argument("--paths", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--assay", required=True, choices=["wgs", "bulk_rna", "sc_rna"])
    ap.add_argument("--catalog-dir", default="output")
    # WGS
    ap.add_argument("--chrom")
    ap.add_argument("--library", default="tumor", choices=["tumor", "normal"])
    ap.add_argument("--depth", type=float, default=30.0)
    # RNA
    ap.add_argument("--chroms")
    ap.add_argument("--reads", type=int, default=20_000_000, help="target reads for the RNA assays")
    ap.add_argument("--cells", type=int, default=4000)
    ap.add_argument("--mean-molecules", type=int, default=8000)
    ap.add_argument("--min-tpm", type=float, default=0.01)
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--seed", type=int, default=9100,
                    help="single-cell assays MUST use 9100 to share the roster; WGS and bulk RNA are "
                         "independent of it")
    ap.add_argument("--jobs", type=int, default=None)
    ap.add_argument("--records-per-chunk", type=int, default=40_000)
    ap.add_argument("--label", default=None)
    ap.add_argument("--max-genes", type=int, default=None)
    a = ap.parse_args()

    design = yaml.safe_load(open(a.design))
    paths = yaml.safe_load(open(a.paths))
    env = build_env(paths, design, a.dataset)
    ds_cfg = design["datasets"][a.dataset]
    tag = f"{a.dataset}_ont_{a.assay}" + (f"_{a.chrom}_{a.library}" if a.assay == "wgs" else "")
    work = os.path.join(a.work, tag)
    os.makedirs(work, exist_ok=True)
    os.makedirs(a.out, exist_ok=True)
    rng = random.Random(a.seed)
    t0 = time.time()
    log = lambda m: print(m, flush=True)

    events = read_events(os.path.join(a.catalog_dir, f"{a.dataset}.snv_indel.tsv"))
    weights = clone_weights(env.clones, ds_cfg["purity"])

    # Model constants: genomic and cDNA ONT are different libraries and different reads.
    if a.assay == "wgs":
        lmean = int(paths.get("ont_wgs_length_mean", 19117))
        lsd = int(paths.get("ont_wgs_length_sd", 15530))
        ident = paths.get("ont_wgs_identity", "98.9,99.6,1.5")
    elif a.assay == "bulk_rna":
        # A separate fit from the single-cell one. The sc model (910/461) is fitted from IPISRC044's 10x
        # 5' single-cell cDNA, which is a 5'-biased, deliberately truncated library -- correct for the sc
        # assay and wrong for this one, which the design calls "cDNA, full length". At 910 bp most
        # transcripts are not covered end to end, so isoform detection, which is the whole reason to have
        # long-read RNA and what the 30 designed splice events are measured by, had nothing to work with.
        lmean = int(paths.get("ont_bulk_rna_length_mean", 2091))
        lsd = int(paths.get("ont_bulk_rna_length_sd", 843))
        # Identity is shared with the single-cell fit: both are cDNA on the same R10.4.1 SUP chemistry,
        # and it is length, not accuracy, that separates a full-length library from a 5'-biased one.
        ident = paths.get("ont_identity", "98.41,99.6,1.5")
    else:
        lmean = int(paths.get("ont_length_mean", 910))
        lsd = int(paths.get("ont_length_sd", 461))
        ident = paths.get("ont_identity", "98.22,99.6,1.5")
    emodel = paths.get("ont_error_model", "nanopore2023")
    qmodel = paths.get("ont_qscore_model", "nanopore2023")
    log(f"  ONT model: length {lmean:,}/{lsd:,}, identity {ident}, {emodel}/{qmodel}")

    fa = os.path.join(work, "src.fa")
    n_src = 0
    truth = []          # (record_name, meta...) for the RNA assays

    pieces = []
    if a.assay == "wgs":
        if not a.chrom:
            raise SystemExit("--assay wgs needs --chrom")
        # The same derived-chromosome construction the PacBio arm uses, shared rather than reimplemented:
        # one copy per (source, haplotype, copy kind), weighted by that source's share of sequenced
        # molecules, so allele fractions follow the clone model by construction.
        drng = random.Random(f"{a.seed}:{a.dataset}:{a.chrom}:{a.library}")
        pieces, plan_rows = build_derived(env, ds_cfg, a.dataset, a.chrom, a.library, a.catalog_dir,
                                          events, a.depth, drng, work, log=log)
        n_src = len(pieces)
        log(f"  {n_src} derived chromosome copy/copies")
    else:
        rb = assemble(env, ds_cfg, paths, a.catalog_dir, a.dataset,
                      set((a.chroms or "").split(",")), weights, events,
                      min_tpm=a.min_tpm, max_genes=a.max_genes, log=log)
        rb.assign_record_ids()
        seq_by_rec = {r["rec"]: r["sequence"] for r in rb.records}
        sz = SizeSelection(paths.get("kinnex_sc_size_selection" if a.assay == "sc_rna"
                                     else "kinnex_size_selection"))
        if a.assay == "sc_rna":
            cell_seed = f"{a.seed}:{a.dataset}"
            roster = CellRoster(a.cells, load_whitelist(paths["single_cell_whitelist"]),
                                cell_seed, env.clones)
            tx = [(r["rec"], r.get("gene", ""), len(r["sequence"]), r["abundance"],
                   r.get("chrom", ""), r.get("pos", 0), r.get("hap", 0),
                   "" if r.get("source") == "reference" else r.get("clone", ""))
                  for r in rb.records if r["abundance"] > 0 and r["sequence"]]
            pool = MoleculePool(roster, tx, env.clones, cell_seed, mean_molecules=a.mean_molecules)
            mols = [(c["barcode"], u, tx[i][0]) for c in roster.cells for i, u in pool.draw(c)]
            log(f"  roster {len(roster.cells):,} cells, pool {len(mols):,} molecules "
                f"(seed {cell_seed})")
        else:
            live = [r for r in rb.records if r["abundance"] > 0 and r["sequence"]]
            w = [r["abundance"] * sz.factor(len(r["sequence"])) for r in live]
            mols = [("", "", r["rec"]) for r in rng.choices(live, weights=w, k=a.reads)]
            log(f"  {len(mols):,} molecules sampled by abundance x size selection")
        if sz.enabled and a.assay == "sc_rna":
            mols = [m for m in mols if rng.random() < min(1.0, sz.factor(len(seq_by_rec[m[2]])))]
            log(f"  size selection kept {len(mols):,}")
        if a.reads < len(mols):
            mols = rng.sample(mols, a.reads)
        with open(fa, "w") as fh:
            for i, (bc, u, rec) in enumerate(mols):
                s = ont_cdna(bc, u, seq_by_rec[rec], rng, single_cell=(a.assay == "sc_rna"))
                if len(s) < 120:
                    continue
                nm = f"m{i:09d}"
                fh.write(f">{nm}\n{s}\n")
                truth.append((nm, bc, u, rec))
                n_src += 1
        log(f"  {n_src:,} cDNA molecules written; one molecule is one read (no MAS array)")
        quantity = "1x"

    br = paths["badread_cmd"]
    n_par = a.jobs or min(48, int(os.environ.get("SLURM_CPUS_PER_TASK", 8)))
    timeout_s = int(paths.get("pacbio_chunk_timeout_s", 3600))
    want = max(lmean * 6, 200_000) if a.assay != "wgs" else lmean
    jobs = []
    if a.assay == "wgs":
        # Each derived copy carries its own coverage, so it is simulated separately, and each copy is split
        # by BASE PAIRS. An earlier version passed one chunk per copy, which left a 249 Mb chromosome whole:
        # 14 chunks over 48 slots, each needing about 4 hours against a 1 hour timeout. That cost 61 of 96
        # tasks. At 20 Mb a chunk, chr1 becomes ~13 chunks per copy and every slot is used.
        chunk_bp = int(paths.get("pacbio_chunk_bp", 20_000_000))
        k = 0
        for tag, pfa, cov in pieces:
            for sub in split_fasta_bp(pfa, chunk_bp, work, f"w_{tag}"):
                out = os.path.join(work, f"ont_{k:05d}.fq")
                args = (f"simulate --reference {sub} --quantity {cov:.4f}x "
                        f"--length {lmean},{lsd} --identity {ident} "
                        f"--error_model {emodel} --qscore_model {qmodel} --seed {a.seed + k} "
                        f"--start_adapter_seq '' --end_adapter_seq '' --junk_reads 0 "
                        f"--random_reads 0 --chimeras 0")
                jobs.append((f"ont/{tag}", out, br.format(args=args)))
                k += 1
    else:
        for j, sub in enumerate(chunk_fasta(fa, a.records_per_chunk, work, "src")):
            out = os.path.join(work, f"ont_{j:05d}.fq")
            args = (f"simulate --reference {sub} --quantity 1x "
                    f"--length {want},{want//10} --identity {ident} "
                    f"--error_model {emodel} --qscore_model {qmodel} --seed {a.seed + j} "
                    f"--start_adapter_seq '' --end_adapter_seq '' --junk_reads 0 --random_reads 0 "
                    f"--chimeras 0")
            jobs.append((f"ont/{j}", out, br.format(args=args)))
    label = a.label or (f"{a.chrom}_{a.library}" if a.assay == "wgs" else "full")
    cat = os.path.join(work, "all.fq")
    # The ONT FASTQ *is* the deliverable, so the simulator's description is stripped and written to
    # a separate map: Badread's header names the source reference, strand and coordinates for a
    # genomic read and the source molecule for an RNA one. See combine_fastq's docstring.
    map_path = os.path.join(a.out, f"{a.dataset}_{label}_ont_{a.assay}_read_map.tsv.gz")

    # Resume. The combined FASTQ is the expensive artefact: the ds-02 bulk RNA run spent 58,506 s in
    # badread and another hour combining, then hit the 24 h wall limit during the final gzip and the whole
    # thing would otherwise be redone. If a complete, frame-aligned all.fq is already sitting in the work
    # directory, the simulation and the combine are both already paid for.
    if os.path.exists(cat) and os.path.getsize(cat) > 1_000_000 and _framed(cat):
        with open(cat, "rb") as fh:
            n_reads = sum(blk.count(b"\n") for blk in iter(lambda: fh.read(1 << 24), b"")) // 4
        log(f"  resuming: {os.path.basename(cat)} is complete ({os.path.getsize(cat) / 1e9:.1f} GB, "
            f"{n_reads:,} reads); skipping simulation and combine")
    else:
        parts = run_chunks(jobs, n_par, timeout_s, log=log)
        n_reads = combine_fastq(parts, cat, log=log, map_path=map_path)

    out_fq = os.path.join(a.out, f"{a.dataset}_{label}_ont_{a.assay}.fastq.gz")
    # pigz, not Python's gzip. Compressing 85 GB single-threaded is what ran out of wall clock.
    t_gz = time.time()
    if shutil.which("pigz"):
        subprocess.run(f"pigz -p {min(16, n_par)} -c {cat} > {out_fq}", shell=True, check=True)
    else:
        with open(cat, "rb") as src, gzip.open(out_fq, "wb") as dst:
            for b in iter(lambda: src.read(1 << 24), b""):
                dst.write(b)
    log(f"  compressed in {time.time() - t_gz:.0f}s")
    os.remove(cat)
    log(f"  {n_reads:,} reads -> {os.path.basename(out_fq)} "
        f"({os.path.getsize(out_fq)/1e9:.1f} GB)")

    if truth:
        tp = os.path.join(a.out, f"{a.dataset}_{label}_ont_{a.assay}_molecules.tsv.gz")
        with gzip.open(tp, "wt") as fh:
            fh.write("molecule\tbarcode\tumi\trecord\n")
            for nm, bc, u, rec in truth:
                fh.write(f"{nm}\t{bc}\t{u}\t{rec}\n")

    meta = {"dataset": a.dataset, "seed": a.seed, "assay": f"ONT {a.assay}",
            "platform": "Oxford Nanopore R10.4.1", "length_mean": lmean, "length_sd": lsd,
            "identity": ident, "error_model": emodel, "sources": n_src, "reads": n_reads,
            "fastq": out_fq, "runtime_s": round(time.time() - t0)}
    if a.assay == "wgs":
        meta.update({"chrom": a.chrom, "library": a.library, "depth": a.depth})
    with open(os.path.join(a.out, f"{a.dataset}_{label}_ont_{a.assay}.json"), "w") as fh:
        json.dump(meta, fh, indent=2); fh.write("\n")
    log(f"  done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
