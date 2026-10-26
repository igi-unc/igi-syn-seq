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
import hashlib
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
from igi_catalog.derived import build_derived, viral_pieces
from igi_catalog.kinnex import SizeSelection, ont_cdna
from igi_catalog.readnames import stream_seed
from igi_catalog import longread as _longread_mod
from igi_catalog.longread import (ONT_BASECALL_GPU, _framed, chunk_fasta, combine_fastq,
                                  run_chunks, split_fasta_bp)
from igi_catalog.molecules import MoleculePool
from igi_catalog.rna_assembly import assemble
from igi_catalog.simulate import clone_weights


def _combine_stamp():
    """A fingerprint of the code that produces all.fq, so a resume cannot cross a code change.

    WHY THIS EXISTS. On 2026-10-12 the ONT header was replaced with a measured one (section 30). The job
    was cancelled and resubmitted to pick the fix up. Two of its 96 tasks had been cancelled mid-flight
    with a complete all.fq already written by the PRE-FIX code, and the resume below reused it -- so
    IGI-SYN-SEQ-01 chr1 and chr8 tumour shipped the invented header the rebuild existed to remove, with a
    fresh mtime on them. Cancelling a job does not unwrite what its tasks already wrote, and a resume that
    asks only "is this file complete?" will happily carry the old behaviour across the fix.

    The hash covers this file and longread.py, which between them decide every byte of all.fq: the
    simulator invocation here and the header rewrite plus read map in combine_fastq there.
    """
    h = hashlib.sha256()
    for mod in (os.path.abspath(__file__), os.path.abspath(_longread_mod.__file__)):
        with open(mod, "rb") as fh:
            h.update(fh.read())
    return h.hexdigest()[:16]


def _resume_block(cat, map_path, log):
    """Why all.fq beside it cannot be resumed from, or None when it can.

    Three things have to hold, and `_framed` plus a size floor is none of them:

      * a completion record exists. `_framed` tests that the file ENDS ON A RECORD BOUNDARY, which is
        exactly where a killed combine leaves it, since it writes whole records. The old check logged
        "all.fq is complete" on the strength of a test that cannot see truncation at all.
      * the record's code stamp matches this code.
      * the read map is still where the record says, because combine_fastq writes the FASTQ and the map
        together. Skipping the combine skipped the map too, so a resumed task delivered a library with no
        truth beside it -- which is what finally surfaced this, as a chr1to6 merge failing on a missing
        shard two stages downstream.
    """
    done = cat + ".done"
    if not os.path.exists(done):
        return "no completion record beside it"
    try:
        rec = json.load(open(done))
    except (ValueError, OSError) as e:
        return f"its completion record is unreadable ({e})"
    stamp = _combine_stamp()
    if rec.get("code") != stamp:
        return f"it was written by code {rec.get('code')}, and this is {stamp}"
    if not os.path.exists(rec.get("map") or ""):
        return f"its read map is gone ({rec.get('map')})"
    with open(cat, "rb") as fh:
        n = sum(blk.count(b"\n") for blk in iter(lambda: fh.read(1 << 24), b"")) // 4
    if n != rec.get("reads"):
        return f"it holds {n:,} reads against the {rec.get('reads'):,} recorded; it is truncated"
    log(f"  resuming: all.fq is complete and current ({os.path.getsize(cat) / 1e9:.1f} GB, "
        f"{n:,} reads, code {stamp}); skipping simulation and combine")
    return None


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
        pieces += viral_pieces(paths.get("viral_fasta"), a.catalog_dir, a.dataset, a.chrom,
                               a.library, env, ds_cfg["purity"], a.depth, work, log=log)
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
                        f"--error_model {emodel} --qscore_model {qmodel} "
                        f"--seed {stream_seed(a.seed, a.dataset, 'ont', a.assay) + k} "
                        f"--start_adapter_seq '' --end_adapter_seq '' --junk_reads 0 "
                        f"--random_reads 0 --chimeras 0")
                jobs.append((f"ont/{tag}", out, br.format(args=args)))
                k += 1
    else:
        for j, sub in enumerate(chunk_fasta(fa, a.records_per_chunk, work, "src")):
            out = os.path.join(work, f"ont_{j:05d}.fq")
            args = (f"simulate --reference {sub} --quantity 1x "
                    f"--length {want},{want//10} --identity {ident} "
                    f"--error_model {emodel} --qscore_model {qmodel} "
                    f"--seed {stream_seed(a.seed, a.dataset, 'ont', a.assay) + j} "
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
    # thing would otherwise be redone. But resuming is only safe when the file on disk was produced by
    # THIS code and is provably whole -- see _resume_block for the two deliverables that proved otherwise.
    n_reads = None
    if os.path.exists(cat) and os.path.getsize(cat) > 1_000_000 and _framed(cat):
        why = _resume_block(cat, map_path, log)
        if why is None:
            with open(cat, "rb") as fh:
                n_reads = sum(blk.count(b"\n") for blk in iter(lambda: fh.read(1 << 24), b"")) // 4
        else:
            # Loud, because this is a silent-wrong-data path whenever it is quiet: the file looks finished.
            log(f"  NOT resuming from all.fq: {why}; re-simulating from scratch")
            os.remove(cat)
            if os.path.exists(cat + ".done"):
                os.remove(cat + ".done")
    if n_reads is None:
        parts = run_chunks(jobs, n_par, timeout_s, log=log)
        # Real ONT run metadata, in the layout MEASURED from 511,401 real PromethION R10.4.1 reads --
        # see ONT_HEADER_FIELDS in longread.py for the measurement and for what the first version of
        # this got wrong by inventing the shape instead of reading one off disk.
        #
        # Every value is derived from the dataset and assay, so it is stable across reruns and distinct
        # per library. `runid` is a UUID because that is what the real field holds; it is built from the
        # scoped seed rather than uuid4() so a rebuild reproduces it.
        #
        # sample_id is the SYNTHETIC dataset name, never an upstream donor accession. The ds-02 germline
        # baseline carries its real donor's cohort accession as its VCF sample name, and this field ends
        # up in every read of every delivered file.
        rs = stream_seed(a.seed, a.dataset, "ont_runid", a.assay)
        h = f"{rs:08x}{stream_seed(a.seed, a.dataset, 'ont_runid2', a.assay):08x}"
        run_meta = {
            "runid": f"{h[:8]}-{h[8:12]}-{h[12:16]}-{h[:4]}-{h[8:]}{h[:4]}",
            "flow_cell_id": f"PAW{stream_seed(a.seed, a.dataset, 'ont_fc', a.assay) % 100000:05d}",
            "basecall_gpu": ONT_BASECALL_GPU,
            "protocol_group_id": f"{a.dataset}_ont_{a.assay}",
            "sample_id": a.dataset,
            "basecall_model_version_id": "dna_r10.4.1_e8.2_400bps_sup@v4.2.0",
            # A fixed, deterministic run start. Not the wall clock: two builds of the same library would
            # then differ in every read header, and the release is meant to be reproducible from a seed.
            "start_epoch": 1760000000 + (stream_seed(a.seed, a.dataset, "ont_t0", a.assay) % 86400),
        }
        n_reads = combine_fastq(parts, cat, log=log, map_path=map_path, run_meta=run_meta)
        # The completion record, written only after combine_fastq returns -- so its presence is the one
        # trustworthy statement that all.fq is whole. It carries the code stamp and the map path so a
        # later resume can check that neither the code nor the map moved underneath it.
        with open(cat + ".done", "w") as fh:
            json.dump({"reads": n_reads, "code": _combine_stamp(), "map": map_path}, fh)

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
    if os.path.exists(cat + ".done"):
        os.remove(cat + ".done")   # a record must never outlive the file it describes
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
    else:
        # The RNA assays are given a read COUNT, and it belongs in the record beside the delivered one.
        # Without it a release-level yield check has nothing to compare against and can only report the
        # number it found -- which is the shape of check that let the 10x GEX library ship missing every
        # transcript under 450 nt. The WGS assay is given a depth instead, and acceptance measures that.
        meta["reads_target"] = a.reads
    with open(os.path.join(a.out, f"{a.dataset}_{label}_ont_{a.assay}.json"), "w") as fh:
        json.dump(meta, fh, indent=2); fh.write("\n")
    log(f"  done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
