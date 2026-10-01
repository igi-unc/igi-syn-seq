#!/usr/bin/env python3
"""Build the Kinnex bulk RNA library (PacBio MAS-seq, 8-mer array) for one dataset.

Delivers both halves of what a real run produces: the pre-skera `hifi_reads.bam`, where each read is one
whole array, and the post-skera `segmented.bam`, produced by running the real `skera split`. Shipping both
is the point -- the segmented BAM is what LENS consumes, and the pre-skera BAM is what lets skera itself be
validated.

The chain:

    transcript records (shared with bulk RNA)  ->  molecules sampled by abundance x size selection
      ->  arrays: A0 S0 A1 S1 ... S7 A8, reverse-complemented about half the time
      ->  badread, one full-length read per array        ->  pre-skera BAM (zm, np, rq tags)
      ->  real `skera split` with the recovered adapters ->  segmented BAM

Why one read per array works: badread emits a whole reference sequence when the requested read length
exceeds it. Verified on 30 synthetic arrays -- all 32 reads started at position 0 and ran to the end -- so
an array plays the part of a ZMW insert and badread plays the part of ccs.

    python3 run_kinnex_bulk.py --design design.yaml --paths paths.yaml --dataset IGI-SYN-SEQ-01 \
        --chroms chr1,...,chrY --segments 15000000 --out /path/out --work /scratch
"""
import argparse
import json
import math
import os
import random
import subprocess
import time

import yaml

from igi_catalog.designer import build_env
from igi_catalog.genome_build import read_events
from igi_catalog.kinnex import MasArrays, SizeSelection
from igi_catalog.longread import chunk_fasta, combine_fastq, run_chunks
from igi_catalog.pacbio_bam import source_of, write_hifi_bam
from igi_catalog.rna_assembly import assemble
from igi_catalog.simulate import clone_weights


def write_arrays(mas, molecules, rng, path, start=0):
    """Write arrays to FASTA and return the membership map {array_name: [record ids]}."""
    members = {}
    with open(path, "w") as fh:
        for name, seq, mols in mas.build(molecules, rng, start_index=start):
            members[name] = [m["rec"] for m in mols]
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    return members


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", default="design.yaml")
    ap.add_argument("--paths", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--chroms", required=True)
    ap.add_argument("--catalog-dir", default="output")
    ap.add_argument("--segments", type=int, default=15_000_000,
                    help="target segmented reads, i.e. molecules. One molecule yields one segment")
    ap.add_argument("--min-tpm", type=float, default=0.01)
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--jobs", type=int, default=None)
    ap.add_argument("--arrays-per-chunk", type=int, default=2000)
    ap.add_argument("--fill-passes", type=int, default=3,
                    help="badread samples arrays at random, so some get no read and others two. Each extra "
                         "pass re-simulates only the arrays still missing. At 1x quantity a pass covers "
                         "1-1/e = 63%% of what is left, so coverage runs 63%%, 86.5%%, 95%%, 98.2%% for a "
                         "cumulative 1.00x, 1.37x, 1.50x, 1.55x of compute. Three passes is the knee. "
                         "Asking for 2x on the first pass instead reaches 99.3%% but costs 2.2x, so the "
                         "fills are the cheaper way to the same place")
    ap.add_argument("--label", default=None)
    ap.add_argument("--max-genes", type=int, default=None)
    a = ap.parse_args()

    chroms = set(a.chroms.split(","))
    design = yaml.safe_load(open(a.design))
    paths = yaml.safe_load(open(a.paths))
    env = build_env(paths, design, a.dataset)
    ds_cfg = design["datasets"][a.dataset]
    work = os.path.join(a.work, f"{a.dataset}_kinnex_bulk")
    os.makedirs(work, exist_ok=True)
    os.makedirs(a.out, exist_ok=True)
    rng = random.Random(a.seed)
    t0 = time.time()

    events = read_events(os.path.join(a.catalog_dir, f"{a.dataset}.snv_indel.tsv"))
    weights = clone_weights(env.clones, ds_cfg["purity"])
    rb = assemble(env, ds_cfg, paths, a.catalog_dir, a.dataset, chroms, weights, events,
                  min_tpm=a.min_tpm, max_genes=a.max_genes, log=lambda m: print(m, flush=True))
    rb.assign_record_ids()

    ad = paths["kinnex_adapters_fasta"]
    pr = paths["kinnex_mas_profile"]
    sz = SizeSelection(paths.get("kinnex_size_selection"))
    mas = MasArrays(ad, pr, size_selection=sz)
    print(f"  MAS profile: {mas.n_seg} segments/array, {len(mas.adapters)} adapters, "
          f"{mas.p_reverse_expected:.0%} reverse-oriented expected, size selection "
          f"{'on' if sz.enabled else 'OFF (no curve configured)'}", flush=True)

    molecules = mas.sample_molecules(rb.records, a.segments, rng)
    mean_len = sum(len(m["sequence"]) for m in molecules) / len(molecules)
    print(f"  {len(molecules):,} molecules sampled, mean cDNA length {mean_len:,.0f} bp", flush=True)

    # Arrays, chunked by record count: a chunk boundary inside an array would hand the simulator a
    # fragment as though it were a whole insert.
    fa_all = os.path.join(work, "arrays.fa")
    members = write_arrays(mas, molecules, rng, fa_all)
    print(f"  {len(members):,} arrays, mean {len(molecules)/len(members):.2f} segments each", flush=True)

    br = paths["badread_cmd"]
    n_par = a.jobs or min(8, int(os.environ.get("SLURM_CPUS_PER_TASK", 8)))
    timeout_s = int(paths.get("pacbio_chunk_timeout_s", 3600))
    # The requested length must exceed the longest array so every read spans its whole array.
    want_len = int(mas.n_seg * max(mean_len, 1000) * 3) + 50_000

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
        return run_chunks(jobs, n_par, timeout_s, log=lambda m: print(m, flush=True))

    parts = simulate(fa_all, "br0", a.seed + 1000)
    covered = set()
    for p in parts:
        with open(p) as fh:
            for line in fh:
                if line.startswith("@"):
                    s = source_of(line.split(None, 1)[1] if " " in line else "")
                    if s:
                        covered.add(s)

    # Fill passes: re-simulate only the arrays that got no read.
    for p_i in range(a.fill_passes):
        missing = [n for n in members if n not in covered]
        if not missing:
            break
        frac = len(missing) / len(members)
        print(f"  fill pass {p_i + 1}: {len(missing):,} array(s) with no read ({frac:.2%})", flush=True)
        sub_fa = os.path.join(work, f"refill{p_i}.fa")
        keep, w = set(missing), None
        with open(fa_all) as src, open(sub_fa, "w") as dst:
            for line in src:
                if line.startswith(">"):
                    w = line[1:].split()[0] in keep
                if w:
                    dst.write(line)
        more = simulate(sub_fa, f"br{p_i + 1}", a.seed + 2000 + p_i * 500)
        parts += more
        for p in more:
            with open(p) as fh:
                for line in fh:
                    if line.startswith("@"):
                        s = source_of(line.split(None, 1)[1] if " " in line else "")
                        if s:
                            covered.add(s)
    lost = [n for n in members if n not in covered]
    if lost:
        print(f"  {len(lost):,} array(s) never produced a read ({len(lost)/len(members):.2%}); their "
              f"molecules are recorded as not sequenced", flush=True)

    label = a.label or "full"
    cat = os.path.join(work, "all.fq")
    combine_fastq(parts, cat, log=lambda m: print(m, flush=True))

    pre_bam = os.path.join(a.out, f"{a.dataset}_{label}_kinnex_bulk_hifi_reads.bam")
    n_pre, by_src = write_hifi_bam([cat], pre_bam, one_per_source=True,
                                   np_passes=int(round(float(
                                       json.load(open(pr)).get("passes_np_mean") or 8))),
                                   sample=a.dataset, library=f"{a.dataset}_kinnex_bulk")
    print(f"  pre-skera: {n_pre:,} reads (one per array) -> {os.path.basename(pre_bam)}", flush=True)

    seg_bam = os.path.join(a.out, f"{a.dataset}_{label}_kinnex_bulk_segmented.bam")
    sk = paths["skera_cmd"]
    subprocess.run(sk.format(args=f"split {pre_bam} {ad} {seg_bam}"), shell=True, check=True)
    sam = paths["samtools_cmd"]
    n_seg = int(subprocess.run(sam.format(args=f"view -c {seg_bam}"), shell=True,
                               capture_output=True, text=True).stdout.strip() or 0)
    print(f"  post-skera: {n_seg:,} segments -> {os.path.basename(seg_bam)}", flush=True)

    # Truth: which molecules each array carried, and which arrays reached a read.
    mpath = os.path.join(a.out, f"{a.dataset}_{label}_kinnex_bulk_arrays.tsv.gz")
    import gzip
    with gzip.open(mpath, "wt") as fh:
        fh.write("array\tread_name\tsegment_index\trecord\n")
        for arr in sorted(members):
            rn = by_src.get(arr, "")
            for i, rec in enumerate(members[arr]):
                fh.write(f"{arr}\t{rn}\t{i}\t{rec}\n")
    rb.write_manifest(os.path.join(a.out, f"{a.dataset}_{label}_kinnex_bulk_transcripts.tsv"))

    meta = {"dataset": a.dataset, "seed": a.seed, "assay": "Kinnex bulk RNA",
            "platform": "PacBio Revio", "array": f"MAS {mas.n_seg}-mer",
            "min_tpm": a.min_tpm, "segments_target": a.segments,
            "molecules": len(molecules), "arrays": len(members),
            "arrays_with_a_read": len(covered), "arrays_lost": len(lost),
            "mean_cdna_length": round(mean_len, 1),
            "size_selection": "measured" if sz.enabled else "flat",
            "reads_pre_skera": n_pre, "segments_post_skera": n_seg,
            "pre_skera_bam": pre_bam, "segmented_bam": seg_bam, "array_map": mpath,
            "runtime_s": round(time.time() - t0)}
    with open(os.path.join(a.out, f"{a.dataset}_{label}_kinnex_bulk.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
        fh.write("\n")
    print(f"  done in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
