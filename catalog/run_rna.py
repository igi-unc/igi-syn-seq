#!/usr/bin/env python3
"""Build the bulk RNA library for one dataset.

Depth per transcript follows its abundance, both haplotypes of every clone are built so allele-specific
expression and copy-number dosage appear in the reads, and the designed classes (fusions, ERVs, splice
isoforms, cancer-testis antigens and expressed viruses) are included and inherited by descendant clones.

    python3 run_rna.py --design design.yaml --paths paths.yaml --dataset IGI-SYN-SEQ-01 \
        --chroms chr1,chr2,chr3,chr4,chr5,chr6 --pairs 20000000 --out /path/out --work /scratch
"""
import argparse
import csv
import json
import os
import random
import subprocess
import time

import yaml

from igi_catalog.designer import build_env
from igi_catalog.genome_build import read_events
from igi_catalog.pipeline import quality_model
from igi_catalog.readnames import shuffle_and_rename
from igi_catalog.rna_assembly import assemble, load_rows
from igi_catalog.simulate import clone_weights


def _chrkey(c):
    t = c[3:]
    return (0, int(t)) if t.isdigit() else (1, t)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", required=True)
    ap.add_argument("--paths", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--chroms", required=True)
    ap.add_argument("--catalog-dir", default="output")
    ap.add_argument("--pairs", type=int, default=20_000_000)
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--bins", type=int, default=60)
    ap.add_argument("--min-tpm", type=float, default=0.01,
                    help="expression floor. LENS filters on the 50th percentile of the sample's own "
                         "non-zero TPM, so a high floor removes the low tail and lifts that bar: 0.5 "
                         "keeps 27%% of non-zero transcripts and puts p50 at 7.8x the real value. The "
                         "default is low enough that truncation comes from sequencing depth instead.")
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--max-genes", type=int, default=None)
    ap.add_argument("--label", default=None, help="slice name used in output filenames")
    a = ap.parse_args()

    chroms = set(a.chroms.split(","))
    design = yaml.safe_load(open(a.design))
    paths = yaml.safe_load(open(a.paths))
    env = build_env(paths, design, a.dataset)
    ds_cfg = design["datasets"][a.dataset]
    work = os.path.join(a.work, a.dataset)
    os.makedirs(work, exist_ok=True)
    os.makedirs(a.out, exist_ok=True)
    rng = random.Random(a.seed)

    t0 = time.time()
    events = read_events(os.path.join(a.catalog_dir, f"{a.dataset}.snv_indel.tsv"))
    weights = clone_weights(env.clones, ds_cfg["purity"])
    rb = assemble(env, ds_cfg, paths, a.catalog_dir, a.dataset, chroms, weights, events,
                  min_tpm=a.min_tpm, max_genes=a.max_genes,
                  log=lambda m: print(m, flush=True))
    print(f"  records assembled in {time.time() - t0:.0f}s", flush=True)

    rb.assign_record_ids()
    plan, realised = rb.coverage_plan(a.pairs, a.read_len, n_bins=a.bins)
    print(f"  {len(plan)} abundance bins, coverage {min(c for _i, c, _r in plan):.3f}-"
          f"{max(c for _i, c, _r in plan):.1f}x, ~{realised:,} pairs", flush=True)

    qual = quality_model(paths)
    print(f"  quality model: {qual}", flush=True)
    pieces = []
    for n, (idx, cov, recs) in enumerate(plan):
        fa = rb.write_bin(recs, os.path.join(work, f"bin{idx:03d}.fa"))
        pre = os.path.join(work, f"bin{idx:03d}_")
        cmd = paths["art_cmd"].format(args=(f"{qual} -i {fa} -p -l {a.read_len} -f {cov:.4f} "
                                            f"-m 250 -s 50 -rs {a.seed + n} -na -o {pre}"))
        subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        if os.path.exists(f"{pre}1.fq"):
            pieces.append((f"{pre}1.fq", f"{pre}2.fq"))
        os.remove(fa)

    # concatenate, then shuffle and rename in one disk-based pass so peak memory does not scale with
    # the library, and write the map from each read name back to the record it came from
    # the slice label comes from --chroms, not a hardcoded "chr1to6": a chr4 run was being written under
    # a chr1to6 name, which would have put a development slice in a release's filenames
    label = a.label or ("chr1to6" if chroms == {f"chr{i}" for i in range(1, 7)}
                        else "wg" if len(chroms) >= 23 else "_".join(sorted(chroms, key=_chrkey)))
    r1 = os.path.join(a.out, f"{a.dataset}_{label}_rna_R1.fastq.gz")
    r2 = os.path.join(a.out, f"{a.dataset}_{label}_rna_R2.fastq.gz")
    rmap = os.path.join(a.out, f"{a.dataset}_{label}_rna_readmap.tsv.gz")
    cat1 = os.path.join(work, "all_1.fq")
    cat2 = os.path.join(work, "all_2.fq")
    for idx, out in ((0, cat1), (1, cat2)):
        with open(out, "w") as fh:
            subprocess.run(["cat"] + [p[idx] for p in pieces], stdout=fh, check=True)
    n_reads = shuffle_and_rename(cat1, cat2, r1, r2, rmap, work, seed=a.seed)
    print(f"  {n_reads:,} pairs written, names Illumina-style, map in {os.path.basename(rmap)}", flush=True)

    rb.attach_counts(rmap, a.read_len)
    rb.write_manifest(os.path.join(a.out, f"{a.dataset}_{label}_rna_transcripts.tsv"))
    # Every builder records its seed; see the note in run_wes.py.
    meta = {"dataset": a.dataset, "seed": a.seed, "assay": "bulk RNA", "min_tpm": a.min_tpm,
            "chroms": sorted(chroms), "pairs_target": a.pairs,
            "pairs_planned": realised, "pairs_written": n_reads, "records": len(rb.records),
            "bins": len(plan), "quality_model": qual, "read_map": rmap,
            "r1": r1, "r2": r2, "runtime_s": round(time.time() - t0)}
    with open(os.path.join(a.out, f"{a.dataset}_{label}_rna.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    for p in pieces:
        for f in p:
            if os.path.exists(f):
                os.remove(f)
    for f in (cat1, cat2):
        if os.path.exists(f):
            os.remove(f)
    print(f"  wrote {r1} ({os.path.getsize(r1) / 1e6:.0f} MB)", flush=True)


if __name__ == "__main__":
    main()
