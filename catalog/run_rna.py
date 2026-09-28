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
from igi_catalog import fusion_core
from igi_catalog.readnames import shuffle_and_rename
from igi_catalog.rna import RnaBuilder
from igi_catalog.simulate import clone_weights
from igi_catalog.transcriptome import TranscriptomeBuilder


def load_rows(path):
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def viral_records(paths, viruses, env):
    """Expressed viral transcripts, taken from the viral reference by accession."""
    out = []
    fa = paths.get("viral_fasta")
    if not fa or not os.path.exists(fa):
        return out
    import pysam
    ref = pysam.FastaFile(fa)
    names = {r.split("|")[-2].strip() if "|" in r else r: r for r in ref.references}
    for v in viruses:
        if str(v.get("expressed", "")).lower() in ("", "false", "none"):
            continue
        acc = v.get("accession", "")
        key = next((r for r in ref.references if acc.split(".")[0] in r), None)
        if key is None:
            continue
        seq = ref.fetch(key)
        out.append({"id": f"{v['event_id']}|{v['virus']}", "source": "virus", "gene": v["virus"],
                    "sequence": seq, "tpm": 80.0 if v.get("expressed") == "True" else 8.0,
                    "clone": v.get("clone", "T"), "hap": 0,
                    "in_normal": bool(float(v.get("normal_trace_copies") or 0) > 0),
                    "normal_tpm": 0.5 if float(v.get("normal_trace_copies") or 0) > 0 else 0.0})
    return out


def cta_records(expressed, env):
    """Cancer-testis antigen transcripts at their designed expression tier."""
    tier_tpm = {"T0": 0.0, "T1": 1.5, "T10": 15.0, "T100": 120.0, "T1000": 600.0}
    out = []
    for e in expressed:
        if e.get("class") != "cta":
            continue
        t = env.tx.get(e.get("transcript", ""))
        if t is None:
            continue
        tpm = tier_tpm.get(e.get("target_expression_tier", "T10"), 15.0)
        if tpm <= 0:
            continue
        from igi_catalog.transcriptome import exon_sequence
        out.append({"id": f"{e['event_id']}|{e['gene']}", "source": "cta", "gene": e["gene"],
                    "sequence": exon_sequence(env.genome, t), "tpm": tpm,
                    "clone": e.get("clone", "T"), "hap": int(e.get("haplotype", 0) or 0)})
    return out


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
    ap.add_argument("--min-tpm", type=float, default=0.5)
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--max-genes", type=int, default=None)
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
    snv_rows = load_rows(os.path.join(a.catalog_dir, f"{a.dataset}.snv_indel.tsv"))
    events = read_events(os.path.join(a.catalog_dir, f"{a.dataset}.snv_indel.tsv"))
    fusions = load_rows(os.path.join(a.catalog_dir, f"{a.dataset}.fusions.tsv"))
    expressed = load_rows(os.path.join(a.catalog_dir, f"{a.dataset}.expressed.tsv"))
    viruses = load_rows(os.path.join(a.catalog_dir, f"{a.dataset}.viruses.tsv"))

    weights = clone_weights(env.clones, ds_cfg["purity"])
    rb = RnaBuilder(env, ds_cfg, weights, min_tpm=a.min_tpm)
    rb.index_events(snv_rows)

    # every transcript of an expressed gene, not just one per gene
    tx = [t for t in env.tx.values() if t.chrom in chroms and t.exons
          and t.gene_type in ("protein_coding", "lncRNA")]
    if a.max_genes:
        tx = tx[:a.max_genes]
    rb.add_transcripts(tx, events)
    print(f"[{a.dataset}] {len(tx)} transcripts considered -> {len(rb.records)} records "
          f"({time.time() - t0:.0f}s)", flush=True)

    # designed classes
    tb = TranscriptomeBuilder(env, a.dataset, env.clones, env.expr)
    for f in fusions:
        if f.get("chrom_5p") in chroms:
            tb.add_fusion_transcript(f, fusion_core)
    for e in expressed:
        if e.get("chrom") not in chroms:
            continue
        if e["class"] == "erv":
            tb.add_erv_transcript(e)
        elif e["class"] == "splice":
            tb.add_splice_isoform(e)
    designed = [dict(r) for r in tb.records]
    designed += cta_records(expressed, env)
    designed += viral_records(paths, viruses, env)
    rb.add_designed(designed)
    print(f"  designed classes: {len(designed)} source transcripts -> {len(rb.records)} total records",
          flush=True)

    rb.assign_record_ids()
    plan, realised = rb.coverage_plan(a.pairs, a.read_len, n_bins=a.bins)
    print(f"  {len(plan)} abundance bins, coverage {min(c for _i, c, _r in plan):.3f}-"
          f"{max(c for _i, c, _r in plan):.1f}x, ~{realised:,} pairs", flush=True)

    pieces = []
    for n, (idx, cov, recs) in enumerate(plan):
        fa = rb.write_bin(recs, os.path.join(work, f"bin{idx:03d}.fa"))
        pre = os.path.join(work, f"bin{idx:03d}_")
        cmd = paths["art_cmd"].format(args=(f"-ss HS25 -i {fa} -p -l {a.read_len} -f {cov:.4f} "
                                            f"-m 250 -s 50 -rs {a.seed + n} -na -o {pre}"))
        subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        if os.path.exists(f"{pre}1.fq"):
            pieces.append((f"{pre}1.fq", f"{pre}2.fq"))
        os.remove(fa)

    # concatenate, then shuffle and rename in one disk-based pass so peak memory does not scale with
    # the library, and write the map from each read name back to the record it came from
    r1 = os.path.join(a.out, f"{a.dataset}_chr1to6_rna_R1.fastq.gz")
    r2 = os.path.join(a.out, f"{a.dataset}_chr1to6_rna_R2.fastq.gz")
    rmap = os.path.join(a.out, f"{a.dataset}_chr1to6_rna_readmap.tsv.gz")
    cat1 = os.path.join(work, "all_1.fq")
    cat2 = os.path.join(work, "all_2.fq")
    for idx, out in ((0, cat1), (1, cat2)):
        with open(out, "w") as fh:
            subprocess.run(["cat"] + [p[idx] for p in pieces], stdout=fh, check=True)
    n_reads = shuffle_and_rename(cat1, cat2, r1, r2, rmap, work, seed=a.seed)
    print(f"  {n_reads:,} pairs written, names Illumina-style, map in {os.path.basename(rmap)}", flush=True)

    rb.attach_counts(rmap, a.read_len)
    rb.write_manifest(os.path.join(a.out, f"{a.dataset}_chr1to6_rna_transcripts.tsv"))
    meta = {"dataset": a.dataset, "chroms": sorted(chroms), "pairs_target": a.pairs,
            "pairs_planned": realised, "pairs_written": n_reads, "records": len(rb.records),
            "bins": len(plan), "read_map": rmap,
            "r1": r1, "r2": r2, "runtime_s": round(time.time() - t0)}
    with open(os.path.join(a.out, f"{a.dataset}_chr1to6_rna.json"), "w") as fh:
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
