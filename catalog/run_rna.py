#!/usr/bin/env python3
"""Build the bulk RNA library for one dataset over a chromosome set.

Transcripts are rebuilt per clone so germline and somatic variants appear in RNA, designed events override
the baseline expression of their gene, and molecules are drawn per clone in proportion to that clone's
share of the tumour plus the normal fraction. Fusion, ERV and splice-isoform transcripts are added as
their own records so they can be recovered by a caller that assembles or aligns to them.

    python3 run_rna.py --design design.yaml --paths paths.yaml --dataset IGI-SYN-SEQ-01 \
        --chroms chr1,chr2,chr3,chr4,chr5,chr6 --pairs 20000000 --out /path/out
"""
import argparse
import csv
import json
import os
import subprocess
import time

import yaml

from igi_catalog.designer import build_env
from igi_catalog.genome_build import read_events, germline_edits, somatic_edits, EditSet
from igi_catalog import fusion_core
from igi_catalog.simulate import clone_weights
from igi_catalog.transcriptome import TranscriptomeBuilder, exon_sequence


def load_rows(path):
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", required=True)
    ap.add_argument("--paths", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--chroms", required=True)
    ap.add_argument("--catalog-dir", default="output")
    ap.add_argument("--pairs", type=int, default=20_000_000)
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--max-genes", type=int, default=None)
    a = ap.parse_args()

    chroms = set(a.chroms.split(","))
    design = yaml.safe_load(open(a.design))
    paths = yaml.safe_load(open(a.paths))
    env = build_env(paths, design, a.dataset)
    purity = design["datasets"][a.dataset]["purity"]
    os.makedirs(a.out, exist_ok=True)
    os.makedirs(a.work, exist_ok=True)

    t0 = time.time()
    events = read_events(os.path.join(a.catalog_dir, f"{a.dataset}.snv_indel.tsv"))
    fusions = load_rows(os.path.join(a.catalog_dir, f"{a.dataset}.fusions.tsv"))
    expressed = load_rows(os.path.join(a.catalog_dir, f"{a.dataset}.expressed.tsv"))

    tb = TranscriptomeBuilder(env, a.dataset, env.clones, env.expr)
    weights = clone_weights(env.clones, purity)

    # reference transcripts, one record per clone, carrying that clone's variants
    genes = [t for t in env.rep.values() if t.chrom in chroms]
    if a.max_genes:
        genes = genes[:a.max_genes]
    for t in genes:
        base = env.expr.transcript(t.tid, t.gene_id)
        if base < 0.05:
            continue
        for src in weights:
            clone = "T" if src == "NORMAL" else src
            hap = 0
            es = EditSet(t.chrom)
            g = germline_edits(env.germline, t.chrom, hap, t.start, t.end)
            es.edits += g.edits
            if src != "NORMAL":
                s = somatic_edits(events.get(t.chrom, []), env.clones, t.chrom, hap, clone, t.start, t.end)
                es.edits += s.edits
            seq = exon_sequence(env.genome, t, es)
            tb.records.append({"id": f"{t.tid}|{src}", "source": "reference", "gene": t.gene_name,
                               "chrom": t.chrom, "hap": hap, "clone": src, "sequence": seq, "tpm": base})

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
    print(f"[{a.dataset}] {len(tb.records)} transcript records in {time.time() - t0:.0f}s", flush=True)

    # one FASTA per clone source; ART is run at a coverage that yields that source's share of the pairs
    pieces = []
    for i, (src, w) in enumerate(sorted(weights.items())):
        recs = [r for r in tb.records if r["clone"] in (src, "T" if src == "NORMAL" else src) and r["sequence"]]
        if src == "NORMAL":
            recs = [r for r in tb.records if r["clone"] == "NORMAL" and r["sequence"]]
        if not recs:
            continue
        fa = os.path.join(a.work, f"{a.dataset}_{src}_rna.fa")
        total_bp = 0
        with open(fa, "w") as fh:
            for r in recs:
                tpm = float(r["tpm"] or 0)
                if tpm < 0.05:
                    continue
                fh.write(f">{r['id']} tpm={tpm} source={r['source']}\n")
                for j in range(0, len(r["sequence"]), 60):
                    fh.write(r["sequence"][j:j + 60] + "\n")
                total_bp += len(r["sequence"])
        if total_bp == 0:
            continue
        want_pairs = a.pairs * w
        cov = max(1.0, want_pairs * 2 * a.read_len / total_bp)
        pre = os.path.join(a.work, f"{a.dataset}_{src}_rna_")
        cmd = paths["art_cmd"].format(args=(f"-ss HS25 -i {fa} -p -l {a.read_len} -f {cov:.3f} "
                                            f"-m 250 -s 50 -rs {a.seed + i} -na -o {pre}"))
        subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        pieces.append((f"{pre}1.fq", f"{pre}2.fq", src, round(cov, 2)))
        print(f"    {src}: {len(recs)} transcripts, {total_bp / 1e6:.1f} Mb, coverage {cov:.2f}", flush=True)

    r1 = os.path.join(a.out, f"{a.dataset}_chr1to6_rna_R1.fastq.gz")
    r2 = os.path.join(a.out, f"{a.dataset}_chr1to6_rna_R2.fastq.gz")
    for idx, out in ((0, r1), (1, r2)):
        parts = [p[idx] for p in pieces if os.path.exists(p[idx])]
        with open(out, "wb") as fh:
            cat = subprocess.Popen(["cat"] + parts, stdout=subprocess.PIPE)
            gz = subprocess.Popen(["gzip", "-c"], stdin=cat.stdout, stdout=fh)
            cat.stdout.close()
            if gz.wait() or cat.wait():
                raise RuntimeError(f"failed writing {out}")
    tb.write_table(os.path.join(a.out, f"{a.dataset}_chr1to6_rna_transcripts.tsv"))
    meta = {"dataset": a.dataset, "chroms": sorted(chroms), "pairs_target": a.pairs,
            "records": len(tb.records), "sources": {s: c for _x, _y, s, c in pieces},
            "r1": r1, "r2": r2, "runtime_s": round(time.time() - t0)}
    with open(os.path.join(a.out, f"{a.dataset}_chr1to6_rna.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    for p in pieces:
        for f in p[:2]:
            if os.path.exists(f):
                os.remove(f)
    print(f"  wrote {r1} ({os.path.getsize(r1) / 1e6:.0f} MB)", flush=True)


if __name__ == "__main__":
    main()
