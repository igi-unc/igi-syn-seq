"""Assemble the per-clone, per-haplotype transcript records that every RNA assay is built from.

This was inline in `run_rna.py`. Kinnex bulk, Kinnex single cell and the ONT RNA assays need exactly the
same record set -- the same reference transcripts, the same designed fusions, ERVs, splice isoforms,
cancer-testis antigens and expressed viruses, reconciled the same way -- so it lives here and each driver
calls it rather than keeping its own copy that would drift.

What a record is: one transcript, on one clone, on one haplotype, with that clone's and haplotype's germline
and somatic edits applied, plus an `abundance`. Everything downstream differs only in how abundance is
turned into reads: short-read RNA spends it as coverage, and the long-read assays spend it as a molecule
count.
"""
import csv
import os

from .rna import RnaBuilder
from .transcriptome import TranscriptomeBuilder
from . import fusion_core


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
        # the transcript is passed through rather than a fixed sequence, so the builder can apply this
        # clone and haplotype's germline and somatic edits: a CTA built from bare reference exons would
        # carry none of the variants designed into it
        out.append({"id": f"{e['event_id']}|{e['gene']}", "source": "cta", "gene": e["gene"],
                    "transcript": t, "tpm": tpm,
                    "clone": e.get("clone", "T"), "hap": int(e.get("haplotype", 0) or 0),
                    # de-repression is epigenetic and acts on both alleles, so a CTA is expressed from
                    # both haplotypes; emitting it from one made every germline het inside it homozygous
                    "haplotypes": [0, 1], "chrom": t.chrom, "pos": (t.start + t.end) // 2,
                    "reconcile": "replace",
                    "normal_gene_tpm": float(e.get("normal_tissue_tpm_p95") or 0)})
    return out


def assemble(env, ds_cfg, paths, catalog_dir, dataset, chroms, weights, events,
             min_tpm=0.01, max_genes=None, log=print):
    """Build the RnaBuilder for one dataset and return it with all records in place.

    `events` is the parsed somatic event table (from `genome_build.read_events`); `weights` the clone
    weights from `simulate.clone_weights`. Both are passed in rather than derived here so a caller that
    already has them does not build them twice.
    """
    snv_rows = load_rows(os.path.join(catalog_dir, f"{dataset}.snv_indel.tsv"))
    fusions = load_rows(os.path.join(catalog_dir, f"{dataset}.fusions.tsv"))
    expressed = load_rows(os.path.join(catalog_dir, f"{dataset}.expressed.tsv"))
    viruses = load_rows(os.path.join(catalog_dir, f"{dataset}.viruses.tsv"))

    rb = RnaBuilder(env, ds_cfg, weights, min_tpm=min_tpm)
    rb.index_events(snv_rows)

    # every transcript of an expressed gene, not just one per gene
    tx = [t for t in env.tx.values() if t.chrom in chroms and t.exons
          and t.gene_type in ("protein_coding", "lncRNA")]
    if max_genes:
        tx = tx[:max_genes]
    rb.add_transcripts(tx, events)
    log(f"[{dataset}] {len(tx)} transcripts considered -> {len(rb.records)} records")

    tb = TranscriptomeBuilder(env, dataset, env.clones, env.expr)
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
    rb.add_designed(designed, events)
    rb.reconcile_designed(log=log)
    log(f"  designed classes: {len(designed)} source transcripts -> {len(rb.records)} total records")
    return rb
