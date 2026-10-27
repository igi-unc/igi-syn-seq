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

from .genome import normalize_bases
from .rna import RnaBuilder
from .transcriptome import TranscriptomeBuilder
from . import fusion_core


def load_rows(path):
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _viral_cds(gff_path, accession):
    """(gene, start1, end1, strand) for every annotated CDS of one viral accession."""
    import re
    out = []
    if not gff_path or not os.path.exists(gff_path):
        return out
    acc = accession.split(".")[0]
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS" or acc not in f[0]:
                continue
            attrs = f[8]
            m = re.search(r"Parent=([^;.]+)", attrs) or re.search(r"product=([^;]+)", attrs)
            gene = (m.group(1) if m else "CDS").replace("%20", "_")
            try:
                out.append((gene, int(f[3]), int(f[4]), f[6] or "+"))
            except ValueError:
                continue
    return out


def viral_records(paths, viruses, env):
    """Expressed viral TRANSCRIPTS, one per annotated CDS, taken from the viral reference.

    This used to emit the whole viral genome as a single "transcript" -- HPV16 7,905 bp, EBV 171,823 bp.
    Bulk and single-cell draw from the same record set, but the single-cell size-selection curve is built
    around ~1 kb 10x cDNA and returns 0 above its range, so every viral molecule was discarded before it
    became a read: ds-01's single-cell long-read libraries carried no HPV16 or EBV at all while the ONT
    bulk library from the same design held 1,572 HPV16 reads. The single-cell lane was right -- it was
    rejecting an object that was never shaped like a transcript.

    Widening the size curve would have been the wrong fix, because it would admit a 171 kb molecule into
    a 10x library, which no real run produces. Emitting the annotated CDS features instead gives both
    lanes viral content at realistic lengths: HPV16 E6 is 477 bp and E7 297 bp, well inside the window.

    The designed TPM is split across a virus's CDS features rather than inflated, so total viral
    expression is unchanged and only its shape differs. That split is uniform, which is a simplification:
    in a real integrated HPV the early genes dominate and the late ones are largely silent.
    """
    out = []
    fa = paths.get("viral_fasta")
    gff = paths.get("viral_cds_gff")
    if not fa or not os.path.exists(fa):
        return out
    import pysam
    from .genome import revcomp
    ref = pysam.FastaFile(fa)
    for v in viruses:
        if str(v.get("expressed", "")).lower() in ("", "false", "none"):
            continue
        acc = v.get("accession", "")
        key = next((r for r in ref.references if acc.split(".")[0] in r), None)
        if key is None:
            continue
        genome = normalize_bases(ref.fetch(key).upper())
        tpm = 80.0 if v.get("expressed") == "True" else 8.0
        trace = float(v.get("normal_trace_copies") or 0)
        cds = _viral_cds(gff, acc)
        if not cds:
            # No annotation for this accession: fall back to the genome rather than emitting nothing,
            # and say so, because a silent empty result here is how the original defect read.
            print(f"  [viral] no CDS annotation for {acc}; falling back to the whole genome "
                  f"({len(genome):,} bp) -- single-cell size selection will drop it", flush=True)
            cds = [(v.get("virus", "genome"), 1, len(genome), "+")]
        share = tpm / len(cds)
        for gene, start1, end1, strand in cds:
            seq = genome[start1 - 1:end1]
            if strand == "-":
                seq = revcomp(seq)
            if len(seq) < 100:
                continue
            # EBV's internal repeat carries BWRF1 eight times, and a couple of genes are annotated as
            # two CDS blocks, so the gene name alone is not unique. The start coordinate makes the record
            # id unique without hiding that these are genuinely separate copies.
            out.append({"id": f"{v['event_id']}|{v['virus']}|{gene}@{start1}", "source": "virus",
                        "gene": f"{v['virus']}_{gene}", "sequence": seq, "tpm": share,
                        "clone": v.get("clone", "T"), "hap": 0,
                        "in_normal": bool(trace > 0),
                        "normal_tpm": (0.5 / len(cds)) if trace > 0 else 0.0})
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
