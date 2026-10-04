"""Per-clone, per-haplotype derived chromosomes: the input every long-read WGS assay simulates from.

Lifted out of run_pacbio_wgs.py so the ONT WGS arm uses the same construction rather than a second
implementation of it. The logic is unchanged, including the parts that are load-bearing:

- A copy is built per (source, haplotype, copy kind) from `source_plan`, and weighted by that source's
  share of sequenced molecules, so allele fractions follow the clone model by construction.
- `plan_chromosome` returns applied AND dropped SVs. The dropped list is truth, not noise: an SV that
  could not be placed -- overlapping another, or too close to a chromosome end -- is genuinely absent from
  the reads, and a truth table that still claimed it would be wrong.
- `pre_cna` copies are realised with only pre-CNA events, which is what makes a variant that predates an
  amplification appear on every copy and one that postdates it appear on one.

run_pacbio_wgs.py still carries its own copy of this loop. It has already produced a validated 96-chromosome
release, so it was left untouched rather than refactored mid-flight; unifying the two is worth doing with an
equivalence test, the way rna_assembly was.
"""
import csv
import os

from . import rearrange
from .pipeline import depth_denominator, source_plan


def sv_rows(catalog_dir, dataset, chrom):
    """Designed structural variants on this chromosome, as rearrange.plan_chromosome expects them."""
    path = os.path.join(catalog_dir, f"{dataset}.svs.tsv")
    if not os.path.exists(path):
        return []
    out = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("chrom") != chrom or r.get("svtype") == "TRA":
                continue        # translocations need both partners; handled separately
            try:
                out.append({"event_id": r["event_id"], "svtype": r["svtype"],
                            "start": int(r["start"]), "end": int(r["end"]),
                            "clone": r.get("clone") or "T",
                            "haplotype": int(r.get("haplotype") or 0),
                            "inserted_sequence": r.get("inserted_sequence") or ""})
            except (ValueError, KeyError):
                continue
    return out


def mei_sequence_factory(rng):
    """Sequence for a mobile-element insertion. The catalog records a family and a length, not bases."""
    def make(event):
        return "".join(rng.choice("ACGT") for _ in range(300))
    return make


def build_derived(env, dcfg, dataset, chrom, library, catalog_dir, events, depth, rng, work,
                  log=print):
    """Write one FASTA per derived chromosome copy. Returns [(tag, fasta_path, coverage)] and plan rows."""
    purity = dcfg["purity"]
    chrom_len = env.genome.lengths[chrom]
    denom = depth_denominator(env.clones, purity)
    tumor = library == "tumor"
    svs = sv_rows(catalog_dir, dataset, chrom) if tumor else []
    pieces, plan_rows = [], []
    for src, hap, kind, w in source_plan(env.clones, purity, chrom, chrom_len // 2, tumor=tumor):
        clone = "T" if (src == "NORMAL" or kind == "pre_cna") else src
        pre_only = kind == "pre_cna"
        mine = ([s for s in svs
                 if s["haplotype"] == hap and env.clones.is_descendant(clone, s["clone"])]
                if src != "NORMAL" else [])
        segs, applied, dropped = rearrange.plan_chromosome(chrom, chrom_len, mine)
        if dropped:
            log(f"    {src} hap{hap} {kind}: {len(dropped)} SV(s) not placeable")
        seq = rearrange.realise(env.genome, env.germline, events, env.clones, clone, hap, segs,
                                mei_sequence=mei_sequence_factory(rng), only_pre_cna=pre_only)
        if len(seq) < 10_000:
            continue
        tag = f"{src}_hap{hap}_{kind}"
        fa = rearrange.write_derived(os.path.join(work, f"{tag}.fa"), f"{chrom}_{tag}", seq)
        cov = max(0.01, depth * w / denom)
        plan_rows.append({"source": src, "haplotype": hap, "copy_kind": kind, "weight": round(w, 5),
                          "coverage": round(cov, 4), "derived_length": len(seq),
                          "reference_length": chrom_len, "segments": len(segs),
                          "svs_applied": applied,
                          "svs_dropped": [(d[-1].get("event_id") if isinstance(d, tuple) else str(d))
                                          for d in dropped],
                          "fasta": fa})
        pieces.append((tag, fa, cov))
        log(f"    {tag}: {len(segs)} segments, derived {len(seq):,} bp "
            f"(reference {chrom_len:,}), {len(applied)}/{len(mine)} SVs applied, coverage {cov:.3f}x")
    return pieces, plan_rows
