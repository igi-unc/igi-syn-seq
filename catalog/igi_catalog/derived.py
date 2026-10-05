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



MIN_SLICE_BP = 20_000


def cn_slices(cmap, boundaries, min_bp=MIN_SLICE_BP):
    """Cut a derived chromosome into pieces of uniform reference copy number.

    `cmap` is `rearrange.coordinate_map(segs)`: (derived_start, derived_end, chrom, ref_start, ref_end,
    strand, kind, event), 1-based inclusive. `boundaries` are 1-based reference coordinates where copy
    number can change, from `Clones.cn_boundaries`.

    Returns [(derived_start, derived_end, ref_probe)], where `ref_probe` is a reference position inside
    the piece that `source_plan` can be asked about.

    Why this exists: the long-read builders evaluated `source_plan` ONCE per chromosome, at
    `chrom_len // 2`, and simulated the whole derived chromosome at that single coverage. So a copy-number
    event was represented only if it happened to span the chromosome midpoint, and none of the designed
    focal events do. PTEN_homdel, RB1_homdel and MYC_amp were absent from PacBio and ONT WGS entirely --
    not under-represented, absent -- because every base of the chromosome was given the midpoint's copy
    number. That is a stronger version of the 5 Mb window defect on the Illumina side.

    Pieces shorter than `min_bp` are folded into the previous one rather than emitted: badread needs a
    reference comfortably longer than a read, and a 2 kb slice at 30x is a rounding error against a 0.3 Mb
    event. Inserted sequence with no reference span (an MEI, a viral integration) inherits the copy number
    of the segment it sits in, which is what it would have in a real genome.
    """
    cuts = sorted(set(boundaries))
    out = []
    for d_start, d_end, _chrom, r_start, r_end, strand, _kind, _event in cmap:
        n = d_end - d_start + 1
        if r_start is None or r_end is None or r_end < r_start:
            # no reference span: inherit the preceding piece, or probe the following segment
            if out:
                out[-1] = (out[-1][0], d_end, out[-1][2])
            else:
                out.append((d_start, d_end, r_start or 1))
            continue
        lo, hi = min(r_start, r_end), max(r_start, r_end)
        inner = [c for c in cuts if lo < c <= hi]
        marks = [lo] + inner + [hi + 1]
        for i in range(len(marks) - 1):
            a, b = marks[i], marks[i + 1] - 1
            # map the reference sub-span back onto derived coordinates, honouring strand
            off_a, off_b = a - lo, b - lo
            if strand == "-":
                off_a, off_b = hi - b, hi - a
            ds, de = d_start + off_a, d_start + off_b
            ds, de = max(d_start, ds), min(d_end, de)
            if de < ds:
                continue
            probe = (a + b) // 2
            if out and de - ds + 1 < min_bp and out[-1][2] == probe:
                out[-1] = (out[-1][0], de, out[-1][2])
            elif out and de - ds + 1 < min_bp:
                out[-1] = (out[-1][0], de, out[-1][2])
            else:
                out.append((ds, de, probe))
        assert n >= 0
    return out


def slice_weight(env, purity, chrom, probe, tumor, src, hap, kind):
    """The weight `source_plan` gives this one (source, haplotype, copy kind) at a reference position.

    Zero -- the key absent from the plan -- means this copy does not exist there, which is what a deletion
    IS. The slice then gets no reads, and that is the signal a caller is meant to see.
    """
    for s, h, k, w in source_plan(env.clones, purity, chrom, probe, tumor=tumor):
        if (s, h, k) == (src, hap, kind):
            return w
    return 0.0


def build_derived(env, dcfg, dataset, chrom, library, catalog_dir, events, depth, rng, work,
                  log=print):
    """Write one FASTA per derived chromosome copy. Returns [(tag, fasta_path, coverage)] and plan rows."""
    purity = dcfg["purity"]
    chrom_len = env.genome.lengths[chrom]
    denom = depth_denominator(env.clones, purity)
    tumor = library == "tumor"
    svs = sv_rows(catalog_dir, dataset, chrom) if tumor else []
    pieces, plan_rows = [], []
    boundaries = [c for c in env.clones.cn_boundaries(chrom) if 1 < c <= chrom_len] if tumor else []
    # The plan is still enumerated at the midpoint, but only to discover WHICH copies exist on this
    # chromosome. Each copy's coverage is then computed per copy-number-uniform slice, not once.
    for src, hap, kind, _w_mid in source_plan(env.clones, purity, chrom, chrom_len // 2, tumor=tumor):
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
        cmap = rearrange.coordinate_map(segs)
        cmap_path = rearrange.write_coordinate_map(
            os.path.join(work, f"{tag}.coordmap.tsv"),
            [(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]) for r in cmap])
        slices = cn_slices(cmap, boundaries) if boundaries else [(1, len(seq), chrom_len // 2)]
        sub_rows, n_emitted, n_empty = [], 0, 0
        for si, (ds, de, probe) in enumerate(slices):
            w = slice_weight(env, purity, chrom, probe, tumor, src, hap, kind)
            if w <= 0:
                # this copy is absent here: a deletion, and the absence of reads is the point
                n_empty += 1
                sub_rows.append({"slice": si, "derived_start": ds, "derived_end": de,
                                 "ref_probe": probe, "weight": 0.0, "coverage": 0.0})
                continue
            piece_seq = seq[ds - 1:de]
            if len(piece_seq) < 10_000:
                continue
            stag = tag if len(slices) == 1 else f"{tag}_s{si:03d}"
            fa = rearrange.write_derived(os.path.join(work, f"{stag}.fa"),
                                         f"{chrom}_{stag}", piece_seq)
            cov = max(0.01, depth * w / denom)
            pieces.append((stag, fa, cov))
            n_emitted += 1
            sub_rows.append({"slice": si, "derived_start": ds, "derived_end": de,
                             "ref_probe": probe, "weight": round(w, 5), "coverage": round(cov, 4),
                             "length": len(piece_seq), "fasta": fa})
        plan_rows.append({"source": src, "haplotype": hap, "copy_kind": kind,
                          "derived_length": len(seq), "reference_length": chrom_len,
                          "segments": len(segs), "svs_applied": applied,
                          "svs_offered": [s["event_id"] for s in mine],
                          "svs_dropped": [(d[-1].get("event_id") if isinstance(d, tuple) else str(d))
                                          for d in dropped],
                          "coordinate_map": cmap_path, "cn_slices": sub_rows})
        covs = [r["coverage"] for r in sub_rows if r["coverage"]]
        log(f"    {tag}: {len(segs)} segments, derived {len(seq):,} bp "
            f"(reference {chrom_len:,}), {len(applied)}/{len(mine)} SVs applied; "
            f"{len(slices)} cn slice(s), {n_emitted} simulated, {n_empty} with no copies, "
            f"coverage {min(covs) if covs else 0:.3f}-{max(covs) if covs else 0:.3f}x")
    return pieces, plan_rows
