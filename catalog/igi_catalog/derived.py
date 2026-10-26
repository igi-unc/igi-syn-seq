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

from . import mei, rearrange
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


def tra_rows(catalog_dir, dataset):
    """Every designed translocation, with both endpoints. Not filtered by chromosome.

    `sv_rows` skips TRA because `plan_chromosome` lays out one chromosome and a translocation spans two.
    They were then skipped by everything: `plan_translocation` had no caller anywhere in the repository,
    so the derived chromosomes -- and therefore both long-read libraries -- contained no trace of any
    designed translocation, while the short-read path carried them as junction contigs the whole time.

    Both endpoints are returned because both derivatives must be built, and because whether a given copy
    can carry a translocation at all depends on what else that copy carries. See `placeable_tra`.
    """
    path = os.path.join(catalog_dir, f"{dataset}.svs.tsv")
    if not os.path.exists(path):
        return []
    out = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("svtype") != "TRA":
                continue
            try:
                a, b = r.get("chrom"), r.get("end_chrom") or ""
                if not a or not b or a == b:
                    continue
                out.append({"event_id": r["event_id"], "chrom": a, "start": int(r["start"]),
                            "end_chrom": b, "end": int(r["end"]),
                            "clone": r.get("clone") or "T",
                            "haplotype": int(r.get("haplotype") or 0)})
            except (ValueError, KeyError):
                continue
    out.sort(key=lambda t: (t["chrom"], t["start"]))
    return out


def placeable_tra(my_tras, chrom):
    """The translocation whose derivative begins on `chrom` for one copy, and the ones that cannot.

    Returns (chosen, cut, partner, partner_start, skipped).

    A copy may carry translocations that chain -- ds-02 has three meeting on chr17 hap0 of clone T, two
    with chr17 as the 5' partner and one with chr17 as the 3'. A chain makes derivatives that are not
    simply "this chromosome cut once and continued into one partner", and building one cut anyway would
    double-count the bases the other cut moved elsewhere.

    So a translocation is placed only when BOTH its chromosomes carry exactly one translocation on this
    copy. That is the condition under which the two owning tasks independently reach the same two
    derivatives, and every base of both chromosomes appears exactly once across them. Everything else is
    returned as skipped and logged; those events are still present in the short-read libraries as junction
    contigs, which is where they always were.
    """
    n = {}
    for t in my_tras:
        n[t["chrom"]] = n.get(t["chrom"], 0) + 1
        n[t["end_chrom"]] = n.get(t["end_chrom"], 0) + 1
    chosen = cut = partner = partner_start = None
    skipped = []
    for t in my_tras:
        if chrom not in (t["chrom"], t["end_chrom"]):
            continue
        if n[t["chrom"]] != 1 or n[t["end_chrom"]] != 1:
            skipped.append(t)
            continue
        if chosen is None:
            chosen = t
            if t["chrom"] == chrom:
                cut, partner, partner_start = t["start"], t["end_chrom"], t["end"] + 1
            else:
                cut, partner, partner_start = t["end"], t["chrom"], t["start"] + 1
    return chosen, cut, partner, partner_start, skipped


def mei_sequence_factory(genome, svs, dataset):
    """Bases for each designed insertion, looked up by the event id the plan carries.

    `realise` hands this an event id, not a row, so the rows are indexed here.

    Each event draws from its OWN seeded stream. That is not tidiness: a truncal insertion is carried by
    every descendant copy, and with one shared stream each copy drew a different sequence for the same
    event -- so the clones disagreed about what the truth table called a single insertion, and the reads
    supporting it never stacked. Seeding per event makes a shared site produce identical bases on every
    copy that carries it, whatever else that copy carries.
    """
    import random
    by_id = {s["event_id"]: s for s in svs if s.get("event_id")}

    def make(event):
        row = by_id.get(event)
        if not row:
            return ""
        erng = random.Random(f"{dataset}:{event}:mei")
        try:
            pos = int(row.get("start") or 0)
        except (TypeError, ValueError):
            return ""
        return mei.insertion(genome, row.get("inserted_sequence", ""), row.get("chrom", ""), pos, erng)

    return make



MIN_SLICE_BP = 20_000


def cn_slices(cmap, cuts_for, min_bp=MIN_SLICE_BP):
    """Cut a derived chromosome into pieces of uniform reference copy number.

    `cmap` is `rearrange.coordinate_map(segs)`: (derived_start, derived_end, chrom, ref_start, ref_end,
    strand, kind, event), 1-based inclusive. `cuts_for(chrom)` returns the 1-based reference coordinates
    where copy number can change on that chromosome, from `Clones.cn_boundaries`.

    Returns [(derived_start, derived_end, ref_chrom, ref_probe)], where `ref_probe` is a reference
    position inside the piece that `source_plan` can be asked about on `ref_chrom`.

    The chromosome travels with the slice because a derived chromosome is no longer guaranteed to be made
    of one chromosome's segments: a translocation continues into its partner, and that tail must be
    weighted by the partner's copy number at the partner's coordinates. Asking this chromosome's clone
    model about a position on another chromosome returns a plausible, wrong answer.

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
    cache = {}

    def cuts(chrom):
        if chrom not in cache:
            cache[chrom] = sorted(set(cuts_for(chrom)))
        return cache[chrom]

    out = []
    for d_start, d_end, chrom, r_start, r_end, strand, _kind, _event in cmap:
        n = d_end - d_start + 1
        if _kind == "insert" or r_start is None or r_end is None or r_end < r_start:
            # no reference span: inherit the preceding piece, or probe the following segment
            if out:
                out[-1] = (out[-1][0], d_end, out[-1][2], out[-1][3])
            else:
                out.append((d_start, d_end, chrom, r_start or 1))
            continue
        lo, hi = min(r_start, r_end), max(r_start, r_end)
        inner = [c for c in cuts(chrom) if lo < c <= hi]
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
            # a short piece folds into the previous one only when that piece is on the same chromosome;
            # folding across a translocation breakpoint would hand the partner's bases this side's weight
            if out and de - ds + 1 < min_bp and out[-1][2] == chrom:
                out[-1] = (out[-1][0], de, out[-1][2], out[-1][3])
            else:
                out.append((ds, de, chrom, probe))
        assert n >= 0

    # Coalesce neighbours that would be simulated identically. `source_plan` depends only on the
    # chromosome and which inter-cut band the position falls in, so two adjacent pieces in the same band
    # of the same chromosome get the same weight and the same coverage; emitting them separately costs a
    # FASTA and a simulator invocation each and changes nothing. Without this, every segment boundary --
    # and chr5's chromothripsis region alone has 60 -- became its own slice on a chromosome that has no
    # copy-number change at all, where the builder used to emit exactly one.
    import bisect

    def band(chrom, probe):
        return bisect.bisect_right(cuts(chrom), probe)

    merged = []
    for ds, de, chrom, probe in out:
        if merged and merged[-1][2] == chrom and band(chrom, merged[-1][3]) == band(chrom, probe) \
                and merged[-1][1] + 1 == ds:
            merged[-1] = (merged[-1][0], de, chrom, merged[-1][3])
        else:
            merged.append((ds, de, chrom, probe))
    return merged


def slice_weight(env, purity, chrom, probe, tumor, src, hap, kind):
    """The weight `source_plan` gives this one (source, haplotype, copy kind) at a reference position.

    Zero -- the key absent from the plan -- means this copy does not exist there, which is what a deletion
    IS. The slice then gets no reads, and that is the signal a caller is meant to see.
    """
    for s, h, k, w in source_plan(env.clones, purity, chrom, probe, tumor=tumor):
        if (s, h, k) == (src, hap, kind):
            return w
    return 0.0


def viral_pieces(viral_fasta, catalog_dir, dataset, chrom, library, env, purity, depth, work,
                 log=print):
    """Viral genomes as extra simulation pieces for a long-read library.

    The long-read path builds from derived chromosomes and never read `viruses.tsv` at all, so PacBio
    HiFi and ONT WGS carried no viral sequence of any kind. Pinned to one chromosome task for the same
    reason the short-read path pins it: a virus has no human coordinate, so nothing else emits it once.

    This gives the viral genomes their designed dosage. It does NOT put the integrated virus inside the
    derived chromosome, so a long read still never crosses a host-virus boundary; integration-site
    evidence remains short-read only. Doing that properly means an insert segment carrying the viral
    sequence, which is the same machinery an MEI needs and is best done with it.
    """
    from .junctions import VIRAL_PAYLOAD_CHROM, viral_payloads
    if chrom != VIRAL_PAYLOAD_CHROM:
        return []
    path = os.path.join(catalog_dir, f"{dataset}.viruses.tsv")
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    tumor = library == "tumor"
    denom = depth_denominator(env.clones, purity, tumor=tumor)
    out = []
    for pl in viral_payloads(viral_fasta, rows, tumor=tumor):
        w = (purity * pl["ccf"] * pl["copies"]) if pl["tumor_weighted"] else pl["copies"]
        if w <= 0:
            continue
        tag = f"virus_{pl['slug']}"
        fa = rearrange.write_derived(os.path.join(work, f"{tag}.fa"), tag, pl["sequence"])
        cov = max(0.01, depth * w / denom)
        out.append((tag, fa, cov))
        log(f"    viral payload: {pl['desc']}, {len(pl['sequence']):,} bp at {cov:.2f}x")
    return out


def build_derived(env, dcfg, dataset, chrom, library, catalog_dir, events, depth, rng, work,
                  log=print):
    """Write one FASTA per derived chromosome copy. Returns [(tag, fasta_path, coverage)] and plan rows."""
    purity = dcfg["purity"]
    chrom_len = env.genome.lengths[chrom]
    tumor = library == "tumor"
    # tumor first: the denominator depends on it. A normal library is diploid and is normalised by 2.0
    # whatever the tumour's ploidy -- see depth_denominator. This is the ONT WGS and HiFi path.
    denom = depth_denominator(env.clones, purity, tumor=tumor)
    svs = sv_rows(catalog_dir, dataset, chrom) if tumor else []
    tras = tra_rows(catalog_dir, dataset) if tumor else []
    pieces, plan_rows = [], []

    def cuts_for(c):
        if not tumor:
            return []
        return [b for b in env.clones.cn_boundaries(c) if 1 < b <= env.genome.lengths[c]]
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
        # A translocation belongs to the copy that carries it, and only to the copy that carries
        # everything -- a pre-CNA copy predates it, exactly as it does for a junction contig on the
        # short-read path, where `kind != "all"` is excluded.
        my_tras = [t for t in tras
                   if t["haplotype"] == hap and env.clones.is_descendant(clone, t["clone"])] \
            if (not pre_only and src != "NORMAL") else []
        chosen, cut, partner, p_start, skipped = placeable_tra(my_tras, chrom)
        for t in skipped:
            log(f"    {src} hap{hap} {kind}: {t['event_id']} chains with another translocation on "
                f"this copy and is not placed in the derived chromosome")
        if chosen is not None and partner in env.genome.lengths:
            segs, ok = rearrange.plan_with_translocation(
                segs, chrom, cut, partner, p_start, env.genome.lengths[partner],
                chosen["event_id"])
            if ok:
                applied.append(chosen["event_id"])
                log(f"    {src} hap{hap} {kind}: der({chrom}) cut at {cut:,} and continued into "
                    f"{partner}:{p_start:,} ({chosen['event_id']})")
            else:
                log(f"    {src} hap{hap} {kind}: {chosen['event_id']} breakpoint falls inside a "
                    f"rearranged segment and is not placed")
        seq = rearrange.realise(env.genome, env.germline, events, env.clones, clone, hap, segs,
                                mei_sequence=mei_sequence_factory(env.genome, svs, dataset),
                                only_pre_cna=pre_only)
        if len(seq) < 10_000:
            continue
        tag = f"{src}_hap{hap}_{kind}"
        cmap = rearrange.coordinate_map(segs)
        cmap_path = rearrange.write_coordinate_map(
            os.path.join(work, f"{tag}.coordmap.tsv"),
            [(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]) for r in cmap])
        slices = (cn_slices(cmap, cuts_for) if tumor
                  else [(1, len(seq), chrom, chrom_len // 2)])
        sub_rows, n_emitted, n_empty = [], 0, 0
        for si, (ds, de, s_chrom, probe) in enumerate(slices):
            w = slice_weight(env, purity, s_chrom, probe, tumor, src, hap, kind)
            if w <= 0:
                # this copy is absent here: a deletion, and the absence of reads is the point
                n_empty += 1
                sub_rows.append({"slice": si, "derived_start": ds, "derived_end": de,
                                 "ref_chrom": s_chrom, "ref_probe": probe,
                                 "weight": 0.0, "coverage": 0.0})
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
                             "ref_chrom": s_chrom, "ref_probe": probe,
                             "weight": round(w, 5), "coverage": round(cov, 4),
                             "length": len(piece_seq), "fasta": fa})
        plan_rows.append({"source": src, "haplotype": hap, "copy_kind": kind,
                          "derived_length": len(seq), "reference_length": chrom_len,
                          "segments": len(segs), "svs_applied": applied,
                          "svs_offered": [s["event_id"] for s in mine],
                          "translocations_offered": [t["event_id"] for t in my_tras],
                          "translocations_skipped": [t["event_id"] for t in skipped],
                          "svs_dropped": [(d[-1].get("event_id") if isinstance(d, tuple) else str(d))
                                          for d in dropped],
                          "coordinate_map": cmap_path, "cn_slices": sub_rows})
        covs = [r["coverage"] for r in sub_rows if r["coverage"]]
        log(f"    {tag}: {len(segs)} segments, derived {len(seq):,} bp "
            f"(reference {chrom_len:,}), {len(applied)}/{len(mine)} SVs applied; "
            f"{len(slices)} cn slice(s), {n_emitted} simulated, {n_empty} with no copies, "
            f"coverage {min(covs) if covs else 0:.3f}-{max(covs) if covs else 0:.3f}x")
    return pieces, plan_rows
