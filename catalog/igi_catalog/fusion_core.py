"""Fusion transcript construction: exon-level joins, frame arithmetic, and junction peptide windows."""
from .annotation import translate
from .genome import revcomp


def cds_segments(t):
    """CDS segments in transcript 5'->3' order."""
    return t.cds if t.strand == "+" else t.cds[::-1]


def coding_exon_bounds(t):
    """For each CDS segment, the (start, end) of the exon containing it, in transcript order."""
    out = []
    for s, e in cds_segments(t):
        ex = next(((a, b) for a, b in t.exons if a <= s and e <= b), (s, e))
        out.append(ex)
    return out


def cds_prefix_len(t, n_segs):
    """Coding length contributed by the first n_segs CDS segments."""
    segs = cds_segments(t)
    return sum(e - s + 1 for s, e in segs[:n_segs])


def segment_seq(genome, t, s, e):
    seq = genome.seq(t.chrom, s - 1, e)
    return seq if t.strand == "+" else revcomp(seq)


def five_prime_cds(genome, t, n_segs):
    """Coding sequence of the first n_segs CDS segments of the 5' partner."""
    return "".join(segment_seq(genome, t, s, e) for s, e in cds_segments(t)[:n_segs])


def three_prime_cds(genome, t, from_seg):
    """Coding sequence from CDS segment index `from_seg` to the end of the 3' partner."""
    return "".join(segment_seq(genome, t, s, e) for s, e in cds_segments(t)[from_seg:])


def intron_breakpoint(t, seg_index, side, rng, near_exon=False):
    """A genomic position in the intron flanking a CDS segment.

    side='after': in the intron following segment seg_index (5' partner's break).
    side='before': in the intron preceding segment seg_index (3' partner's break).
    near_exon: place the break within 20-45 bp of the exon boundary, so it falls inside a capture
    region's flank and the fusion is detectable in exome data as well as in WGS.
    Returns (chrom, pos1) or None when the required intron does not exist.
    """
    ex = coding_exon_bounds(t)
    if t.strand == "+":
        if side == "after":
            if seg_index >= len(ex) - 1:
                return None
            lo, hi = ex[seg_index][1] + 1, ex[seg_index + 1][0] - 1
        else:
            if seg_index <= 0:
                return None
            lo, hi = ex[seg_index - 1][1] + 1, ex[seg_index][0] - 1
    else:
        if side == "after":
            if seg_index >= len(ex) - 1:
                return None
            lo, hi = ex[seg_index + 1][1] + 1, ex[seg_index][0] - 1
        else:
            if seg_index <= 0:
                return None
            lo, hi = ex[seg_index][1] + 1, ex[seg_index - 1][0] - 1
    if hi - lo < 40:
        return None
    if near_exon:
        # the intron boundary adjacent to the retained exon: lo abuts it for 'after', hi for 'before'
        if (t.strand == "+") == (side == "after"):
            a, b = lo, min(hi, lo + 45)
        else:
            a, b = max(lo, hi - 45), hi
        if b - a >= 25:
            return t.chrom, rng.randrange(a + 20, b)
    return t.chrom, rng.randrange(lo + 20, hi - 20)


def build_fusion(genome, t5, n5, t3, from3):
    """Join the first n5 CDS segments of t5 to CDS segments from3.. of t3.

    Returns a dict with the fusion CDS, protein, whether the join is in frame, the junction
    codon index, and the peptide window spanning the junction.
    """
    left = five_prime_cds(genome, t5, n5)
    right = three_prime_cds(genome, t3, from3)
    in_frame = (len(left) % 3) == 0
    fusion_cds = left + right
    prot = translate(fusion_cds)
    stop = prot.find("*")
    prot_trunc = prot[:stop] if stop >= 0 else prot
    j_aa = len(left) // 3
    window = prot_trunc[max(0, j_aa - 10):j_aa + 11]
    return {
        "cds_5p_len": len(left),
        "cds_len": len(fusion_cds),
        "in_frame": in_frame,
        "junction_codon": j_aa,
        "protein": prot_trunc,
        "protein_len": len(prot_trunc),
        "junction_window": window,
        "stop_before_junction": stop >= 0 and stop < j_aa,
    }


def junction_diff(genome, t5, n5, t3, from3, flank=21):
    """Text lines showing the mRNA sequence either side of the fusion junction against each parent."""
    left = five_prime_cds(genome, t5, n5)
    right = three_prime_cds(genome, t3, from3)
    l_tail = left[-flank:] if len(left) >= flank else left
    r_head = right[:flank]
    # parent continuations for contrast
    segs5 = cds_segments(t5)
    nxt5 = segment_seq(genome, t5, *segs5[n5])[:flank] if n5 < len(segs5) else ""
    prv3 = ""
    segs3 = cds_segments(t3)
    if from3 > 0:
        prv3 = segment_seq(genome, t3, *segs3[from3 - 1])[-flank:]
    return [
        f"  5' parent  ...{l_tail} | {nxt5}...   ({t5.gene_name} continues)",
        f"  3' parent  ...{prv3} | {r_head}...   ({t3.gene_name} continues)",
        f"  fusion     ...{l_tail} | {r_head}...",
        f"             {' ' * (len(l_tail) + 3)}^ junction",
    ]
