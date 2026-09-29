"""Derived chromosome sequences for a clone and haplotype.

Junction contigs are enough for short reads, which only ever see one breakpoint at a time. A long read
crosses a whole rearrangement, so PacBio libraries need the derived sequence itself: a deletion actually
missing, a duplication actually repeated, an inversion actually reverse-complemented, and a translocation
actually continuing into its partner chromosome.

A derived chromosome is planned as an ordered list of reference segments with orientation, then realised
by fetching each segment with the clone's germline and somatic edits already applied. The plan doubles as
the coordinate map from derived position back to reference position, which is what lets the truth bundle
say where a simulated read really came from.
"""
from dataclasses import dataclass


@dataclass
class Segment:
    """One piece of a derived chromosome, in reference coordinates."""
    chrom: str
    start: int          # 1-based inclusive
    end: int            # 1-based inclusive
    strand: str = "+"   # '-' means the piece is reverse complemented
    kind: str = "ref"   # ref, dup, inv, insert
    event: str = ""     # the designed event that produced it

    def length(self):
        return max(0, self.end - self.start + 1)


def _usable(svs, chrom_len, margin=1000):
    """Designed SVs on one chromosome, sorted and stripped of overlaps.

    Two rearrangements whose spans overlap cannot both be laid out on one derived chromosome, so the
    later one is dropped rather than producing a sequence that silently contradicts the truth table.
    """
    rows = []
    for s in svs:
        try:
            start, end = int(s["start"]), int(s["end"])
        except (KeyError, TypeError, ValueError):
            continue
        if end < start:
            start, end = end, start
        if start < margin or end > chrom_len - margin:
            continue
        rows.append((start, end, s))
    rows.sort()
    kept, last_end = [], 0
    dropped = []
    for start, end, s in rows:
        if start <= last_end:
            dropped.append(s.get("event_id", ""))
            continue
        kept.append((start, end, s))
        last_end = end
    return kept, dropped


def plan_chromosome(chrom, chrom_len, svs):
    """Ordered segments composing the derived chromosome, plus the events applied and dropped.

    Intra-chromosomal events only. Translocations are handled by `plan_translocation`, since they produce
    derivative chromosomes rather than modifying one in place.
    """
    kept, dropped = _usable([s for s in svs if s.get("svtype") in ("DEL", "DUP", "INV", "INS")], chrom_len)
    segs = []
    applied = []
    cursor = 1
    for start, end, s in kept:
        t = s["svtype"]
        if start > cursor:
            segs.append(Segment(chrom, cursor, start - 1))
        if t == "DEL":
            pass                                   # the span is simply absent
        elif t == "DUP":
            segs.append(Segment(chrom, start, end, "+", "ref", s.get("event_id", "")))
            segs.append(Segment(chrom, start, end, "+", "dup", s.get("event_id", "")))
        elif t == "INV":
            segs.append(Segment(chrom, start, end, "-", "inv", s.get("event_id", "")))
        elif t == "INS":
            segs.append(Segment(chrom, start, end, "+", "ref", s.get("event_id", "")))
            segs.append(Segment(chrom, start, start, "+", "insert", s.get("event_id", "")))
        applied.append(s.get("event_id", ""))
        cursor = end + 1
    if cursor <= chrom_len:
        segs.append(Segment(chrom, cursor, chrom_len))
    return segs, applied, dropped


def plan_translocation(chrom_a, len_a, pos_a, chrom_b, len_b, pos_b, event=""):
    """The two derivative chromosomes of a reciprocal translocation.

    der(A) is A up to its breakpoint joined to B beyond its own; der(B) is the reciprocal.
    """
    der_a = [Segment(chrom_a, 1, pos_a, "+", "ref", event),
             Segment(chrom_b, pos_b + 1, len_b, "+", "ref", event)]
    der_b = [Segment(chrom_b, 1, pos_b, "+", "ref", event),
             Segment(chrom_a, pos_a + 1, len_a, "+", "ref", event)]
    return der_a, der_b


def coordinate_map(segs):
    """(derived_start, derived_end, chrom, ref_start, ref_end, strand, kind, event) per segment.

    Derived coordinates are 1-based inclusive and run over the concatenated sequence.
    """
    out = []
    pos = 1
    for s in segs:
        n = s.length()
        if n <= 0:
            continue
        out.append((pos, pos + n - 1, s.chrom, s.start, s.end, s.strand, s.kind, s.event))
        pos += n
    return out


def derived_length(segs):
    return sum(s.length() for s in segs)


def realise(genome, germline, events_by_chrom, clones, clone, hap, segs,
            mei_sequence=None, only_pre_cna=False, chunk=2_000_000):
    """Build the derived sequence for a plan, with germline and somatic edits applied.

    Segments are fetched in chunks so a whole chromosome never has to be held twice, and an inverted
    segment is reverse complemented after its edits are applied, not before, so edit coordinates stay in
    reference space throughout.
    """
    from .genome_build import haplotype_sequence
    from .genome import revcomp

    pieces = []
    for s in segs:
        if s.kind == "insert":
            pieces.append(mei_sequence(s.event) if mei_sequence else "")
            continue
        parts = []
        pos = s.start
        while pos <= s.end:
            hi = min(s.end, pos + chunk - 1)
            seq, _stats = haplotype_sequence(genome, germline, events_by_chrom.get(s.chrom, []),
                                             clones, s.chrom, hap, clone, pos, hi,
                                             strict=False, only_pre_cna=only_pre_cna)
            parts.append(seq)
            pos = hi + 1
        piece = "".join(parts)
        pieces.append(revcomp(piece) if s.strand == "-" else piece)
    return "".join(pieces)


def write_derived(path, name, sequence, line=60):
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(sequence), line):
            fh.write(sequence[i:i + line] + "\n")
    return path


def write_coordinate_map(path, rows, header=True):
    with open(path, "w") as fh:
        if header:
            fh.write("derived_start\tderived_end\tchrom\tref_start\tref_end\tstrand\tkind\tevent\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    return path
