"""Stage 2: build per-clone, per-haplotype genome sequences.

A clone's genome is the reference with, in order:
  1. the individual's phased germline variants for that haplotype,
  2. the designed somatic small variants assigned to that haplotype and to a clone in this lineage,
  3. copy-number state, expressed as how many times each haplotype's sequence is emitted per region.

Edits are applied in one left-to-right pass over reference coordinates, and every applied edit is
checked against the sequence it replaces, so a silent mismatch cannot pass.
"""
import csv
from collections import defaultdict


class EditSet:
    """Reference-coordinate substitutions for one chromosome."""

    def __init__(self, chrom):
        self.chrom = chrom
        self.edits = []          # (pos1, ref, alt, label)

    def add(self, pos1, ref, alt, label):
        self.edits.append((int(pos1), ref, alt, label))

    def apply(self, seq, offset1=1, strict=True):
        """Apply to `seq`, whose first base is reference position `offset1`.

        Single left-to-right pass: reference stretches and alternate alleles are appended to a list and
        joined once, so cost is linear in sequence length rather than quadratic in the number of edits.
        Returns (sequence, n_applied, [rejected labels]).
        """
        pieces = []
        cursor = 0
        applied = 0
        rejected = []
        last_end = -1
        for pos1, ref, alt, label in sorted(self.edits, key=lambda e: (e[0], len(e[1]))):
            off = pos1 - offset1
            if off < 0 or off + len(ref) > len(seq):
                rejected.append(f"{label}:out_of_range")
                continue
            if off < last_end:
                rejected.append(f"{label}:overlapping_edit")
                continue
            if seq[off:off + len(ref)].upper() != ref.upper():
                rejected.append(f"{label}:ref_mismatch")
                if strict:
                    raise ValueError(
                        f"{self.chrom}:{pos1} {label}: expected {ref!r}, found "
                        f"{seq[off:off + len(ref)].upper()!r}")
                continue
            pieces.append(seq[cursor:off])
            pieces.append(alt)
            cursor = off + len(ref)
            last_end = cursor
            applied += 1
        pieces.append(seq[cursor:])
        return "".join(pieces), applied, rejected


def germline_edits(germline, chrom, hap, start1, end1):
    """Phased germline edits for one haplotype over a region."""
    es = EditSet(chrom)
    for pos, ref, alt, gt, _phased in germline.variants(chrom, start1 - 1, end1):
        allele = gt[hap] if len(gt) == 2 else gt[0]
        if allele == 1:
            es.add(pos, ref, alt, "germline")
    return es


def somatic_edits(events, clones, chrom, hap, clone, start1=None, end1=None, only_pre_cna=False):
    """Designed somatic small-variant edits visible in `clone` on haplotype `hap`.

    An event is present in a clone when that clone is a descendant of the clone the event arose in.
    """
    es = EditSet(chrom)
    for e in events:
        if e["chrom"] != chrom or int(e["haplotype"]) != hap:
            continue
        if not clones.is_descendant(clone, e["clone"]):
            continue
        if start1 is not None and (e["pos"] < start1 or e["pos"] + len(e["ref"]) - 1 > end1):
            continue          # outside the region being built, not a failure
        if only_pre_cna and e.get("timing") != "pre_cna":
            continue          # this copy predates the event, so it does not carry it
        es.add(e["pos"], e["ref"], e["alt"], e["event_id"])
    return es


def read_events(tsv):
    """Designed SNV/indel events from a catalog table, keyed by chromosome."""
    by_chrom = defaultdict(list)
    with open(tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("class") not in ("snv", "indel"):
                continue
            by_chrom[row["chrom"]].append({
                "event_id": row["event_id"], "chrom": row["chrom"], "pos": int(row["pos"]),
                "ref": row["ref"], "alt": row["alt"], "haplotype": int(row["haplotype"]),
                "clone": row["clone"], "timing": row.get("timing", "pre_cna"),
            })
    return by_chrom


def haplotype_sequence(genome, germline, events, clones, chrom, hap, clone, start1=1, end1=None, strict=True,
                       only_pre_cna=False):
    """Reference sequence for one chromosome region with germline and somatic edits applied.

    Returns (sequence, stats dict).
    """
    end1 = end1 or genome.lengths[chrom]
    seq = genome.seq(chrom, start1 - 1, end1)
    g = germline_edits(germline, chrom, hap, start1, end1)
    s = somatic_edits(events, clones, chrom, hap, clone, start1, end1, only_pre_cna)
    # germline and somatic edits are both in reference coordinates, so they must be applied in a single
    # pass; applying one set first would shift the coordinates the other set refers to
    merged = EditSet(chrom)
    merged.edits = s.edits + g.edits
    somatic_ids = {e[3] for e in s.edits}
    seq, n_applied, rejected = merged.apply(seq, start1, strict=False)
    rej_som = [r for r in rejected if r.split(":")[0] in somatic_ids]
    if rej_som and strict:
        raise ValueError(f"{chrom}: designed somatic edits could not be applied: {rej_som[:5]}")
    n_som = len(somatic_ids) - len(rej_som)
    return seq, {
        "chrom": chrom, "hap": hap, "clone": clone, "length": len(seq),
        "somatic_applied": n_som, "somatic_rejected": rej_som,
        "germline_applied": n_applied - n_som,
        "germline_rejected": len(rejected) - len(rej_som),
        "reject_reasons": _reasons(rejected),
    }


def _reasons(rejected):
    out = {}
    for r in rejected:
        out[r.rsplit(":", 1)[-1]] = out.get(r.rsplit(":", 1)[-1], 0) + 1
    return out
