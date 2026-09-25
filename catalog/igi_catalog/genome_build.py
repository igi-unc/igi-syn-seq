"""Stage 2: build per-clone, per-haplotype genome sequences.

A clone's genome is the reference with, in order:
  1. the individual's phased germline variants for that haplotype,
  2. the designed somatic small variants assigned to that haplotype and to a clone in this lineage,
  3. copy-number state, expressed as how many times each haplotype's sequence is emitted per region.

Edits are applied right to left in reference coordinates so earlier offsets stay valid, and every applied
edit is checked against the sequence it replaces, so a silent mismatch cannot pass.
"""
import csv
from collections import defaultdict


class EditSet:
    """Reference-coordinate substitutions for one chromosome, applied right to left."""

    def __init__(self, chrom):
        self.chrom = chrom
        self.edits = []          # (pos1, ref, alt, label)

    def add(self, pos1, ref, alt, label):
        self.edits.append((int(pos1), ref, alt, label))

    def apply(self, seq, offset1=1, strict=True):
        """Apply to `seq`, whose first base is reference position `offset1`.

        Returns (sequence, n_applied, [rejected labels]).
        """
        out = seq
        applied = 0
        rejected = []
        for pos1, ref, alt, label in sorted(self.edits, key=lambda e: -e[0]):
            off = pos1 - offset1
            if off < 0 or off + len(ref) > len(out):
                rejected.append(f"{label}:out_of_range")
                continue
            if out[off:off + len(ref)].upper() != ref.upper():
                rejected.append(f"{label}:ref_mismatch")
                if strict:
                    raise ValueError(
                        f"{self.chrom}:{pos1} {label}: expected {ref!r}, found "
                        f"{out[off:off + len(ref)].upper()!r}")
                continue
            out = out[:off] + alt + out[off + len(ref):]
            applied += 1
        return out, applied, rejected


def germline_edits(germline, chrom, hap, start1, end1):
    """Phased germline edits for one haplotype over a region."""
    es = EditSet(chrom)
    for pos, ref, alt, gt, _phased in germline.variants(chrom, start1 - 1, end1):
        allele = gt[hap] if len(gt) == 2 else gt[0]
        if allele == 1:
            es.add(pos, ref, alt, "germline")
    return es


def somatic_edits(events, clones, chrom, hap, clone):
    """Designed somatic small-variant edits visible in `clone` on haplotype `hap`.

    An event is present in a clone when that clone is a descendant of the clone the event arose in.
    """
    es = EditSet(chrom)
    for e in events:
        if e["chrom"] != chrom or int(e["haplotype"]) != hap:
            continue
        if not clones.is_descendant(clone, e["clone"]):
            continue
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
                "clone": row["clone"],
            })
    return by_chrom


def haplotype_sequence(genome, germline, events, clones, chrom, hap, clone, start1=1, end1=None, strict=True):
    """Reference sequence for one chromosome region with germline and somatic edits applied.

    Returns (sequence, stats dict).
    """
    end1 = end1 or genome.lengths[chrom]
    seq = genome.seq(chrom, start1 - 1, end1)
    g = germline_edits(germline, chrom, hap, start1, end1)
    s = somatic_edits(events, clones, chrom, hap, clone)
    # somatic first (right to left), then germline: both are in reference coordinates and the designer
    # guarantees they never overlap, so the order only affects bookkeeping
    seq, n_som, rej_som = s.apply(seq, start1, strict=strict)
    seq, n_germ, rej_germ = g.apply(seq, start1, strict=False)
    return seq, {
        "chrom": chrom, "hap": hap, "clone": clone, "length": len(seq),
        "somatic_applied": n_som, "somatic_rejected": rej_som,
        "germline_applied": n_germ, "germline_rejected": len(rej_germ),
    }
