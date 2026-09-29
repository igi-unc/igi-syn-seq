"""Junction contigs for structural variants, fusions and viral integrations.

A rearranged clone genome would need chromosome-level assembly. What a caller actually needs in order to
find these events is reads that cross the new joins, so each event contributes a short contig made of the
two flanks it brings together. Reads simulated from that contig align as split or discordant pairs at the
breakpoint, which is the evidence a structural-variant or fusion caller looks for.

The wild-type sequence across each breakpoint is reduced by the same weight on the same clone and
haplotype, so the fraction of reads supporting the junction matches the event's expected allele fraction
rather than being added on top of an unchanged background.

What this does not give: correct copy number across a rearranged genome. Deletions and duplications also
carry a depth change, which is applied separately through the copy-number profile of the affected
intervals. Anything requiring the rearranged sequence itself, such as a caller that assembles across a
whole inversion, is out of scope and documented as such.
"""
from .genome import revcomp

FLANK = 600


def _seq(genome, chrom, start1, end1):
    if chrom not in genome.lengths:
        return ""
    s = max(1, start1)
    e = min(end1, genome.lengths[chrom])
    return genome.seq(chrom, s - 1, e) if e > s else ""


def junction_sequence(genome, chrom_a, pos_a, chrom_b, pos_b, orient_a="+", orient_b="+", flank=FLANK):
    """Sequence spanning a join between two breakpoints.

    orient '+' takes the flank ending at the breakpoint; '-' takes the flank beginning at it, reverse
    complemented, which is how an inversion or a translocation to the opposite strand presents.
    """
    if orient_a == "+":
        left = _seq(genome, chrom_a, pos_a - flank + 1, pos_a)
    else:
        left = revcomp(_seq(genome, chrom_a, pos_a, pos_a + flank - 1))
    if orient_b == "+":
        right = _seq(genome, chrom_b, pos_b, pos_b + flank - 1)
    else:
        right = revcomp(_seq(genome, chrom_b, pos_b - flank + 1, pos_b))
    return left + right


def sv_junctions(genome, sv_rows, flank=FLANK):
    """One junction record per structural variant that creates a new join."""
    out = []
    for r in sv_rows:
        t = r.get("svtype")
        chrom = r.get("chrom")
        try:
            start = int(r["start"])
            end = int(r["end"]) if r.get("end") not in ("", None) else start
        except (KeyError, ValueError):
            continue
        end_chrom = r.get("end_chrom") or chrom
        if t == "DEL":
            seq = junction_sequence(genome, chrom, start, chrom, end, "+", "+", flank)
            desc = "deletion join"
        elif t == "DUP":
            # a tandem duplication reads as the end of the duplicated span joined back to its start
            seq = junction_sequence(genome, chrom, end, chrom, start, "+", "+", flank)
            desc = "tandem duplication join"
        elif t == "INV":
            seq = junction_sequence(genome, chrom, start, chrom, end, "+", "-", flank)
            desc = "inversion join"
        elif t == "TRA":
            seq = junction_sequence(genome, chrom, start, end_chrom, int(r["end"]), "+", "+", flank)
            desc = "translocation join"
        elif t == "INS":
            # a mobile element insertion: the two host-element boundaries, with the element unknown to the
            # reference, so only the host side of each is emitted
            seq = junction_sequence(genome, chrom, start, chrom, start + 1, "+", "+", flank)
            desc = "insertion breakpoint"
        else:
            continue
        if len(seq) < 2 * flank * 0.5:
            continue
        out.append({
            "id": r.get("event_id", ""), "kind": "sv", "svtype": t, "desc": desc,
            "chrom": chrom, "pos": start, "end_chrom": end_chrom, "end": end,
            "clone": r.get("clone", "T"), "hap": int(r.get("haplotype", 0) or 0),
            "sequence": seq,
        })
    return out


def fusion_junctions(genome, fusion_rows, flank=FLANK):
    """One junction record per fusion that has a genomic breakpoint.

    A read-through fusion has no DNA breakpoint and contributes nothing here; it exists only in RNA.
    """
    out = []
    for r in fusion_rows:
        bp5, bp3 = r.get("breakpoint_5p"), r.get("breakpoint_3p")
        if not bp5 or not bp3:
            continue
        try:
            p5, p3 = int(bp5), int(bp3)
        except ValueError:
            continue
        o5 = "+" if r.get("strand_5p", "+") == "+" else "-"
        o3 = "+" if r.get("strand_3p", "+") == "+" else "-"
        seq = junction_sequence(genome, r["chrom_5p"], p5, r["chrom_3p"], p3, o5, o3, flank)
        if len(seq) < flank:
            continue
        out.append({
            "id": r.get("event_id", ""), "kind": "fusion", "svtype": r.get("mechanism", ""),
            "desc": f"{r.get('gene_5p')}-{r.get('gene_3p')} breakpoint",
            "chrom": r["chrom_5p"], "pos": p5, "end_chrom": r["chrom_3p"], "end": p3,
            "clone": r.get("clone", "T"), "hap": int(r.get("haplotype", 0) or 0),
            "sequence": seq,
        })
    return out


def viral_junctions(genome, viral_fasta, virus_rows, flank=FLANK):
    """Host-virus junctions for integrated viruses, plus the episomal genomes themselves."""
    import os
    out = []
    if not viral_fasta or not os.path.exists(viral_fasta):
        return out
    import pysam
    ref = pysam.FastaFile(viral_fasta)
    for r in virus_rows:
        acc = (r.get("accession") or "").split(".")[0]
        contig = next((c for c in ref.references if acc and acc in c), None)
        if contig is None:
            continue
        vseq = ref.fetch(contig)
        if r.get("subclass") == "integrated" and r.get("integration_pos"):
            pos = int(r["integration_pos"])
            chrom = r["chrom"]
            left = _seq(genome, chrom, pos - flank + 1, pos)
            # both sides of the insertion: host into virus, and virus back into host
            out.append({
                "id": r["event_id"] + ":5p", "kind": "virus", "svtype": "integration",
                "desc": f"{r.get('virus')} host junction 5'", "chrom": chrom, "pos": pos,
                "end_chrom": contig, "end": 1, "clone": r.get("clone", "T"),
                "hap": int(r.get("haplotype", 0) or 0), "sequence": left + vseq[:flank],
            })
            out.append({
                "id": r["event_id"] + ":3p", "kind": "virus", "svtype": "integration",
                "desc": f"{r.get('virus')} host junction 3'", "chrom": chrom, "pos": pos,
                "end_chrom": contig, "end": len(vseq), "clone": r.get("clone", "T"),
                "hap": int(r.get("haplotype", 0) or 0),
                "sequence": vseq[-flank:] + _seq(genome, chrom, pos + 1, pos + flank),
            })
        else:
            out.append({
                "id": r["event_id"], "kind": "virus", "svtype": "episome",
                "desc": f"{r.get('virus')} episome", "chrom": contig, "pos": 1,
                "end_chrom": contig, "end": len(vseq), "clone": r.get("clone", "T"),
                "hap": 0, "sequence": vseq,
            })
    return out
