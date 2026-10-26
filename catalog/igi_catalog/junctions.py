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
import re

from .genome import normalize_bases, revcomp
from .mei import insertion as mei_insertion

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


def _mei_rng(dataset, event_id):
    """The SAME stream `derived.mei_sequence_factory` uses, so short and long reads agree on the element.

    If these drifted, the exome and WGS contigs would carry a different insertion from the one in the
    derived chromosome, and a caller comparing assays would see two different events at one site.
    """
    import random
    return random.Random(f"{dataset}:{event_id}:mei")


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
            # A mobile-element insertion has TWO boundaries and neither is a join between two reference
            # positions. This used to emit `junction_sequence(chrom, start, chrom, start + 1)` -- two
            # adjacent bases of the reference -- so the contig was reference sequence, the reads from it
            # aligned perfectly, and the insertion left no soft clip and no improper pair anywhere. The
            # element itself has to be in the contig, or there is nothing for a caller to find.
            element = mei_insertion(genome, r.get("inserted_sequence", ""), chrom, start,
                                    _mei_rng(r.get("dataset", ""), r.get("event_id", "")))
            if not element:
                continue
            host5 = _seq(genome, chrom, start - flank + 1, start)
            host3 = _seq(genome, chrom, start + 1, start + flank)
            common = {"kind": "sv", "svtype": t, "chrom": chrom, "pos": start,
                      "end_chrom": chrom, "end": end,
                      "clone": r.get("clone", "T"), "hap": int(r.get("haplotype", 0) or 0)}
            for side, seq_i, desc_i in (
                    ("5p", host5 + element[:flank], "insertion 5' junction"),
                    ("3p", element[-flank:] + host3, "insertion 3' junction")):
                if len(seq_i) < flank:
                    continue
                out.append(dict(common, id=f"{r.get('event_id', '')}:{side}",
                                desc=desc_i, sequence=seq_i))
            continue
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
        vseq = normalize_bases(ref.fetch(contig).upper())
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
        # An episome is deliberately NOT emitted here. It has no human coordinate, so it was given the
        # viral contig in `chrom` -- and the per-task junction filter and `WesBuilder.junctions_in` are
        # both keyed by human chromosome, so it was looked up under a chromosome that is never queried
        # and silently dropped from every library. It is material, not a breakpoint: see
        # `viral_payloads`. Leaving a second mechanism here that claims to emit it is how the first one
        # went unnoticed.
    return out


# The chromosome task that owns the viral payload. A virus is not on a human chromosome, so nothing
# selects it the way a breakpoint is selected by its position; without pinning it to exactly one task it
# would be emitted by all 24 and the delivered copy number would be 24x the design.
VIRAL_PAYLOAD_CHROM = "chr1"


def viral_payloads(viral_fasta, virus_rows, tumor=True):
    """The viral genomes themselves, as material rather than as breakpoint evidence.

    `viral_junctions` emits the host-virus joins, and those are right: they are what an integration-site
    caller reads, and they are spliced into a host interval in place of an equal number of bases so the
    host's own coverage is unchanged. But a junction is 2*FLANK of sequence, so the integrated virus
    reached the library as 1.2 kb of its two ends and nothing else -- HPV16 is 7.9 kb, so E6 and E7 had
    no DNA coverage at any depth while the truth table claimed 3 copies per cell.

    An episome is worse. It has no human coordinate at all, so `viral_junctions` puts the viral contig in
    `chrom`, and both the per-task junction filter and `WesBuilder.junctions_in` -- which is keyed by
    (human chromosome, clone, haplotype) -- look for it under a chromosome that is never queried. EBV at
    5 copies per cell and HPV18 at 2 measured 0.07 % breadth and 0.0016x depth in the delivered ds-02
    tumour WGS: absent, not thin.

    Both are the same shape of problem: the virus's own genome is extra material, not a replacement for
    host sequence, so it belongs in its own source at its own weight rather than on the junction path.

    Weight is in the units `source_plan` uses, where one haplotype copy of a clone is `purity * frac`:
    `copies_per_cell` copies carried by a clone at `ccf` is `purity * ccf * copies_per_cell`. A normal
    library is flat and diploid, so a trace carry-over is `normal_trace_copies` against that plan.
    """
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
        vseq = normalize_bases(ref.fetch(contig).upper())
        if not vseq:
            continue
        try:
            copies = float(r.get("copies_per_cell") or 0)
            trace = float(r.get("normal_trace_copies") or 0)
            ccf = float(r.get("ccf") or 1.0)
        except ValueError:
            continue
        n = copies if tumor else trace
        if n <= 0:
            continue
        # A filesystem- and shell-safe name. The viral reference uses NCBI headers --
        # `gi|310698439|ref|NC_001526.2|` -- and the contig name reaches a FASTA path and a chunk tag
        # that are interpolated into the simulator's command line. A pipe there is a shell pipe: badread
        # produced nothing for all three viruses, including the 171 kb EBV, so it was the name and not
        # the length. Prefer the accession, which is already a clean identifier.
        slug = re.sub(r"[^A-Za-z0-9._-]", "_", (r.get("accession") or contig)).strip("_") or "virus"
        out.append({
            "slug": slug,
            "id": r["event_id"] + ":genome", "kind": "virus",
            "svtype": r.get("subclass") or "episome",
            "desc": f"{r.get('virus')} genome at {n:g} copies/cell",
            "contig": contig, "accession": r.get("accession", ""),
            "clone": r.get("clone", "T"), "ccf": ccf,
            "copies": n, "tumor_weighted": bool(tumor),
            "sequence": vseq,
        })
    return out
