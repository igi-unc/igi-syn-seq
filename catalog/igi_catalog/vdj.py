"""T-cell receptor V(D)J transcript assembly from real GRCh38 segments.

A synthetic TCR library is only useful if a V(D)J caller can recover the clonotype from it, and that needs
genuine V, J and C framework sequence -- a caller aligns reads to germline segments and infers the junction.
Invented framework would make every read unassignable and the assay would test nothing.

The segments come from the GTF and the genome directly rather than through the annotation cache, because
that cache filters TR gene types out: of 202 TR segments present in GENCODE v37, exactly one survives into
`env.tx`, and it is mistyped as lncRNA. 67 TRBV, 14 TRBJ, 2 TRBC, 54 TRAV, 60 TRAJ and 1 TRAC are available
when read straight from the GTF.

A transcript is assembled as V + CDR3 + J + C. The roster supplies the CDR3 nucleotide sequence per
clonotype, so cells sharing a clonotype share a junction, which is what clonal expansion means and what a
caller is being asked to detect.
"""
import collections
import gzip
import re

COMP = str.maketrans("ACGTN", "TGCAN")
_NAME = re.compile(r'gene_name "(TR[AB][VDJC][^"]*)"')


def load_segments(gtf_path, genome):
    """{segment_name: sequence} for every TR V, D, J and C segment in the annotation."""
    exons = collections.defaultdict(set)
    op = gzip.open if str(gtf_path).endswith(".gz") else open
    with op(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            m = _NAME.search(f[8])
            if m:
                exons[m.group(1)].add((f[0], int(f[3]) - 1, int(f[4]), f[6]))
    out = {}
    for name, ex in exons.items():
        ex = sorted(ex)
        s = "".join(genome.seq(c, a, b) for c, a, b, _ in ex)
        if ex[0][3] == "-":
            s = s.translate(COMP)[::-1]
        if len(s) >= 30:
            out[name] = s
    return out


def pick(segments, prefix, rng):
    """One segment of a family, excluding pseudogene-looking names."""
    names = sorted(n for n in segments if n.startswith(prefix) and "OR" not in n)
    return rng.choice(names) if names else None


def assemble(segments, chain, cdr3_nt, rng):
    """A full-length TCR transcript: V + CDR3 + J + C, with the segment names used.

    Returns (sequence, {"v":..., "j":..., "c":...}) or (None, None) if the chain's segments are missing.
    """
    v = pick(segments, "TRBV" if chain == "TRB" else "TRAV", rng)
    j = pick(segments, "TRBJ" if chain == "TRB" else "TRAJ", rng)
    c = "TRBC2" if chain == "TRB" else "TRAC"
    if not (v and j and c in segments):
        return None, None
    # The V segment's 3' end and the J segment's 5' end are replaced by the junction, which is what the
    # CDR3 is: trimming a little of each is what makes the junction look like a real recombination rather
    # than a clean concatenation a caller would never see.
    vt = segments[v][:-rng.randint(0, 6)] if len(segments[v]) > 20 else segments[v]
    jt = segments[j][rng.randint(0, 4):]
    return vt + cdr3_nt + jt + segments[c], {"v": v, "j": j, "c": c}


# One codon per amino acid, chosen by human usage frequency. Reverse translation is deterministic given the
# peptide, which is what matters: cells of one clonotype must share a nucleotide junction, because that is
# what clonal expansion is and what a caller is asked to detect. Drawing codons per cell would give every
# cell of a clone a different junction and destroy the thing being tested.
_CODON = {
    "A": "GCC", "R": "CGG", "N": "AAC", "D": "GAC", "C": "TGC", "Q": "CAG", "E": "GAG",
    "G": "GGC", "H": "CAC", "I": "ATC", "L": "CTG", "K": "AAG", "M": "ATG", "F": "TTC",
    "P": "CCC", "S": "AGC", "T": "ACC", "W": "TGG", "Y": "TAC", "V": "GTG",
}


def cdr3_nucleotide(peptide):
    """Reverse-translate a CDR3 peptide. Deterministic, so one clonotype has one junction.

    The roster's CDR3s are synthetic peptides rather than real recombinants -- `cells.CellRoster._cdr3`
    says so -- which means this library exercises clonotype clustering and expansion, not germline V-J
    assignment. A caller asked to name the true V gene would be being asked something the data cannot
    answer, and that limit is inherited here rather than papered over.
    """
    return "".join(_CODON.get(aa, "NNN") for aa in (peptide or "").upper())
