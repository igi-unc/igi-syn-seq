"""T-cell receptor V(D)J recombination and transcript assembly from real GRCh38 germline segments.

A synthetic TCR library is only useful if a V(D)J caller can recover the clonotype from it, and that needs
genuine V, J and C framework sequence -- a caller aligns reads to germline segments and infers the
junction. Invented framework would make every read unassignable and the assay would test nothing.

The segments come from the GTF and the genome directly rather than through the annotation cache, because
that cache filters TR gene types out: of 202 TR segments present in GENCODE v37, exactly one survives into
`env.tx`, and it is mistyped as lncRNA.

**The junction is germline-derived, not invented.** CDR3s used to be random peptides -- the literal string
`CAS` or `CAV`, then uniform draws over all twenty amino acids, then `F`. Two things were wrong with that.
It put an internal cysteine in roughly half of them, `CASCAWCQESIIYPATQVF` being a delivered example, which
no selected repertoire contains. And the junction had no relationship to the V and J genes the truth table
named, so a caller that correctly recovered the real V gene from the framework would disagree with the
truth table, and the library would score its user wrong.

A real CDR3 runs from the conserved cysteine at IMGT position 104, which the V segment contributes, to the
conserved phenylalanine of the J segment's `FGXG` motif at position 118. `resources/vdj_germline.json`
holds, for every functional V and J gene, exactly those nucleotides plus the framework either side; it is
built from GENCODE over GRCh38 by `jobs/measure/fit_vdj_germline.py`. `recombine` then does what V(D)J
recombination does: trim each end, add an N region, and keep the result in frame.
"""
import collections
import gzip
import json
import os
import re

COMP = str.maketrans("ACGTN", "TGCAN")
_NAME = re.compile(r'gene_name "(TR[AB][VDJC][^"]*)"')

_B = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_AA = {a + b + c: _AAS[i * 16 + j * 4 + k]
            for i, a in enumerate(_B) for j, b in enumerate(_B) for k, c in enumerate(_B)}

# CDR3 amino-acid length, inclusive of the anchor C and F. Real human repertoires peak at 14 for beta and
# 13 for alpha; the ranges here are the bulk of the observed distribution rather than its extremes.
CDR3_LEN = {"TRB": (11, 19, 14), "TRA": (10, 17, 13)}

_GERM_PATH = os.path.join(os.path.dirname(__file__), "..", "resources", "vdj_germline.json")
_GERM = {}


def translate(nt):
    return "".join(CODON_AA.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def germline():
    """The germline junction anchors and framework, keyed by segment name."""
    if not _GERM:
        with open(os.path.abspath(_GERM_PATH)) as fh:
            _GERM.update(json.load(fh))
    return _GERM


def genes(chain, part):
    """Functional gene names for a chain, e.g. genes("TRB", "v")."""
    g = germline()
    return sorted(n for n in g[part] if n.startswith(chain[:3] + part.upper()))


def load_segments(gtf_path, genome):
    """{segment_name: sequence} for every TR V, D, J and C segment in the annotation.

    Kept for the callers that want whole segments; `germline()` is what the junction is built from.
    """
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
            out[name] = s.upper()
    return out


def _n_region(rng, n_nt):
    """Non-templated nucleotides, as TdT adds them: GC-rich, and free of stops and cysteines.

    Rejecting stops is mandatory -- a stop in the CDR3 is a non-productive contig and Cell Ranger will say
    so. Rejecting cysteines is selection: a repertoire that has been through the thymus carries internal
    cysteines at a percent or two, not at the 3 % two of sixty-four codons would give by chance.
    """
    if n_nt <= 0:
        return ""
    for _ in range(200):
        s = "".join(rng.choices("ACGT", weights=[0.2, 0.3, 0.3, 0.2], k=n_nt))
        if n_nt % 3 == 0 and any(CODON_AA.get(s[i:i + 3]) in ("*", "C") for i in range(0, n_nt, 3)):
            continue
        return s
    return s


def recombine(chain, rng, v=None, j=None):
    """One V(D)J recombinant.

    Returns {"v", "j", "c", "cdr3_nt", "cdr3_aa", "v_trim", "j_trim", "n_len"}. The CDR3 keeps the V's
    conserved cysteine and the J's conserved phenylalanine, so the result is productive and anchored, and
    its length falls in the real range for the chain.
    """
    g = germline()
    v = v or rng.choice(genes(chain, "v"))
    j = j or rng.choice(genes(chain, "j"))
    c = "TRBC2" if chain == "TRB" else "TRAC"
    vn, jn = g["v"][v]["cdr3_nt"], g["j"][j]["cdr3_nt"]

    lo, hi, mode = CDR3_LEN[chain]
    # trimming is in whole codons here so the anchors and the frame are never disturbed; real trimming is
    # per nucleotide and is absorbed by the N region, which is what makes junction length vary at all
    v_max = max(0, len(vn) // 3 - 1)          # never trim the cysteine
    j_max = max(0, len(jn) // 3 - 1)          # never trim the phenylalanine
    for _ in range(100):
        tv = rng.randint(0, min(2, v_max))
        tj = rng.randint(0, min(3, j_max))
        want = max(lo, min(hi, int(round(rng.gauss(mode, 2.0)))))
        keep = len(vn) // 3 - tv + len(jn) // 3 - tj
        n_aa = want - keep
        if 0 <= n_aa <= 8:
            break
    else:
        tv = tj = 0
        n_aa = max(0, mode - (len(vn) // 3 + len(jn) // 3))
    v_part = vn[:len(vn) - tv * 3]
    j_part = jn[tj * 3:]
    cdr3_nt = v_part + _n_region(rng, n_aa * 3) + j_part
    return {"v": v, "j": j, "c": c, "cdr3_nt": cdr3_nt, "cdr3_aa": translate(cdr3_nt),
            "v_trim": tv, "j_trim": tj, "n_len": n_aa * 3}


def assemble(rec, chain=None):
    """A full-length TCR transcript: V framework + CDR3 + J framework + C.

    The CDR3 already carries the V's conserved cysteine and the J's phenylalanine, so the framework either
    side is spliced on exactly where the junction took it from and the transcript reads through in frame.
    """
    g = germline()
    v, j, c = rec["v"], rec["j"], rec["c"]
    if v not in g["v"] or j not in g["j"] or c not in g["c"]:
        return None
    return g["v"][v]["framework"] + rec["cdr3_nt"] + g["j"][j]["framework"] + g["c"][c]["framework"]
