"""Real mobile-element sequence for a designed insertion.

The catalog records an insertion as a family and a length -- `<L1HS:6033>` -- and for a long time nothing
turned that token into bases. The short-read path emitted the insertion junction as plain reference
(`junction_sequence(chrom, start, chrom, start + 1)`, two adjacent bases, so a contig identical to the
reference) and the long-read path inserted 300 random bases whatever the family said. The consequence was
measured on the delivered libraries: at 14 of 16 planted sites across both datasets there were zero
soft-clipped reads and zero improper pairs at ~390x, and GRIDSS and xTea were right to call nothing.
There was nothing to find.

Elements are taken from the reference rather than from a vendored consensus. A real genomic copy carries
the divergence, the 5' truncation and the internal structure a caller actually sees, it needs no new
download or licence, and it cannot drift from the assembly the reads are simulated against. The three
donors below are the highest-scoring full-length copy of each family on a primary chromosome, forward
strand, selected from the UCSC RepeatMasker annotation (`rmsk_bed`) with:

    zcat hg38.rmsk.bed.gz | awk -F'\t' '$1 ~ /^chr([0-9]+|X)$/ && $6=="+" {split($4,a,"|"); ...}'

Their 5' ends are canonical -- L1HS opens AGCTCCGGTCTACAGCTCCCAGCGTGAGCGACGCAGAAGA, AluY opens
GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTT, SVA_E opens with the (CCCTCT)n hexamer -- which is the check
that they are what the annotation says.

A placed insertion is built the way a real retrotransposition product looks:

  * **5' truncation.** Most L1 insertions are 5' truncated, so the body is the 3'-MOST bases of the donor.
    The 3' end is the part that always arrives, because reverse transcription starts there.
  * **A poly-A tail**, drawn per insertion. The donors carry their own, which is stripped first so the
    tail length is set here rather than inherited twice.
  * **A target-site duplication.** The host bases immediately 5' of the insertion point are repeated on
    its 3' side. The first copy is already in the host sequence, so only the second is emitted. This is
    the signature MELT and xTea key on, and without it an insertion looks like a plain novel-sequence SV.
  * **A strand**, since roughly half of insertions are antisense. The body and tail are reverse
    complemented together, so an antisense insertion presents a poly-T at its 5' junction, which is what
    a caller sees in real data.

The realised length is EXACTLY the catalogued length. That is deliberate: `rearrange.plan_chromosome`
sizes the insert segment from the same token without calling this module, so the coordinate map and the
derived sequence cannot disagree. Before this, the insert segment was `Segment(chrom, start, start)` --
length 1 -- while `realise` appended 300 bases, so every segment after an insertion was already mapped
299 bp out of true. At 6 kb that error would have been 6 kb.
"""
from .genome import normalize_bases, revcomp

# family -> (chrom, start1, end1), 1-based inclusive
DONORS = {
    "L1HS":  ("chr5", 152_076_869, 152_082_891),      # 6,023 bp, RepeatMasker score 29,598
    "AluY":  ("chr12", 62_838_640, 62_838_950),       #   311 bp, score 3,139
    "SVA_E": ("chr19", 52_810_861, 52_812_782),       # 1,922 bp, score 9,133
}

TSD_RANGE = (7, 20)        # target-site duplications are short and tight; MELT's own window is similar
POLYA_RANGE = (20, 60)
MIN_BODY = 100             # below this it is not an element any more, it is noise

_cache = {}


def donor_body(genome, family):
    """The donor element with its own poly-A stripped, cached per family."""
    if family in _cache:
        return _cache[family]
    if family not in DONORS:
        return ""
    chrom, start1, end1 = DONORS[family]
    if chrom not in genome.lengths:
        return ""
    seq = normalize_bases(genome.seq(chrom, start1 - 1, end1).upper())
    seq = seq.rstrip("A")              # the tail is drawn per insertion, not inherited
    _cache[family] = seq
    return seq


def family_of(token):
    """`<L1HS:6033>` -> ('L1HS', 6033). Returns (None, 0) for anything else."""
    if not token or not token.startswith("<") or not token.endswith(">") or ":" not in token:
        return None, 0
    fam, _, n = token[1:-1].partition(":")
    try:
        return fam, int(n)
    except ValueError:
        return fam, 0


def insertion(genome, token, host_chrom, host_pos, rng):
    """The bases a designed insertion contributes, exactly `total` long.

    `host_pos` is the 1-based reference position the element lands after, so the target-site duplication
    is the host bases ENDING there -- the copy the host already has.
    """
    family, total = family_of(token)
    if not family or total <= 0:
        return ""
    body_src = donor_body(genome, family)
    if not body_src:
        return ""

    n_tsd = rng.randint(*TSD_RANGE)
    n_tsd = min(n_tsd, max(0, total - MIN_BODY))
    n_polya = rng.randint(*POLYA_RANGE)
    n_body = total - n_tsd - n_polya
    if n_body > len(body_src):
        # the donor is shorter than the design asks for: give the remainder to the tail, which is where
        # the length variation of a real insertion mostly lives anyway
        n_polya += n_body - len(body_src)
        n_body = len(body_src)
    elif n_body < MIN_BODY:
        n_polya = max(0, total - n_tsd - MIN_BODY)
        n_body = total - n_tsd - n_polya

    body = body_src[-n_body:] if n_body else ""          # 5' truncated: the 3' end always arrives
    tsd = ""
    if n_tsd and host_chrom in genome.lengths:
        lo = max(0, host_pos - n_tsd)
        tsd = normalize_bases(genome.seq(host_chrom, lo, host_pos).upper())[:n_tsd]
    tsd = (tsd + "A" * n_tsd)[:n_tsd]                    # near a contig edge, pad rather than shorten

    element = body + "A" * n_polya
    if rng.random() < 0.5:
        element = revcomp(element)                       # antisense: a poly-T at the 5' junction
    out = element + tsd
    if len(out) > total:
        out = out[:total]
    elif len(out) < total:
        out += "A" * (total - len(out))
    return out
