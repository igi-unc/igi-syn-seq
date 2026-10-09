"""Reference genome access and sequence-context helpers."""
import pysam

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def revcomp(s): return s.translate(COMP)[::-1]

# GRCh38 and the viral reference both contain IUPAC ambiguity codes -- chr1 carries one M and one R, chr3
# a B, an R, two W and three Y, and virus_unmasked.02ec8.fa carries 497 of them. Delivered reads must not,
# for three reasons:
#
#  - No instrument emits them. A real sequencer reading an ambiguous locus reports a definite base or N,
#    because the ambiguity is in the reference's knowledge, not in the molecule.
#  - Consumers refuse them, and not only the one we predicted. This comment first said "Cell Ranger
#    rejects any base outside ACGTN ... an ambiguity code reaching a 10x arm would reproduce that
#    failure exactly" -- true, but it understated the blast radius, and the arm that actually broke was
#    the exome. OptiType's razers3 aborts while loading its first chunk with
#        seqan::ParseError: Unexpected character 'Y' found.
#    because SeqAn's reader accepts only A, C, G, T and N. It failed all seven retries while every
#    aligner in the same run read the file without complaint. The rule is therefore not "normalise for
#    the strict consumers we can name" but "a delivered base is a base", because a reference dataset
#    cannot know what will read it.
#  - revcomp() above maps only ACGTN, so an IUPAC base survives reverse-complementing UNCHANGED. That is
#    latent rather than observed: the paths that would hit it are minus-strand transcripts
#    (transcriptome.py), inverted SV segments (rearrange.py) and fusion partners. The one strand-flipped
#    case actually found in delivered reads -- chr3's single B appearing as V on the minus strand -- was
#    Badread complementing correctly on its own. Normalising here removes the class either way, because
#    revcomp then only ever sees ACGTN.
#
# They become N rather than a definite base drawn from the code's constituents. A definite base would be
# more realistic, but it would differ from the reference at a position no truth table records, and a
# caller would report it as a variant absent from the truth VCF -- the same mistake as writing `*` into
# the sequence, in miniature. Callers skip N.
IUPAC_CODES = "RYSWKMBDHVryswkmbdhv"
TO_N = str.maketrans(IUPAC_CODES, "N" * len(IUPAC_CODES))


def normalize_bases(seq):
    """Replace IUPAC ambiguity codes with N. Anything outside ACGTN after this is a defect, not a base."""
    return seq.translate(TO_N)

class Genome:
    def __init__(self, fasta):
        self.fa = pysam.FastaFile(fasta)
        self.lengths = dict(zip(self.fa.references, self.fa.lengths))
    def seq(self, chrom, start, end):
        """0-based half-open, uppercase, ambiguity codes normalised to N (see normalize_bases)."""
        return normalize_bases(
            self.fa.fetch(chrom, max(0, start), min(end, self.lengths[chrom])).upper())
    def homopolymer_run(self, chrom, pos0):
        """Length of the longest homopolymer run touching 0-based position pos0 (run of the base at pos0 or its neighbours)."""
        s = self.seq(chrom, pos0 - 12, pos0 + 13); c = 12
        best = 1
        for anchor in (c - 1, c, c + 1):
            if anchor < 0 or anchor >= len(s): continue
            b = s[anchor]; i = anchor; j = anchor
            while i > 0 and s[i - 1] == b: i -= 1
            while j < len(s) - 1 and s[j + 1] == b: j += 1
            best = max(best, j - i + 1)
        return best
    def gc(self, chrom, start, end):
        s = self.seq(chrom, start, end); n = len(s)
        return (s.count("G") + s.count("C")) / n if n else 0.0


def left_align(genome, chrom, pos1, ref, alt):
    """Shift an indel to its leftmost equivalent representation, as `bcftools norm` would.

    Inside a repeat the same deletion has several equally valid representations, and a caller's VCF is
    normalised while a truth table built from the designer's own coordinates is not. An unnormalised
    truth row does not merely fail to match: it is scored as a false negative for the event that is
    really there and a false positive for the call that found it.
    """
    while len(ref) > 1 and len(alt) > 1:
        if ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
        elif ref[0] == alt[0]:
            ref, alt, pos1 = ref[1:], alt[1:], pos1 + 1
        else:
            break
    while pos1 > 1 and (len(ref) == 1 or len(alt) == 1) and ref[-1] == alt[-1]:
        prev = genome.seq(chrom, pos1 - 2, pos1 - 1)
        ref, alt, pos1 = prev + ref[:-1], prev + alt[:-1], pos1 - 1
    return pos1, ref, alt
