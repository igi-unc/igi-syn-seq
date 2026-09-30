"""Reference genome access and sequence-context helpers."""
import pysam

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def revcomp(s): return s.translate(COMP)[::-1]

class Genome:
    def __init__(self, fasta):
        self.fa = pysam.FastaFile(fasta)
        self.lengths = dict(zip(self.fa.references, self.fa.lengths))
    def seq(self, chrom, start, end):
        """0-based half-open, uppercase."""
        return self.fa.fetch(chrom, max(0, start), min(end, self.lengths[chrom])).upper()
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
