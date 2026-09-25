"""Phased germline access: haplotype sequences for a window and nearby heterozygous sites."""
import pysam
from .genome import revcomp

class Germline:
    def __init__(self, vcf, sample=None, sex="male"):
        self.vcf = pysam.VariantFile(vcf); self.sample = sample or self.vcf.header.samples[0]; self.sex = sex
    def variants(self, chrom, start, end):
        """Yield (pos1, ref, alt, gt_tuple, phased) for records overlapping [start,end) 0-based; biallelic GT only."""
        try: it = self.vcf.fetch(chrom, start, end)
        except ValueError: return
        for rec in it:
            s = rec.samples[self.sample]; gt = s.get("GT")
            if gt is None or None in gt or rec.alts is None: continue
            if any(a > 1 for a in gt): continue  # multiallelic genotypes are skipped for window construction
            yield rec.pos, rec.ref, rec.alts[0], tuple(gt), bool(s.phased) or len(gt) == 1
    def overlaps_variant(self, chrom, pos1, ref_len=1):
        """Whether any germline record overlaps the span [pos1, pos1 + ref_len). A somatic allele that
        overlaps a germline variant cannot be expressed unambiguously in reference coordinates, so the
        designer excludes such sites."""
        start, end = pos1 - 1, pos1 - 1 + ref_len
        for pos, r, _a, _gt, _ph in self.variants(chrom, max(0, start - 60), end + 60):
            if pos - 1 < end and pos - 1 + len(r) > start:
                return True
        return False

    def hets_near(self, chrom, pos1, window=30):
        return [v for v in self.variants(chrom, pos1 - 1 - window, pos1 + window) if len(v[3]) == 2 and v[3][0] != v[3][1] and v[0] != pos1]
    def haplotype_seqs(self, genome, chrom, start1, end1, extra=None):
        """Two haplotype sequences (hap0, hap1) for the 1-based inclusive window with phased germline variants applied.
        extra: optional (pos1, ref, alt, hap_index) somatic change applied on one haplotype -> returns (h0, h1, m0, m1)."""
        ref = genome.seq(chrom, start1 - 1, end1)
        edits = {0: [], 1: []}
        for pos, r, a, gt, phased in self.variants(chrom, start1 - 1, end1):
            if pos < start1 or pos + len(r) - 1 > end1: continue
            haps = (gt[0], gt[1]) if len(gt) == 2 else (gt[0], gt[0])
            for h, allele in enumerate(haps):
                if allele == 1: edits[h].append((pos - start1, r, a))
        def apply(seq, ed, strict_from=None):
            """Apply reference-coordinate edits right to left. Germline records that do not match are
            skipped (overlapping or inconsistent calls); a somatic edit that does not match is an error,
            because silently dropping it would produce a truth row describing a change that is not there."""
            out = seq
            for off, r, a in sorted(ed, key=lambda x: -x[0]):
                if out[off:off + len(r)] != r:
                    if strict_from is not None and (off, r, a) == strict_from:
                        raise ValueError(
                            f"somatic edit {r}>{a} at window offset {off} does not match the haplotype "
                            f"sequence {out[off:off + len(r)]!r}; the site overlaps a germline variant")
                    continue
                out = out[:off] + a + out[off + len(r):]
            return out
        h0, h1 = apply(ref, edits[0]), apply(ref, edits[1])
        if extra is None: return h0, h1
        pos, r, a, hi = extra
        som = (pos - start1, r, a)
        m = [h0, h1]
        m[hi] = apply(ref, edits[hi] + [som], strict_from=som)
        return h0, h1, m[0], m[1]
