"""Check the declared allele-to-haplotype mapping against the germline VCF's own phase.

`design.yaml`'s `hla_haplotypes` says which HLA allele sits on which haplotype. Nothing derives it: a
SNV-level phased VCF does not name HLA alleles, so the arrangement is declared from linkage
disequilibrium and the declaration carried a comment saying it was unverified. The CNA that models HLA
LOH removes a haplotype *by index*, and `hla_loss.py` turns that index back into allele names using the
declaration -- so if the declaration is inverted relative to the VCF, every HLA LOH truth row names the
wrong alleles while the reads are perfectly correct. That is exactly what happened to IGI-SYN-SEQ-01:
measured here, hap1 carries A*01:01/B*35:08/C*04:01 and the design declared hap0. LOHHLA, reading the
reads, called those three lost at ~0.3x of normal coverage, matching purity 0.70; the truth bundle named
the other three. See `docs/build-reference.md` section 34.

The check needs no external database. GRCh38's analysis set already ships 525 `HLA-*` allele contigs, so
each declared two-field allele has a full-length sequence in the same FASTA the build already opens. The
test is k-mer containment: reconstruct both haplotypes over the gene window from the VCF, then ask what
fraction of each allele contig's k-mers appear in each haplotype. The haplotype genuinely carrying an
allele contains far more of it. Containment rather than Jaccard because a reconstructed haplotype is
reference-plus-variants and so is longer and noisier than the clean allele sequence, and both
orientations are tried because the class I genes are not all on the same strand as the contigs.

Why this is not merely a restatement of the declaration: the two inputs are independent. The declaration
comes from LD reasoning in `design.yaml`; the measurement comes from phase in the VCF plus allele
sequence in the reference. Nothing in the pipeline connected them until this module.
"""
import os

K = 31

# GRCh38 class I gene windows, matching the ones run_germline_spikein.py asserts over.
CLASS_I_WINDOWS = {"HLA-A": ("chr6", 29942532, 29945870),
                   "HLA-B": ("chr6", 31353875, 31357188),
                   "HLA-C": ("chr6", 31268749, 31272130)}

# A margin below this is reported as indeterminate rather than as an inversion. ds-01's margins are
# 0.30-0.46; ds-02's, which is the less divergent genotype and whose MHC comes from short-read calls
# rather than an assembly projection, are 0.07-0.15. Anything under 0.02 is noise either way.
MIN_MARGIN = 0.02


def _kmers(s):
    s = s.upper()
    return {s[i:i + K] for i in range(len(s) - K + 1) if "N" not in s[i:i + K]}


def _revcomp(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def _containment(hap, allele):
    """Fraction of the allele's k-mers present in the haplotype, best of the two orientations."""
    a = _kmers(allele)
    if not a:
        return 0.0
    return len(a & (_kmers(hap) | _kmers(_revcomp(hap)))) / len(a)


def _contig_for(fasta, allele):
    """The longest non-null GRCh38 contig for a two-field allele, or None.

    Allele names in the design are two-field (`HLA-B*35:08`); the contigs carry the full field path
    (`HLA-B*35:08:01`), so the match is a prefix up to the next colon. Null and low-expression alleles
    (the `N`, `L`, `Q`, `S` suffixes) are skipped when anything else is available, because their
    sequences carry the very defect that makes them null.
    """
    names = [r for r in fasta.references if r.startswith(allele + ":")]
    if not names:
        return None
    plain = [n for n in names if not n[-1] in "NLQS"] or names
    return max(plain, key=lambda n: fasta.get_reference_length(n))


def verify(fasta, germline, by_hap):
    """Measure which haplotype each declared class I allele sits on.

    `fasta`: an open pysam.FastaFile on the reference. `germline`: a Germline. `by_hap`: {hap: [alleles]}
    as `HlaLoss` parses it. Returns (rows, problems). A row is one allele's measurement; `problems` holds
    the clear inversions, which the caller should treat as fatal.
    """
    class _G:
        def seq(self, chrom, s0, e1):
            return fasta.fetch(chrom, s0, e1).upper()

    declared = {}
    for hap, alleles in by_hap.items():
        for a in alleles:
            declared.setdefault(a, set()).add(int(hap))

    rows, problems = [], []
    for locus, (chrom, start1, end1) in sorted(CLASS_I_WINDOWS.items()):
        here = sorted(a for a in declared if a.split("*")[0] == locus)
        # A homozygous locus carries no phase information: the same allele is on both haplotypes, so
        # there is nothing to invert and nothing to check. Recorded, not silently dropped.
        if len(here) < 2:
            rows.append({"locus": locus, "allele": here[0] if here else None, "status": "homozygous",
                         "declared": None, "measured": None, "hap0": None, "hap1": None})
            continue
        h0, h1 = germline.haplotype_seqs(_G(), chrom, start1, end1)
        for allele in here:
            contig = _contig_for(fasta, allele)
            if contig is None:
                rows.append({"locus": locus, "allele": allele, "status": "no reference contig",
                             "declared": sorted(declared[allele]), "measured": None,
                             "hap0": None, "hap1": None})
                continue
            seq = fasta.fetch(contig).upper()
            c0, c1 = _containment(h0, seq), _containment(h1, seq)
            decl = sorted(declared[allele])
            if c0 == c1 and h0 == h1:
                # Not indeterminate -- the two haplotypes are the same sequence, so the locus carries no
                # allelic difference for anything to find. That is the IGI-SYN-SEQ-02 defect of section
                # 33: a locus declared heterozygous whose germline VCF has no variants in it, leaving
                # both haplotypes as bare reference. run_germline_spikein.py refuses it earlier and more
                # specifically, but a second reading costs nothing and this one would otherwise pass it.
                status, measured = "NO ALLELIC DIFFERENCE", None
            elif abs(c0 - c1) < MIN_MARGIN:
                status, measured = "indeterminate", None
            else:
                measured = 0 if c0 > c1 else 1
                status = "ok" if measured in decl else "INVERTED"
            rows.append({"locus": locus, "allele": allele, "contig": contig, "status": status,
                         "declared": decl, "measured": measured,
                         "hap0": round(c0, 4), "hap1": round(c1, 4)})
            if status == "INVERTED":
                problems.append(f"{allele} is declared on hap{'/'.join(str(d) for d in decl)} but its "
                                f"sequence is on hap{measured} (containment hap0={c0:.4f} hap1={c1:.4f} "
                                f"against {contig}); the HLA LOH truth would name the wrong alleles")
            elif status == "NO ALLELIC DIFFERENCE":
                # one fact about the locus, not one per allele, so report it once
                msg = (f"{locus} is declared heterozygous but both haplotypes are the same sequence over "
                       f"{chrom}:{start1}-{end1}, so the reads there carry neither "
                       f"{' nor '.join(here)}; the germline VCF has no variants in it")
                if msg not in problems:
                    problems.append(msg)
    return rows, problems


def verify_or_die(reference_fasta, germline, by_hap, dataset, log=print):
    """Run `verify` and raise on a clear inversion. Prints every row, including the ones it cannot judge.

    An indeterminate locus does not fail the build. Failing on one would block a genotype that is simply
    less divergent -- ds-02's margins are a fifth of ds-01's -- and the honest outcome there is a loud
    "not checked", not a refusal. A clear inversion is different: it means the truth bundle and the reads
    disagree, and that has shipped once already.
    """
    try:
        import pysam
        fasta = pysam.FastaFile(reference_fasta)
    except Exception as e:                                  # noqa: BLE001 - any failure here is the same
        log(f"[hla] {dataset}: NOT CHECKED -- could not open {reference_fasta} ({e})")
        return []
    if not any(r.startswith("HLA-") for r in fasta.references):
        log(f"[hla] {dataset}: NOT CHECKED -- {os.path.basename(reference_fasta)} carries no HLA-* "
            f"contigs, so the declared haplotype assignment cannot be measured against it")
        return []
    rows, problems = verify(fasta, germline, by_hap)
    for r in rows:
        if r["status"] == "homozygous":
            log(f"[hla] {dataset}: {r['locus']} homozygous ({r['allele']}), carries no phase information")
        elif r["measured"] is None:
            log(f"[hla] {dataset}: {r['locus']} {r['allele']} {r['status']}"
                + (f" (hap0={r['hap0']} hap1={r['hap1']})" if r["hap0"] is not None else ""))
        else:
            log(f"[hla] {dataset}: {r['allele']} on hap{r['measured']} "
                f"(declared hap{'/'.join(str(d) for d in r['declared'])}, "
                f"hap0={r['hap0']} hap1={r['hap1']})  {r['status']}")
    if problems:
        raise SystemExit(f"[hla] FAIL: {dataset}'s hla_haplotypes disagrees with the phase in its "
                         "germline VCF:\n  - " + "\n  - ".join(problems))
    return rows
