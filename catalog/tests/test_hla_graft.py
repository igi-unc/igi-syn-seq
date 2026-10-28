"""Regression tests for the class I allele graft.

The graft turns an alignment of an allele contig into phased germline records, and every way that can go
wrong produces a plausible-looking variant list rather than an error:

  * coordinate base  -- the cs tag walks the target in 0-based offsets while VCF counts from 1. The first
                        implementation mixed the two and placed every record one base off; the records
                        rebuilt nothing, which is what `_realize` refuses.
  * orientation      -- HLA-B and HLA-C are on the minus strand, so their contigs align reversed. minimap2
                        reports cs against the forward target, but assuming the opposite also yields a
                        full-looking list of substitutions.
  * anchoring        -- writing every record with a left anchor base turns each SNV into `GA>GG`, which is
                        legal and unreadable; writing none turns a deletion into an empty ALT.

These need no cluster, no reference and no network: they run the parser and the record builder over a
hand-made reference.
"""
import os
import sys

try:
    import pytest
except ImportError:
    pytest = None

CATALOG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, CATALOG)

from igi_catalog import hla_graft                                      # noqa: E402


class FakeFasta:
    """Just enough of pysam.FastaFile for `_records`: one chromosome, 1-based coordinates outside."""

    def __init__(self, seq, offset=1):
        self.seq, self.offset = seq, offset

    def fetch(self, _chrom, s0, e1):
        return self.seq[s0 - (self.offset - 1):e1 - (self.offset - 1)]


def paf(cs, tstart, qlen=10, qs=0, qe=10, strand="+", nmatch=10, blk=10):
    """A PAF row with only the fields `_cs_edits` reads."""
    return ["q", str(qlen), str(qs), str(qe), strand, "t", "0", str(tstart), "0",
            str(nmatch), str(blk), "60", f"cs:Z:{cs}"]


def test_cs_edits_are_zero_based_against_a_one_based_window():
    # window starts at chr6:101, alignment starts 5 bases in, then 3 matches before a substitution:
    # the substituted base is the 9th of the window, chr6:109, which is 0-based 108.
    edits = hla_graft._cs_edits(paf(":3*ag", tstart=5), win_start1=101)
    assert edits == [(108, "A", "G")]


def test_cs_edits_walk_insertions_without_consuming_reference():
    edits = hla_graft._cs_edits(paf(":2+tt:2*ca", tstart=0), win_start1=1)
    assert edits == [(2, "", "TT"), (4, "C", "A")]


def test_cs_edits_consume_reference_for_a_deletion():
    edits = hla_graft._cs_edits(paf(":2-acg:1*tc", tstart=0), win_start1=1)
    assert edits == [(2, "ACG", ""), (6, "T", "C")]


def test_realize_rebuilds_the_query_and_refuses_a_misplaced_edit():
    assert hla_graft._realize("AAACAAA", 0, [(3, "C", "G")]) == "AAAGAAA"
    try:
        hla_graft._realize("AAACAAA", 0, [(2, "C", "G")])
    except SystemExit as e:
        assert "does not match the reference" in str(e)
    else:                                                              # pragma: no cover
        raise AssertionError("a misplaced edit must not be applied silently")


def _records(ref, by_hap, span=None):
    fa = FakeFasta(ref)
    span = span or (1, len(ref))
    return hla_graft._records(fa, "chr6", by_hap, span)


def test_a_substitution_carries_no_anchor_base():
    recs = _records("AAACAAA", {0: [(3, "C", "G")], 1: []})
    assert len(recs) == 1
    r = recs[0]
    assert (r["pos"], r["ref"], r["alts"], r["gt"]) == (4, "C", ["G"], (1, 0))


def test_a_deletion_takes_an_anchor_base():
    recs = _records("AAACGTAAA", {0: [(3, "CGT", "")], 1: []})
    r, = recs
    assert (r["pos"], r["ref"], r["alts"], r["gt"]) == (3, "ACGT", ["A"], (1, 0))


def test_an_insertion_takes_an_anchor_base_so_the_alleles_share_a_first_base():
    recs = _records("AAACAAA", {0: [(3, "", "TT")], 1: []})
    r, = recs
    assert r["ref"][0] == r["alts"][0][0], "an indel's alleles must share their first base"
    assert (r["pos"], r["ref"], r["alts"], r["gt"]) == (3, "AC", ["ATTC"], (1, 0))


def test_both_haplotypes_differing_at_one_site_make_one_biallelic_record():
    recs = _records("AAACAAA", {0: [(3, "C", "G")], 1: [(3, "C", "T")]})
    r, = recs
    assert (r["pos"], r["ref"], r["alts"], r["gt"]) == (4, "C", ["G", "T"], (1, 2))


def test_the_same_change_on_both_haplotypes_is_homozygous():
    recs = _records("AAACAAA", {0: [(3, "C", "G")], 1: [(3, "C", "G")]})
    r, = recs
    assert (r["alts"], r["gt"]) == (["G"], (1, 1))


def test_overlapping_changes_merge_into_one_record_rather_than_two_at_one_position():
    # hap0 deletes CGT, hap1 substitutes inside it. Emitted separately these would be two records at
    # overlapping spans whose REFs disagree, and `Germline.apply` would skip one of them.
    recs = _records("AAACGTAAA", {0: [(3, "CGT", "")], 1: [(4, "G", "A")]})
    assert len(recs) == 1
    r = recs[0]
    assert (r["pos"], r["ref"], r["alts"], r["gt"]) == (3, "ACGT", ["A", "ACAT"], (1, 2))


def test_a_site_where_both_haplotypes_match_the_reference_is_not_written():
    assert _records("AAACAAA", {0: [], 1: []}) == []


def test_adjacent_substitutions_become_one_mnp_without_an_anchor():
    recs = _records("AAACGAAA", {0: [(3, "C", "T"), (4, "G", "A")], 1: []})
    r, = recs
    assert (r["pos"], r["ref"], r["alts"]) == (4, "CG", ["TA"])


def test_a_length_change_at_the_very_start_of_the_span_is_refused_not_anchored_outside_it():
    try:
        _records("CGTAAA", {0: [(0, "CGT", "")], 1: []}, span=(1, 6))
    except SystemExit as e:
        assert "no anchor base" in str(e)
    else:                                                              # pragma: no cover
        raise AssertionError("anchoring outside the grafted span would write a record the graft "
                             "does not own")


if __name__ == "__main__":                                             # pragma: no cover
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"{len(fns)} passed")
