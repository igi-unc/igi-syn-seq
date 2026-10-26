"""Regression tests for truth classes that were described but not present in the reads.

Each of these exists because the catalog promised something the libraries did not contain, and nothing
in the suite asked. They are cheap, need no cluster and no data, and run either as
`python3 catalog/tests/test_absent_classes.py` or under pytest.

  * translocations  -- `plan_translocation` had no caller anywhere in the repository, and `sv_rows` skips
                       TRA, so 8 designed translocations per dataset reached the short-read libraries as
                       junction contigs and the long-read libraries not at all.
  * episomes        -- a viral episome has no human coordinate, so `viral_junctions` put the viral contig
                       in `chrom`; both the per-task junction filter and `junctions_in`, which is keyed by
                       (human chromosome, clone, haplotype), then looked for it under a chromosome that is
                       never queried. EBV at 5 copies/cell measured 0.07 % breadth in delivered ds-02.
  * integrated body -- an integration emitted only its two host junctions, so 1.2 kb of a 7.9 kb HPV16
                       genome reached the reads and E6/E7 had no coverage at any depth.
"""
import os
import sys
import tempfile

try:
    import pytest
except ImportError:
    pytest = None

CATALOG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, CATALOG)

from igi_catalog import rearrange                                      # noqa: E402
from igi_catalog.derived import cn_slices, placeable_tra               # noqa: E402
from igi_catalog.junctions import viral_payloads, FLANK                # noqa: E402

Segment = rearrange.Segment


# --- translocations -------------------------------------------------------------------------------

def test_translocation_derivatives_partition_both_chromosomes():
    """der(A) and der(B) together contain every base of A and of B exactly once.

    This is the invariant that lets each derivative be emitted by the task owning its 5' chromosome
    without double-counting: A's head and B's tail here, B's head and A's tail there.
    """
    len_a, len_b, pos_a, pos_b = 1000, 800, 400, 300
    plan_a, ok_a = rearrange.plan_with_translocation(
        [Segment("chrA", 1, len_a)], "chrA", pos_a, "chrB", pos_b + 1, len_b, "EV1")
    plan_b, ok_b = rearrange.plan_with_translocation(
        [Segment("chrB", 1, len_b)], "chrB", pos_b, "chrA", pos_a + 1, len_a, "EV1")
    assert ok_a and ok_b

    covered = {}
    for seg in plan_a + plan_b:
        for base in range(seg.start, seg.end + 1):
            key = (seg.chrom, base)
            covered[key] = covered.get(key, 0) + 1

    assert set(covered.values()) == {1}, "a base is emitted twice or not at all"
    assert len([k for k in covered if k[0] == "chrA"]) == len_a
    assert len([k for k in covered if k[0] == "chrB"]) == len_b


def test_translocation_continues_into_its_partner():
    plan, ok = rearrange.plan_with_translocation(
        [Segment("chr1", 1, 1000)], "chr1", 400, "chr18", 301, 800, "EV1")
    assert ok
    assert [(s.chrom, s.start, s.end) for s in plan] == [("chr1", 1, 400), ("chr18", 301, 800)]


def test_translocation_refuses_to_clip_a_rearranged_segment():
    """A breakpoint inside an inversion or a duplicate is left unplaced rather than silently clipped."""
    segs = [Segment("chr1", 1, 200), Segment("chr1", 201, 600, "-", "inv", "SV9")]
    plan, ok = rearrange.plan_with_translocation(segs, "chr1", 400, "chr18", 301, 800, "EV1")
    assert not ok
    assert plan == segs, "the original plan must survive a refusal unchanged"


def test_translocation_keeps_intra_chromosomal_events_before_the_breakpoint():
    segs = [Segment("chr1", 1, 100), Segment("chr1", 101, 200, "-", "inv", "SV9"),
            Segment("chr1", 201, 1000)]
    plan, ok = rearrange.plan_with_translocation(segs, "chr1", 500, "chr18", 301, 800, "EV1")
    assert ok
    kinds = [(s.chrom, s.kind) for s in plan]
    assert ("chr1", "inv") in kinds, "an inversion before the cut must be retained"
    assert plan[-1].chrom == "chr18"


def test_both_partners_claim_the_same_translocation():
    """Each chromosome builds the derivative that starts with itself, so both get built.

    Returning only the events whose catalogued `chrom` matched built der(A) and never der(B): B's task
    still produced an intact B, so B's tail was emitted twice and A's was lost.
    """
    tra = [{"event_id": "EV1", "chrom": "chr1", "start": 400,
            "end_chrom": "chr18", "end": 300, "clone": "T", "haplotype": 0}]
    a = placeable_tra(tra, "chr1")
    b = placeable_tra(tra, "chr18")
    assert a[0]["event_id"] == "EV1" and b[0]["event_id"] == "EV1"
    assert (a[1], a[2], a[3]) == (400, "chr18", 301), "der(chr1) = chr1 head + chr18 tail"
    assert (b[1], b[2], b[3]) == (300, "chr1", 401), "der(chr18) = chr18 head + chr1 tail"


def test_chained_translocations_are_skipped_not_half_built():
    """Two translocations meeting on one chromosome of one copy make a chain, not two simple cuts."""
    tras = [{"event_id": "EV1", "chrom": "chr17", "start": 2_416_216,
             "end_chrom": "chr8", "end": 133_969_364, "clone": "T", "haplotype": 0},
            {"event_id": "EV2", "chrom": "chr17", "start": 41_895_180,
             "end_chrom": "chr7", "end": 7_496_665, "clone": "T", "haplotype": 0}]
    chosen, _cut, _p, _ps, skipped = placeable_tra(tras, "chr17")
    assert chosen is None
    assert {t["event_id"] for t in skipped} == {"EV1", "EV2"}
    # and the partners must agree, or one side would build a derivative the other does not
    assert placeable_tra(tras, "chr8")[0] is None
    assert placeable_tra(tras, "chr7")[0] is None


def test_an_unrelated_translocation_on_the_same_copy_still_places():
    tras = [{"event_id": "EV1", "chrom": "chr17", "start": 100,
             "end_chrom": "chr8", "end": 200, "clone": "T", "haplotype": 0},
            {"event_id": "EV2", "chrom": "chr17", "start": 300,
             "end_chrom": "chr7", "end": 400, "clone": "T", "haplotype": 0},
            {"event_id": "EV3", "chrom": "chr1", "start": 500,
             "end_chrom": "chr2", "end": 600, "clone": "T", "haplotype": 0}]
    assert placeable_tra(tras, "chr1")[0]["event_id"] == "EV3"
    assert placeable_tra(tras, "chr2")[0]["event_id"] == "EV3"


# --- copy number across a translocation -----------------------------------------------------------

def test_cn_slices_carry_their_own_chromosome():
    """A partner tail must be weighted at the partner's coordinates, not this chromosome's.

    cn_slices used to return (start, end, probe) and the caller asked `source_plan` about `chrom`, so a
    translocated tail would have been given the copy number of a position on the wrong chromosome.
    """
    cmap = [(1, 400, "chr1", 1, 400, "+", "ref", "EV1"),
            (401, 900, "chr18", 301, 800, "+", "ref", "EV1")]
    cuts = {"chr1": [200], "chr18": []}
    out = cn_slices(cmap, lambda c: cuts.get(c, []), min_bp=1)
    chroms = {s[2] for s in out}
    assert chroms == {"chr1", "chr18"}
    for ds, de, chrom, probe in out:
        lo, hi = (1, 400) if chrom == "chr1" else (301, 800)
        assert lo <= probe <= hi, f"{chrom} probe {probe} is outside that chromosome's segment"


def test_cn_slices_split_on_the_right_chromosomes_boundaries():
    cmap = [(1, 400, "chr1", 1, 400, "+", "ref", "EV1"),
            (401, 900, "chr18", 301, 800, "+", "ref", "EV1")]
    out = cn_slices(cmap, lambda c: [200] if c == "chr1" else [], min_bp=1)
    assert len([s for s in out if s[2] == "chr1"]) == 2, "chr1's boundary must cut chr1"
    assert len([s for s in out if s[2] == "chr18"]) == 1, "chr1's boundary must not cut chr18"


def test_a_chromosome_with_no_copy_number_change_is_one_slice():
    """Segment boundaries are not copy-number boundaries, and must not become simulator invocations.

    The builder used to emit exactly one slice for a chromosome with no CNA. Keying slices off the
    coordinate map without coalescing would have given chr5's chromothripsis region 60-odd slices, each
    with its own FASTA and its own simulator run, all at identical coverage.
    """
    cmap = [(1, 1000, "chr1", 1, 1000, "+", "ref", ""),
            (1001, 1500, "chr1", 2000, 2499, "+", "dup", "SV1"),
            (1501, 2000, "chr1", 2500, 2999, "-", "inv", "SV2")]
    out = cn_slices(cmap, lambda c: [], min_bp=1)
    assert len(out) == 1, f"expected one slice with no cuts, got {len(out)}"
    assert out[0][0] == 1 and out[0][1] == 2000


def test_coalescing_never_merges_across_a_translocation():
    cmap = [(1, 400, "chr1", 1, 400, "+", "ref", "EV1"),
            (401, 900, "chr18", 301, 800, "+", "ref", "EV1")]
    out = cn_slices(cmap, lambda c: [], min_bp=1)
    assert len(out) == 2
    assert [s[2] for s in out] == ["chr1", "chr18"]


# --- viral payloads -------------------------------------------------------------------------------

def _viral_fasta(tmp, seqs):
    import pysam
    path = os.path.join(tmp, "viral.fa")
    with open(path, "w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    pysam.faidx(path)
    return path


EBV = "ACGT" * 5000          # 20 kb stand-in
HPV = "ACGT" * 2000          # 8 kb stand-in

ROWS = [
    {"event_id": "VIR-1", "subclass": "integrated", "virus": "HPV16", "accession": "NC_000001.1",
     "copies_per_cell": "3", "normal_trace_copies": "0", "ccf": "1.0", "clone": "T",
     "chrom": "chr8", "integration_pos": "127400000", "haplotype": "0"},
    {"event_id": "VIR-2", "subclass": "episomal", "virus": "EBV", "accession": "NC_000002.1",
     "copies_per_cell": "5", "normal_trace_copies": "0.05", "ccf": "1.0", "clone": "T",
     "chrom": "", "integration_pos": "", "haplotype": "0"},
]


def test_every_designed_virus_produces_a_payload():
    with tempfile.TemporaryDirectory() as tmp:
        fa = _viral_fasta(tmp, {"NC_000001.1": HPV, "NC_000002.1": EBV})
        out = viral_payloads(fa, ROWS, tumor=True)
        assert {p["id"] for p in out} == {"VIR-1:genome", "VIR-2:genome"}, \
            "an episome has no human locus and was dropped by the junction path entirely"


def test_payload_is_the_whole_genome_not_two_flanks():
    with tempfile.TemporaryDirectory() as tmp:
        fa = _viral_fasta(tmp, {"NC_000001.1": HPV, "NC_000002.1": EBV})
        by_id = {p["id"]: p for p in viral_payloads(fa, ROWS, tumor=True)}
        assert len(by_id["VIR-1:genome"]["sequence"]) == len(HPV)
        assert len(by_id["VIR-1:genome"]["sequence"]) > 2 * FLANK, \
            "the integrated virus reached the reads as its two 600 bp ends and nothing between"
        assert len(by_id["VIR-2:genome"]["sequence"]) == len(EBV)


def test_payload_name_is_safe_on_a_command_line():
    """The viral reference uses NCBI headers; the name reaches a path the simulator is run with.

    `gi|310698439|ref|NC_001526.2|` became a FASTA filename and a chunk tag that were interpolated into
    badread's command line, where a pipe is a shell pipe. All three viruses produced nothing -- the
    171 kb EBV included, so it was never about the reference being shorter than a read.
    """
    import re as _re
    with tempfile.TemporaryDirectory() as tmp:
        fa = _viral_fasta(tmp, {"gi|1|ref|NC_000001.1|": HPV, "gi|2|ref|NC_000002.1|": EBV})
        rows = [dict(r) for r in ROWS]
        out = viral_payloads(fa, rows, tumor=True)
        assert out, "the NCBI-style contig must still be matched by accession"
        for pl in out:
            assert "slug" in pl, "a payload must carry a name safe to put in a path"
            assert _re.fullmatch(r"[A-Za-z0-9._-]+", pl["slug"]), f"unsafe slug: {pl['slug']!r}"
            for bad in "|;&$ ><*?()":
                assert bad not in pl["slug"]


def test_payload_weight_follows_purity_ccf_and_copy_number():
    """Weight is in `source_plan` units, where one haplotype copy of a clone is purity * frac."""
    with tempfile.TemporaryDirectory() as tmp:
        fa = _viral_fasta(tmp, {"NC_000001.1": HPV, "NC_000002.1": EBV})
        by_id = {p["id"]: p for p in viral_payloads(fa, ROWS, tumor=True)}
        purity = 0.45
        w = purity * by_id["VIR-2:genome"]["ccf"] * by_id["VIR-2:genome"]["copies"]
        assert abs(w - 0.45 * 1.0 * 5) < 1e-9
        assert by_id["VIR-2:genome"]["tumor_weighted"] is True


def test_normal_library_gets_only_the_trace_carry_over():
    with tempfile.TemporaryDirectory() as tmp:
        fa = _viral_fasta(tmp, {"NC_000001.1": HPV, "NC_000002.1": EBV})
        out = viral_payloads(fa, ROWS, tumor=False)
        assert {p["id"] for p in out} == {"VIR-2:genome"}, \
            "only the virus with normal_trace_copies > 0 belongs in a normal library"
        assert out[0]["copies"] == 0.05
        assert out[0]["tumor_weighted"] is False


def test_missing_viral_reference_is_not_a_silent_empty_library():
    """A missing reference still returns nothing -- but the caller logs it, so record the contract."""
    assert viral_payloads(None, ROWS, tumor=True) == []
    assert viral_payloads("/nonexistent/viral.fa", ROWS, tumor=True) == []




# --- mobile element insertions --------------------------------------------------------------------

class _FakeGenome:
    """A genome whose donor loci hold recognisable elements, so tests need no reference FASTA."""
    def __init__(self):
        from igi_catalog import mei
        self.lengths = {}
        self._seqs = {}
        for fam, (chrom, start1, end1) in mei.DONORS.items():
            body = {"L1HS": "GACT", "AluY": "GGCC", "SVA_E": "CCCT"}[fam] * 4000
            self.lengths[chrom] = end1 + 100_000
            self._seqs[chrom] = ("T" * (start1 - 1)) + body[:end1 - start1 + 1] + ("A" * 30) + ("G" * 100_000)
        self.lengths["chr9"] = 500_000
        self._seqs["chr9"] = "CATG" * 125_000

    def seq(self, chrom, start, end):
        return self._seqs.get(chrom, "")[start:end]


def test_insertion_is_exactly_the_catalogued_length():
    """The plan sizes the insert from the token without calling the builder; they must agree exactly."""
    from igi_catalog import mei
    import random
    g = _FakeGenome()
    for token in ("<L1HS:6033>", "<AluY:352>", "<SVA_E:2013>", "<AluY:300>"):
        _fam, total = mei.family_of(token)
        got = mei.insertion(g, token, "chr9", 200_000, random.Random(1))
        assert len(got) == total, f"{token}: realised {len(got)}, catalogued {total}"


def test_insertion_carries_element_sequence_not_reference():
    from igi_catalog import mei
    import random
    g = _FakeGenome()
    host = g.seq("chr9", 199_000, 201_000)
    got = mei.insertion(g, "<L1HS:6033>", "chr9", 200_000, random.Random(7))
    assert got not in host and host[:200] not in got, "the insertion must not be host reference"
    assert "GACT" in got or "AGTC" in got, "the L1HS donor body should be present on one strand or the other"


def test_insertion_ends_in_a_target_site_duplication():
    from igi_catalog import mei
    import random
    g = _FakeGenome()
    pos = 200_000
    got = mei.insertion(g, "<AluY:352>", "chr9", pos, random.Random(3))
    host_before = g.seq("chr9", pos - 40, pos)
    assert any(got.endswith(host_before[-n:]) for n in range(mei.TSD_RANGE[0], mei.TSD_RANGE[1] + 1)), \
        "the element must be followed by a copy of the host bases it landed on"


def test_the_same_event_gives_the_same_element_on_every_copy():
    """A truncal insertion is carried by every descendant clone and must be identical on each.

    With one shared RNG each copy drew a different sequence for the same event, so the reads supporting
    a single catalogued insertion never stacked.
    """
    from igi_catalog.derived import mei_sequence_factory
    g = _FakeGenome()
    svs = [{"event_id": "EV1", "chrom": "chr9", "start": 200_000,
            "inserted_sequence": "<L1HS:6033>", "svtype": "INS"}]
    a = mei_sequence_factory(g, svs, "DS")("EV1")
    b = mei_sequence_factory(g, svs, "DS")("EV1")
    assert a == b and len(a) == 6033


def test_short_and_long_read_paths_agree_on_the_element():
    """The exome/WGS contig and the derived chromosome must carry the SAME insertion."""
    from igi_catalog.derived import mei_sequence_factory
    from igi_catalog.junctions import sv_junctions, FLANK
    g = _FakeGenome()
    row = {"event_id": "EV1", "dataset": "DS", "svtype": "INS", "chrom": "chr9",
           "start": 200_000, "end": 200_400, "inserted_sequence": "<L1HS:6033>",
           "clone": "T", "haplotype": "0"}
    long_read = mei_sequence_factory(g, [row], "DS")("EV1")
    jns = {j["id"]: j for j in sv_junctions(g, [row], flank=FLANK)}
    assert set(jns) == {"EV1:5p", "EV1:3p"}, "an insertion has two boundaries, not one"
    assert jns["EV1:5p"]["sequence"].endswith(long_read[:FLANK])
    assert jns["EV1:3p"]["sequence"].startswith(long_read[-FLANK:])


def test_insertion_junction_is_not_plain_reference():
    """The regression that made all 16 planted MEIs invisible."""
    from igi_catalog.junctions import sv_junctions, FLANK
    g = _FakeGenome()
    row = {"event_id": "EV1", "dataset": "DS", "svtype": "INS", "chrom": "chr9",
           "start": 200_000, "end": 200_400, "inserted_sequence": "<AluY:352>",
           "clone": "T", "haplotype": "0"}
    ref = g.seq("chr9", 200_000 - FLANK, 200_000 + FLANK)
    for j in sv_junctions(g, [row], flank=FLANK):
        assert j["sequence"] not in ref, f"{j['id']} is reference sequence; a caller sees nothing"


def test_insert_segment_knows_its_own_derived_length():
    """`Segment.length()` sized an insertion at 1 base while `realise` appended the whole element."""
    segs, applied, _dropped = rearrange.plan_chromosome(
        "chr9", 500_000,
        [{"event_id": "EV1", "svtype": "INS", "start": 200_000, "end": 200_400,
          "inserted_sequence": "<L1HS:6033>", "clone": "T", "haplotype": 0}])
    ins = [s for s in segs if s.kind == "insert"]
    assert len(ins) == 1 and ins[0].length() == 6033
    cmap = rearrange.coordinate_map(segs)
    ins_rows = [r for r in cmap if r[6] == "insert"]
    assert ins_rows[0][1] - ins_rows[0][0] + 1 == 6033


def test_inserted_bases_are_not_dropped_from_the_slices():
    """cn_slices must cover every derived base; an insert has start == end and no reference span."""
    segs, _a, _d = rearrange.plan_chromosome(
        "chr9", 500_000,
        [{"event_id": "EV1", "svtype": "INS", "start": 200_000, "end": 200_400,
          "inserted_sequence": "<L1HS:6033>", "clone": "T", "haplotype": 0}])
    cmap = rearrange.coordinate_map(segs)
    total = sum(s.length() for s in segs)
    out = cn_slices(cmap, lambda c: [], min_bp=1)
    covered = sum(de - ds + 1 for ds, de, _c, _p in out)
    assert covered == total, f"slices cover {covered} of {total} derived bases"


def _main():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"  PASS  {fn.__name__}")
        except Exception as exc:                       # noqa: BLE001
            failed += 1
            print(f"  FAIL  {fn.__name__}: {exc}")
    print(f"{len(fns) - failed}/{len(fns)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(_main())
