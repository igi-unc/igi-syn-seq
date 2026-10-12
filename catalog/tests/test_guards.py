"""Regression tests for the delivery guards added on 2026-10-08.

Each of these exists because a check that did not do what it claimed let a defect through, or because a
setting's effect was not what its value suggested. They are cheap, they need no cluster and no data, and
they run either as `python3 catalog/tests/test_guards.py` or under pytest.

  * check_alphabet          -- the release had checks on length, accuracy, names, truth-map resolution and
                               format, and none that asked whether the delivered bases were bases
                               (docs/build-reference.md §15). Its verdict now also depends on the consumer
                               (§17.3), which is a rule worth pinning down.
  * sort_buffer_mb          -- a --mem request never bounded GNU sort's buffer (§16).
  * normalize_bases         -- reference ambiguity codes reached delivered reads (§17).
  * rescale_quality         -- the quality correction is an integer decibel shift, and asking for 5.15x
                               delivered 3.70x (§14.3).
"""
import glob
import gzip
import json
import importlib.util
import math
import os
import shutil
import subprocess
import sys
try:
    import pytest
except ImportError:
    pytest = None
import tempfile
import time
from array import array

HERE = os.path.dirname(os.path.abspath(__file__))
CATALOG = os.path.dirname(HERE)
sys.path.insert(0, CATALOG)

from igi_catalog.genome import normalize_bases, revcomp          # noqa: E402
from igi_catalog.pacbio_bam import rescale_quality               # noqa: E402
from igi_catalog.readnames import sort_buffer_mb                 # noqa: E402


def _qa():
    """qa_release.py is a script, not a package module, so load it by path."""
    spec = importlib.util.spec_from_file_location("qa", os.path.join(CATALOG, "qa_release.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------------------------- the alphabet check

def _write_fq(path, seqs):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt") as fh:
        for i, s in enumerate(seqs):
            fh.write(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n")


# These currency checks compare mtimes, and the tests for them must not depend on the filesystem's
# clock agreeing with time.time(). TMPDIR here is NFS scratch: a file created by the client is stamped
# by the SERVER, so `os.utime(out, (time.time(), ...))` can land BEHIND an input written a moment
# earlier and invert the very ordering under test. Three tests failed that way once the skew grew.
# Pinning both sides to an explicit base removes the clock from the test entirely.
def _pin(path, offset_s, base=None):
    """Set path's mtime to base+offset_s (base defaults to a fixed recent epoch)."""
    b = base if base is not None else 1_760_000_000
    os.utime(path, (b + offset_s, b + offset_s))
    return b + offset_s

def test_check_alphabet_verdicts():
    """Any non-ACGTN character fails, in every arm -- the consumer enumeration was wrong."""
    qa = _qa()
    root = tempfile.mkdtemp(prefix="igi_alphabet_")
    try:
        _write_fq(f"{root}/rna/DS/DS_full_rna_R1.fastq.gz", ["ACGTACGTAC"] * 5)
        _write_fq(f"{root}/rna/DS/DS_bad_rna_R1.fastq.gz", ["ACGT*CGTAC"])
        _write_fq(f"{root}/wes/DS/DS_chr1_tumor_R1.fastq.gz", ["ACGTACGTAC"] * 4 + ["ACGTMCGTAC"])
        _write_fq(f"{root}/ont_wgs/DS/DS_chr3_tumor_ont_wgs.fastq.gz", ["ACGTYCGTAC", "ACGTBCGTAC"])
        _write_fq(f"{root}/tenx_gex/DS/DS-GEX_S1_L001_R1_001.fastq.gz",
                  ["ACGTACGTAC"] * 3 + ["ACGTYCGTAC"])

        # ont_deliverables resolves the real release layout; point it at the fixture instead.
        qa.ont_deliverables = lambda rel: sorted(glob.glob(f"{rel}/ont_wgs/*/*_ont_wgs.fastq.gz"))
        results = []
        qa.check_alphabet(root, results)
        got = {r["file"]: r["verdict"].split(" ")[0] for r in results}

        assert got["DS_full_rna_R1.fastq.gz"] == "PASS", got
        assert got["DS_bad_rna_R1.fastq.gz"] == "FAIL", "a `*` in the sequence must fail"
        # These two were WARN until razers3 aborted on a 'Y' in the normal exome. An ambiguity code is
        # not tolerable anywhere: we do not get to decide what reads a reference dataset.
        assert got["DS_chr1_tumor_R1.fastq.gz"] == "FAIL", "IUPAC in the exome aborts razers3"
        assert got["DS_chr3_tumor_ont_wgs.fastq.gz"] == "FAIL", "IUPAC must fail in every arm"
        assert got["DS-GEX_S1_L001_R1_001.fastq.gz"] == "FAIL", "Cell Ranger refuses the run"
        assert sum(1 for r in results if r["verdict"].startswith("FAIL")) == 4, results
        assert not any(r["verdict"].startswith("WARN") for r in results), \
            "there is no WARN tier any more; a tolerated defect is how this reached a user"
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_alphabet_check_is_complete_not_sampled():
    """The offence must be found when it sits past the sample window, which is where it really sits.

    The exome carries ambiguity codes at about one read per million -- chr10 normal held exactly one `R`
    in 1,860,880 reads -- so the 100,000-read sample the check started with reported `none` on a library
    that aborts razers3. This is that failure in miniature: one offending read behind 150,000 clean ones.
    """
    qa = _qa()
    root = tempfile.mkdtemp(prefix="igi_deep_")
    try:
        _write_fq(f"{root}/wes/DS/DS_chr1_tumor_R1.fastq.gz",
                  ["ACGTACGTAC"] * 150000 + ["ACGTYCGTAC"])
        qa.ont_deliverables = lambda rel: []
        f = f"{root}/wes/DS/DS_chr1_tumor_R1.fastq.gz"

        n, counts, empty = qa._scan_alphabet(f, complete=True)
        assert (n, counts) == (150001, {"Y": 1}), (n, counts)

        ns, cs, _ = qa._scan_alphabet(f, complete=False, limit=400000)
        assert (ns, cs) == (100000, {}), "the sample is expected to miss it -- that is the point"

        results = []
        qa.check_alphabet(root, results)
        assert results[0]["verdict"].startswith("FAIL"), results
        assert "complete" in results[0]["value"], "the report must say the scan was complete"
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_alphabet_check_reads_both_mates():
    """R2 was never scanned, and on a 10x arm R2 is the only mate Cell Ranger aligns."""
    qa = _qa()
    root = tempfile.mkdtemp(prefix="igi_mates_")
    try:
        _write_fq(f"{root}/tenx_gex/DS/DS-GEX_S1_L001_R1_001.fastq.gz", ["ACGTACGTAC"] * 4)
        _write_fq(f"{root}/tenx_gex/DS/DS-GEX_S1_L001_R2_001.fastq.gz", ["ACGTACGTAC", "ACGTYCGTAC"])
        qa.ont_deliverables = lambda rel: []
        results = []
        qa.check_alphabet(root, results)
        got = {r["file"]: r["verdict"].split(" ")[0] for r in results}
        assert "DS-GEX_S1_L001_R2_001.fastq.gz" in got, f"R2 was not examined at all: {got}"
        assert got["DS-GEX_S1_L001_R2_001.fastq.gz"] == "FAIL", got
        assert got["DS-GEX_S1_L001_R1_001.fastq.gz"] == "PASS", got
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_sam_placeholder_is_not_a_symbolic_allele():
    """SEQ == "*" means "no sequence stored", not a base. Counting it as a `*` was a false positive.

    ds-02's chr1to6 Kinnex BAM reported `* x7` with zero `*` embedded in any sequence -- ds-02's VCF is
    SHAPEIT5-normalised and contains no `*` at all, which is what made the report suspicious. The two need
    different fixes, so they get different verdicts: a sequence-less record is a real defect in a
    deliverable, but it is not an alphabet defect.
    """
    qa = _qa()
    lines = "ACGT\n*\nACGTY\n*\nACG*T\n"

    def tally(ph):
        r = subprocess.run(["awk", "-v", f"ph={ph}", qa._TALLY],
                           input=lines, capture_output=True, text=True)
        out = {}
        for line in r.stdout.splitlines():
            f = line.split()
            out[f[0]] = out.get(f[0], []) + [f[1:]]
        return out

    bam = tally(1)                                   # BAM: `*` alone is the placeholder
    assert bam["EMPTY"][0][0] == "2", bam
    chars = {f[0]: int(f[1]) for f in bam.get("CHAR", [])}
    assert chars == {"Y": 1, "*": 1}, f"the embedded * must still count: {chars}"

    fq = tally(0)                                    # FASTQ: there is no placeholder, so all three count
    assert fq["EMPTY"][0][0] == "0", fq
    chars = {f[0]: int(f[1]) for f in fq.get("CHAR", [])}
    assert chars == {"Y": 1, "*": 3}, chars


def test_alphabet_verdict_cause_naming():
    """The verdict names the cause, because the two causes need different fixes."""
    qa = _qa()
    assert qa._alphabet_verdict([]) == "PASS"
    assert "symbolic allele" in qa._alphabet_verdict(["*"], "/r/rna/DS/x_R1.fastq.gz")
    v = qa._alphabet_verdict(["Y"], "/r/ont_wgs/DS/x.fastq.gz")
    assert "IUPAC" in v and "razers3" in v, "the verdict must name the consumer that actually aborted"
    # Mixed: the symbolic allele dominates, because it is a generator defect and the other is not.
    assert qa._alphabet_verdict(["*", "Y"], "/r/ont_wgs/DS/x.fastq.gz").startswith("FAIL")
    assert qa._alphabet_verdict(["Y"], "/r/tenx_tcr/DS/x_R1_001.fastq.gz").startswith("FAIL")


def test_currency_catches_a_deliverable_older_than_its_inputs():
    """A stage that never ran leaves a stale deliverable, and no per-stage mtime guard can see that.

    This is the shape of the real failure: nine arms regenerated and reporting COMPLETED, no merge run
    afterwards, so every library a consumer could open predated the rebuild by up to nine days.
    """
    qa = _qa()
    root = tempfile.mkdtemp(prefix="igi_currency_")
    try:
        os.makedirs(f"{root}/wes/DS"); os.makedirs(f"{root}/merged/DS")
        # The merged library is written first, then the per-chromosome reads are regenerated.
        open(f"{root}/merged/DS/DS_normal_R1.fastq.gz", "w").close()
        open(f"{root}/wes/DS/DS_chr1_normal_R1.fastq.gz", "w").close()
        _pin(f"{root}/merged/DS/DS_normal_R1.fastq.gz", -9 * 86400)
        _pin(f"{root}/wes/DS/DS_chr1_normal_R1.fastq.gz", 0)

        results = []
        qa.check_currency(root, results)
        assert len(results) == 1, results
        assert results[0]["verdict"].startswith("FAIL"), results
        assert "stale by 21" in results[0]["verdict"] or "stale by 2" in results[0]["verdict"], \
            f"the lag should be about 216 h: {results[0]['verdict']}"
        assert "wes/" in results[0]["verdict"], "it must name the input directory to re-run from"

        # Re-running the merge clears it.
        _pin(f"{root}/merged/DS/DS_normal_R1.fastq.gz", 60)
        results = []
        qa.check_currency(root, results)
        assert results[0]["verdict"] == "PASS", results
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_subset_currency_catches_a_chr1to6_older_than_its_source():
    """A chr1to6 subset filtered from a full deliverable that was later regenerated is stale.

    This is the shape of the real failure, four times over on 2026-10-10: the full was rebuilt after a
    generator fix, the subset was not, and the subset went on carrying a defect its own source no longer
    had -- `*` characters in ds-01's rna, ont_bulk_rna and ont_sc_rna, and an unstranded R2 in both
    tenx_gex arms. check_currency cannot see any of it, because these subsets live beside their source
    rather than in a directory assembled from per-chromosome parts.
    """
    qa = _qa()
    root = tempfile.mkdtemp(prefix="igi_subcur_")
    try:
        os.makedirs(f"{root}/rna/DS")
        sub = f"{root}/rna/DS/DS_chr1to6_rna_R1.fastq.gz"
        src = f"{root}/rna/DS/DS_full_rna_R1.fastq.gz"
        # The subset is filtered first, then the full is regenerated by the fix.
        open(sub, "w").close()
        open(src, "w").close()
        _pin(sub, -16 * 3600)
        _pin(src, 0)

        results = []
        qa.check_subset_currency(root, results)
        assert len(results) == 1, results
        assert results[0]["verdict"].startswith("FAIL"), results
        assert "90_subset_chr1to6" in results[0]["verdict"], \
            "it must name the stage to re-run"
        assert "rna" in results[0]["verdict"], "it must name the arm"

        # Re-subsetting clears it.
        _pin(sub, 60)
        results = []
        qa.check_subset_currency(root, results)
        assert results[0]["verdict"] == "PASS", results

        # An intermediate or index touched after the subset must NOT fail it: only the reads the
        # subset is actually filtered from count, or every .pbi rebuild would manufacture a failure.
        for nuisance in ("DS_full_rna_transcripts.tsv", "DS_full_kinnex_bulk_hifi_reads.bam",
                         "DS_full_rna_R1.fastq.gz.pbi"):
            open(f"{root}/rna/DS/{nuisance}", "w").close()
            _pin(f"{root}/rna/DS/{nuisance}", 3600)  # later than the subset, and must not matter
        results = []
        qa.check_subset_currency(root, results)
        assert results[0]["verdict"] == "PASS", \
            f"a later index/intermediate is not staleness: {results}"

        # And an arm with no subset at all (tenx_tcr) asserts nothing rather than failing.
        os.makedirs(f"{root}/tenx_gex/DS")
        open(f"{root}/tenx_gex/DS/DS-GEX_S1_L001_R2_001.fastq.gz", "w").close()
        results = []
        qa.check_subset_currency(root, results)
        assert all("tenx_gex" not in r["file"] for r in results), \
            f"no subset present means nothing to assert: {results}"
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_ont_read_map_must_be_current_not_merely_present():
    """A truth map older than the reads beside it is a wrong answer key, and presence does not catch it.

    Shipped on 2026-10-10: ont_wgs_merged held reads merged that morning next to a read map from Oct 5,
    because the merge concatenated the FASTQ and never rebuilt the map. The check of the day asked only
    whether the file existed and exceeded 1 kB, which a two-day-old map satisfies. Separately,
    ont_wgs_chr1to6 shipped with NO map at all and no ONT check looked at it, because ont_deliverables
    omitted that release.
    """
    qa = _qa()
    root = tempfile.mkdtemp(prefix="igi_ontmap_")
    try:
        os.makedirs(f"{root}/ont_wgs_merged/DS")
        fq = f"{root}/ont_wgs_merged/DS/DS_tumor_ont_wgs.fastq.gz"
        mp = f"{root}/ont_wgs_merged/DS/DS_tumor_ont_wgs_read_map.tsv.gz"
        # A readable FASTQ head is not needed for the map check, but the header check iterates first.
        import gzip as _gz
        with _gz.open(fq, "wt") as fh:
            for i in range(4):
                fh.write(f"@read{i}\nACGT\n+\nIIII\n")
        with open(mp, "wb") as fh:
            fh.write(b"x" * 2000)

        def verdicts():
            res = []
            qa.check_ont_format(root, res)
            return {r["check"]: r["verdict"] for r in res}

        _pin(mp, -48 * 3600)
        _pin(fq, 0)
        v = verdicts()
        assert v["ont read map present and current"].startswith("FAIL"), v
        assert "not from the same run" in v["ont read map present and current"], v

        # Re-merging the map clears it.
        _pin(mp, 60)
        assert verdicts()["ont read map present and current"] == "PASS", verdicts()

        # CONTENT BEATS MTIME. The build records its read count, so comparing it to the map's row
        # count says directly whether the map describes these reads -- which mtimes cannot, because
        # they cannot tell a rewritten envelope from a stale map. ds-02's ONT bulk RNA map is 8.8 h
        # older than its reads and is CORRECT: 20,070,655 rows for 20,070,655 reads, every sampled
        # name present. The mtime rule called that a failure; three ONT RNA arms failed the same way.
        # Enough rows that the file clears the 1 kB truncation guard, which fires before this test.
        import uuid as _uuid
        names = [_uuid.uuid4().hex for _ in range(2000)]
        with _gz.open(mp, "wt") as fh:
            fh.write("read_name\tsource\n")
            for nm in names:
                fh.write(f"{nm}\tchr1_T_hap0_all_s000_chunk0,+strand,1-2000\n")
        json.dump({"reads": len(names), "runtime_s": 10},
                  open(fq.replace(".fastq.gz", ".json"), "w"))
        _pin(mp, -48 * 3600)          # far older than the reads, and still correct
        _pin(fq, 0)
        assert verdicts()["ont read map present and current"] == "PASS", \
            f"a map whose rows match the read count describes those reads: {verdicts()}"

        # A row count that disagrees is the real defect, whatever the mtimes say.
        json.dump({"reads": len(names) + 7, "runtime_s": 10},
                  open(fq.replace(".fastq.gz", ".json"), "w"))
        _pin(mp, 60)                  # NEWER than the reads, so mtime alone would pass it
        v = verdicts()["ont read map present and current"]
        assert v.startswith("FAIL") and "does not describe these reads" in v, v

        # With no recorded count there is nothing to compare, so mtime ordering is the only signal
        # left and still applies.
        os.remove(fq.replace(".fastq.gz", ".json"))
        _pin(mp, -48 * 3600)
        _pin(fq, 0)
        v = verdicts()["ont read map present and current"]
        assert v.startswith("FAIL"), v

        # Absent is still a failure, and distinguishable from stale.
        os.remove(mp)
        assert "absent" in verdicts()["ont read map present and current"], verdicts()

        # And the chr1to6 ONT release must be in scope at all.
        os.makedirs(f"{root}/ont_wgs_chr1to6/DS")
        sub = f"{root}/ont_wgs_chr1to6/DS/DS_chr1to6_tumor_ont_wgs.fastq.gz"
        with _gz.open(sub, "wt") as fh:
            fh.write("@r\nACGT\n+\nIIII\n")
        assert any("chr1to6" in f for f in qa.ont_deliverables(root)), \
            "ont_wgs_chr1to6 is shipped and must be checked"
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_drop_empty_records_removes_sequenceless_segments():
    """A skera segment of length zero must not be delivered, and the header must survive the rewrite.

    skera emits these when two adapter alignments overlap -- the release scan found 11 in ds-01's
    chr1to6 Kinnex library and 7 in ds-02's. The @RG/PU assertion is not incidental: skera names each
    segment from PU, and a previous defect produced `/1/ccs/17_2371` across every delivered Kinnex
    segment when PU was missing, so a filter that quietly dropped the header would reintroduce it.
    """
    pysam = pytest.importorskip("pysam") if pytest else __import__("pysam")
    qa_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    sys.path.insert(0, qa_dir)
    from igi_catalog.pacbio_bam import drop_empty_records

    root = tempfile.mkdtemp(prefix="igi_empty_")
    try:
        sam = os.path.join(root, "mini.sam")
        with open(sam, "w") as fh:
            fh.write("@HD\tVN:1.6\tSO:unknown\n")
            fh.write("@RG\tID:synthetic\tPL:PACBIO\tPU:m84000_260101_000000_s0\n")
            for name, seq in (("1/ccs/0_10", "ACGTACGTAC"), ("1/ccs/27_27", "*"),
                              ("1/ccs/44_52", "ACGTACGT"), ("2/ccs/0_0", "*"),
                              ("2/ccs/17_23", "ACGTAC")):
                q = "*" if seq == "*" else "I" * len(seq)
                fh.write(f"m84000_260101_000000_s0/{name}\t4\t*\t0\t255\t*\t*\t0\t0\t"
                         f"{seq}\t{q}\tRG:Z:synthetic\n")
        bam = os.path.join(root, "mini.bam")
        with pysam.AlignmentFile(sam, "r", check_sq=False) as src:
            with pysam.AlignmentFile(bam, "wb", header=src.header) as out:
                for r in src.fetch(until_eof=True):
                    out.write(r)

        assert drop_empty_records(bam) == 2

        with pysam.AlignmentFile(bam, "rb", check_sq=False) as f:
            recs = list(f.fetch(until_eof=True))
            rg = f.header.to_dict().get("RG")
        assert [r.query_name.split("/")[-1] for r in recs] == ["0_10", "44_52", "17_23"], \
            "surviving records must keep their order"
        assert all(r.query_sequence for r in recs)
        assert rg and rg[0].get("PU") == "m84000_260101_000000_s0", \
            "the rewrite must preserve @RG/PU, which skera's segment naming depends on"

        # A BAM with nothing to drop must be LEFT ALONE, not rewritten with identical content.
        # kinnex_sc measured zero empties on both datasets, so the unconditional rewrite copied a
        # 21 GB BAM and rebuilt a 105 MB .pbi to remove nothing. Byte-identity is the assertion:
        # mtime alone would also catch it, but content equality says the file was never touched.
        before = open(bam, "rb").read()
        mtime_before = os.path.getmtime(bam)
        assert drop_empty_records(bam) == 0, "a clean BAM has nothing to drop"
        assert open(bam, "rb").read() == before, \
            "a BAM with no empty records must not be rewritten"
        assert os.path.getmtime(bam) == mtime_before, \
            "leaving the mtime alone keeps check_subset_currency honest"
        assert not os.path.exists(bam + ".tmp"), "no temporary file should be left behind"

        # A rewrite INVALIDATES the .pbi: dropping a record shifts every later file offset, so the
        # index resolves to the wrong records rather than merely looking old. skera writes that .pbi
        # and neither kinnex builder rebuilds it, so the function must not rewrite silently.
        sam2 = os.path.join(root, "mini2.sam")
        with open(sam, "r") as src, open(sam2, "w") as dst:
            dst.write(src.read())
        bam2 = os.path.join(root, "mini2.bam")
        with pysam.AlignmentFile(sam2, "r", check_sq=False) as src:
            with pysam.AlignmentFile(bam2, "wb", header=src.header) as out:
                for r in src.fetch(until_eof=True):
                    out.write(r)
        open(bam2 + ".pbi", "w").close()
        untouched = open(bam2, "rb").read()

        try:
            drop_empty_records(bam2)
        except RuntimeError as exc:
            assert ".pbi" in str(exc) and "offset" in str(exc), str(exc)
            assert "Nothing has been modified" in str(exc), str(exc)
            # and it must be TRUE: the empties are still there for the next call to find
            assert open(bam2, "rb").read() == untouched, \
                "the refusal must come before the rewrite, not after it"
        else:
            raise AssertionError("rewriting next to a .pbi without a pbindex command must raise")

        # Given the command, it is run so the index is rebuilt.
        marker = os.path.join(root, "pbindex.ran")
        n = drop_empty_records(bam2, "sh -c 'echo {args} > %s'" % marker)
        assert n == 2, n
        assert os.path.exists(marker), "pbindex_cmd must be invoked after a rewrite"
        assert bam2 in open(marker).read(), "pbindex must be pointed at the rewritten BAM"
    finally:
        shutil.rmtree(root, ignore_errors=True)


# --------------------------------------------------------------------------- 10x R2 orientation

def _tenx():
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from igi_catalog import tenx
    return tenx


def test_art_strand_recovery_and_orientation():
    """R2 must come out antisense whichever strand ART drew, with the quality reversed alongside."""
    t = _tenx()
    root = tempfile.mkdtemp(prefix="igi_art_")
    try:
        aln = os.path.join(root, "c.aln")
        with open(aln, "w") as fh:
            fh.write("##ART_Illumina\tread_length\t90\n@CM\tart_illumina\n@SQ\tm0\t450\n")
            fh.write(">m0\tm0-1\t12\t+\nACGT\nACGT\n")
            fh.write(">m7\tm7-1\t99\t-\nTTGC\nTTGC\n")
        strands = t.art_read_strands(aln)
        assert strands == {"m0-1": "+", "m7-1": "-"}, strands

        # '+' is the sense window, so it is flipped; '-' is already antisense and is left alone.
        assert t.orient_antisense("AACCGGTT", "ABCDEFGH", "+") == ("AACCGGTT"[::-1].translate(
            str.maketrans("ACGT", "TGCA")), "HGFEDCBA")
        assert t.orient_antisense("AACCGGTT", "ABCDEFGH", "-") == ("AACCGGTT", "ABCDEFGH")

        # A missing ALN must not silently leave the library unstranded: default '+' means "flip", so an
        # absent strand record produces a consistently oriented library rather than a random one.
        assert t.art_read_strands(os.path.join(root, "nope.aln")) == {}
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_art_chunk_index_survives_a_skipped_record():
    """The read id carries the chunk index, so a skipped short window cannot shift the truth map.

    Both builders write the ART input as `>m{i}` over the chunk but SKIP any window shorter than the
    read length, then paired the k-th FASTQ read with chunk[k]. One skip put every later read on the
    wrong barcode, UMI and record.
    """
    t = _tenx()
    assert t.art_chunk_index("m0-1") == 0
    assert t.art_chunk_index("m12-1") == 12
    assert t.art_chunk_index("m7") == 7
    assert t.art_chunk_index("notaread") is None
    # chunk of 3 with the middle window too short: ART emits m0 and m2, never m1.
    emitted = ["m0-1", "m2-1"]
    assert [t.art_chunk_index(r) for r in emitted] == [0, 2], \
        "positional pairing would have mapped these to chunk[0] and chunk[1]"


def test_r2_strand_check_detects_an_unstranded_library():
    """The QA check must fail a library whose molecule-mates disagree, and pass one that agrees."""
    qa = _qa()
    t = _tenx()
    root = tempfile.mkdtemp(prefix="igi_strand_")
    try:
        import gzip as gz
        # A RANDOM window, not a repeat. My first fixture used "ACGT"*40 + "GATC"*40, and both of those
        # repeats are their own reverse complement, so every seed matched in both orientations and the
        # check reported PASS on a deliberately unstranded library. The fixture was wrong, not the check.
        import random as _r
        _rng = _r.Random(20261009)
        window = "".join(_rng.choice("ACGT") for _ in range(200))
        a = window[10:100]
        b = window[30:120]

        def write(kind, reads):
            d = f"{root}/tenx_gex/DS"
            os.makedirs(d, exist_ok=True)
            with gz.open(f"{d}/DS-GEX_S1_L001_R2_001.fastq.gz", "wt") as fh:
                for i, r in enumerate(reads):
                    fh.write(f"@r{i} 2:N:0:1\n{r}\n+\n{'I' * len(r)}\n")
            with gz.open(f"{d}/DS_full_gex_molecules.tsv.gz", "wt") as fh:
                fh.write("read\tbarcode\tumi\trecord\n")
                for i in range(len(reads)):
                    fh.write(f"r{i}\tBC\tUMI\tENST1\n")

        # Consistently oriented: every read antisense to the window.
        write("good", [t.orient_antisense(x, "I" * len(x), "+")[0] for x in (a, b)] * 60)
        res = []
        qa.check_tenx_r2_strand(root, res, limit=1000, min_pairs=10)
        assert res and res[0]["verdict"] == "PASS", res

        # Unstranded: half the reads left sense.
        mixed = []
        for i in range(60):
            mixed.append(a if i % 2 else t.orient_antisense(a, "I" * len(a), "+")[0])
            mixed.append(b)
        write("bad", mixed)
        res = []
        qa.check_tenx_r2_strand(root, res, limit=1000, min_pairs=10)
        assert res and res[0]["verdict"].startswith("FAIL"), res
        assert "unstranded" in res[0]["verdict"]
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_tcr_direction_check_catches_a_wholesale_flip():
    """A library flipped the wrong way passes the consistency check, so direction needs its own test."""
    qa = _qa()
    root = tempfile.mkdtemp(prefix="igi_dir_")
    try:
        import gzip as gz
        import random as _r
        rng = _r.Random(20261009)
        cdr3 = "".join(rng.choice("ACGT") for _ in range(45))
        flank = "".join(rng.choice("ACGT") for _ in range(30))
        sense_read = flank + cdr3[:40] + flank[:20]
        anti_read = qa.revcomp(sense_read)

        def write(reads):
            d = f"{root}/tenx_tcr/DS"
            os.makedirs(d, exist_ok=True)
            with gz.open(f"{d}/DS-TCR_S1_L001_R2_001.fastq.gz", "wt") as fh:
                for i, r in enumerate(reads):
                    fh.write(f"@r{i} 2:N:0:1\n{r}\n+\n{'I' * len(r)}\n")
            with gz.open(f"{d}/DS_full_tcr_molecules.tsv.gz", "wt") as fh:
                fh.write("read\tbarcode\tumi\tchain\tcdr3_aa\tcdr3_nt\tv\tj\tc\n")
                for i in range(len(reads)):
                    fh.write(f"r{i}\tBC\tUMI\tTRB\tCASS\t{cdr3}\tV\tJ\tC\n")

        write([anti_read] * 300)
        res = []
        qa.check_tcr_r2_direction(root, res, min_reads=10)
        assert res and res[0]["verdict"] == "PASS", res

        # Wholesale flip: internally consistent, but the wrong way round.
        write([sense_read] * 300)
        res = []
        qa.check_tcr_r2_direction(root, res, min_reads=10)
        assert res and res[0]["verdict"].startswith("FAIL"), res
        assert "SENSE" in res[0]["verdict"]
    finally:
        shutil.rmtree(root, ignore_errors=True)


# --------------------------------------------------------------------------- alignment QC

def test_alignment_check_judges_and_never_stays_quiet():
    """Bands are enforced, an unstranded library fails, and a MISSING report is a failure not a skip."""
    qa = _qa()
    root = tempfile.mkdtemp(prefix="igi_aln_")
    try:
        # No reports at all: the job not having run must be visible. A quiet check is how nine
        # COMPLETED arms sat behind stale deliverables.
        res = []
        qa.check_alignment(root, res)
        assert len(res) == 1 and res[0]["verdict"].startswith("FAIL"), res
        assert "has not run" in res[0]["verdict"]

        d = f"{root}/align_qc"
        os.makedirs(d)
        good = {"label": "ds01_wes_normal", "kind": "wes", "mapped_fraction": 1.0,
                "properly_paired": 0.99982, "mean_identity": 0.997307,
                "softclip_fraction": 0.00026, "plus_strand_fraction": 0.5}
        json.dump(good, open(f"{d}/a.json", "w"))
        res = []
        qa.check_alignment(root, res)
        assert all(r["verdict"] == "PASS" for r in res), [r for r in res if r["verdict"] != "PASS"]

        # The real defect: R2 drawn from either strand. This is the number Cell Ranger refused.
        # Asserted on TCR, not GEX, and the distinction is the point. TCR reads come from TRA/TRB,
        # whose transcripts sit on one genomic strand, so a stranded library measures ~0.000 and 0.5
        # means unstranded. GEX reads come from genes on BOTH strands, so ~0.5 is what a correctly
        # stranded library measures too -- this test previously used GEX and asserted a failure that
        # the rebuilt, verified-stranded libraries then produced (0.50219 and 0.50826).
        json.dump({"label": "ds01_tcr", "kind": "tenx_tcr", "mapped_fraction": 0.93,
                   "mean_identity": 0.9987, "softclip_fraction": 0.13,
                   "plus_strand_fraction": 0.5},
                  open(f"{d}/b.json", "w"))
        res = []
        qa.check_alignment(root, res)
        v = [r for r in res if r["file"] == "ds01_tcr" and r["check"] == "strand balance"]
        assert v and v[0]["verdict"].startswith("FAIL"), res
        assert "unstranded" in v[0]["verdict"]

        # And GEX must NOT be judged on that figure at all: no band, so no strand row for it.
        json.dump({"label": "ds02_gex", "kind": "tenx_gex", "mapped_fraction": 0.99976,
                   "mean_identity": 0.99495, "softclip_fraction": 0.0854,
                   "plus_strand_fraction": 0.50826},
                  open(f"{d}/d.json", "w"))
        res = []
        qa.check_alignment(root, res)
        assert not [r for r in res if r["file"] == "ds02_gex" and r["check"] == "strand balance"], \
            "a genomic plus fraction carries no strandedness signal for GEX; it must not be banded"
        assert all(r["verdict"] == "PASS" for r in res if r["file"] == "ds02_gex"), \
            [r for r in res if r["file"] == "ds02_gex" and r["verdict"] != "PASS"]

        # An arm with no band must fail rather than pass silently.
        json.dump({"label": "mystery", "kind": "something_new", "mapped_fraction": 1.0},
                  open(f"{d}/c.json", "w"))
        res = []
        qa.check_alignment(root, res)
        v = [r for r in res if r["file"] == "mystery"]
        assert v and v[0]["verdict"].startswith("FAIL"), res
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_consumer_preflight_separates_a_parse_error_from_a_setup_failure():
    """A tool rejecting the data is a release defect; a tool unhappy with its own setup is not."""
    qa = _qa()
    root = tempfile.mkdtemp(prefix="igi_pf_")
    try:
        res = []
        qa.check_consumer_preflight(root, res)
        assert len(res) == 1 and "has not run" in res[0]["verdict"], res

        d = f"{root}/preflight"
        os.makedirs(d)
        json.dump({"label": "ds02_wes_normal", "tool": "razers3", "fraction": 0.05,
                   "exit_code": 0, "ok": True, "parse_error": False},
                  open(f"{d}/a.json", "w"))
        json.dump({"label": "ds01_wes_normal", "tool": "razers3", "fraction": 0.05,
                   "exit_code": 1, "ok": False, "parse_error": True,
                   "signature": "Unexpected character"}, open(f"{d}/b.json", "w"))
        json.dump({"label": "ds01_rna", "tool": "salmon", "fraction": 0.05,
                   "exit_code": 2, "ok": False, "parse_error": False}, open(f"{d}/c.json", "w"))
        res = []
        qa.check_consumer_preflight(root, res)
        by = {r["file"]: r["verdict"] for r in res}
        assert by["ds02_wes_normal"] == "PASS", by
        assert "rejected the data" in by["ds01_wes_normal"], by
        assert "before blaming the reads" in by["ds01_rna"], by
        assert sum(1 for v in by.values() if v.startswith("FAIL")) == 2, by
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_consumer_preflight_names_the_failing_stage_of_a_chain():
    """Bulk RNA is preflighted as a chain, so the verdict has to say WHICH tool spoke.

    LENS quantifies bulk RNA with STAR followed by `salmon quant -t <fa> -l a -a <bam>`, so the
    preflight row runs both. "star_salmon exited 1" sends you to the wrong tool half the time, and in
    the case that prompted this it sent me nowhere at all: the row that was failing was failing
    because the old check ran salmon in a mode LENS never calls, against a path that was not a salmon
    index. The verdict named the data; the fault was in the harness.
    """
    qa = _qa()
    root = tempfile.mkdtemp(prefix="igi_pfc_")
    try:
        d = f"{root}/preflight"
        os.makedirs(d)
        # STAR fine, salmon not -- the shape the old wording could not express
        json.dump({"label": "ds01_rna", "tool": "star_salmon", "fraction": 0.05,
                   "exit_code": 1, "ok": False, "parse_error": False, "stage": "salmon_aln",
                   "stages": [{"stage": "star", "exit_code": 0, "ok": True},
                              {"stage": "salmon_aln", "exit_code": 1, "ok": False}]},
                  open(f"{d}/rna1.json", "w"))
        json.dump({"label": "ds02_rna", "tool": "star_salmon", "fraction": 0.05,
                   "exit_code": 0, "ok": True, "parse_error": False, "stage": "salmon_aln",
                   "stages": [{"stage": "star", "exit_code": 0, "ok": True},
                              {"stage": "salmon_aln", "exit_code": 0, "ok": True}]},
                  open(f"{d}/rna2.json", "w"))
        # a one-stage row must not pick up the chain wording
        json.dump({"label": "ds01_hifi", "tool": "pbmm2", "fraction": 0.05,
                   "exit_code": 1, "ok": False, "parse_error": False, "stage": "pbmm2",
                   "stages": [{"stage": "pbmm2", "exit_code": 1, "ok": False}]},
                  open(f"{d}/hifi.json", "w"))
        res = []
        qa.check_consumer_preflight(root, res)
        by = {r["file"]: r for r in res}

        assert by["ds01_rna"]["verdict"].startswith("FAIL"), by["ds01_rna"]
        assert "stage salmon_aln" in by["ds01_rna"]["verdict"], by["ds01_rna"]
        assert "1/2 stages" in by["ds01_rna"]["value"], by["ds01_rna"]
        assert by["ds02_rna"]["verdict"] == "PASS", by["ds02_rna"]
        assert "2/2 stages" in by["ds02_rna"]["value"], by["ds02_rna"]
        # the single-stage row names the tool and says nothing about stages
        assert by["ds01_hifi"]["verdict"].startswith("FAIL (pbmm2 exited 1"), by["ds01_hifi"]
        assert "stage" not in by["ds01_hifi"]["value"], by["ds01_hifi"]
    finally:
        shutil.rmtree(root, ignore_errors=True)


# --------------------------------------------------------------------------- record structure

def _fq_records(path, recs):
    import gzip as gz
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gz.open(path, "wt") as fh:
        for name, seq, qual in recs:
            fh.write(f"@{name} 1:N:0:1\n{seq}\n+\n{qual}\n")


def test_read_structure_catches_malformed_records():
    """A quality string of the wrong length, a Q0 base, or a non-Phred character must each fail."""
    qa = _qa()
    root = tempfile.mkdtemp(prefix="igi_struct_")
    try:
        qa.ont_deliverables = lambda rel: []
        good = [("r%d" % i, "ACGTACGTAC", "IIIIIIIIII") for i in range(20)]

        _fq_records(f"{root}/wes/DS/DS_a_R1.fastq.gz", good)
        res = []
        qa.check_read_structure(root, res)
        assert all(r["verdict"] == "PASS" for r in res), [r for r in res if r["verdict"] != "PASS"]

        shutil.rmtree(f"{root}/wes")
        _fq_records(f"{root}/wes/DS/DS_b_R1.fastq.gz", good + [("bad", "ACGTACGTAC", "III")])
        res = []
        qa.check_read_structure(root, res)
        v = {r["check"]: r["verdict"] for r in res}
        assert v["sequence and quality lengths agree"].startswith("FAIL"), v

        shutil.rmtree(f"{root}/wes")
        _fq_records(f"{root}/wes/DS/DS_c_R1.fastq.gz", good + [("q0", "ACGTACGTAC", "IIIII!IIII")])
        res = []
        qa.check_read_structure(root, res)
        v = {r["check"]: r["verdict"] for r in res}
        assert v["no Phred 0 bases"].startswith("FAIL"), v

        # chr(127) is above '~' (Q93) and so outside Phred+33. Note what CANNOT be used here: a tab or
        # a space inside the quality string is consumed as a field separator by the scan, so it is
        # reported as a length mismatch rather than as a bad character. Still a FAIL, under a different
        # name -- but a fixture built from a tab tests the wrong path, which is what it did at first.
        shutil.rmtree(f"{root}/wes")
        _fq_records(f"{root}/wes/DS/DS_d_R1.fastq.gz",
                    good + [("bq", "ACGTACGTAC", "IIIII" + chr(127) + "IIII")])
        res = []
        qa.check_read_structure(root, res)
        v = {r["check"]: r["verdict"] for r in res}
        assert v["quality string is Phred+33"].startswith("FAIL"), v
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_read_structure_flags_implausible_gc():
    """GC is a tripwire for gross corruption: a poly-A library is not a human library."""
    qa = _qa()
    root = tempfile.mkdtemp(prefix="igi_gc_")
    try:
        qa.ont_deliverables = lambda rel: []
        _fq_records(f"{root}/wes/DS/DS_gc_R1.fastq.gz",
                    [("r%d" % i, "GCGCGCGCGC", "IIIIIIIIII") for i in range(20)])
        res = []
        qa.check_read_structure(root, res)
        v = {r["check"]: r["verdict"] for r in res}
        assert v["GC fraction is plausible"].startswith("FAIL"), v
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_mate_pairing_catches_desynchronised_mates():
    """Both mates can be individually perfect and still be out of step."""
    qa = _qa()
    root = tempfile.mkdtemp(prefix="igi_mates2_")
    try:
        names = ["r%04d" % i for i in range(200)]
        mk = lambda ns: [(n, "ACGTACGTAC", "IIIIIIIIII") for n in ns]

        _fq_records(f"{root}/wgs_merged/DS/DS_x_R1.fastq.gz", mk(names))
        _fq_records(f"{root}/wgs_merged/DS/DS_x_R2.fastq.gz", mk(names))
        res = []
        qa.check_mate_pairing(root, res)
        assert res and res[0]["verdict"] == "PASS", res

        # one record dropped from R2: every pair after it is shifted
        _fq_records(f"{root}/wgs_merged/DS/DS_x_R2.fastq.gz", mk(names[:100] + names[101:]))
        res = []
        qa.check_mate_pairing(root, res)
        assert res[0]["verdict"].startswith("FAIL"), res
        assert "out of step" in res[0]["verdict"] or "one mate only" in res[0]["verdict"]

        # R2 missing entirely
        os.remove(f"{root}/wgs_merged/DS/DS_x_R2.fastq.gz")
        res = []
        qa.check_mate_pairing(root, res)
        assert res[0]["verdict"].startswith("FAIL") and "unpaired" in res[0]["verdict"], res
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_pbi_must_be_present_and_newer_than_its_bam():
    """A stale .pbi stops pbmm2 as firmly as a bad read, and it has been missing before."""
    qa = _qa()
    root = tempfile.mkdtemp(prefix="igi_pbi_")
    try:
        d = f"{root}/pacbio_merged/DS"
        os.makedirs(d)
        bam = f"{d}/DS_tumor_hifi.bam"
        open(bam, "w").write("x")
        res = []
        qa.check_bam_integrity(root, res, samtools="true")       # stub: quickcheck always succeeds
        v = {r["check"]: r["verdict"] for r in res}
        assert v[".pbi present and current"].startswith("FAIL"), v

        open(bam + ".pbi", "w").write("y")
        _pin(bam, 0)
        _pin(bam + ".pbi", -3600)
        res = []
        qa.check_bam_integrity(root, res, samtools="true")
        v = {r["check"]: r["verdict"] for r in res}
        assert "stale" in v[".pbi present and current"], v

        _pin(bam + ".pbi", 60)
        res = []
        qa.check_bam_integrity(root, res, samtools="true")
        v = {r["check"]: r["verdict"] for r in res}
        assert v[".pbi present and current"] == "PASS", v
    finally:
        shutil.rmtree(root, ignore_errors=True)


# --------------------------------------------------------------------------- junction placement

def _pipeline():
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from igi_catalog import pipeline
    return pipeline


def test_every_junction_in_an_interval_is_handed_out():
    """Two breakpoints in one 5 Mb interval must both be placed, and each only once.

    Returning a single junction lost five designed events across the two datasets -- each absent from
    the library while the truth table claimed it -- because WGS intervals are 5 Mb wide.
    """
    pl = _pipeline()

    class Fake:
        junctions_in = pl.WesBuilder.junctions_in
        def __init__(self):
            self._jn_placed = set()
            self.jn_by_site = {("chr13", "T", 1): [
                (58218657, {"id": "SV-0040", "sequence": "A" * 1200}),
                (61002389, {"id": "SV-0036", "sequence": "C" * 1200}),
                (99000000, {"id": "SV-0099", "sequence": "G" * 1200}),
            ]}

    f = Fake()
    got = f.junctions_in("chr13", 57700001, 62700000, "T", 1)
    assert [j["id"] for _p, j in got] == ["SV-0040", "SV-0036"], got

    # Handed out once: the off-target bands overlap, so a second pass must not double-count.
    assert f.junctions_in("chr13", 57700001, 62700000, "T", 1) == []
    # The out-of-interval one is still available.
    assert [j["id"] for _p, j in f.junctions_in("chr13", 97000001, 102000000, "T", 1)] == ["SV-0099"]
    # Wrong clone or haplotype gets nothing.
    assert f.junctions_in("chr13", 57700001, 62700000, "A", 1) == []


def test_splice_conserves_interval_length():
    """A copy must contribute the same number of bases whether or not it carries a rearrangement.

    The builder used to replace the whole interval with one 1.2 kb contig. Every one of ds-01's 10 chr6
    T_LOH sites sat inside such an interval, and T_LOH was the one tier whose allele fraction failed
    acceptance -- 24 % deviation against a 12 % bound.
    """
    pl = _pipeline()
    start1 = 1_000_000
    seq = "ACGT" * 50_000                      # 200 kb interval
    jns = [(1_050_000, {"id": "J1", "sequence": "A" * 1200}),
           (1_120_000, {"id": "J2", "sequence": "C" * 1200})]

    out, extra = pl.splice_junctions(seq, start1, jns)
    assert extra == ["A" * 1200, "C" * 1200]
    assert len(out) + sum(len(e) for e in extra) == len(seq), \
        "total bases emitted must equal the wild-type interval"
    assert len(out) == len(seq) - 2400

    # No junctions: untouched, and nothing extra.
    out2, extra2 = pl.splice_junctions(seq, start1, [])
    assert out2 == seq and extra2 == []

    # A breakpoint past the end of the edited sequence must not raise or corrupt: germline indels mean
    # the edited interval is not exactly end1-start1+1 long.
    out3, extra3 = pl.splice_junctions("ACGT" * 10, start1, [(start1 + 10_000, {"id": "J", "sequence": "T" * 8})])
    assert len(extra3) == 1 and len(out3) <= 40


# --------------------------------------------------------------------------- the sort buffer

def test_sort_buffer_follows_the_allocation():
    """GNU sort reads physical memory, not the cgroup limit, so the buffer must come from the job."""
    assert sort_buffer_mb({"SLURM_MEM_PER_NODE": "16384"}) == 4096
    assert sort_buffer_mb({"SLURM_MEM_PER_CPU": "2048", "SLURM_CPUS_ON_NODE": "8"}) == 4096
    # Clamped above, so a huge allocation cannot hand sort the whole node.
    assert sort_buffer_mb({"SLURM_MEM_PER_NODE": "983040"}) == 8192
    # Floored below, and outside SLURM, so the step is slow rather than killed.
    assert sort_buffer_mb({"SLURM_MEM_PER_NODE": "512"}) == 256
    assert sort_buffer_mb({}) == 256
    assert sort_buffer_mb({"SLURM_MEM_PER_NODE": "unlimited"}) == 256


# --------------------------------------------------------------------------- base normalisation

def test_normalize_bases_maps_every_ambiguity_code():
    assert normalize_bases("ACGTN") == "ACGTN", "definite bases are untouched"
    assert normalize_bases("RYSWKMBDHV") == "N" * 10
    assert normalize_bases("acgtrykm") == "acgtNNNN"


def test_revcomp_is_safe_once_normalised():
    """revcomp's table maps only ACGTN, so it must never see an ambiguity code."""
    assert revcomp("ACGTRY") == "YRACGT", "documents the untreated behaviour: R and Y pass uncomplemented"
    assert revcomp(normalize_bases("ACGTRY")) == "NNACGT"


# --------------------------------------------------------------------------- the quality rescale

def test_rescale_quality_is_an_exact_integer_decibel_shift():
    """Phred is integer and so is every input, so a constant scale shifts every base by the same dB."""
    q = array("B", range(1, 94))
    for scale, shift in ((10.0, 10), (5.15, 7), (2.0, 3)):
        out = rescale_quality(q, scale)
        assert [max(1, min(93, v - shift)) for v in q] == list(out), scale
    assert list(rescale_quality(q, 1.0)) == list(q), "a scale of 1 is a no-op"


def test_rescale_quality_matches_the_per_base_form():
    """The refactor to a single shift must not change any output; only the constant did."""
    def per_base(qual, scale):
        out = array("B", bytes(len(qual)))
        for i, v in enumerate(qual):
            e = min(0.75, (10.0 ** (-v / 10.0)) * scale)
            out[i] = max(1, min(93, int(round(-10.0 * math.log10(e)))))
        return out

    q = array("B", range(0, 94))
    for scale in (2.0, 3.0, 5.15, 9.0, 10.0, 11.0):
        assert list(per_base(q, scale)) == list(rescale_quality(q, scale)), scale


def test_rescale_quality_saturates_on_the_low_quality_tail():
    """Why the scale is 10.0 and not the measured ratio 5.15: the correction cannot exceed p = 1."""
    # One very bad base carries most of the error mass, and it cannot get five times worse.
    q = array("B", [4] + [40] * 999)
    before = sum(10.0 ** (-v / 10.0) for v in q)
    after = sum(10.0 ** (-v / 10.0) for v in rescale_quality(q, 5.15))
    assert after / before < 5.15, "a multiplicative rescale saturates; 5.15 delivered ~3.7x in practice"


def test_class_presence_can_say_no():
    """The class-presence summary wrote numbers and exited 0 whatever they were.

    It is the only statement in the release that every designed expression class reached the RNA
    library -- the exome acceptance cannot test ERVs, splice isoforms or CTAs at all -- and for its
    whole life it could not fail. An empty manifest gave an empty "rna" block, an empty verdict and a
    job that printed DONE. The three conditions below are the ones that mean the path is broken; how
    many records cleared the expression floor stays a number, because a thin record is not a defect.
    """
    cp = os.path.join(CATALOG, "class_presence.py")
    root = tempfile.mkdtemp(prefix="igi_cp_")
    try:
        ds = "DS"
        cd = os.path.join(root, "output")
        os.makedirs(cd)
        def tsv(name, header, rows):
            with open(os.path.join(cd, name), "w") as fh:
                fh.write("\t".join(header) + "\n")
                for r in rows:
                    fh.write("\t".join(str(x) for x in r) + "\n")
        tsv(f"{ds}.expressed.tsv", ["event_id", "class", "chrom"],
            [(f"e{i}", "erv", "chr1") for i in range(3)] + [(f"c{i}", "cta", "chr1") for i in range(2)])
        tsv(f"{ds}.fusions.tsv", ["event_id", "chrom_5p"], [("f0", "chr1")])

        man = os.path.join(root, "manifest.tsv")
        def manifest(records):
            with open(man, "w") as fh:
                fh.write("transcript_id\tsource\tobserved_pairs\n")
                for i, (src, pairs) in enumerate(records):
                    fh.write(f"t{i}\t{src}\t{pairs}\n")

        def run():
            r = subprocess.run([sys.executable, cp, "--dataset", ds, "--catalog-dir", cd,
                                "--rna-manifest", man, "--out", os.path.join(root, "out.json")],
                               capture_output=True, text=True)
            return r.returncode, (r.stdout + r.stderr)

        full = [("reference", 10), ("erv", 4), ("cta", 7), ("fusion", 2)]
        manifest(full)
        rc, out = run()
        assert rc == 0, out
        assert json.load(open(os.path.join(root, "out.json")))["verdict"][
            "rna_classes_with_every_record_read"] == ["cta", "erv", "fusion", "reference"], out

        # a reference transcript below the expression floor is a number, not a failure
        manifest([("reference", 0)] + full[1:])
        rc, out = run()
        assert rc == 0, out
        assert "failed" not in out, out

        manifest([])
        rc, out = run()
        assert rc == 2 and "yielded no records" in out, (rc, out)

        manifest([r for r in full if r[0] != "erv"])
        rc, out = run()
        assert rc == 2 and "erv: 3 designed event(s), no record in the manifest" in out, (rc, out)

        manifest([(src, 0 if src == "erv" else pairs) for src, pairs in full])
        rc, out = run()
        assert rc == 2 and "erv: 1 record(s) in the manifest, none produced a read" in out, (rc, out)

        # the report is still written when the verdict is no, because the numbers are the diagnosis
        assert json.load(open(os.path.join(root, "out.json")))["rna"]["erv"]["records_with_reads"] == 0
    finally:
        shutil.rmtree(root, ignore_errors=True)


def _truth_af():
    spec = importlib.util.spec_from_file_location("truth_af", os.path.join(CATALOG, "truth_af.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_truth_af_reads_the_clone_tree_and_the_copy_kind():
    """Which copies carry a variant is a property of the tree and of when the variant arose.

    Hard-coding "the tumour populations" is right for the three truncal tiers and wrong for the three
    subclonal ones, and treating every copy of a haplotype as carrying the variant is wrong for anything
    that arose after the copy-number change -- which is exactly what the catalog's multiplicity column
    records. Both rules have to come out of the data, not out of an assumption about ds-02.
    """
    ta = _truth_af()
    clones = {"T": {"parent": None}, "A": {"parent": "T"}, "A1": {"parent": "A"}, "B": {"parent": "T"}}
    desc = ta.descendants(clones)
    assert desc["T"] == {"T", "A", "A1", "B"}, desc
    assert desc["A"] == {"A", "A1"}, desc
    assert desc["B"] == {"B"}, desc
    assert desc["A1"] == {"A1"}, desc

    site = {"clone": "A", "haplotype": "1", "timing": "pre_cna"}
    def rec(source, hap, kind):
        return {"source": source, "haplotype": hap, "copy_kind": kind}
    assert ta.carries(rec("A", "1", "all"), site, desc)
    assert ta.carries(rec("A1", "1", "pre_cna"), site, desc), "A1 descends from A"
    assert not ta.carries(rec("T", "1", "all"), site, desc), "the trunk predates a subclonal variant"
    assert not ta.carries(rec("B", "1", "all"), site, desc), "B is a sibling, not a descendant"
    assert not ta.carries(rec("A", "0", "all"), site, desc), "wrong haplotype"
    assert not ta.carries(rec("NORMAL", "1", "all"), site, desc)
    post = dict(site, timing="post_cna")
    assert ta.carries(rec("A", "1", "all"), post, desc)
    assert not ta.carries(rec("A", "1", "pre_cna"), post, desc), \
        "a post-CNA variant is absent from the copy that predates the CNA"


def test_truth_af_separates_a_generator_defect_from_an_alignment_loss():
    """The check that can tell a release defect from a consumer's recovery loss.

    Two WGS checks failed on ds-02 by comparing the design against a pileup of an alignment, which is
    three things at once. Measured from the truth map the same reads were within 0.2 % of design. This
    test pins the arithmetic that made that statement: reads from the copies carrying the alt over reads
    from every copy of the interval, and a failure when THAT disagrees with the design.
    """
    ta = _truth_af()
    root = tempfile.mkdtemp(prefix="igi_taf_")
    try:
        ds = "DS"
        os.makedirs(f"{root}/output")
        with open(f"{root}/design.yaml", "w") as fh:
            fh.write("datasets:\n  DS:\n    clones:\n"
                     "      T:  {parent: null, ccf: 1.0}\n"
                     "      A:  {parent: T, ccf: 0.4}\n"
                     "      A1: {parent: A, ccf: 0.12}\n"
                     "      B:  {parent: T, ccf: 0.25}\n")
        cols = ["event_id", "chrom", "pos", "ref", "alt", "haplotype", "clone", "timing",
                "clonality_tier", "expected_vaf_tumor", "multiplicity", "tumor_cn_at_locus"]
        # purity 0.45 at cn 4: two mutant copies of four give 0.9/2.9, one gives 0.45/2.9
        with open(f"{root}/output/{ds}.snv_indel.tsv", "w") as fh:
            fh.write("\t".join(cols) + "\n")
            fh.write("\t".join(["v1", "chr1", "500", "A", "G", "0", "T", "pre_cna",
                                 "T_het", "0.3103", "2", "4.0"]) + "\n")
            fh.write("\t".join(["v2", "chr1", "600", "C", "T", "1", "T", "post_cna",
                                 "T_het", "0.1552", "1", "4.0"]) + "\n")
        recs, reads = [], []
        plan = {}   # record -> read count
        rid = 0
        for source in ("T", "A", "A1", "B"):
            for hap in ("0", "1"):
                for kind in ("all", "pre_cna"):
                    rid += 1
                    name = f"e{rid:04d}"
                    recs.append((name, "chr1", 1, 1000, source, hap, kind))
                    plan[name] = 450 if source == "T" else 0
        for hap in ("0", "1"):
            rid += 1
            name = f"e{rid:04d}"
            recs.append((name, "chr1", 1, 1000, "NORMAL", hap, "all"))
            plan[name] = 550

        def write(plan):
            with open(f"{root}/records.tsv", "w") as fh:
                fh.write("record\tchrom\tstart\tend\tsource\thaplotype\tcopy_kind\tinterval_set\n")
                for r in recs:
                    fh.write("\t".join(str(x) for x in r) + "\ton_target\n")
            with open(f"{root}/readmap.tsv", "w") as fh:
                fh.write("read_name\tsource_record\n")
                for name, n in plan.items():
                    for i in range(n):
                        fh.write(f"r{name}_{i}\t{name}\n")

        def run(extra=()):
            r = subprocess.run([sys.executable, os.path.join(CATALOG, "truth_af.py"),
                                "--dataset", ds, "--catalog-dir", f"{root}/output",
                                "--design", f"{root}/design.yaml",
                                "--records", f"{root}/records.tsv",
                                "--readmap", f"{root}/readmap.tsv",
                                "--report", f"{root}/r.json", *extra],
                               capture_output=True, text=True)
            return r.returncode, r.stdout + r.stderr

        write(plan)
        # the awk fast path and the Python fallback must agree; the WGS readmap is a billion lines and
        # only one of these two is fast enough to be worth running on it
        want = {name: [("chr1:500", True)] for name in plan}
        fast = ta.record_counts(f"{root}/readmap.tsv", want, prefer_awk=True)
        slow = ta.record_counts(f"{root}/readmap.tsv", want, prefer_awk=False)
        assert dict(fast) == {k: v for k, v in slow.items() if v}, (fast, slow)
        assert sum(fast.values()) == sum(plan.values()) > 0, fast

        rc, out = run()
        assert rc == 0, out
        rep = json.load(open(f"{root}/r.json"))
        got = {s["site"]: s for s in rep["sites"]}
        # v1 sits on both copies of haplotype 0; v2 only on the copy that postdates the CNA
        assert got["chr1:500"]["alt_reads"] == 900 and got["chr1:500"]["reads"] == 2900, rep
        assert abs(got["chr1:500"]["truth"] - 0.3103) < 0.001, rep
        assert got["chr1:600"]["alt_reads"] == 450, rep
        assert abs(got["chr1:600"]["truth"] - 0.1552) < 0.001, rep
        assert rep["tiers"]["T_het"]["sites"] == 2 and not rep["failures"], rep

        # the generator writes one copy where the design says two: a release defect, and this is the
        # check that calls it one
        halved = dict(plan)
        halved["e0002"] = 0
        write(halved)
        rc, out = run()
        assert rc == 1 and "generator and the design disagree" in out and "T_het" in out, (rc, out)

        # a readmap that names none of the records measured nothing, and must not pass
        write({k: 0 for k in plan})
        rc, out = run()
        assert rc == 2 and "nothing was measured" in out, (rc, out)

        # the counts cache must be rebuilt when it is older than the readmap it came from. A stale cache
        # is the §20.2 defect in miniature: an answer that is current-looking and wrong, and here it
        # would report the previous build's allele fractions against this build's design.
        write(plan)
        cache = f"{root}/counts.tsv"
        rc, out = run(("--counts", cache))
        assert rc == 0 and os.path.exists(cache), (rc, out)
        first = open(cache).read()
        rc, out = run(("--counts", cache))      # reused, not rebuilt
        assert rc == 0 and open(cache).read() == first, (rc, out)
        assert json.load(open(f"{root}/r.json"))["tiers"]["T_het"]["ratio"] == 1.0

        os.utime(cache, (1, 1))                 # older than the readmap
        halved = dict(plan); halved["e0002"] = 0
        write(halved)
        rc, out = run(("--counts", cache))
        assert rc == 1, "a stale cache was reused and hid a generator defect\n" + out
        assert os.path.getmtime(cache) > 1, "the stale cache was not rebuilt"

        os.utime(cache, (1, 1))                 # and a stale cache must not hide a correct build either
        write(plan)
        rc, out = run(("--counts", cache))
        assert rc == 0, out
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_normal_library_depth_ignores_tumour_ploidy():
    """A normal library is diploid, so its depth denominator is 2.0 whatever the tumour's ploidy.

    The normal plan is flat -- `source_plan(..., tumor=False)` returns two haplotypes at weight 1.0 each,
    and its own comment says "the normal library is diploid everywhere". Dividing that by the TUMOUR
    mixture denominator made every IGI-SYN-SEQ-02 normal library 69 % of its designed depth: with WGD and
    purity 0.45 the denominator is 0.45*1*4 + 0.55*2 = 2.9, so 2/2.9 = 0.690. Normal WES came out at 69x
    against a designed 100x, and normal WGS, HiFi and ONT at 20.7x against 30x.

    IGI-SYN-SEQ-01 cannot show this: no WGD and purity 0.70 give 0.70*1*2 + 0.30*2 = 2.0 exactly, the
    right answer by coincidence. A test covering only ds-01 would pass on the defect, so both datasets
    are asserted and the WGD one carries the point.
    """
    import yaml
    from igi_catalog.clones import CloneModel
    from igi_catalog.pipeline import depth_denominator, source_plan

    root = tempfile.mkdtemp(prefix="igi_normdepth_")
    try:
        # Every arm the design's CNA list names, so CloneModel can resolve "chr17:p" and friends. The
        # coordinates are stand-ins: this test is arithmetic on the clone tree and never reads a base.
        arms = os.path.join(root, "arms.bed")
        with open(arms, "w") as fh:
            for i in list(range(1, 23)) + ["X", "Y"]:
                c = f"chr{i}"
                fh.write(f"{c}\t0\t100000000\tp\n{c}\t100000000\t250000000\tq\n")
        design = yaml.safe_load(open(os.path.join(CATALOG, "design.yaml")))

        seen_wgd = False
        for ds, cfg in design["datasets"].items():
            cm = CloneModel(cfg, arms)
            seen_wgd = seen_wgd or cm.wgd
            den_n = depth_denominator(cm, cm.purity, tumor=False)
            assert abs(den_n - 2.0) < 1e-9, (
                f"{ds}: normal depth denominator is {den_n:.4f}, not 2.0; every normal library it "
                f"builds lands at {2.0 / den_n:.3f} of its designed depth")

            plan = source_plan(cm, cm.purity, "chr1", 50_000_000, tumor=False)
            assert sum(w for *_, w in plan) == 2.0, plan
            for _src, _hap, _kind, w in plan:
                cov = 100.0 * w / den_n
                assert abs(cov - 50.0) < 1e-9, (
                    f"{ds}: a normal haplotype plans {cov:.2f}x of a requested 100x, not 50x")

            # and the tumour denominator must NOT have been flattened by the fix
            den_t = depth_denominator(cm, cm.purity, tumor=True)
            want = cm.purity * sum(cm.excl.values()) * sum(cm.base) + (1.0 - cm.purity) * 2.0
            assert abs(den_t - want) < 1e-9, (
                f"{ds}: tumour denominator {den_t:.4f} != {want:.4f}; the normal fix reached the tumour")

        assert seen_wgd, "no WGD dataset in the design, so this test cannot see the defect it exists for"
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_vdj_reads_reach_the_constant_region():
    """A TCR amplicon must cover J and the proximal constant region, not just the V framework.

    The TCR builder drew reads through `five_prime_window`, GEX's 450 bp 5' constant. A V(D)J library is
    a different enrichment -- nested PCR off TRAC/TRBC, amplifying V through C -- and 450 bp stops at the
    V-J/C boundary: C begins at a median 404 nt on TRA and 434 on TRB, so the window reached 46 and 16 nt
    into constant regions of 978 and 758 nt. Cell Ranger called a constant gene on 89 % of TRA contigs
    and 17 % of TRB, and paired only 38-54 % of cells.

    Both chains are asserted, and TRB carries the point: its longer junction is what pushed C out of
    reach, so a test on TRA alone would have been nearly satisfied by the defect.
    """
    import random
    from igi_catalog.vdj import assemble, germline, recombine, vdj_amplicon
    from igi_catalog.tenx import five_prime_window

    g = germline()
    rng = random.Random(20261012)
    for chain in ("TRA", "TRB"):
        c_cov, old_cov, n = [], [], 0
        for _ in range(200):
            rec = recombine(chain, rng)
            seq = assemble(rec)
            if not seq:
                continue
            n += 1
            c_start = len(seq) - len(g["c"][rec["c"]]["framework"])
            c_cov.append(len(vdj_amplicon(seq, rec)) - c_start)
            old_cov.append(len(five_prime_window(seq)) - c_start)
        assert n > 100, f"{chain}: only {n} recombinants assembled"
        med = sorted(c_cov)[len(c_cov) // 2]
        assert med >= 150, (
            f"{chain}: the amplicon covers a median {med} nt of the constant region; Cell Ranger needs "
            f"enough C to call the constant gene and annotate the chain as productive")
        assert min(c_cov) > 0, f"{chain}: some amplicons do not reach the constant region at all"
        # the amplicon must stay a V(D)J amplicon, not become the whole transcript
        assert max(len(vdj_amplicon(assemble(recombine(chain, rng)) or "", {"c": rec["c"]}) or "x")
                   for _ in range(20)) < 1200, f"{chain}: amplicon is implausibly long"
        # and the regression it exists for must be visible: the old window did not reach C
        old_med = sorted(old_cov)[len(old_cov) // 2]
        assert old_med < med, (
            f"{chain}: the 5' GEX window covered {old_med} nt of C and the amplicon {med}; this test "
            f"cannot tell the fix from the defect")


def test_kinnex_bulk_segment_carries_its_primers():
    """A bulk Kinnex segment must start and end in primer, not in cDNA.

    The builder emitted bare cDNA at both ends. Nothing downstream could orient or trim it, and LENS's
    bulk path goes straight to pbmm2 with no lima, so pigeon called 34 % of reads antisense and 51 % of
    isoforms novel-not-in-catalog -- against 0.01 % and 89 % full-splice-match for the single-cell
    library, which carries its 10x structure and does pass through lima.

    The constants are measured off the real HG002 bulk `segmented.bam`, so this test also pins the
    recovered layout: a later edit that changes a primer silently would change what lima can demultiplex.
    """
    import random
    from igi_catalog.kinnex import (BULK_ISOSEQ_5P, BULK_POLYA, BULK_SAMPLE_BARCODE, BULK_SMART_RC,
                                   BULK_TRUSEQ_R1, kinnex_bulk_segment)

    rng = random.Random(7)
    cdna = "".join(rng.choice("ACGT") for _ in range(1200))
    polya = []
    for _ in range(300):
        seg = kinnex_bulk_segment(cdna, rng)
        assert seg.startswith(BULK_TRUSEQ_R1 + BULK_SAMPLE_BARCODE + BULK_ISOSEQ_5P), seg[:60]
        assert seg.endswith(BULK_SMART_RC), seg[-40:]
        assert cdna in seg, "the cDNA did not survive into the segment"
        body = seg[len(BULK_TRUSEQ_R1 + BULK_SAMPLE_BARCODE + BULK_ISOSEQ_5P):-len(BULK_SMART_RC)]
        assert body.startswith(cdna), "something was inserted between the 5' primer and the cDNA"
        tail = body[len(cdna):]
        assert set(tail) == {"A"}, f"the polyA run is not pure A: {tail!r}"
        polya.append(len(tail))
    lo, mode, hi = BULK_POLYA
    assert lo <= min(polya) and max(polya) <= hi, (min(polya), max(polya), BULK_POLYA)
    med = sorted(polya)[len(polya) // 2]
    assert abs(med - mode) <= 2, f"polyA median {med} against a measured {mode}"

    # the outer primers must be the pair the consumer preflight demultiplexes with
    assert BULK_TRUSEQ_R1 == "CTACACGACGCTCTTCCGATCT"
    assert BULK_SMART_RC == "GTACTCTGCGTTGATACCACTGCTT"


def test_interleave_mixes_chromosomes_and_keeps_mates():
    """A merged library's head must look like the whole library, and R1 record n must stay mated to R2 n.

    The merge stages used `cat`, so every delivered library was chromosome-block ordered: ds-01 tumour
    WES began with chr10 and its first 3.4 M pairs were chr10 and chr11 only. Reads within a chromosome
    are shuffled, so no aligner and no acceptance result depended on it -- but the consumer preflight
    pushes "the first 5 %" through each tool, and on the exome that slice held no chr6 reads at all. The
    razers3 and OptiType preflight, whose entire purpose is the MHC region on chr6, tested a sample with
    none of it in.

    The inputs here are deliberately lopsided (200 against 50 against 10), because equal inputs would
    make a naive round robin look correct.
    """
    root = tempfile.mkdtemp(prefix="igi_interleave_")
    try:
        sizes = {"chrA": 50, "chrB": 200, "chrC": 10}
        for name, n in sizes.items():
            for m in (1, 2):
                with gzip.open(os.path.join(root, f"{name}_R{m}.fastq.gz"), "wt") as fh:
                    for i in range(n):
                        fh.write(f"@{name}_{i}/{m}\nACGT\n+\nIIII\n")
        names = sorted(sizes)
        rc = subprocess.run(
            [sys.executable, os.path.join(CATALOG, "interleave_fastq.py"),
             "--in-r1", ",".join(os.path.join(root, f"{n}_R1.fastq.gz") for n in names),
             "--in-r2", ",".join(os.path.join(root, f"{n}_R2.fastq.gz") for n in names),
             "--out-r1", os.path.join(root, "m_R1.fastq.gz"),
             "--out-r2", os.path.join(root, "m_R2.fastq.gz")],
            capture_output=True, text=True)
        assert rc.returncode == 0, rc.stdout + rc.stderr

        def heads(path):
            with gzip.open(path, "rt") as fh:
                return [l.strip() for i, l in enumerate(fh) if i % 4 == 0]
        h1, h2 = heads(os.path.join(root, "m_R1.fastq.gz")), heads(os.path.join(root, "m_R2.fastq.gz"))
        assert len(h1) == sum(sizes.values()) == len(h2), (len(h1), len(h2))
        assert [x.rsplit("/", 1)[0] for x in h1] == [x.rsplit("/", 1)[0] for x in h2], \
            "the mates desynced, which breaks every paired consumer"

        src = [x.split("_")[0][1:] for x in h1]
        head = src[:max(1, len(src) // 20)]
        assert set(head) == set(sizes), (
            f"the first 5 % holds only {sorted(set(head))}; this is the defect, not the fix")
        # every decile must carry the library's composition, not just the first slice
        for d in range(10):
            seg = src[d * len(src) // 10:(d + 1) * len(src) // 10]
            for name, n in sizes.items():
                want = n / sum(sizes.values()) * len(seg)
                got = seg.count(name)
                assert abs(got - want) <= 1.5, (d, name, got, want)
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_simulator_streams_differ_per_library():
    """The four single-cell assays share a cell roster and must NOT share a read-error stream.

    `--seed 9100` is mandatory and documented in three job headers: the roster and molecule pool are
    shared, and that is what makes a cell the same cell in 10x GEX, 10x TCR, Kinnex sc and ONT sc. The
    same integer was then passed to the read simulator, so the i-th GEX R2 quality string came out
    byte-identical to the i-th TCR R2 string and to IGI-SYN-SEQ-02's -- across all four libraries and
    both datasets. Nothing downstream breaks; it is a tell, and it defeats treating the assays as
    independent measurements.

    So the roster seed stays shared and only the simulator stream is scoped.
    """
    from igi_catalog.readnames import stream_seed

    assays = ("tenx_gex", "tenx_tcr", "kinnex_sc", "ont")
    datasets = ("IGI-SYN-SEQ-01", "IGI-SYN-SEQ-02")
    seen = {}
    for ds in datasets:
        for assay in assays:
            v = stream_seed(9100, ds, assay)
            assert v not in seen, (
                f"{ds}/{assay} and {seen[v]} draw the same simulator stream ({v})")
            seen[v] = f"{ds}/{assay}"
            assert 0 <= v < 2_000_000_000, v
    # stable across processes: hash() is salted per process and would make a build unreproducible
    assert stream_seed(9100, "IGI-SYN-SEQ-01", "tenx_gex") == 754194991, \
        "stream_seed changed value; every single-cell library's error stream moves with it"


def _main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"  ok    {t.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"  FAIL  {t.__name__}: {e}")
    print(f"\n{len(tests) - failed}/{len(tests)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(_main())
