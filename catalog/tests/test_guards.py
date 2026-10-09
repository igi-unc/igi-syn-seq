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
        old = time.time() - 9 * 86400
        os.utime(f"{root}/merged/DS/DS_normal_R1.fastq.gz", (old, old))
        open(f"{root}/wes/DS/DS_chr1_normal_R1.fastq.gz", "w").close()

        results = []
        qa.check_currency(root, results)
        assert len(results) == 1, results
        assert results[0]["verdict"].startswith("FAIL"), results
        assert "stale by 21" in results[0]["verdict"] or "stale by 2" in results[0]["verdict"], \
            f"the lag should be about 216 h: {results[0]['verdict']}"
        assert "wes/" in results[0]["verdict"], "it must name the input directory to re-run from"

        # Re-running the merge clears it.
        now = time.time()
        os.utime(f"{root}/merged/DS/DS_normal_R1.fastq.gz", (now, now))
        results = []
        qa.check_currency(root, results)
        assert results[0]["verdict"] == "PASS", results
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
