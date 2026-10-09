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
import sys
import tempfile
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
    """A symbolic allele always fails; an ambiguity code fails only where the consumer is strict."""
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
        qa.check_alphabet(root, results, limit=400000)
        got = {r["file"]: r["verdict"].split(" ")[0] for r in results}

        assert got["DS_full_rna_R1.fastq.gz"] == "PASS", got
        assert got["DS_bad_rna_R1.fastq.gz"] == "FAIL", "a `*` in the sequence must fail"
        assert got["DS_chr1_tumor_R1.fastq.gz"] == "WARN", "IUPAC in a bwa-consumed arm warns"
        assert got["DS_chr3_tumor_ont_wgs.fastq.gz"] == "WARN", "IUPAC in a minimap2-consumed arm warns"
        assert got["DS-GEX_S1_L001_R1_001.fastq.gz"] == "FAIL", \
            "IUPAC in a 10x arm must FAIL: Cell Ranger refuses the run on one character"
        # A WARN must not inflate the release's failure count, or the distinction is pointless.
        assert sum(1 for r in results if r["verdict"].startswith("FAIL")) == 2, results
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_alphabet_verdict_cause_naming():
    """The verdict names the cause, because the two causes need different fixes."""
    qa = _qa()
    assert qa._alphabet_verdict([]) == "PASS"
    assert "symbolic allele" in qa._alphabet_verdict(["*"], "/r/rna/DS/x_R1.fastq.gz")
    assert "IUPAC" in qa._alphabet_verdict(["Y"], "/r/ont_wgs/DS/x.fastq.gz")
    # Mixed: the symbolic allele dominates and the arm no longer matters.
    assert qa._alphabet_verdict(["*", "Y"], "/r/ont_wgs/DS/x.fastq.gz").startswith("FAIL")
    assert qa._alphabet_verdict(["Y"], "/r/tenx_tcr/DS/x_R1_001.fastq.gz").startswith("FAIL")


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
