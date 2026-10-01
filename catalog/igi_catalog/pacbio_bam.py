"""Write a PacBio-style unaligned HiFi BAM that the real PacBio tools will accept.

`samtools import` produces a valid unaligned BAM, but not one skera will read: it refuses with
"Bam record missing the read quality tag". PacBio tools expect the per-read tags a Revio run carries, so
they are written here.

    zm  the ZMW: one per insert, and for Kinnex one per array
    np  number of passes behind the consensus
    rq  predicted read accuracy

`rq` is computed from the simulated base qualities rather than asserted from a constant, so it describes the
read that was actually produced. It is a mean error PROBABILITY converted once at the end, not a mean of
Phred scores: Phred is a logarithm, and averaging it then converting reports 98.92 % where the true answer is
92.15 %.
"""
import math
import os
import re

import pysam

SRC = re.compile(r"([A-Za-z0-9_.|-]+),([+-])strand,(\d+)-(\d+)")


def mean_accuracy(qual):
    """Mean accuracy over a quality string, averaging error probabilities."""
    if not qual:
        return 0.0
    e = sum(10 ** (-q / 10.0) for q in qual) / len(qual)
    return 1.0 - e


def source_of(description):
    """The reference sequence a Badread read came from, parsed out of its description."""
    m = SRC.search(description or "")
    return m.group(1) if m else None


def write_hifi_bam(fastq_paths, out_bam, movie="m84000_260101_000000_s0", one_per_source=False,
                   np_passes=8, sample="SAMPLE", library=None):
    """Convert Badread FASTQ(s) into a HiFi-style unaligned BAM.

    Returns (n_written, {source: read_name}). With `one_per_source` only the first read from each reference
    sequence is kept, which is what a Kinnex array needs: one ZMW yields one HiFi read, and Badread's
    sampling can hand the same array two reads and another none.
    """
    header = {
        "HD": {"VN": "1.6", "SO": "unknown"},
        "RG": [{"ID": "synthetic", "PL": "PACBIO", "SM": sample,
                "LB": library or sample,
                "DS": "READTYPE=CCS;BINDINGKIT=102-739-100;SEQUENCINGKIT=102-118-800;"
                      "BASECALLERVERSION=5.0;FRAMERATEHZ=100.000000"}],
        "PG": [{"ID": "igi-syn-seq", "PN": "igi-syn-seq", "DS": "synthetic HiFi reads"}],
    }
    seen = {}
    n = 0
    zmw = 0
    with pysam.AlignmentFile(out_bam, "wb", header=header) as out:
        for fq in fastq_paths:
            if not fq or not os.path.exists(fq):
                continue
            with pysam.FastxFile(fq) as fh:
                for rec in fh:
                    src = source_of(rec.comment)
                    if one_per_source:
                        if src is None or src in seen:
                            continue
                    zmw += 1
                    name = f"{movie}/{zmw}/ccs"
                    if src is not None:
                        seen.setdefault(src, name)
                    a = pysam.AlignedSegment()
                    a.query_name = name
                    a.query_sequence = rec.sequence
                    a.flag = 4
                    a.query_qualities = pysam.qualitystring_to_array(rec.quality)
                    rq = mean_accuracy(a.query_qualities)
                    a.set_tags([("zm", zmw, "i"), ("np", int(np_passes), "i"),
                                ("rq", round(rq, 6), "f"), ("RG", "synthetic", "Z")])
                    out.write(a)
                    n += 1
    return n, seen
