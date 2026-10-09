"""Write a PacBio-style unaligned HiFi BAM that the real PacBio tools will accept.

`samtools import` produces a valid unaligned BAM, but not one any PacBio tool will read: it carries
eleven fields and no tags at all. PacBio tools expect what a Revio run carries, so it is written here.

    zm  the ZMW: one per insert, and for Kinnex one per array
    np  number of passes behind the consensus -- FABRICATED. Badread has no pass structure, so
        there are no passes behind these reads; np is drawn from the real HG002 distribution so a
        tool filtering on pass count sees a realistic spread rather than one constant, and it is
        not evidence about consensus depth. ec inherits the same caveat.
    ec  effective coverage, a float that tracks np
    rq  predicted read accuracy
    qs  query start and qe query end, which skera uses to name segments

Three details were each paid for by a defect in delivered data:

- **`rq` is computed from the simulated base qualities** rather than asserted from a constant, so it
  describes the read that was actually produced. It is a mean error PROBABILITY converted once at the end,
  not a mean of Phred scores: Phred is a logarithm, and averaging it then converting reports 98.92 % where
  the true answer is 92.15 %.

- **`np` is drawn per read, not fixed.** Every synthetic read carried `np:i:8` while `rq` varied per read,
  so a tool filtering on pass count saw no variation at all across a library. np now comes from the
  empirical CDF measured off the real HG002 BAMs (`resources/pacbio_np_model.json`) and ec from np times
  the measured ec/np ratio, which is how the two relate in real data.

- **`@RG` must carry `PU`.** skera names each segment `<movie>/<zmw>/ccs/<qs>_<qe>` and takes the movie
  from the read group's PU field, not from the read name it was handed. With no PU it emitted
  `/1/ccs/17_2371` -- a leading slash where the movie belongs -- across every delivered Kinnex segment,
  which is not a name any PacBio tool will parse.
"""
import bisect
import hashlib
import json
import math
import os
import re
import subprocess
from array import array

import pysam

SRC = re.compile(r"([A-Za-z0-9_.|-]+),([+-])strand,(\d+)-(\d+)")

# Revio ZMW hole numbers are large and sparsely spread rather than counted from one; the real movie
# measured for the np model spans 201,329,666 - 266,145,206. Synthetic ZMWs are drawn from that range so
# that a hole number looks like one, while staying unique and reproducible (see zmw_for).
ZMW_LO, ZMW_HI = 201_329_666, 266_145_206

_MODEL_PATH = os.path.join(os.path.dirname(__file__), "..", "resources", "pacbio_np_model.json")
_MODELS = {}


def np_model(kind="hifi_wgs"):
    """The empirical np CDF and ec/np ratio measured from real HG002 data."""
    if not _MODELS:
        with open(os.path.abspath(_MODEL_PATH)) as fh:
            _MODELS.update(json.load(fh))
    m = _MODELS[kind]
    return [p[0] for p in m["np_cdf"]], [p[1] for p in m["np_cdf"]], m["ec_over_np"]


def _u(seed, *parts):
    """A uniform in [0,1) from a name. Deterministic and independent of read order, so a rerun or a
    reordered merge gives a read the same np it had before."""
    h = hashlib.blake2b(":".join([str(seed)] + [str(p) for p in parts]).encode(), digest_size=8)
    return int.from_bytes(h.digest(), "big") / 2 ** 64


def draw_np(seed, key, kind="hifi_wgs"):
    """(np, ec) for one read, from the empirical distribution."""
    vals, cdf, ratio = np_model(kind)
    n = vals[min(bisect.bisect_left(cdf, _u(seed, "np", key)), len(vals) - 1)]
    # ec scatters about np*ratio in real data; a tenth of a pass of jitter reproduces that without
    # pretending to a precision the measurement does not support
    ec = n * ratio * (1.0 + 0.04 * (_u(seed, "ec", key) - 0.5))
    return int(n), round(ec, 5)


def zmw_for(index, seed=0):
    """A plausible, unique, reproducible ZMW hole number for the index-th read of a library.

    Uniqueness matters more than realism here: two reads sharing a hole number would make the library
    self-inconsistent, and pbindex would index one of them away. The index is therefore mapped into the
    real hole-number range by a stride that is coprime with the range, which is a permutation.
    """
    span = ZMW_HI - ZMW_LO
    return ZMW_LO + (index * 7_919 + int(seed) * 104_729) % span


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


def rg_id(movie, read_type="CCS"):
    """The read-group ID, which must be EIGHT HEXADECIMAL DIGITS.

    pbbam parses an @RG ID as a hex number -- real ones look like `b0776b05` or `b0776b05/0--0` -- and
    `ID:synthetic` made pbindex abort with "ERROR: stoul" and leave behind a 65-byte index reporting zero
    reads. The same ID was in every delivered Kinnex BAM. Eight hex digits off the movie name reproduce
    pbbam's own scheme closely enough that the whole toolchain accepts it; verified against pbindex,
    pbindexdump, extracthifi and zmwfilter.
    """
    return hashlib.blake2b(f"{movie}//{read_type}".encode(), digest_size=4).hexdigest()


def header(movie, sample, library=None, kind="hifi_wgs"):
    """A Revio-shaped BAM header. PU is the movie name and is what skera reads to name segments."""
    return {
        "HD": {"VN": "1.6", "SO": "unknown", "pb": "5.0.0"},
        "RG": [{"ID": rg_id(movie), "PL": "PACBIO", "PM": "REVIO", "PU": movie,
                "SM": sample, "LB": library or sample,
                "DS": "READTYPE=CCS;BINDINGKIT=102-739-100;SEQUENCINGKIT=102-118-800;"
                      "BASECALLERVERSION=5.0;FRAMERATEHZ=100.000000"}],
        "PG": [{"ID": "igi-syn-seq", "PN": "igi-syn-seq", "DS": f"synthetic {kind} reads"}],
    }


def movie_for(sample, library=None):
    """A stable Revio-style movie name per library, so merging two libraries cannot collide ZMWs.

    Every field after the instrument must be DECIMAL. A real name is m<instrument>_<YYMMDD>_<HHMMSS>_s<N>
    and pbindex parses those fields as numbers: a hex digest here put a letter in the date field and
    pbindex died with "ERROR: stoul", leaving a 65-byte index reporting zero reads.
    """
    h = int(hashlib.blake2b(f"{sample}/{library or ''}".encode(), digest_size=8).hexdigest(), 16)
    mon, day = h % 12 + 1, (h // 12) % 28 + 1
    hh, mm, ss = (h // 400) % 24, (h // 10_000) % 60, (h // 1_000_000) % 60
    return f"m84000_26{mon:02d}{day:02d}_{hh:02d}{mm:02d}{ss:02d}_s{(h // 7) % 4 + 1}"


def _emit(out, name, seq, qual, zmw, np_passes, ec, rq, qs=8, rg="synthetic"):
    a = pysam.AlignedSegment()
    a.query_name = name
    a.query_sequence = seq
    a.flag = 4
    a.query_qualities = pysam.qualitystring_to_array(qual) if isinstance(qual, str) else qual
    # qs/qe bracket the insert inside the ZMW read; real Revio HiFi has qs=8 and qe=qs+len for every
    # read of the 300,000 measured, and skera uses them to name segments
    a.set_tags([("zm", int(zmw), "i"), ("np", int(np_passes), "i"), ("ec", float(ec), "f"),
                ("rq", round(float(rq), 6), "f"), ("qs", int(qs), "i"),
                ("qe", int(qs) + len(seq), "i"), ("RG", rg, "Z")])
    out.write(a)


def rescale_quality(qual, scale):
    """Make the quality string describe the error the reads actually carry.

    Badread builds sequence from an error model and qualities from a separate qscore model, and for the
    PacBio model they disagree: measured against the alignment, the reads carry about 5.15x the error their
    quality strings claim, at EVERY identity request (jobs/81_calibrate_identity.sbatch). The ratio is
    constant because the qscore model is keyed to the REQUESTED identity, not the realised one, so raising
    --identity improves the sequence and leaves the claim exactly as optimistic as before. That is the F24
    defect itself -- reads asserting an accuracy they do not have -- and unlike the error floor it is ours
    to fix, because we write the BAM.

    `scale` multiplies each base's error PROBABILITY, which preserves the qscore model's per-position
    structure. Two things stop that from being a plain multiply, and both had to be measured rather than
    assumed (see docs/build-reference.md §14.3):

    1. It SATURATES. 45.8 % of the claimed error mass sits on the 0.08 % of bases at Q3-Q7, whose
       probability cannot be multiplied by five without exceeding 1. Setting `scale` to the measured 5.15
       therefore raised the aggregate by only 3.70x, leaving reads still 1.39x optimistic. The scale that
       lands on the measured true error is 10.0, solved numerically against both the genomic and the cDNA
       quality distributions, which agree.
    2. Its effect is QUANTIZED. BAM Phred is integer and every input here is an integer, so a constant
       scale shifts every base by the same whole number of decibels: the achievable aggregate moves in ~1 dB
       steps, and no scale lands exactly on the target. 10.0 is chosen as the smallest step that is
       PESSIMISTIC (claimed error 1.05x true, rather than 0.93x one step below) because a read understating
       its own accuracy is harmless where overstating it is the defect being fixed -- the same rule applied
       to the ONT arms. It is also exactly 10 dB, so it needs no rounding at all; the alternative plateau
       boundary sits at scale 8.9125, where a float a few thousandths either way flips every base by a
       whole Phred.

    So the correction is applied as the integer decibel shift it actually is, computed once, rather than
    per base where it would sit on that knife edge. `rq` then follows from the rescaled array.

    Phred is clamped to [1, 93]: 0 means "no quality available" in SAM and 93 is the encodable maximum. The
    floor also caps implied error at Q1 (p = 0.79), which is the right ceiling for a base that is wrong.
    """
    if not scale or scale == 1.0:
        return qual
    shift = int(round(10.0 * math.log10(scale)))
    if shift == 0:
        return qual
    out = array("B", bytes(len(qual)))
    for i, q in enumerate(qual):
        out[i] = max(1, min(93, q - shift))
    return out


def write_hifi_bam(fastq_paths, out_bam, movie=None, one_per_source=False,
                   np_passes=None, sample="SAMPLE", library=None, kind="hifi_wgs", seed=0,
                   index_offset=0, qual_error_scale=1.0):
    """Convert Badread FASTQ(s) into a HiFi-style unaligned BAM.

    Returns (n_written, {source: read_name}). With `one_per_source` only the first read from each
    reference sequence is kept, which is what a Kinnex array needs: one ZMW yields one HiFi read, and
    Badread's sampling can hand the same array two reads and another none.

    `np_passes` fixes the pass count for every read; leaving it None draws per read from the measured
    distribution, which is what real data looks like and what callers filtering on np need.

    `qual_error_scale` calibrates the quality strings against the alignment -- see `rescale_quality`.
    """
    movie = movie or movie_for(sample, library)
    rg = rg_id(movie)
    seen = {}
    n = 0
    with pysam.AlignmentFile(out_bam, "wb", header=header(movie, sample, library, kind)) as out:
        for fq in fastq_paths:
            if not fq or not os.path.exists(fq):
                continue
            with pysam.FastxFile(fq) as fh:
                for rec in fh:
                    src = source_of(rec.comment)
                    if one_per_source:
                        if src is None or src in seen:
                            continue
                    idx = index_offset + n
                    # For Kinnex the ZMW is the array index, counted from 1, because skera's segment
                    # names and the di/dl tags are read per ZMW; for WGS it is a real-looking hole number.
                    zmw = idx + 1 if kind == "kinnex" else zmw_for(idx, seed)
                    name = f"{movie}/{zmw}/ccs"
                    if src is not None:
                        seen.setdefault(src, name)
                    qual = rescale_quality(pysam.qualitystring_to_array(rec.quality),
                                           qual_error_scale)
                    rq = mean_accuracy(qual)
                    if np_passes is None:
                        npass, ec = draw_np(seed, f"{library or sample}:{idx}", kind)
                    else:
                        npass, ec = int(np_passes), round(int(np_passes) * np_model(kind)[2], 5)
                    _emit(out, name, rec.sequence, qual, zmw, npass, ec, rq, rg=rg)
                    n += 1
    return n, seen


def movie_in(bam):
    """The movie name already used by a BAM's read names, or None. Preserving it matters: the Kinnex truth
    tables key arrays by ZMW, so a repair must not renumber or rename anything a truth file points at."""
    with pysam.AlignmentFile(bam, "rb", check_sq=False) as fh:
        for rec in fh.fetch(until_eof=True):
            parts = (rec.query_name or "").split("/")
            return parts[0] or None
    return None


def repair_bam(in_bam, out_bam, sample, library=None, kind="hifi_wgs", seed=0, movie=None,
               preserve_names=False, progress=None):
    """Rewrite an existing unaligned BAM into a Revio-shaped one.

    Two different repairs go through here.

    `preserve_names=False` is the HiFi WGS case. Those BAMs came from `samtools import`, which keeps
    sequence and qualities and discards everything else: eleven fields, no tags, no read group, and
    Badread's UUIDs for read names. Sequence and qualities are the expensive part and they are intact, so
    the libraries are repaired from the BAM rather than re-simulated -- one streaming pass instead of 96
    badread jobs -- and the reads are renamed and tagged on the way through.

    `preserve_names=True` is the Kinnex case. There the read names and ZMWs are already right and the
    truth tables point at them, so only the per-read np/ec and the header's PU are rebuilt. Renaming would
    silently detach every array from its truth row.
    """
    movie = movie or (preserve_names and movie_in(in_bam)) or movie_for(sample, library)
    rg = rg_id(movie)
    n = 0
    with pysam.AlignmentFile(in_bam, "rb", check_sq=False, threads=4) as src:
        with pysam.AlignmentFile(out_bam, "wb", header=header(movie, sample, library, kind),
                                 threads=6) as out:
            for rec in src.fetch(until_eof=True):
                if preserve_names:
                    name = rec.query_name
                    zmw = rec.get_tag("zm") if rec.has_tag("zm") else n + 1
                    qs = rec.get_tag("qs") if rec.has_tag("qs") else 8
                else:
                    zmw = n + 1 if kind == "kinnex" else zmw_for(n, seed)
                    name, qs = f"{movie}/{zmw}/ccs", 8
                npass, ec = draw_np(seed, f"{library or sample}:{n}", kind)
                _emit(out, name, rec.query_sequence, rec.query_qualities, zmw, npass, ec,
                      mean_accuracy(rec.query_qualities), qs=qs, rg=rg)
                n += 1
                if progress and n % 1_000_000 == 0:
                    progress(n)
    return n


def pbindex(bam, cmd=None):
    """Build the .pbi index the design calls for. Returns the index path, or None if pbindex is absent.

    A failed run still leaves a stub index behind -- 65 bytes reporting zero reads -- which looks like a
    .pbi to anything that only checks for the file. It is removed rather than shipped.
    """
    if not cmd:
        return None
    pbi = bam + ".pbi"
    try:
        subprocess.run(cmd.format(args=bam), shell=True, check=True)
    except subprocess.CalledProcessError:
        if os.path.exists(pbi):
            os.remove(pbi)
        raise
    if os.path.exists(pbi) and os.path.getsize(pbi) > 128:
        return pbi
    if os.path.exists(pbi):
        os.remove(pbi)
    raise RuntimeError(f"pbindex produced an empty index for {bam}")

def drop_empty_records(bam):
    """Remove records with no sequence from a skera-split BAM, in place. Returns the number dropped.

    `skera split` can emit a ZERO-LENGTH segment, and it is not our arrays that cause it. The release scan
    found 11 such records in ds-01's chr1to6 Kinnex library and 7 in ds-02's, written out as SAM's `SEQ=*`
    placeholder -- which is not a base, and which an earlier version of check_alphabet mis-reported as a
    symbolic allele. Mechanism, from the segment names skera assigns (`<movie>/<zmw>/ccs/<qs>_<qe>`):

        161745/ccs/5867_12020    len 6153
        161745/ccs/12033_12033   len 0
        161745/ccs/12050_13683   len 1633

    The adapters are 17 bp. Between the two real segments lie 30 bp carrying TWO adapter alignments -- one
    error-shortened to 13 bp (12020-12033) and one exact (12033-12050) -- which overlap in the underlying
    sequence, so the gap between their inner ends is empty and skera emits a segment of length zero.
    KinnexArray.build lays out `adapter + molecule + adapter + molecule ... + adapter` and molecule
    sampling drops any record without a sequence, so two adapters are never adjacent in what we hand
    skera; the double match is skera's aligner on a noisy adapter copy. It is plausibly inflated by our
    reads sitting at 1.72x real HiFi error (see section 14), which gives an adapter more room to match
    twice than a real one has.

    Dropped rather than kept, for two reasons:

      - A delivered read with no bases is invalid whoever produced it. SeqAn-based readers abort on
        malformed input rather than skipping it, which is the razers3 lesson (section 17.6).
      - It corrupts the truth-map join. The arrays map records `segment_index` over OUR molecule order,
        so a consumer pairing skera's emitted segments to it positionally is shifted by every phantom
        segment after the first. Dropping them makes that join correct, so this is not cosmetic.

    At 7 in 5,460,256 the volume is immaterial; the malformedness is not.
    """
    import pysam
    tmp = bam + ".tmp"
    dropped = 0
    with pysam.AlignmentFile(bam, "rb", check_sq=False) as src:
        with pysam.AlignmentFile(tmp, "wb", header=src.header) as out:
            for rec in src.fetch(until_eof=True):
                if rec.query_sequence is None or len(rec.query_sequence) == 0:
                    dropped += 1
                    continue
                out.write(rec)
    os.replace(tmp, bam)
    return dropped
