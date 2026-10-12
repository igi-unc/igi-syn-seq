"""Running a long-read simulator over many chunks, and combining the result safely.

Extracted from run_pacbio_wgs.py so the Kinnex builders inherit these details instead of rediscovering
them. Each line of it was paid for:

- **Chunks must run concurrently, not merely exist.** `subprocess.run` blocks, so splitting the work and
  then looping over it is serialisation with extra steps: 5 chunks an hour. A bounded pool does ~98.
- **No pipe to gzip.** Piping badread into gzip left six shells hung at zero CPU: gzip held the pipe open,
  sh waited on gzip, and each stuck shell held a pool slot until the map deadlocked. The compression was
  wasted anyway, since the deliverable is a BAM.
- **A timeout per chunk, and retries.** Container starts contend, and sixteen at once on one node was
  enough that some failed to start at all. One bad start should cost a retry, not the run.
- **Combining is not `cat`.** Six of 54 badread chunks ended without a trailing newline, so concatenation
  glued each one's last line onto the next chunk's first. Every file was individually well formed and the
  combination was not. Framing is enforced and then asserted.
"""
import concurrent.futures as cf
import gzip
import os
import random
import subprocess
import time
import zlib


def _framed(path, tail_bytes=1 << 16, max_tail=1 << 24):
    """Whether a FASTQ ends on a record boundary. Reads only the tail, so it costs nothing on a 170 MB file.

    A complete FASTQ has a line count divisible by four and ends with a newline. Checking the tail catches
    the truncation a killed simulator leaves behind; a full line count would be exact but would mean
    reading every chunk twice.

    THE WINDOW MUST BE SIZED FROM THE DATA. A fixed 64 KiB tail rejected seven of 96 ONT WGS tasks, and
    the files were perfectly good: ONT reads run to tens of kilobases -- the genomic length model has a
    standard deviation of 15.5 kb and a measured maximum of 551 kb -- so one record is a header plus two
    lines of equal length plus a '+', and at a 32,725 bp read that is 65,605 bytes against a 65,536 byte
    window. Sixty-nine bytes over, and the first line in the window was a fragment of a sequence rather
    than a header, so the check declared a complete file truncated. badread had reached 100.0 %, written
    117 MB and exited cleanly; the retry then re-ran it four more times to the same byte and the task
    failed. The window therefore grows until it holds at least five lines or the whole file.
    """
    try:
        size = os.path.getsize(path)
        want = tail_bytes
        while True:
            with open(path, "rb") as fh:
                fh.seek(max(0, size - want))
                tail = fh.read()
            if not tail.endswith(b"\n"):
                return False
            lines = tail.split(b"\n")[:-1]
            # Enough lines to see a whole record AND know the first one is not a fragment, or the window
            # already covers the file, in which case the first line really is the first line.
            if len(lines) >= 5 or want >= size:
                break
            if want >= max_tail:
                # a single record larger than 16 MiB is not a FASTQ this pipeline produces
                return False
            want = min(want * 8, max_tail)
        if len(lines) < 4:
            return True
        h, seq, plus, qual = lines[-4:]
        return h.startswith(b"@") and plus.startswith(b"+") and len(seq) == len(qual)
    except OSError:
        return False


def run_chunks(jobs, n_parallel, timeout_s, attempts=5, log=print):
    """Run (label, out_path, command) jobs concurrently. Returns the list of outputs that produced data.

    Raises if any job produced nothing after `attempts` tries: a silently missing chunk is a silently
    truncated library.
    """
    t0 = time.time()
    log(f"  {len(jobs)} chunk(s) over {n_parallel} concurrent slot(s)")

    def run_one(spec, attempt=1):
        label, out, cmd = spec
        # Stagger the launches. The failure this retries is container start contention -- 48 simultaneous
        # singularity invocations on one node -- so retrying immediately retries into the same storm. The
        # jitter is derived from the output path so it is reproducible, and it grows with the attempt.
        time.sleep((zlib.crc32(out.encode()) % 1000) / 1000.0 * 2.0 * attempt)
        err_path = out + f".err{attempt}"
        try:
            with open(out, "wb") as fh, open(err_path, "wb") as eh:
                rc = subprocess.run(cmd, shell=True, stdout=fh, stderr=eh,
                                    timeout=timeout_s).returncode
        except subprocess.TimeoutExpired:
            rc = -1
        # Judge on the OUTPUT, not only the return code. A timeout kills the simulator mid-record, so a
        # truncated chunk can be large and still useless: the ONT bulk ds-01 run left 500 chunks totalling
        # 82 GB of which 19 in 20 had a line count that was not a multiple of four. Size alone would have
        # accepted them, and concatenating them would have produced a corrupt library that looked right.
        ok = rc == 0 and os.path.exists(out) and os.path.getsize(out) > 100 and _framed(out)
        if ok:
            for a in range(1, attempts + 1):
                p = out + f".err{a}"
                if os.path.exists(p):
                    os.remove(p)
        if not ok and attempt < attempts:
            return run_one(spec, attempt + 1)
        return label, out, ok, attempt

    done, failed = [], []
    with cf.ThreadPoolExecutor(max_workers=n_parallel) as ex:
        for label, out, ok, tries in ex.map(run_one, jobs):
            if ok:
                done.append(out)
                if tries > 1:
                    log(f"    {label}: succeeded on attempt {tries}")
            else:
                failed.append((label, out))
    if failed:
        # Report WHY, not just that. stderr used to go to DEVNULL, so a reproducible failure left nothing
        # to diagnose: seven ONT tasks each lost one chunk of eighty and the only evidence was the chunk's
        # label, which several chunks share. The tail of the last attempt's stderr is kept instead.
        detail = []
        for label, out in failed[:4]:
            tail = ""
            for a in range(attempts, 0, -1):
                p = out + f".err{a}"
                if os.path.exists(p) and os.path.getsize(p):
                    with open(p, "rb") as fh:
                        fh.seek(max(0, os.path.getsize(p) - 4096))
                        lines = [x for x in fh.read().decode(errors="replace").replace("\r", "\n").split("\n")
                                 if x.strip() and "Simulating" not in x]
                    tail = " | ".join(lines[-3:])[:300]
                    break
            detail.append(f"{label} [{os.path.basename(out)}]: {tail or 'no stderr captured'}")
        raise RuntimeError(f"the simulator produced nothing for {len(failed)} chunk(s) after {attempts} "
                           f"attempts:\n    " + "\n    ".join(detail))
    log(f"  simulation done in {time.time() - t0:.0f}s, {len(done)} chunk(s)")
    return done


# ONT header structure, MEASURED from 511,401 real PromethION R10.4.1 reads -- the IPISRC044 T1 library
# in the site's own data-sources, basecalled with dorado sup@v4.3.0. Written out because the first
# version of this invented a plausible header instead of reading one that was already on disk, which is
# the mistake §26 exists to record. A real header:
#
#   @f870f5e6-... runid=aa773787-6b3e-4ae6-b0e4-99917e0808ff ch=1462
#   start_time=2025-07-17T16:35:15.369691-07:00 flow_cell_id=PBE81376
#   basecall_gpu=NVIDIA_A100_80GB_PCIe protocol_group_id=PM0349_20250717
#   sample_id=PROM0236_... parent_read_id=f870f5e6-... basecall_model_version_id=dna_r10.4.1_...
#
# What the measurement settled, and what the invented version had wrong:
#
#   field order   fixed, and this is it: 511,400 of 511,401 reads. ONT_HEADER_FIELDS records it.
#   read=         DOES NOT EXIST. Older guppy emitted it; this basecaller does not. It was invented.
#   runid         ONE UUID for the whole run, not 16 hex characters.
#   ch            15 to 2997, 2,664 distinct over 511k reads -- a ~3,000-position flow cell with some
#                 channels inactive, not a counter and not 2,675. ONT_CHANNELS/ONT_ACTIVE_CHANNELS.
#   start_time    1-SECOND resolution with a CONSTANT sub-second part (.369691 on every one of the
#                 511,401 reads: the run's start phase), and 835.6 reads/s over 612 distinct seconds.
#   parent_read_id  equals the read id for 493,479 of 511,401 reads (96.5 %); the other 3.5 % are reads
#                 split from a longer one and keep the parent's id. Not modelled: always equal here, and
#                 this comment is the record of what that leaves out.
ONT_HEADER_FIELDS = ("runid", "ch", "start_time", "flow_cell_id", "basecall_gpu",
                     "protocol_group_id", "sample_id", "parent_read_id",
                     "basecall_model_version_id")
ONT_CHANNELS = 3000          # flow cell positions; measured max channel 2997
ONT_ACTIVE_CHANNELS = 2664   # distinct channels seen in the real run; the rest were inactive pores
ONT_READS_PER_S = 836        # measured 835.6 reads/s
ONT_SUBSECOND = 369691       # the constant microsecond part of every real start_time
ONT_BASECALL_GPU = "NVIDIA_A100_80GB_PCIe"
# The offset is arbitrary and says nothing about anyone: the real run's -07:00 is its sequencing site's,
# and a synthetic dataset should not carry it.
ONT_UTC_OFFSET = "-04:00"

# Simulator truth that must never appear in a delivered header. Badread writes all of it; the strings are
# listed rather than tested as "the header contains a space", which is what the first check did -- and
# every real ONT header contains spaces, so that test would have failed this release on all 20 files
# while passing a header that leaked coordinates and no run metadata at all.
ONT_SIMULATOR_FIELDS = ("length=", "error-free_length=", "read_identity=", "strand,")


def _ont_start_time(t0_epoch, second_index):
    """One real ONT `start_time`: 1-second resolution, constant sub-second part, fixed offset.

    Both properties are measured, not stylistic: every one of 511,401 real reads carried `.369691`, so
    the fractional part is the run's start phase rather than a per-read time, and emitting a varying one
    would be a detail no basecaller produces.
    """
    t = time.gmtime(t0_epoch + second_index)
    return (f"{time.strftime('%Y-%m-%dT%H:%M:%S', t)}.{ONT_SUBSECOND:06d}{ONT_UTC_OFFSET}").encode()


def combine_fastq(parts, out_path, log=print, map_path=None, keep_description=False,
                  run_meta=None):
    """Concatenate FASTQ chunks with record framing enforced, then assert it. Returns the read count.

    By default the simulator's description is STRIPPED from every header and, if `map_path` is given,
    written to a separate gzipped map instead.

    Badread puts the answer in the header. A genomic read carries the reference it came from, the strand
    and the coordinates:

        @a1887f88-... chr10_T_hap0_all,-strand,89396886-89421712 length=24764 error-free_length=24792 read_identity=99.593%

    which names the chromosome, the clone, the haplotype and the copy, and gives the true position; an RNA
    read carries the source molecule id, which joins straight to `*_molecules.tsv.gz` and so hands over the
    cell, the clone and the transcript. `error-free_length` and `read_identity` do not exist in real data
    at all. A benchmark whose reads carry their own truth is not a benchmark, and the Illumina path was
    fixed for exactly this (see `readnames.shuffle_and_rename`) while the long-read path was not.

    The map keeps the link for whoever needs it, in the one place a tool under test will not read.

    `run_meta` puts REAL ONT run metadata back in the stripped header. Stripping the description left a
    bare UUID, which no basecaller produces: a real header carries runid, read, ch, start_time,
    flow_cell_id, protocol_group_id and sample_id, and tools and users alike expect them. Supply a dict
    of fixed fields and `read=` and `ch=` are filled in per record.

    `sample_id` must be the SYNTHETIC dataset name. The ds-02 germline baseline is called from a real
    donor and the baseline VCF still carries that donor's cohort accession as its sample name; putting
    the same string in 390 million read headers would propagate an upstream identifier into every
    delivered file, which is the opposite of what a synthetic release is for.
    """
    n_fixed = 0
    n_reads = 0
    if run_meta:
        missing = [k for k in ONT_HEADER_FIELDS
                   if k not in ("ch", "start_time", "parent_read_id") and k not in run_meta]
        if missing:
            raise ValueError(f"run_meta is missing {missing}; the delivered header would differ from "
                             f"the measured ONT layout in ONT_HEADER_FIELDS")
        runid = str(run_meta["runid"]).encode()
        # The active channel subset, so the library shows ~2,664 distinct channels over a 3,000-position
        # flow cell rather than all 3,000 -- the real run's inactive pores. Derived from the runid so it
        # is stable across reruns and differs between libraries.
        r = random.Random(f"ch:{run_meta['runid']}")
        channels = r.sample(range(1, ONT_CHANNELS + 1), ONT_ACTIVE_CHANNELS)
        t0 = run_meta.get("start_epoch") or 0
        ts_bucket = -1
        ts = b""
        # A cheap LCG rather than channels[n % len(channels)]. Indexing by the read number cycles, so
        # channel order repeats exactly every 2,664 reads and the library would carry a periodicity no
        # real run has. One multiply-add per read is free at this scale, and seeding from the runid keeps
        # a rebuild reproducible.
        lcg = (int(run_meta["runid"].replace("-", "")[:8], 16) or 1) & 0x7FFFFFFF
        mid = (f" flow_cell_id={run_meta['flow_cell_id']}"
               f" basecall_gpu={run_meta['basecall_gpu']}"
               f" protocol_group_id={run_meta['protocol_group_id']}"
               f" sample_id={run_meta['sample_id']}"
               f" parent_read_id=").encode()
        tail = f" basecall_model_version_id={run_meta['basecall_model_version_id']}".encode()
    gm = gzip.open(map_path, "wt") if map_path else None
    try:
        if gm:
            gm.write("read_name\tsource\n")
        with open(out_path, "wb") as out:
            for part in parts:
                tail_nl = True
                with open(part, "rb") as fh:
                    for i, line in enumerate(fh):
                        if i % 4 == 0 and not keep_description:
                            raw = line.rstrip(b"\n")
                            sp = raw.split(b" ", 1)
                            if gm and len(sp) > 1:
                                gm.write(f"{sp[0][1:].decode()}\t{sp[1].decode()}\n")
                            if run_meta:
                                rid = sp[0][1:]
                                lcg = (lcg * 1103515245 + 12345) & 0x7FFFFFFF
                                ch = channels[lcg % len(channels)]
                                if n_reads // ONT_READS_PER_S != ts_bucket:
                                    ts_bucket = n_reads // ONT_READS_PER_S
                                    ts = _ont_start_time(t0, ts_bucket)
                                # Assembled from precomputed halves in ONT_HEADER_FIELDS order rather
                                # than by joining a dict: the order is a measured property of the real
                                # header, and this runs once per read over hundreds of millions.
                                line = (sp[0] + b" runid=" + runid + b" ch=" + str(ch).encode()
                                        + b" start_time=" + ts + mid + rid + tail + b"\n")
                            else:
                                line = sp[0] + b"\n"
                        if i % 4 == 0:
                            n_reads += 1
                        tail_nl = line.endswith(b"\n")
                        if not tail_nl:
                            line += b"\n"
                            n_fixed += 1
                        out.write(line)
    finally:
        if gm:
            gm.close()
    with open(out_path, "rb") as fh:
        n_lines = sum(1 for _ in fh)
    if n_lines % 4:
        raise RuntimeError(f"combined FASTQ has {n_lines} lines, not a multiple of 4; record framing is "
                           f"broken and any BAM built from it would be corrupt")
    if n_lines // 4 != n_reads:
        raise RuntimeError(f"counted {n_reads} headers but wrote {n_lines // 4} records")
    log(f"  combined {len(parts)} chunk(s), {n_lines // 4:,} reads"
        + (f", {n_fixed} needed a trailing newline" if n_fixed else "")
        + ("" if keep_description else ", headers stripped")
        + (f", map -> {os.path.basename(map_path)}" if map_path else ""))
    return n_lines // 4


def chunk_fasta(path, chunk_records, work, tag):
    """Split a FASTA into chunks of at most `chunk_records` sequences. Returns the chunk paths.

    Chunking is by record count rather than by base count because a Kinnex array is an indivisible unit:
    a chunk boundary inside an array would hand the simulator half an array as if it were a whole one.
    """
    out, fh, n, idx = [], None, 0, 0
    for line in open(path):
        if line.startswith(">"):
            if fh is None or n >= chunk_records:
                if fh:
                    fh.close()
                p = os.path.join(work, f"{tag}_c{idx:04d}.fa")
                fh = open(p, "w")
                out.append(p)
                idx += 1
                n = 0
            n += 1
        if fh:
            fh.write(line)
    if fh:
        fh.close()
    return out


def split_fasta_bp(path, chunk_bp, work, tag, overlap=60_000):
    """Split one long sequence into overlapping chunks of `chunk_bp` bases.

    Distinct from `chunk_fasta`, which splits by RECORD count. Use this when the records are whole
    chromosomes: a 249 Mb chromosome is one record, so record-count chunking leaves it whole, which both
    confines it to a single core and makes the chunk take hours. That is exactly how the ONT WGS arm lost
    61 of 96 tasks -- every chunk was a whole derived chromosome, 14 chunks over 48 slots, and badread
    needed about 4 hours per chunk against a 1 hour timeout.

    The overlap is one maximum read length so molecules spanning a boundary are not lost; only the chunk
    that owns a region emits reads for it, so none is duplicated.
    """
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].strip()
            else:
                seq.append(line.strip())
    s = "".join(seq)
    if len(s) <= chunk_bp:
        return [path]
    out, pos = [], 0
    while pos < len(s):
        end = min(len(s), pos + chunk_bp)
        sub = os.path.join(work, f"{tag}_{pos // chunk_bp:03d}.fa")
        with open(sub, "w") as fh:
            fh.write(f">{name}_chunk{pos // chunk_bp}\n")
            piece = s[pos:min(len(s), end + overlap)]
            for i in range(0, len(piece), 60):
                fh.write(piece[i:i + 60] + "\n")
        out.append(sub)
        pos = end
    return out
