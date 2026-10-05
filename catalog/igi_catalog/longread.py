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
import subprocess
import time


def _framed(path, tail_bytes=1 << 16):
    """Whether a FASTQ ends on a record boundary. Reads only the tail, so it costs nothing on a 170 MB file.

    A complete FASTQ has a line count divisible by four and ends with a newline. Checking the tail catches
    the truncation a killed simulator leaves behind; a full line count would be exact but would mean reading
    every chunk twice.
    """
    try:
        size = os.path.getsize(path)
        with open(path, "rb") as fh:
            fh.seek(max(0, size - tail_bytes))
            tail = fh.read()
        if not tail.endswith(b"\n"):
            return False
        # the last four lines should be header, sequence, plus, quality, with matching lengths
        lines = tail.split(b"\n")[:-1]
        if len(lines) < 4:
            return True
        h, seq, plus, qual = lines[-4:]
        return h.startswith(b"@") and plus.startswith(b"+") and len(seq) == len(qual)
    except OSError:
        return False


def run_chunks(jobs, n_parallel, timeout_s, attempts=3, log=print):
    """Run (label, out_path, command) jobs concurrently. Returns the list of outputs that produced data.

    Raises if any job produced nothing after `attempts` tries: a silently missing chunk is a silently
    truncated library.
    """
    t0 = time.time()
    log(f"  {len(jobs)} chunk(s) over {n_parallel} concurrent slot(s)")

    def run_one(spec, attempt=1):
        label, out, cmd = spec
        try:
            with open(out, "wb") as fh:
                rc = subprocess.run(cmd, shell=True, stdout=fh, stderr=subprocess.DEVNULL,
                                    timeout=timeout_s).returncode
        except subprocess.TimeoutExpired:
            rc = -1
        # Judge on the OUTPUT, not only the return code. A timeout kills the simulator mid-record, so a
        # truncated chunk can be large and still useless: the ONT bulk ds-01 run left 500 chunks totalling
        # 82 GB of which 19 in 20 had a line count that was not a multiple of four. Size alone would have
        # accepted them, and concatenating them would have produced a corrupt library that looked right.
        ok = rc == 0 and os.path.exists(out) and os.path.getsize(out) > 100 and _framed(out)
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
                failed.append(label)
    if failed:
        raise RuntimeError(f"the simulator produced nothing for {len(failed)} chunk(s) after {attempts} "
                           f"attempts: {failed[:5]}")
    log(f"  simulation done in {time.time() - t0:.0f}s, {len(done)} chunk(s)")
    return done


def combine_fastq(parts, out_path, log=print, map_path=None, keep_description=False):
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
    """
    n_fixed = 0
    n_reads = 0
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
                            line = sp[0] + b"\n"
                            if gm and len(sp) > 1:
                                gm.write(f"{sp[0][1:].decode()}\t{sp[1].decode()}\n")
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
