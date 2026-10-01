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
import os
import subprocess
import time


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
        ok = rc == 0 and os.path.exists(out) and os.path.getsize(out) > 100
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


def combine_fastq(parts, out_path, log=print):
    """Concatenate FASTQ chunks with record framing enforced, then assert it. Returns the read count."""
    n_fixed = 0
    with open(out_path, "wb") as out:
        for part in parts:
            with open(part, "rb") as fh:
                data = fh.read()
            if not data:
                continue
            if not data.endswith(b"\n"):
                data += b"\n"
                n_fixed += 1
            out.write(data)
    with open(out_path, "rb") as fh:
        n_lines = sum(1 for _ in fh)
    if n_lines % 4:
        raise RuntimeError(f"combined FASTQ has {n_lines} lines, not a multiple of 4; record framing is "
                           f"broken and any BAM built from it would be corrupt")
    log(f"  combined {len(parts)} chunk(s), {n_lines // 4:,} reads"
        + (f", {n_fixed} needed a trailing newline" if n_fixed else ""))
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
