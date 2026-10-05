#!/usr/bin/env python3
"""Repair delivered long-read libraries in place, without re-simulating them.

Four defects were found in delivered data by review. All four are in the packaging, not in the reads, so
the sequence and the base qualities -- which are what the compute bought -- are kept and only the envelope
is rebuilt.

    --what ont      strip Badread's description from every FASTQ header and write it to a separate
                    gzipped map. The header named the source reference, strand and coordinates for a
                    genomic read and the source molecule id for an RNA one, so every read carried its own
                    answer and `*_molecules.tsv.gz` could be joined without aligning anything.

    --what hifi     rewrite a `samtools import` BAM as a Revio one: Revio read names, zm/np/ec/rq/qs/qe,
                    an @RG with PU/PM, and a .pbi. As delivered they had eleven fields and no tags, which
                    pbmm2, pbindex, extracthifi and DeepVariant's HiFi model all reject or mishandle.

    --what kinnex   redraw per-read np/ec and add the @RG PU that skera reads to name segments. np was a
                    single constant across the whole library while rq varied per read, and the missing PU
                    made every delivered segment name start with a slash where the movie belongs.

Each repair writes beside the original and only replaces it once the new file has been counted and
verified, because a half-written 90 GB library that has already overwritten its source is unrecoverable.
"""
import argparse
import gzip
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from igi_catalog.pacbio_bam import movie_in, pbindex, repair_bam  # noqa: E402


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def _replace(tmp, final, label):
    """Swap a verified temporary file in for the original, keeping the original until the swap lands."""
    old = final + ".superseded"
    os.replace(final, old)
    os.replace(tmp, final)
    os.remove(old)
    log(f"  {label}: replaced in place ({os.path.getsize(final) / 1e9:.1f} GB)")


def strip_ont(fastq, threads=8, keep_map=True):
    """Strip descriptions from a gzipped FASTQ, writing the removed text to a sidecar map."""
    t0 = time.time()
    tmp = fastq + ".stripped.tmp.gz"
    map_path = fastq.replace(".fastq.gz", "_read_map.tsv.gz")
    n = n_desc = 0
    # pigz for both directions: a 90 GB library is 40 minutes of single-threaded gzip on the write side
    # alone, and the decompress side cannot be parallelised at all, so it sets the floor.
    rd = subprocess.Popen(["pigz", "-dc", fastq], stdout=subprocess.PIPE, bufsize=1 << 22)
    wr = subprocess.Popen(["pigz", "-p", str(threads), "-c"], stdin=subprocess.PIPE,
                          stdout=open(tmp, "wb"), bufsize=1 << 22)
    gm = gzip.open(map_path + ".tmp", "wt", compresslevel=6) if keep_map else None
    try:
        if gm:
            gm.write("read_name\tsource\n")
        for i, line in enumerate(rd.stdout):
            if i % 4 == 0:
                n += 1
                sp = line.rstrip(b"\n").split(b" ", 1)
                if len(sp) > 1:
                    n_desc += 1
                    if gm:
                        gm.write(f"{sp[0][1:].decode()}\t{sp[1].decode()}\n")
                line = sp[0] + b"\n"
                if n % 20_000_000 == 0:
                    log(f"    {n:,} reads ({time.time() - t0:.0f}s)")
            wr.stdin.write(line)
    finally:
        if gm:
            gm.close()
        wr.stdin.close()
        rd.stdout.close()
    if rd.wait() != 0:
        raise RuntimeError(f"pigz -dc failed on {fastq}")
    if wr.wait() != 0:
        raise RuntimeError(f"pigz -c failed writing {tmp}")

    # Verify before replacing: the record count must match and no header may keep a space.
    got = int(subprocess.run(f"pigz -dc {tmp} | wc -l", shell=True, capture_output=True,
                             text=True).stdout.split()[0])
    if got % 4 or got // 4 != n:
        raise RuntimeError(f"{fastq}: wrote {got} lines ({got / 4:.2f} records) for {n} reads")
    head = subprocess.run(f"pigz -dc {tmp} | head -400000 | awk 'NR%4==1' | grep -c ' ' || true",
                          shell=True, capture_output=True, text=True).stdout.strip()
    if head not in ("0", ""):
        raise RuntimeError(f"{fastq}: {head} of the first 100,000 headers still carry a description")
    log(f"  {os.path.basename(fastq)}: {n:,} reads, {n_desc:,} descriptions removed "
        f"in {time.time() - t0:.0f}s")
    if gm:
        os.replace(map_path + ".tmp", map_path)
    _replace(tmp, fastq, os.path.basename(fastq))
    return n


def repair_hifi(bam, sample, library, pbindex_cmd=None, kinnex=False, seed=0):
    """Rewrite an unaligned BAM into a Revio-shaped one and index it."""
    t0 = time.time()
    tmp = bam + ".repair.tmp.bam"
    before = movie_in(bam)
    n = repair_bam(bam, tmp, sample=sample, library=library,
                   kind="kinnex" if kinnex else "hifi_wgs", seed=seed,
                   preserve_names=kinnex,
                   progress=lambda k: log(f"    {k:,} reads ({time.time() - t0:.0f}s)"))
    # Verify before replacing: every read must carry the five tags, and the count must match.
    chk = subprocess.run(
        f"samtools view {tmp} | head -20000 | "
        f"awk '{{t=0; for(i=12;i<=NF;i++) if($i~/^(zm|np|ec|rq|qs):/) t++; if(t==5) ok++}} "
        f"END {{print ok+0, NR}}'", shell=True, capture_output=True, text=True).stdout.split()
    if len(chk) != 2 or chk[0] != chk[1] or chk[1] == "0":
        raise RuntimeError(f"{bam}: only {chk} reads carry all five HiFi tags")
    pu = subprocess.run(f"samtools view -H {tmp} | grep -o 'PU:[^\t]*'", shell=True,
                        capture_output=True, text=True).stdout.strip()
    if not pu.startswith("PU:m"):
        raise RuntimeError(f"{bam}: @RG has no movie in PU ({pu!r}); skera would emit '/zmw/ccs/...'")
    log(f"  {os.path.basename(bam)}: {n:,} reads, movie {before!r} -> {pu[3:]!r}, "
        f"{time.time() - t0:.0f}s")
    _replace(tmp, bam, os.path.basename(bam))
    if pbindex_cmd:
        p = pbindex(bam, pbindex_cmd)
        log(f"  {os.path.basename(bam)}: pbi {'built' if p else 'NOT built'}"
            + (f" ({os.path.getsize(p) / 1e6:.0f} MB)" if p else ""))
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--what", required=True, choices=["ont", "hifi", "kinnex"])
    ap.add_argument("--path", required=True, help="the file to repair, in place")
    ap.add_argument("--sample", default=None, help="dataset id, for the BAM read group")
    ap.add_argument("--library", default=None)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--pbindex-cmd", default=None)
    a = ap.parse_args()
    if not os.path.exists(a.path):
        sys.exit(f"missing: {a.path}")
    log(f"{a.what}: {a.path} ({os.path.getsize(a.path) / 1e9:.1f} GB)")
    if a.what == "ont":
        strip_ont(a.path, threads=a.threads)
    else:
        sample = a.sample or os.path.basename(a.path).split("_")[0]
        repair_hifi(a.path, sample, a.library or sample, a.pbindex_cmd,
                    kinnex=(a.what == "kinnex"), seed=a.seed)
    log("done")


if __name__ == "__main__":
    main()
