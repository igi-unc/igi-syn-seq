#!/usr/bin/env python3
"""Merge per-chromosome FASTQ pieces into one library in mixed order, instead of concatenating them.

    interleave_fastq.py --out-r1 merged_R1.fastq.gz --out-r2 merged_R2.fastq.gz \
                        --in-r1 chr1_R1.fastq.gz,chr2_R1.fastq.gz,... \
                        --in-r2 chr1_R2.fastq.gz,chr2_R2.fastq.gz,...

WHY. The merge stages concatenated the per-chromosome arms with `cat`, so every delivered library is
chromosome-block ordered: IGI-SYN-SEQ-01's tumour WES begins with chr10 and its first 3.4 M pairs are
chr10 and chr11 only. Reads WITHIN a chromosome are shuffled, so aligners are unaffected and no
acceptance result depends on this -- but any head-based slice of the file is unrepresentative, and two
things in this repo took such a slice:

  * `consumer_preflight.py` pushes "the first 5 %" of a library through each consumer. On the exome that
    sample contained no chr6 reads at all, so the razers3 and OptiType preflight -- whose entire purpose
    is the MHC region on chr6 -- tested a sample with none of it in.
  * a user doing `head` or `zcat | head` to look at the data, which is the first thing anyone does.

Block order is also a fingerprint: no real sequencer emits reads sorted by reference position, so the
ordering alone says the library is synthetic. That is cosmetic next to the preflight problem.

HOW. A weighted round robin, so the chromosomes exhaust together rather than one at a time. Records are
counted once per input (a decompression pass), then drawn with a largest-remaining-fraction scheduler:
at each step take from whichever input has the greatest share of itself still unread. That needs no sort
and no temp space, holds one record per input in memory, and is deterministic.

COST, MEASURED rather than guessed. 33,191,885 pairs over three chromosomes on 8 cores: 772 s, against
4 s for `cat` of the same inputs. So ~43,000 pairs/s, roughly 190x the cost of concatenation, because
`cat` of gzip members needs no decompression and this needs two passes of it plus a Python loop per
record.

What that means per arm, which matters more than the ratio:

  Illumina WGS tumour   943,646,057 pairs   ~6.1 h, and this is the expensive one
  Illumina WGS normal   ~280 M pairs        ~1.8 h
  exome                 68,114,175 pairs    ~26 min
  ONT WGS               ~5 M single reads   a couple of minutes; the decompression dominates

An earlier note in 50_merge_longread.sbatch put the cost on "the ~180 GB ONT WGS arm", which is the
wrong arm: ONT libraries hold millions of reads, not hundreds of millions. The 24 h limits on the merge
stages hold with room to spare, and the figure above is what to check that against when a depth changes.

MATES. R1 and R2 are advanced in lockstep from the same chosen input, so mate n of the output R1 is
still mate n of the output R2. Cutting them independently is what §20.5's preflight failure was, and the
same mistake here would break every paired consumer rather than one check.
"""
import argparse
import os
import subprocess
import sys


def pigz_reader(path):
    """A line iterator over a gzip file, via pigz when available."""
    exe = "pigz" if _have("pigz") else "gzip"
    p = subprocess.Popen([exe, "-dc", path], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    return p


def _have(prog):
    return subprocess.run(["bash", "-c", f"command -v {prog}"],
                          capture_output=True).returncode == 0


def count_records(path):
    """Records in a FASTQ, counted by lines/4.

    Counted rather than estimated because the scheduler's exhaustion alignment depends on it, and a
    wrong count would pile one chromosome up at the end of the file -- the defect this script exists to
    remove, in a less obvious form.
    """
    exe = "pigz" if _have("pigz") else "gzip"
    r = subprocess.run(f"{exe} -dc {path!r} | wc -l", shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"could not read {path}: {r.stderr.strip()[-200:]}")
    n = int(r.stdout.split()[0])
    if n % 4:
        raise SystemExit(f"{path} holds {n} lines, not a multiple of 4; it is truncated")
    return n // 4


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-r1", required=True, help="comma-separated per-chromosome R1 files")
    ap.add_argument("--in-r2", default="", help="comma-separated R2 files, in the same order")
    ap.add_argument("--out-r1", required=True)
    ap.add_argument("--out-r2", default="")
    ap.add_argument("--level", type=int, default=6, help="output gzip level")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()

    r1 = [p for p in a.in_r1.split(",") if p]
    r2 = [p for p in a.in_r2.split(",") if p] if a.in_r2 else []
    if r2 and len(r2) != len(r1):
        raise SystemExit(f"{len(r1)} R1 files against {len(r2)} R2 files")
    if r2 and not a.out_r2:
        raise SystemExit("--in-r2 given without --out-r2")
    for p in r1 + r2:
        if not os.path.exists(p):
            raise SystemExit(f"missing input: {p}")

    counts = [count_records(p) for p in r1]
    if r2:
        # The mates must agree file by file. A mismatch here means the builder wrote an unequal pair and
        # the interleave would silently pair read n of one chromosome with read n of another.
        for p1, p2, n1 in zip(r1, r2, counts):
            n2 = count_records(p2)
            if n1 != n2:
                raise SystemExit(f"{p1} has {n1} records and {p2} has {n2}; the pair is not mated")
    total = sum(counts)
    if total == 0:
        raise SystemExit("every input is empty, so there is nothing to merge")
    print(f"  {len(r1)} input(s), {total:,} pairs", flush=True)

    procs1 = [pigz_reader(p) for p in r1]
    procs2 = [pigz_reader(p) for p in r2] if r2 else []
    exe = "pigz" if _have("pigz") else "gzip"
    comp = f"{exe} -{a.level}" + (f" -p {a.threads}" if exe == "pigz" else "")
    o1 = subprocess.Popen(f"{comp} > {a.out_r1!r}", shell=True, stdin=subprocess.PIPE, text=True,
                          bufsize=1 << 20)
    o2 = (subprocess.Popen(f"{comp} > {a.out_r2!r}", shell=True, stdin=subprocess.PIPE, text=True,
                           bufsize=1 << 20) if r2 else None)

    remaining = list(counts)
    written = 0
    try:
        while True:
            # largest remaining fraction of its own total; ties break to the lower index, so the result
            # is deterministic and reproducible from the same inputs
            best, best_frac = -1, -1.0
            for i, rem in enumerate(remaining):
                if rem <= 0:
                    continue
                frac = rem / counts[i]
                if frac > best_frac:
                    best, best_frac = i, frac
            if best < 0:
                break
            rec1 = [procs1[best].stdout.readline() for _ in range(4)]
            if not rec1[0]:
                raise SystemExit(f"{r1[best]} ended after {counts[best] - remaining[best]} of "
                                 f"{counts[best]} records; it changed under us")
            o1.stdin.write("".join(rec1))
            if o2:
                rec2 = [procs2[best].stdout.readline() for _ in range(4)]
                if not rec2[0]:
                    raise SystemExit(f"{r2[best]} ended early; it changed under us")
                o2.stdin.write("".join(rec2))
            remaining[best] -= 1
            written += 1
            if written % 20_000_000 == 0:
                print(f"    {written:,} / {total:,}", flush=True)
    finally:
        for p in procs1 + procs2:
            if p.stdout:
                p.stdout.close()
            p.wait()
        for o in (o1, o2):
            if o:
                o.stdin.close()
                o.wait()

    if written != total:
        raise SystemExit(f"wrote {written} pairs of an expected {total}")
    print(f"  interleaved {written:,} pairs", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
