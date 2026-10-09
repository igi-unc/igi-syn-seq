#!/usr/bin/env python3
"""Push the first ~5 % of a delivered library through a tool LENS will actually use, and record whether
the tool choked.

    consumer_preflight.py --tool razers3 --reads a_R1.fastq.gz --out out.json --label x \
                          --cmd '<container prefix> {args}' [--fraction 0.05]

WHY A FRACTION AND NOT THE WHOLE FILE. The question this answers is "does the consumer accept this
input", which is a property of the FORMAT and of the first record that breaks it -- not something that
needs every read. razers3 aborted on the ds-02 exome while loading; it never reached read two million.
Running a whole library through every tool would cost more than generating it.

WHAT THIS DOES NOT ESTABLISH, stated plainly because the distinction has already cost this release once.
A 5 % sample cannot prove the ABSENCE of a rare offender. The character that aborted razers3 was a single
`Y` in 4,858,789,950 bases: a 5 % sample would have missed it 95 % of the time. Absence is established by
the complete alphabet scan (section 17.6.2), which reads every base of every file. The two are
complementary and neither replaces the other --

    complete scan      every base, one cheap property (is it ACGTN)
    this preflight      5 % of the bases, the real tool's entire opinion of them

A tool that chokes on the 5 % would have choked on the whole; a tool that accepts the 5 % has told us the
format is right, not that the library is clean.
"""
import argparse
import json
import os
import subprocess
import sys


def run(cmd, timeout=None):
    return subprocess.run(["bash", "-c", cmd], capture_output=True, text=True, timeout=timeout)


def uncompressed_size(path, samtools=None):
    """Best estimate of the decompressed size, for working out what 5 % means."""
    if path.endswith(".bam"):
        return os.path.getsize(path) * 3        # BAM is already compressed; SAM text is ~3x
    r = run(f"pigz -l {path} 2>/dev/null | awk 'NR==2{{print $2}}'")
    try:
        n = int(r.stdout.strip())
        if n > 0:
            return n                            # pigz -l wraps above 4 GiB, so sanity-check it
    except ValueError:
        pass
    return os.path.getsize(path) * 4            # FASTQ gzips to about a quarter


def sample(path, out, fraction, samtools):
    """First `fraction` of the file, cut at a record boundary."""
    want = max(400_000, int(uncompressed_size(path) * fraction))
    if path.endswith(".bam"):
        # SEQ/QUAL only; enough for any tool that reads unaligned records
        cmd = (f"{samtools.format(args=f'view -h {path}')} | "
               f"awk -v m={want} '/^@/{{print; next}} {{b+=length($0)+1; print; if (b>m) exit}}' "
               f"> {out}")
    else:
        cmd = (f"pigz -dc {path} | "
               f"awk -v m={want} '{{b+=length($0)+1; print; if (b>m && NR%4==0) exit}}' > {out}")
    r = run(cmd)
    return r.returncode == 0 and os.path.getsize(out) > 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tool", required=True)
    ap.add_argument("--reads", required=True, help="one path, or R1,R2")
    ap.add_argument("--out", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--cmd", required=True, help="container prefix with {args}")
    ap.add_argument("--samtools-cmd", default="samtools {args}")
    ap.add_argument("--ref", default="")
    ap.add_argument("--fraction", type=float, default=0.05)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--timeout", type=int, default=7200)
    a = ap.parse_args()

    work = os.path.dirname(os.path.abspath(a.out))
    os.makedirs(work, exist_ok=True)
    paths = a.reads.split(",")
    cut = []
    for i, p in enumerate(paths, 1):
        ext = "sam" if p.endswith(".bam") else "fq"
        q = os.path.join(work, f"pf_{a.label}_{i}.{ext}")
        if not sample(p, q, a.fraction, a.samtools_cmd):
            json.dump({"label": a.label, "tool": a.tool, "ok": False,
                       "detail": f"could not sample {os.path.basename(p)}"}, open(a.out, "w"))
            print(f"could not sample {p}", file=sys.stderr)
            return 1
        cut.append(q)

    # Each tool is invoked the way LENS invokes it; the point is the tool's own parser, not the result.
    if a.tool == "razers3":
        args = (f"-i 95 -m 1 -dr 0 -tc {a.threads} -o {work}/pf_{a.label}.bam "
                f"/usr/local/bin/data/hla_reference_dna.fasta {cut[0]}")
    elif a.tool == "salmon":
        args = (f"quant -i {a.ref} -l A -1 {cut[0]} -2 {cut[1]} -p {a.threads} "
                f"--validateMappings -o {work}/pf_{a.label}_salmon")
    elif a.tool == "pbmm2":
        args = f"align --preset HIFI -j {a.threads} {a.ref} {cut[0]} {work}/pf_{a.label}.bam"
    elif a.tool == "lima":
        args = (f"{cut[0]} {a.ref} {work}/pf_{a.label}_lima.bam --per-read --num-threads {a.threads}")
    elif a.tool == "cellranger_vdj":
        # cellranger insists on its own directory layout and sample naming
        fqdir = os.path.join(work, f"pf_{a.label}_fq")
        os.makedirs(fqdir, exist_ok=True)
        for i, q in enumerate(cut, 1):
            os.replace(q, os.path.join(fqdir, f"{a.label}_S1_L001_R{i}_001.fastq"))
        args = (f"vdj --id=pf_{a.label} --reference={a.ref} --fastqs={fqdir} --sample={a.label} "
                f"--localcores={a.threads} --localmem=16")
    else:
        print(f"unknown tool {a.tool}", file=sys.stderr)
        return 2

    try:
        r = run(a.cmd.format(args=args), timeout=a.timeout)
        rc, err = r.returncode, (r.stderr or "")[-600:]
    except subprocess.TimeoutExpired:
        rc, err = -1, f"timed out after {a.timeout}s"

    # A parse error is the signature that matters: it means the tool rejected the DATA, not the setup.
    choke = [k for k in ("ParseError", "Unexpected character", "not a valid", "malformed",
                         "invalid character", "Invalid FASTQ", "unexpected end")
             if k.lower() in err.lower()]
    out = {"label": a.label, "tool": a.tool, "reads": paths, "fraction": a.fraction,
           "exit_code": rc, "ok": rc == 0, "parse_error": bool(choke),
           "signature": choke[0] if choke else "", "stderr_tail": err.strip()[-300:]}
    json.dump(out, open(a.out, "w"), indent=2)
    print(f"  {a.label:<22} {a.tool:<15} exit={rc} ok={out['ok']} parse_error={out['parse_error']}")
    if choke:
        print(f"    CHOKED: {choke[0]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
