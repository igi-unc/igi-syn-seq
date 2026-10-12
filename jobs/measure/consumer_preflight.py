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


def run(cmd, timeout=None, cwd=None):
    return subprocess.run(["bash", "-c", cmd], capture_output=True, text=True,
                          timeout=timeout, cwd=cwd)


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


def sample(path, out, fraction, samtools, records=None):
    """First `fraction` of the file, cut at a record boundary.

    `records`, when given, takes exactly that many reads instead of a byte budget. MATES MUST BE CUT TO
    THE SAME RECORD COUNT, and cutting them independently is what made this whole stage test itself
    rather than the release: R1 and R2 headers differ slightly in length, so equal byte budgets gave
    2,996,325 and 2,991,981 records, and every paired consumer then ran off the end of the shorter
    file. STAR reported `quality string length is not equal to sequence length` (it had started reading
    a header as a sequence) and Cell Ranger reported `Input FASTQ file ended prematurely`, both of which
    read exactly like malformed delivered data and were neither. The reads themselves had zero
    length mismatches, here and across all 612 files of the complete scan.
    """
    if path.endswith(".bam"):
        want = max(400_000, int(uncompressed_size(path) * fraction))
        # A BAM slice must stay a BAM: pbmm2 refuses `.sam` with "Unknown file suffix" and lima with
        # "[pbbam] dataset I/O ERROR: unsupported extension", so both failed on the container rather
        # than on the reads.
        n = max(1000, int(want / 2000))
        cmd = (f"{samtools.format(args=f'view -h {path}')} | "
               f"awk -v n={n} '/^@/{{print; next}} {{print; if (++c >= n) exit}}' | "
               f"{samtools.format(args=f'view -b -o {out} -')}")
    else:
        if records is None:
            want = max(400_000, int(uncompressed_size(path) * fraction))
            cmd = (f"pigz -dc {path} | "
                   f"awk -v m={want} '{{b+=length($0)+1; print; if (b>m && NR%4==0) exit}}' > {out}")
        else:
            cmd = (f"pigz -dc {path} | head -n {records * 4} > {out}")
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
    ap.add_argument("--ref2", default="",
                    help="second-stage reference (star_salmon: the transcript FASTA)")
    ap.add_argument("--cmd2", default="",
                    help="second-stage container prefix with {args}; defaults to --cmd")
    ap.add_argument("--fraction", type=float, default=0.05)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--timeout", type=int, default=7200)
    a = ap.parse_args()

    work = os.path.dirname(os.path.abspath(a.out))
    os.makedirs(work, exist_ok=True)
    paths = a.reads.split(",")
    cut = []
    # Mates get the SAME record count, taken from the first one, so the streams stay in step.
    n_rec = None
    for i, p in enumerate(paths, 1):
        ext = "bam" if p.endswith(".bam") else "fq"
        q = os.path.join(work, f"pf_{a.label}_{i}.{ext}")
        if not sample(p, q, a.fraction, a.samtools_cmd, records=n_rec if i > 1 else None):
            json.dump({"label": a.label, "tool": a.tool, "ok": False,
                       "detail": f"could not sample {os.path.basename(p)}"}, open(a.out, "w"))
            print(f"could not sample {p}", file=sys.stderr)
            return 1
        if i == 1 and ext == "fq":
            n_rec = sum(1 for _ in open(q)) // 4
        cut.append(q)
    if len(cut) > 1 and n_rec is not None:
        counts = [sum(1 for _ in open(c)) // 4 for c in cut]
        if len(set(counts)) != 1:
            json.dump({"label": a.label, "tool": a.tool, "ok": False,
                       "detail": f"mate slices disagree: {counts}"}, open(a.out, "w"))
            print(f"mate slices disagree: {counts}", file=sys.stderr)
            return 1

    # Each tool is invoked the way LENS invokes it; the point is the tool's own parser, not the result.
    #
    # `stages` is a LIST because some of those invocations are chains. The stages run in order and stop
    # at the first failure, since a later stage's input is an earlier stage's output.
    pfx = f"{work}/pf_{a.label}_"
    if a.tool == "razers3":
        stages = [("razers3", a.cmd,
                   f"-i 95 -m 1 -dr 0 -tc {a.threads} -o {work}/pf_{a.label}.bam "
                   f"/usr/local/bin/data/hla_reference_dna.fasta {cut[0]}")]
    elif a.tool == "star_salmon":
        # THE CHAIN LENS ACTUALLY RUNS ON BULK RNA, and the reason this row has this shape.
        # rna_quant.nf calls gffread_make_tx_fa and then salmon_aln_quant, whose script line is
        # `salmon quant --threads N -t <tx fasta> -l a -a <aln>`. The alignment is a transcriptome
        # BAM, and STAR --quantMode TranscriptomeSAM is what produces it.
        #
        # What was here before ran `salmon quant -i <index> -1 -2`, selective-alignment mode, which
        # LENS never calls -- against a path that is not a salmon index at all. It failed with
        # `versionInfo.json doesn't seem to exist`: a verdict about this harness, sitting in the
        # column where a verdict about the release belongs. Two of fourteen rows were red for that.
        #
        # Both stages in ONE task rather than two array rows, for two reasons: sibling array tasks
        # have no ordering, so a separate salmon row could start before the BAM it reads exists; and
        # STAR is what sets this stage's 48G, so running it once per dataset keeps that to one task.
        stages = [
            ("star", a.cmd,
             f"--genomeDir {a.ref} --readFilesIn {' '.join(cut)} --runThreadN {a.threads} "
             f"--outFileNamePrefix {pfx} --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM"),
            ("salmon_aln", a.cmd2 or a.cmd,
             f"quant --threads {a.threads} -t {a.ref2} -l a "
             f"-a {pfx}Aligned.toTranscriptome.out.bam -o {work}/pf_{a.label}_salmon"),
        ]
    elif a.tool == "pbmm2":
        stages = [("pbmm2",  a.cmd,
                   f"align --preset HIFI -j {a.threads} {a.ref} {cut[0]} {work}/pf_{a.label}.bam")]
    elif a.tool == "lima":
        stages = [("lima", a.cmd,
                   f"{cut[0]} {a.ref} {work}/pf_{a.label}_lima.bam --per-read "
                   f"--num-threads {a.threads}")]
    elif a.tool == "cellranger_vdj":
        # cellranger insists on its own directory layout and sample naming
        fqdir = os.path.join(work, f"pf_{a.label}_fq")
        os.makedirs(fqdir, exist_ok=True)
        for i, q in enumerate(cut, 1):
            os.replace(q, os.path.join(fqdir, f"{a.label}_S1_L001_R{i}_001.fastq"))
        stages = [("cellranger_vdj", a.cmd,
                   f"vdj --id=pf_{a.label} --reference={a.ref} --fastqs={fqdir} "
                   f"--sample={a.label} --localcores={a.threads} --localmem=16")]
    else:
        print(f"unknown tool {a.tool}", file=sys.stderr)
        return 2

    ran, rc, err, stage = [], 0, "", ""
    for name, tmpl, args in stages:
        if not tmpl:
            rc, err, stage = 2, f"no command was given for stage {name}", name
            ran.append({"stage": name, "exit_code": rc, "ok": False, "stderr_tail": err})
            break
        try:
            # cwd=work, because cellranger writes its pipestance tree into the CURRENT directory and
            # the job's cwd is the repo: the first run left pf_ds01_tcr/ and pf_ds02_tcr/ in the
            # checkout.
            r = run(tmpl.format(args=args), timeout=a.timeout, cwd=work)
            # STDOUT AS WELL AS STDERR. STAR prints `EXITING because of FATAL ERROR in reads input`
            # to stdout, so reading stderr alone would classify a STAR rejection of the DATA as a
            # setup failure -- the opposite of the distinction this script exists to make.
            srv, serr = r.returncode, ((r.stdout or "") + (r.stderr or ""))[-2000:]
        except subprocess.TimeoutExpired:
            srv, serr = -1, f"timed out after {a.timeout}s"
        ran.append({"stage": name, "exit_code": srv, "ok": srv == 0,
                    "stderr_tail": serr.strip()[-300:]})
        rc, err, stage = srv, serr, name
        # Stop at the first failure: stage two reads stage one's output, so running it anyway would
        # report a second failure that only restates the first.
        if srv != 0:
            break

    # A parse error is the signature that matters: it means the tool rejected the DATA, not the setup.
    # One signature per tool family, taken from the message each actually emits:
    #   razers3/SeqAn   seqan::ParseError: Unexpected character 'Y' found.
    #   STAR            EXITING because of FATAL ERROR in reads input: ...
    #   salmon          Invalid FASTQ / malformed
    #   cellranger      "not a valid" / unsupported
    choke = [k for k in ("ParseError", "Unexpected character", "FATAL ERROR in reads input",
                         "not a valid", "malformed", "invalid character", "Invalid FASTQ",
                         "unexpected end", "quality string length")
             if k.lower() in err.lower()]
    # YIELD, NOT JUST EXIT CODE, for the demultiplexers. lima exits 0 having found a valid primer pair
    # in one read of a thousand: it writes whatever passed and reports success, so `ok = rc == 0` called
    # a library with no primer structure at all a pass. That is this script's own documented failure mode
    # -- asking whether the tool RAN rather than whether it did anything -- committed a second time in
    # the same file. The same applies to any stage whose job is to recognise structure.
    yield_info = {}
    if rc == 0 and a.tool == "lima":
        lima_bam = f"{work}/pf_{a.label}_lima.bam"
        n_in = n_out = None
        for name, path in (("in", cut[0]), ("out", lima_bam)):
            r = run(a.samtools_cmd.format(args=f"view -c {path}"))
            if r.returncode == 0 and r.stdout.strip().isdigit():
                if name == "in":
                    n_in = int(r.stdout.strip())
                else:
                    n_out = int(r.stdout.strip())
        frac = (n_out / n_in) if (n_in and n_out is not None) else None
        # A real Iso-Seq or 10x library demultiplexes most of its reads; 0.50 is loose enough for a
        # 5 % slice of a mixed-orientation library and nowhere near 1-in-1000.
        FLOOR = 0.50
        yield_info = {"reads_in": n_in, "reads_out": n_out,
                      "fraction_demultiplexed": round(frac, 4) if frac is not None else None,
                      "floor": FLOOR}
        if frac is None:
            rc, err = 1, (err + "\nlima exited 0 but its yield could not be counted, so nothing was "
                                "established about whether it recognised the primer structure")
        elif frac < FLOOR:
            rc, err = 1, (err + f"\nlima exited 0 and demultiplexed {n_out} of {n_in} reads "
                                f"({frac:.1%}), under the {FLOOR:.0%} floor: it ran, and it did not "
                                f"recognise this library's primers")

    out = {"label": a.label, "tool": a.tool, "reads": paths, "fraction": a.fraction,
           "exit_code": rc, "ok": rc == 0, "parse_error": bool(choke), "yield": yield_info,
           "signature": choke[0] if choke else "", "stderr_tail": err.strip()[-300:],
           # Which stage the verdict came from, and what every stage did. For a one-stage tool these
           # restate `tool` and the fields above; for a chain they are the difference between "bulk
           # RNA failed" and "STAR was fine and salmon could not read its own reference".
           "stage": stage, "stages": ran}
    json.dump(out, open(a.out, "w"), indent=2)
    shown = a.tool if len(stages) == 1 else f"{a.tool}/{stage}"
    print(f"  {a.label:<22} {shown:<20} exit={rc} ok={out['ok']} parse_error={out['parse_error']}")
    if choke:
        print(f"    CHOKED: {choke[0]}")
    # EXIT NON-ZERO WHEN THE TOOL DID NOT SUCCEED. This returned 0 unconditionally, so the stage
    # reported success to the driver while ten of fourteen tools had failed -- the driver logs
    # `94_consumer_preflight OK` purely from SLURM exit status, and the verdict lived only inside a
    # JSON nobody had read yet. A verification stage that cannot fail is not a gate.
    #
    # Non-zero for a setup failure too, not just a data rejection. Distinguishing them is
    # check_consumer_preflight's job and it does; the thing that must not happen here is a tool
    # failing for ANY reason and the pipeline calling it OK. The driver treats this stage with
    # `|| say WARN`, so a failure reports loudly without halting the release.
    return 0 if rc == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
