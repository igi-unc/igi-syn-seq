#!/usr/bin/env python3
"""Align a sample of a delivered library and report what a real aligner makes of it.

    align_qc.py --reads a_R1.fastq.gz[,a_R2.fastq.gz] --kind wgs --ref ref.fa --out out.json
                --aligner-cmd '<container prefix> {args}' --samtools-cmd '... {args}' [--max-reads 200000]

WHY THIS EXISTS, AND WHY IT IS DIFFERENT FROM EVERY OTHER CHECK HERE.

Every other check in this release asserts a property someone thought to assert: that the bases are ACGTN,
that the quality string is the same length as the sequence, that mates are in step, that R2 is antisense.
Each was written AFTER a consumer hit the defect it now catches -- the `*` allele, the IUPAC bases, the
unstranded R2. That pattern has an obvious weakness: it cannot catch the defect nobody has thought of.

Alignment is the one check that does not work that way. It does not test a property; it hands the reads
to the same class of tool the consumer will use and asks whether they behave like sequencing reads of
this genome. A read that is reverse complemented when it should not be, built from the wrong reference,
shifted off its coordinates, corrupted mid-record, or paired with the wrong mate does not announce which
invariant it broke -- but it fails to align, or aligns with the wrong orientation, insert size or clipping.

Reported per library, with no verdicts: thresholds belong with the release QA, and this script exists to
MEASURE so the bands can be set from data rather than guessed.

    mapped_fraction          primary alignments / reads examined
    properly_paired          for paired libraries, the samtools flag
    mean_identity            1 - NM / aligned length (M+I+D), over primary alignments
    softclip_fraction        clipped bases / total query bases; a jump means the ends are wrong
    plus_strand_fraction     a DNA library should sit near 0.5; a cDNA library need not
    mean_insert              template length for properly paired reads
    mapq0_fraction           multi-mapping share, which repeats make nonzero even in perfect data

The sample is taken from the head of the file, which is a real limitation and a deliberate one: unlike
the absence of a rare character, every one of these is a BULK property that every read contributes to, so
a head sample estimates it well. It is biased in WHERE the reads come from (the first chromosomes of a
merged library), so a per-chromosome arm should be sampled per chromosome rather than from the merge.
"""
import argparse
import json
import os
import re
import subprocess
import sys

PRESET = {            # kind -> minimap2 preset, or None to use bwa mem
    "wgs": None, "wes": None, "rna": None, "tenx_gex": None, "tenx_tcr": None,
    "ont_wgs": "map-ont", "ont_rna": "splice", "pacbio": "map-hifi", "kinnex": "splice:hq",
}


def run(cmd, **kw):
    return subprocess.run(["bash", "-c", cmd], capture_output=True, text=True, **kw)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reads", required=True, help="one path, or R1,R2 for a paired library")
    ap.add_argument("--kind", required=True, choices=sorted(PRESET))
    ap.add_argument("--ref", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--label", default="")
    ap.add_argument("--aligner-cmd", required=True, help="container prefix with {args}")
    ap.add_argument("--samtools-cmd", required=True)
    ap.add_argument("--max-reads", type=int, default=200000)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()

    paths = a.reads.split(",")
    paired = len(paths) == 2
    lines = a.max_reads * 4
    work = os.path.dirname(os.path.abspath(a.out)) or "."
    os.makedirs(work, exist_ok=True)

    # The sample is materialised rather than piped, because bwa needs to read both mates in step and a
    # process-substitution pair can deadlock on a full pipe buffer when one mate is consumed faster.
    cut = []
    for i, p in enumerate(paths, 1):
        q = os.path.join(work, f"sample_{a.label or 'lib'}_{i}.fq")
        if p.endswith(".bam"):
            r = run(f"{a.samtools_cmd.format(args=f'view {p}')} | head -{a.max_reads} | "
                    f"awk -F'\\t' '{{print \"@\"$1\"\\n\"$10\"\\n+\\n\"$11}}' > {q}")
        else:
            r = run(f"pigz -dc {p} | head -{lines} > {q}")
        if r.returncode != 0 and not os.path.getsize(q):
            print(f"could not sample {p}: {r.stderr[:200]}", file=sys.stderr)
            sys.exit(1)
        cut.append(q)

    preset = PRESET[a.kind]
    sam = os.path.join(work, f"aln_{a.label or 'lib'}.sam")
    if preset is None:
        args = f"mem -t {a.threads} {a.ref} " + " ".join(cut)
    else:
        args = f"-ax {preset} -t {a.threads} --secondary=no {a.ref} " + " ".join(cut)
    r = run(f"{a.aligner_cmd.format(args=args)} > {sam} 2>/dev/null")
    if r.returncode != 0:
        print(f"aligner failed: {r.stderr[:300]}", file=sys.stderr)
        sys.exit(2)

    # One awk pass over the SAM: flags, NM, CIGAR and TLEN are all there.
    prog = r'''
      $1 ~ /^@/ { next }
      { n++; flag = $2
        if (and(flag, 256) || and(flag, 2048)) { next }      # secondary/supplementary
        prim++
        if (and(flag, 4)) { unmapped++; next }
        mapped++
        if ($5 == 0) mapq0++
        if (!and(flag, 16)) plus++
        if (and(flag, 1)) { pairs++; if (and(flag, 2)) proper++ }
        if (and(flag, 2) && $9 > 0) { ins += $9; insn++ }
        nm = -1
        for (i = 12; i <= NF; i++) if ($i ~ /^NM:i:/) { split($i, t, ":"); nm = t[3] }
        alen = 0; clip = 0; qlen = 0
        c = $6; while (match(c, /^[0-9]+[MIDNSHP=X]/)) {
          L = substr(c, RSTART, RLENGTH - 1) + 0; op = substr(c, RSTART + RLENGTH - 1, 1)
          if (op ~ /[MID=X]/) alen += L
          if (op ~ /[MIS=X]/) qlen += L
          if (op ~ /[SH]/) clip += L
          c = substr(c, RLENGTH + 1) }
        if (nm >= 0 && alen > 0) { nmsum += nm; alensum += alen }
        clipsum += clip; qlensum += qlen + 0 }
      END { printf "%d %d %d %d %d %d %d %d %d %d %d %d\n",
              prim+0, mapped+0, unmapped+0, mapq0+0, plus+0, pairs+0, proper+0,
              nmsum+0, alensum+0, clipsum+0, qlensum+0, (insn ? int(ins/insn) : 0) }
    '''
    r = run(f"awk '{prog}' {sam}")
    f = [int(x) for x in r.stdout.split()]
    (prim, mapped, unmapped, mapq0, plus, pairs, proper, nmsum, alensum,
     clipsum, qlensum, mean_insert) = f

    out = {
        "label": a.label, "kind": a.kind, "reads": paths, "paired": paired,
        "primary": prim, "mapped": mapped, "unmapped": unmapped,
        "mapped_fraction": round(mapped / prim, 5) if prim else None,
        "properly_paired": round(proper / pairs, 5) if pairs else None,
        "mean_identity": round(1 - nmsum / alensum, 6) if alensum else None,
        "softclip_fraction": round(clipsum / qlensum, 5) if qlensum else None,
        "plus_strand_fraction": round(plus / mapped, 5) if mapped else None,
        "mapq0_fraction": round(mapq0 / mapped, 5) if mapped else None,
        "mean_insert": mean_insert or None,
    }
    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=2); fh.write("\n")
    for k, v in out.items():
        if k not in ("reads",):
            print(f"  {k:<22} {v}")
    for p in cut + [sam]:
        if os.path.exists(p):
            os.remove(p)


if __name__ == "__main__":
    main()
