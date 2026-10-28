"""Give a dataset's class I loci the sequence of the alleles it declares.

`design.yaml` names six class I alleles per dataset, and `hla_loss.py` builds the HLA LOH truth from them.
Nothing made the READS carry them. A haplotype is GRCh38 plus the donor's germline calls, so what the
reads hold at HLA-A, -B and -C is whatever the donor's variant caller could recover there -- which, for a
short-read donor, is not an HLA allele.

IGI-SYN-SEQ-02 is that case. Its baseline is IPISRC044's short-read WGS, and the MHC recovery
(`03_mhc_recover.sbatch`) fixed the DEPTH problem -- HLA-B went from 7-11x to 48.9x and gained 93 PASS
records, 64 of them heterozygous -- but depth is not sequence. Scored by 31-mer containment against
GRCh38's own allele contigs, the two ds-02 haplotypes reach 0.84/0.77 at HLA-B and 0.73/0.68 at HLA-C:
they lean toward the declared alleles, so the phase and the declaration are right, but the sequence is a
mosaic that is no allele at all. Every typer then picks a different nearest neighbour -- OptiType called
C*07:27/C*12:03 on the normal and C*07:05/C*12:03 on the tumour, seq2hla called C*01:02/C*07:02 on the
RNA -- and a consensus over disagreeing typers invents alleles the patient does not carry. IGI-SYN-SEQ-01
does not have the problem, because HG002's baseline is assembly-derived and reconstructs its alleles
exactly (1.000/1.000 at HLA-B).

The fix needs no new data. GRCh38's analysis set already ships 525 `HLA-*` allele contigs -- the same
FASTA the build opens for everything else, and the one `hla_phase_check.py` already reads to MEASURE
containment. This module WRITES with them: each declared allele's contig is aligned back to its own class
I window and the differences become ordinary phased germline records, so a haplotype carries the allele's
real sequence and a typer has something exact to find.

What this costs in honesty, stated plainly: ds-02's MHC stops being IPISRC044's genome over three 3.4 kb
windows and becomes synthetic. That is a real departure from "GRCh38 plus a real donor's calls" and it is
deliberate -- the alternative is a diploid-assembly baseline for IPISRC044, which does not exist. ds-01
earns its alleles rather than being handed them, and is left alone: the graft is opt-in per dataset via
`hla_graft` in `design.yaml`, not applied wherever containment happens to be low.

Outside the grafted spans the donor's own calls are untouched, so the rest of the MHC keeps its
provenance.

Every step is checked rather than assumed. The alignment must cover the whole contig; the records derived
from it must REBUILD that contig when applied to the reference span, which is what catches a cs tag read
in the wrong orientation; and the finished VCF must put each allele on the haplotype the design declares,
which `hla_phase_check.verify` measures independently.
"""
import os
import re
import subprocess
import tempfile

from .hla_phase_check import _contig_for, CLASS_I_WINDOWS

# How far outside the gene window to look for the contig. The contigs run a little past the gene at both
# ends -- HLA-B*08:01:01 aligns to chr6:31,354,120-31,357,442 against a window ending at 31,357,188 -- so a
# window with no padding would clip the alignment and leave the ends of the allele unbuilt.
PAD = 2000

# An allele contig against its own locus aligns at 0.96-0.99; the next locus over would not align at all.
# These reject a contig that matched the wrong window, not a divergent allele.
MIN_IDENTITY = 0.90
MIN_COVERED = 0.97

# Donor records this far before a grafted span are dropped too, so an indel that starts outside the span
# and reaches into it cannot survive and overlap a grafted record. 60 bp is the margin
# `Germline.overlaps_variant` already uses for the same question.
EDGE = 60


def _run(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"[hla-graft] command failed: {cmd}\n{r.stderr.strip()[-400:]}")
    return r.stdout


def _revcomp(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def _best_paf(paf):
    """The alignment with the most matching bases, or None if the contig did not align."""
    rows = [l.split("\t") for l in paf.strip().splitlines() if l.strip()]
    return max(rows, key=lambda f: int(f[9])) if rows else None


def _cs_edits(f, win_start1):
    """Differences from one PAF row as (start0_in_chr6, ref, alt), reference-oriented.

    Starts are 0-based, matching `_realize` and `_records`; `win_start1` is the 1-based chr6 coordinate of
    the window's first base.

    minimap2 walks the cs tag along the target, and for a reverse-strand hit the query bases it reports
    are already reverse-complemented, so nothing here has to know the strand. `_realize` proves that
    rather than trusting it.
    """
    cs = next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), None)
    if cs is None:
        raise SystemExit("[hla-graft] minimap2 returned no cs tag; --cs was not honoured")
    pos = (win_start1 - 1) + int(f[7])                # f[7] is a 0-based target offset into the window
    out = []
    for op, val in re.findall(r"([:*+-])([0-9]+|[A-Za-z]+)", cs):
        if op == ":":
            pos += int(val)
        elif op == "*":
            out.append((pos, val[0].upper(), val[1].upper()))
            pos += 1
        elif op == "-":                               # present in the reference, absent from the allele
            out.append((pos, val.upper(), ""))
            pos += len(val)
        elif op == "+":                               # present in the allele, absent from the reference
            out.append((pos, "", val.upper()))
    return out


def _realize(ref_seq, span_start0, edits):
    """`ref_seq` with `edits` applied, right to left so earlier offsets stay valid."""
    out = ref_seq
    for start0, r, a in sorted(edits, key=lambda e: -e[0]):
        off = start0 - span_start0
        if out[off:off + len(r)] != r:
            raise SystemExit(f"[hla-graft] edit {r or '-'}>{a or '-'} at chr6:{start0 + 1} does not match "
                             f"the reference {out[off:off + len(r)]!r}; the alignment was misread")
        out = out[:off] + a + out[off + len(r):]
    return out


def _align_allele(fasta, mm2, workdir, locus, allele, chrom, lo, hi, log):
    """Align the allele's own GRCh38 contig to its class I window and return (contig, edits, span)."""
    contig = _contig_for(fasta, allele)
    if contig is None:
        raise SystemExit(f"[hla-graft] {allele} has no contig in the reference; the graft cannot build a "
                         f"sequence for an allele the reference does not carry")
    win = os.path.join(workdir, f"win_{locus}.fa")
    alf = os.path.join(workdir, f"allele_{locus}_{allele.replace('*', '_').replace(':', '_')}.fa")
    with open(win, "w") as fh:
        fh.write(f">{chrom}:{lo}-{hi}\n{fasta.fetch(chrom, lo - 1, hi).upper()}\n")
    with open(alf, "w") as fh:
        fh.write(f">{contig}\n{fasta.fetch(contig).upper()}\n")
    f = _best_paf(_run(mm2.format(args=f"-cx asm10 --cs -t 2 {win} {alf}")))
    if f is None:
        raise SystemExit(f"[hla-graft] {allele} ({contig}) did not align to {chrom}:{lo}-{hi}")
    qlen, qs, qe, strand, nmatch, blk = int(f[1]), int(f[2]), int(f[3]), f[4], int(f[9]), int(f[10])
    identity, covered = nmatch / blk, (qe - qs) / qlen
    tspan = (lo + int(f[7]), lo + int(f[8]) - 1)      # 1-based inclusive, in chr6 coordinates
    log(f"[hla-graft] {allele:<14} {contig:<22} {strand} chr6:{tspan[0]}-{tspan[1]}  "
        f"identity {identity:.4f}  contig covered {covered:.4f}")
    if identity < MIN_IDENTITY or covered < MIN_COVERED:
        raise SystemExit(f"[hla-graft] {allele} ({contig}) aligned at identity {identity:.4f} over "
                         f"{covered:.4f} of the contig, below {MIN_IDENTITY}/{MIN_COVERED}; this is not "
                         f"the right contig for {chrom}:{lo}-{hi}")
    edits = _cs_edits(f, lo)
    # Prove the edits rebuild the allele. A cs tag read in the wrong orientation, or an off-by-one in the
    # walk, still produces a plausible-looking variant list; it does not survive this.
    ref_span = fasta.fetch(chrom, tspan[0] - 1, tspan[1]).upper()
    want = fasta.fetch(contig).upper()[qs:qe]
    if strand == "-":
        want = _revcomp(want)
    if _realize(ref_span, tspan[0] - 1, edits) != want:
        raise SystemExit(f"[hla-graft] the records derived for {allele} do not rebuild {contig} when "
                         f"applied to {chrom}:{tspan[0]}-{tspan[1]}")
    return contig, edits, tspan


def _records(fasta, chrom, by_hap_edits, span):
    """Phased VCF records over one locus from each haplotype's edit list.

    Edits that touch the same reference bases on either haplotype become ONE record spanning all of them,
    with each haplotype's realization of that span as an allele. Representing them separately would be
    legal VCF but would leave two records at one position whose REFs disagree, and the point of this
    module is to stop producing sequence that tools have to guess at.
    """
    marks = []                                        # (start0, end0) each edit occupies in the reference
    for hap, edits in by_hap_edits.items():
        for start0, r, _a in edits:
            marks.append((start0, start0 + max(len(r), 1)))
    if not marks:
        return []
    groups = []
    for s, e in sorted(marks):
        if groups and s <= groups[-1][1]:
            groups[-1][1] = max(groups[-1][1], e)
        else:
            groups.append([s, e])
    def realize(s, e):
        ref = fasta.fetch(chrom, s, e).upper()
        alleles, gt = [ref], []
        for hap in (0, 1):
            seq = _realize(ref, s, [ed for ed in by_hap_edits[hap] if s <= ed[0] < e])
            if seq not in alleles:
                alleles.append(seq)
            gt.append(alleles.index(seq))
        return ref, alleles, gt

    out = []
    for s, e in groups:
        ref, alleles, gt = realize(s, e)
        # A substitution needs no anchor base and should not carry one: writing every SNV as `GA>GG`
        # would be legal but unreadable, and no germline VCF a caller produces looks like that. Anything
        # that changes length does need one -- a deletion would otherwise have an empty ALT, and an
        # insertion written as `A>TA` shares no first base and would be read as a 1 bp MNP plus an
        # insertion by anything applying the VCF convention.
        if len({len(x) for x in alleles}) != 1 or not ref:
            if s <= span[0] - 1:
                raise SystemExit(f"[hla-graft] a length-changing difference sits at the very start of "
                                 f"the aligned span at chr6:{s + 1}; there is no anchor base for it")
            s -= 1
            ref, alleles, gt = realize(s, e)
        if gt == [0, 0]:
            continue                                  # both haplotypes match the reference here
        out.append({"chrom": chrom, "pos": s + 1, "ref": ref, "alts": alleles[1:], "gt": tuple(gt)})
    return out


def plan(fasta, dcfg, paths, workdir, log=print):
    """Every record the graft would write, plus the spans it replaces. Reads nothing, writes nothing."""
    mm2 = paths.get("minimap2_cmd")
    if not mm2:
        raise SystemExit("[hla-graft] paths.yaml has no `minimap2_cmd`; the graft aligns the declared "
                         "alleles' reference contigs and cannot run without it")
    by_hap = {int(k): list(v) for k, v in (dcfg.get("hla_haplotypes") or {}).items()}
    if set(by_hap) != {0, 1}:
        raise SystemExit("[hla-graft] `hla_haplotypes` must declare both haplotypes")
    records, spans, report = [], [], []
    for locus, (chrom, lo, hi) in sorted(CLASS_I_WINDOWS.items()):
        chosen = {h: next((a for a in by_hap[h] if a.split("*")[0] == locus), None) for h in (0, 1)}
        if chosen[0] is None or chosen[1] is None:
            raise SystemExit(f"[hla-graft] {locus} is not declared on both haplotypes")
        edits, tspans, contigs = {}, [], {}
        for hap in (0, 1):
            contig, ed, ts = _align_allele(fasta, mm2, workdir, locus, chosen[hap],
                                           chrom, lo - PAD, hi + PAD, log)
            edits[hap], contigs[hap] = ed, contig
            tspans.append(ts)
        # Only the span both alleles cover is grafted: outside it one haplotype would be the allele and the
        # other bare reference, which is the mosaic this module exists to remove.
        span = (max(t[0] for t in tspans), min(t[1] for t in tspans))
        edits = {h: [e for e in edits[h] if span[0] - 1 <= e[0] and e[0] + len(e[1]) <= span[1]]
                 for h in (0, 1)}
        recs = _records(fasta, chrom, edits, span)
        n_het = sum(1 for r in recs if r["gt"][0] != r["gt"][1])
        log(f"[hla-graft] {locus}: {chosen[0]} | {chosen[1]} over chr6:{span[0]}-{span[1]} "
            f"-> {len(recs)} records, {n_het} heterozygous")
        records += recs
        spans.append((chrom, span[0], span[1]))
        report.append({"locus": locus, "hap0": chosen[0], "hap1": chosen[1],
                       "contig0": contigs[0], "contig1": contigs[1],
                       "chrom": chrom, "start": span[0], "end": span[1],
                       "records": len(recs), "heterozygous": n_het})
    records.sort(key=lambda r: (r["chrom"], r["pos"]))
    return records, spans, report


def write_vcf(base_vcf, records, path):
    """The graft's records alone, carrying the baseline's header and sample."""
    import pysam
    src = pysam.VariantFile(base_vcf)
    hdr = src.header.copy()
    hdr.add_line('##INFO=<ID=IGIHLA,Number=0,Type=Flag,Description="IGI-SYN-SEQ class I allele graft">')
    out = pysam.VariantFile(path, "wz", header=hdr)
    sample = src.header.samples[0]
    for r in records:
        rec = out.new_record(contig=r["chrom"], start=r["pos"] - 1, alleles=(r["ref"], *r["alts"]))
        rec.info["IGIHLA"] = True
        rec.samples[sample]["GT"] = r["gt"]
        rec.samples[sample].phased = True
        out.write(rec)
    out.close()
    src.close()
    pysam.tabix_index(path, preset="vcf", force=True)
    return path


def apply_graft(base_vcf, fasta, dcfg, paths, out_vcf, workdir, dataset, log=print):
    """Write `base_vcf` with its class I windows replaced by the declared alleles. Returns the report."""
    with tempfile.TemporaryDirectory(dir=workdir) as td:
        records, spans, report = plan(fasta, dcfg, paths, td, log=log)
        graft = write_vcf(base_vcf, records, os.path.join(td, "graft.vcf.gz"))
        # ONE leading caret for the whole list: `-t ^a,^b` excludes only `a` and reads `^b` as a contig
        # name, which left 89 donor records sitting inside the HLA-B span. The span check below caught it.
        drop = "^" + ",".join(f"{c}:{s - EDGE}-{e}" for c, s, e in spans)
        kept = os.path.join(td, "kept.vcf.gz")
        _run(f"bcftools view -t {drop} -Oz -o {kept} {base_vcf} && bcftools index -t -f {kept}")
        n_base = int(_run(f"bcftools index -n {base_vcf}") or 0)
        n_kept = int(_run(f"bcftools index -n {kept}") or 0)
        _run(f"bcftools concat -a -Oz -o {out_vcf} {kept} {graft} && bcftools index -t -f {out_vcf}")
        n_out = int(_run(f"bcftools index -n {out_vcf}") or 0)
        if n_out != n_kept + len(records):
            raise SystemExit(f"[hla-graft] concat lost records: {n_kept} + {len(records)} != {n_out}")
        # The donor's calls must be gone from the grafted spans: a survivor would overlap a graft record
        # and `Germline.apply` would silently skip one of the two.
        for c, s, e in spans:
            n_here = len(_run(f"bcftools view -H -r {c}:{s}-{e} {out_vcf}").strip().splitlines())
            n_want = sum(1 for r in records if r["chrom"] == c and s <= r["pos"] <= e)
            if n_here != n_want:
                raise SystemExit(f"[hla-graft] {c}:{s}-{e} holds {n_here} records after the graft but the "
                                 f"graft wrote {n_want}; a baseline record survived inside the span")
        log(f"[hla-graft] {n_base:,} baseline records -> {n_kept:,} kept + {len(records)} grafted "
            f"= {n_out:,} in {out_vcf}")
    for r in report:
        r["dataset"] = dataset
    return report
