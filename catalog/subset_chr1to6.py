#!/usr/bin/env python3
"""Derive the chr1to6 subset of the six whole-library arms (design section 12).

    subset_chr1to6.py map    --dataset D --release R --catalog-out O --gtf G --out-dir M
    subset_chr1to6.py filter --arm A --dataset D --release R --map M/<D>.chr1to6_records.tsv.gz

Ten of the fifteen sample types already have a chr1to6 form for free: WES, Illumina WGS, PacBio
HiFi WGS and ONT WGS are simulated one chromosome at a time, so the subset is a selection of the
chr1-chr6 shards, and 10x TCR is kept whole because TRA/TRB/TRG lie outside chr1-6. The six that
do not are the whole-library RNA arms -- bulk RNA, 10x GEX, Kinnex bulk, Kinnex sc, ONT bulk RNA
and ONT sc RNA -- which are each built as one library over the whole transcriptome.

SELECTION IS BY TRUTH MAP, NOT BY ALIGNMENT, which is a deliberate departure from the design's
wording ("derived from the full release by alignment"). Three reasons, in order of weight:

 1. It is the same rule the other ten arms use. A per-chromosome arm's chr1to6 shard contains the
    reads whose SOURCE locus is on chr1-6, because that is what was simulated. Selecting the RNA
    arms by where an aligner happens to put them would mean the subset release used two different
    definitions of "on chr1-6" depending on the assay, and a benchmark that mixes definitions is
    the kind of thing that costs a week to find later.
 2. It is exact. Every one of these libraries ships a read-level truth sidecar -- readmap,
    gex_molecules, arrays, read_map + molecules -- that names the transcript record each read came
    from. There is no MAPQ, no multimapping and no soft-clip judgement call to get wrong.
 3. It is cheap. Aligning these twelve libraries (276 GB of reads) would cost more than the
    libraries cost to simulate; streaming them costs one pass.

What alignment would add is the read whose source is on chr7 but which maps into chr1-6 anyway.
That read is excluded here. It is background, it is a small fraction, and its absence is stated in
the subset README rather than discovered.

Record -> chromosome resolution, per transcript source class:

    reference        the ENST in the record id, through the GENCODE GTF
    splice_isoform   the SPL event in <ds>.expressed.tsv (falls back to its ENST)
    erv, cta         the ERV/CTA event in <ds>.expressed.tsv
    fusion           the FUS event in <ds>.fusions.tsv; KEPT IF EITHER PARTNER is on chr1-6,
                     because a fusion transcript is one molecule and half of it is in the subset
    virus            kept unconditionally -- design section 12 keeps viral-contig reads, and a
                     viral transcript has no human chromosome to test

Kinnex segment -> molecule, which is the one non-obvious mapping. An array is built as
adapter_0 m_0 adapter_1 m_1 ... adapter_k m_k(closing), so `skera split` returns one segment per
bracketed molecule and the segment's molecule index is min(dl, dr) of its two adapter tags. Taking
min rather than dl makes it orientation-agnostic: about half the arrays are sequenced reverse, and
those come back from skera with dl > dr (verified on both libraries, both orientations, and on
short arrays). The single-cell arrays carry sixteen cDNAs but only nine adapters, so the LAST
bracketed segment holds molecules 7..15 concatenated; it is kept if any of them is on chr1-6,
which is why the Kinnex sc subset retains more off-target molecules than the other five arms. That
is the known MAS-16 adapter-layout caveat showing up in the subset, not a new defect, and the
per-arm JSON reports the molecule-level purity so the size of it is on the record.
"""
import argparse
import gzip
from array import array as _array
import io
import json
import os
import shutil
import subprocess
import sys
import time


def arrayi():
    return _array("i")

CHR1TO6 = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6"}
ARMS = ("rna", "tenx_gex", "kinnex_bulk", "kinnex_sc", "ont_bulk_rna", "ont_sc_rna")


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ---------------------------------------------------------------------------
# streaming helpers. The deliverables run to 49 GB gzipped, so every read and write goes through
# pigz rather than the gzip module; python's decompressor is the bottleneck otherwise.

def pigz_in(path, threads=4):
    p = subprocess.Popen(["pigz", "-dc", "-p", str(threads), path],
                         stdout=subprocess.PIPE, bufsize=1 << 22)
    return p, p.stdout


def pigz_out(path, threads=6):
    fh = open(path, "wb")
    p = subprocess.Popen(["pigz", "-p", str(threads)], stdin=subprocess.PIPE,
                         stdout=fh, bufsize=1 << 22)
    return p, fh


def close_in(p):
    p.stdout.close()
    if p.wait() not in (0,):
        raise RuntimeError(f"pigz -dc exited {p.returncode}")


def close_out(p, fh):
    p.stdin.close()
    rc = p.wait()
    fh.close()
    if rc != 0:
        raise RuntimeError(f"pigz exited {rc}")


def text(stream):
    return io.TextIOWrapper(stream, encoding="utf-8", newline="\n")


def rec_num(rec):
    """r000133899 -> 133899. Record ids are dense, so the keep set is a bytearray, not a hash."""
    return int(rec[1:])


# ---------------------------------------------------------------------------
# the record -> chromosome map

def enst_chroms(gtf):
    """transcript_id (versioned) -> chrom, from every transcript line in the GTF."""
    out = {}
    p, raw = pigz_in(gtf, threads=6)
    for line in text(raw):
        if line[0] == "#":
            continue
        f = line.split("\t", 8)
        if len(f) < 9 or f[2] != "transcript":
            continue
        i = f[8].find('transcript_id "')
        if i < 0:
            continue
        j = f[8].index('"', i + 15)
        out[f[8][i + 15:j]] = f[0]
    close_in(p)
    return out


def event_chroms(catalog_out, dataset):
    """event_id -> set of chroms, from the truth tables that carry a chrom column."""
    out = {}

    def add(event, chrom):
        if chrom:
            out.setdefault(event, set()).add(chrom)

    def table(name):
        return os.path.join(catalog_out, f"{dataset}.{name}.tsv")

    for name, cols in (("expressed", ("chrom",)),
                       ("splice_causal", ("chrom",)),
                       ("viruses", ("chrom",)),
                       ("fusions", ("chrom_5p", "chrom_3p"))):
        path = table(name)
        if not os.path.exists(path):
            log(f"  WARNING: {os.path.basename(path)} absent; its events will be unresolved")
            continue
        with open(path) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            idx = {c: hdr.index(c) for c in cols if c in hdr}
            if len(idx) != len(cols):
                raise RuntimeError(f"{path}: expected columns {cols}, header has {hdr[:12]}")
            ev = hdr.index("event_id")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                for c in cols:
                    add(f[ev], f[idx[c]])
    return out


def build_map(a):
    tx_path = os.path.join(a.release, "rna", a.dataset, f"{a.dataset}_full_rna_transcripts.tsv")
    if not os.path.exists(tx_path):
        raise SystemExit(f"no transcript manifest at {tx_path}")
    log(f"GTF transcripts from {os.path.basename(a.gtf)}")
    enst = enst_chroms(a.gtf)
    log(f"  {len(enst):,} transcripts")
    ev = event_chroms(a.catalog_out, a.dataset)
    log(f"  {len(ev):,} designed events with a chromosome")

    os.makedirs(a.out_dir, exist_ok=True)
    out_path = os.path.join(a.out_dir, f"{a.dataset}.chr1to6_records.tsv.gz")
    n = kept = 0
    by_source = {}
    unresolved = {}
    p, fh = pigz_out(out_path, threads=4)
    w = p.stdin
    w.write(b"rec\tsource\tchroms\tkeep\n")
    with open(tx_path) as tfh:
        hdr = tfh.readline().rstrip("\n").split("\t")
        i_rec, i_id, i_src = hdr.index("rec"), hdr.index("id"), hdr.index("source")
        for line in tfh:
            f = line.rstrip("\n").split("\t")
            rec, rid, src = f[i_rec], f[i_id], f[i_src]
            parts = rid.split("|")
            chroms = set()
            if src == "virus":
                chroms = {"viral"}
            elif src == "reference":
                c = enst.get(parts[0])
                if c:
                    chroms = {c}
            else:
                chroms = set(ev.get(parts[0], ()))
                if not chroms and src == "splice_isoform" and len(parts) > 1:
                    c = enst.get(parts[1])
                    if c:
                        chroms = {c}
                if not chroms and src == "erv" and len(parts) > 1:
                    # Hsap38.chrX.83694121.83694447.- carries its own locus
                    bits = parts[1].split(".")
                    if len(bits) > 1 and bits[1].startswith("chr"):
                        chroms = {bits[1]}
            keep = bool(chroms & CHR1TO6) or src == "virus"
            if not chroms:
                unresolved[src] = unresolved.get(src, 0) + 1
            n += 1
            kept += keep
            s = by_source.setdefault(src, [0, 0])
            s[0] += 1
            s[1] += keep
            w.write(f"{rec}\t{src}\t{','.join(sorted(chroms)) or 'NA'}\t{int(keep)}\n".encode())
    close_out(p, fh)

    log(f"{a.dataset}: {kept:,} of {n:,} records on chr1-6 ({kept / n:.1%}) -> "
        f"{os.path.basename(out_path)}")
    for src in sorted(by_source):
        tot, k = by_source[src]
        log(f"    {src:<16} {k:>8,} / {tot:>9,}  ({k / tot:.1%})")
    if unresolved:
        # A record with no chromosome is DROPPED, so an unresolved class would silently shrink the
        # subset. Fail rather than ship that: the only class allowed to be chromosome-less is virus,
        # and it is handled above.
        raise SystemExit(f"{sum(unresolved.values()):,} records could not be resolved to a "
                         f"chromosome: {unresolved}")
    return out_path


def load_keep(path):
    """keep set as a bytearray indexed by record number."""
    p, raw = pigz_in(path, threads=4)
    keep = bytearray(1 << 21)
    hi = 0
    n = 0
    for line in text(raw):
        f = line.rstrip("\n").split("\t")
        if f[0] == "rec":
            continue
        i = rec_num(f[0])
        if i >= len(keep):
            keep.extend(b"\0" * (i + 1 - len(keep)))
        hi = max(hi, i)
        if f[3] == "1":
            keep[i] = 1
            n += 1
    close_in(p)
    log(f"  keep set: {n:,} records of {hi + 1:,}")
    return keep


# ---------------------------------------------------------------------------
# per-arm filters. Every one of these streams the deliverable and its truth sidecar in LOCKSTEP and
# asserts the read names agree, rather than trusting that they do: the sidecars were written in
# emission order, but an arm that was ever resumed or re-merged could have broken that, and a
# silent misalignment here would mislabel every read downstream.

class Lockstep:
    def __init__(self, what):
        self.what = what
        self.n = 0

    def check(self, fastq_name, side_name):
        self.n += 1
        if fastq_name != side_name:
            raise RuntimeError(
                f"{self.what}: read {self.n} is {fastq_name!r} in the reads but {side_name!r} in "
                f"the truth sidecar; they are not in the same order")


def fastq_records(stream):
    """Yield (name_bytes, record_bytes) for a FASTQ stream."""
    rl = stream.readline
    while True:
        h = rl()
        if not h:
            return
        s, p, q = rl(), rl(), rl()
        if not q:
            raise RuntimeError("truncated FASTQ: incomplete 4-line record")
        name = h[1:h.index(b"\n")].split(b" ", 1)[0].split(b"\t", 1)[0]
        yield name, h + s + p + q


def filter_paired(a, keep, src_r1, src_r2, side_path, side_read_col, side_rec_col,
                  out_r1, out_r2, out_side, what):
    """Illumina paired arm: one sidecar row per read pair, naming its transcript record."""
    pi1, r1 = pigz_in(src_r1, threads=4)
    pi2, r2 = pigz_in(src_r2, threads=4)
    pis, rs = pigz_in(side_path, threads=2)
    po1, fh1 = pigz_out(out_r1, threads=5)
    po2, fh2 = pigz_out(out_r2, threads=5)
    pos, fhs = pigz_out(out_side, threads=2)

    hdr = rs.readline()
    pos.stdin.write(hdr)
    cols = hdr.decode().rstrip("\n").split("\t")
    i_read, i_rec = cols.index(side_read_col), cols.index(side_rec_col)

    ls = Lockstep(what)
    n = w = 0
    g1, g2 = fastq_records(r1), fastq_records(r2)
    for sline in rs:
        sf = sline.rstrip(b"\n").split(b"\t")
        name1, rec1 = next(g1)
        name2, rec2 = next(g2)
        if name1 != name2:
            raise RuntimeError(f"{what}: R1 {name1!r} against R2 {name2!r}; mates out of order")
        ls.check(name1, sf[i_read])
        n += 1
        if keep[rec_num(sf[i_rec].decode())]:
            po1.stdin.write(rec1)
            po2.stdin.write(rec2)
            pos.stdin.write(sline)
            w += 1
        if n % 20_000_000 == 0:
            log(f"    {n:,} pairs, {w:,} kept ({w / n:.1%})")
    for g in (g1, g2):
        leftover = next(g, None)
        if leftover is not None:
            raise RuntimeError(f"{what}: reads outlast the truth sidecar at {leftover[0]!r}")

    close_in(pi1); close_in(pi2); close_in(pis)
    close_out(po1, fh1); close_out(po2, fh2); close_out(pos, fhs)
    return n, w


def filter_single(a, keep, src_fq, read_map, molecules, out_fq, out_read_map, out_molecules, what):
    """ONT arm: read_map names the molecule, molecules names the transcript record."""
    mol_keep = bytearray(1 << 21)
    mol_rec = {}
    pim, rm = pigz_in(molecules, threads=4)
    hdr_mol = rm.readline()
    cols = hdr_mol.decode().rstrip("\n").split("\t")
    i_mol, i_rec = cols.index("molecule"), cols.index("record")
    n_mol = 0
    for line in rm:
        f = line.rstrip(b"\n").split(b"\t")
        i = int(f[i_mol][1:])
        if i >= len(mol_keep):
            mol_keep.extend(b"\0" * (i + 1 - len(mol_keep)))
        mol_keep[i] = 1 if keep[rec_num(f[i_rec].decode())] else 2
        n_mol += 1
    close_in(pim)
    log(f"    {n_mol:,} molecules, {sum(1 for b in mol_keep if b == 1):,} on chr1-6")

    pif, rf = pigz_in(src_fq, threads=6)
    pir, rr = pigz_in(read_map, threads=3)
    pof, fhf = pigz_out(out_fq, threads=8)
    por, fhr = pigz_out(out_read_map, threads=2)
    hdr = rr.readline()
    por.stdin.write(hdr)
    cols = hdr.decode().rstrip("\n").split("\t")
    i_read, i_src = cols.index("read_name"), cols.index("source")

    ls = Lockstep(what)
    n = w = 0
    g = fastq_records(rf)
    for line in rr:
        f = line.rstrip(b"\n").split(b"\t")
        name, rec = next(g)
        ls.check(name, f[i_read])
        n += 1
        mol = int(f[i_src].split(b",", 1)[0][1:])
        if mol_keep[mol] == 0:
            raise RuntimeError(f"{what}: read {name!r} names molecule m{mol:09d}, which is not in "
                               f"the molecule table")
        if mol_keep[mol] == 1:
            pof.stdin.write(rec)
            por.stdin.write(line)
            w += 1
        if n % 5_000_000 == 0:
            log(f"    {n:,} reads, {w:,} kept ({w / n:.1%})")
    leftover = next(g, None)
    if leftover is not None:
        raise RuntimeError(f"{what}: reads outlast the read map at {leftover[0]!r}")
    close_in(pif); close_in(pir)
    close_out(pof, fhf); close_out(por, fhr)

    # the molecule table, filtered to the molecules that survived
    pim, rm = pigz_in(molecules, threads=3)
    pom, fhm = pigz_out(out_molecules, threads=2)
    pom.stdin.write(rm.readline())
    kept_mol = 0
    for line in rm:
        f = line.rstrip(b"\n").split(b"\t")
        if mol_keep[int(f[i_mol][1:])] == 1:
            pom.stdin.write(line)
            kept_mol += 1
    close_in(pim); close_out(pom, fhm)
    return n, w, n_mol, kept_mol


def filter_kinnex(a, keep, arm, src_bam, arrays, out_bam, out_arrays, what):
    """Kinnex arm: one BAM record per skera segment; the array table names each molecule.

    Built as adapter_0 m_0 adapter_1 m_1 ... so a segment's molecule index is min(dl, dr) -- see the
    module docstring on orientation and on the single-cell tail segment.
    """
    # zmw -> record numbers, ordered by segment_index. array('i') rather than a list: the single-cell
    # library has 24.7 M molecules over 1.5 M arrays, and boxed ints cost about 900 MB more.
    members = {}
    pia, ra = pigz_in(arrays, threads=4)
    hdr_a = ra.readline()
    cols = hdr_a.decode().rstrip("\n").split("\t")
    i_rn, i_si, i_rec = cols.index("read_name"), cols.index("segment_index"), cols.index("record")
    n_mol = n_unread = 0
    for line in ra:
        f = line.rstrip(b"\n").split(b"\t")
        rn = f[i_rn]
        n_mol += 1
        if not rn:
            n_unread += 1          # an array that never produced a read
            continue
        zmw = int(rn.split(b"/")[1])
        lst = members.get(zmw)
        if lst is None:
            lst = members[zmw] = arrayi()
        si = int(f[i_si])
        while len(lst) <= si:
            lst.append(-1)
        lst[si] = rec_num(f[i_rec].decode())
    close_in(pia)
    log(f"    {n_mol:,} molecules over {len(members):,} sequenced arrays "
        f"({n_unread:,} molecules in arrays that produced no read)")

    sam = a.samtools_cmd
    rp = subprocess.Popen(sam.format(args=f"view -h {src_bam}"), shell=True,
                          stdout=subprocess.PIPE, bufsize=1 << 22)
    wp = subprocess.Popen(sam.format(args=f"view -b -o {out_bam} -"), shell=True,
                          stdin=subprocess.PIPE, bufsize=1 << 22)
    poa, fha = pigz_out(out_arrays, threads=3)
    poa.stdin.write(hdr_a)

    n = w = 0
    mol_seen = mol_kept = n_unattributed = 0
    kept_pairs = set()            # zmw * 32 + molecule index, packed: 4.8 M tuples cost 290 MB
    for line in rp.stdout:
        if line[:1] == b"@":
            wp.stdin.write(line)
            continue
        name = line[:line.index(b"\t")]
        dl = tag_int(line, b"\tdl:i:", name)
        dr = tag_int(line, b"\tdr:i:", name)
        zmw = int(name.split(b"/")[1])
        lst = members.get(zmw)
        if lst is None:
            raise RuntimeError(f"{what}: segment {name!r} is on a ZMW the array table does not list")
        idx = dl if dl < dr else dr
        if idx >= 32:
            # kept_pairs packs (zmw, idx) as zmw * 32 + idx; nine adapters cannot give idx > 7, so this
            # only fires if the array grammar changes, and it must fire rather than collide silently.
            raise RuntimeError(f"{what}: segment {name!r} has adapter index {idx}, over the packing limit")
        n += 1
        if idx >= len(lst):
            # skera read an adapter index past the end of this array's molecule list, which happens on
            # about 0.02 % of segments: HiFi error on a 17 bp adapter occasionally makes adapter j look
            # like j+1 on a short or degraded array. The segment cannot be attributed to a molecule, so
            # it is KEPT and counted. Keeping it adds a little background; dropping it would remove
            # reads we cannot account for, and an unexplained hole in a benchmark library is worse than
            # an unexplained extra read. The count is reported, so the rate stays visible.
            n_unattributed += 1
            wp.stdin.write(line)
            w += 1
            continue
        # the last bracketed segment of a single-cell array holds every remaining molecule
        span = lst[idx:] if idx == min(8, len(lst)) - 1 else [lst[idx]]
        mol_seen += len(span)
        hit = [r for r in span if r >= 0 and keep[r]]
        if hit:
            wp.stdin.write(line)
            w += 1
            mol_kept += len(span)
            kept_pairs.add(zmw * 32 + idx)
        if n % 2_000_000 == 0:
            log(f"    {n:,} segments, {w:,} kept ({w / n:.1%})")
    rp.stdout.close()
    if rp.wait() != 0:
        raise RuntimeError(f"samtools view -h {src_bam} exited {rp.returncode}")
    wp.stdin.close()
    if wp.wait() != 0:
        raise RuntimeError(f"samtools view -b exited {wp.returncode}")

    # The array table for the subset: the molecules carried by the segments that survived. A molecule
    # that rode along in a kept tail segment is listed, because it IS in the reads -- the alternative
    # would be a truth table that disagrees with the library it describes.
    pia, ra = pigz_in(arrays, threads=3)
    ra.readline()
    on_target = off_target = 0
    for line in ra:
        f = line.rstrip(b"\n").split(b"\t")
        rn = f[i_rn]
        if not rn:
            continue
        zmw = int(rn.split(b"/")[1])
        si = int(f[i_si])
        lst = members.get(zmw)
        bracket = min(8, len(lst)) - 1
        seg = si if si <= bracket else bracket
        if zmw * 32 + seg in kept_pairs:
            poa.stdin.write(line)
            if keep[rec_num(f[i_rec].decode())]:
                on_target += 1
            else:
                off_target += 1
    close_in(pia)
    close_out(poa, fha)
    if n_unattributed:
        log(f"    {n_unattributed:,} segments ({n_unattributed / n:.3%}) had an adapter index past the "
            f"end of their array and were kept unattributed")
    return n, w, on_target, off_target, n_unattributed


def tag_int(line, tag, name):
    i = line.rfind(tag)
    if i < 0:
        raise RuntimeError(f"segment {name!r} carries no {tag.strip().decode()} tag; skera is "
                           f"expected to write dl and dr on every segment")
    j = i + len(tag)
    k = j
    while line[k] not in (9, 10):
        k += 1
    return int(line[j:k])


def pbindex(a, bam):
    if not a.pbindex_cmd:
        log(f"    no pbindex command configured; {os.path.basename(bam)}.pbi NOT written")
        return False
    subprocess.run(a.pbindex_cmd.format(args=bam), shell=True, check=True)
    if not os.path.exists(bam + ".pbi"):
        raise RuntimeError(f"pbindex reported success but {bam}.pbi is absent")
    return True


# ---------------------------------------------------------------------------

def run_filter(a):
    keep = load_keep(a.map)
    d, arm = a.dataset, a.arm
    src = os.path.join(a.release, arm, d)
    out = src if a.out_dir is None else os.path.join(a.out_dir, arm, d)
    os.makedirs(out, exist_ok=True)
    f = lambda n: os.path.join(src, n)
    o = lambda n: os.path.join(out, n)
    meta = {"dataset": d, "arm": arm, "release": "chr1to6",
            "selection": "truth map (source transcript locus), not alignment",
            "chroms": sorted(CHR1TO6), "viral_kept": True}
    t0 = time.time()

    if arm == "rna":
        n, w = filter_paired(
            a, keep,
            f(f"{d}_full_rna_R1.fastq.gz"), f(f"{d}_full_rna_R2.fastq.gz"),
            f(f"{d}_full_rna_readmap.tsv.gz"), "read_name", "source_record",
            o(f"{d}_chr1to6_rna_R1.fastq.gz"), o(f"{d}_chr1to6_rna_R2.fastq.gz"),
            o(f"{d}_chr1to6_rna_readmap.tsv.gz"), "bulk RNA")
        meta.update(pairs_full=n, pairs_written=w, fraction=round(w / n, 4))
        subset_transcripts(f(f"{d}_full_rna_transcripts.tsv"),
                           o(f"{d}_chr1to6_rna_transcripts.tsv"), keep)

    elif arm == "tenx_gex":
        n, w = filter_paired(
            a, keep,
            f(f"{d}-GEX_S1_L001_R1_001.fastq.gz"), f(f"{d}-GEX_S1_L001_R2_001.fastq.gz"),
            f(f"{d}_full_gex_molecules.tsv.gz"), "read", "record",
            o(f"{d}-chr1to6-GEX_S1_L001_R1_001.fastq.gz"),
            o(f"{d}-chr1to6-GEX_S1_L001_R2_001.fastq.gz"),
            o(f"{d}_chr1to6_gex_molecules.tsv.gz"), "10x GEX")
        meta.update(read_pairs_full=n, read_pairs_written=w, fraction=round(w / n, 4))
        # Every barcode keeps its chr1-6 molecules, so the cell table carries over unchanged; the
        # UMI depth per cell drops with the fraction above, which is what design section 12 predicts.
        for side in ("cells", "clonotypes"):
            p = f(f"{d}_full_{side}.tsv")
            if os.path.exists(p):
                shutil.copyfile(p, o(f"{d}_chr1to6_{side}.tsv"))
        meta["barcodes_retained"] = barcodes(o(f"{d}_chr1to6_gex_molecules.tsv.gz"))

    elif arm in ("ont_bulk_rna", "ont_sc_rna"):
        n, w, nm, km = filter_single(
            a, keep,
            f(f"{d}_full_{arm}.fastq.gz"),
            f(f"{d}_full_{arm}_read_map.tsv.gz"),
            f(f"{d}_full_{arm}_molecules.tsv.gz"),
            o(f"{d}_chr1to6_{arm}.fastq.gz"),
            o(f"{d}_chr1to6_{arm}_read_map.tsv.gz"),
            o(f"{d}_chr1to6_{arm}_molecules.tsv.gz"), arm)
        meta.update(reads_full=n, reads_written=w, fraction=round(w / n, 4),
                    molecules_full=nm, molecules_written=km)
        if arm == "ont_sc_rna":
            meta["barcodes_retained"] = barcodes(o(f"{d}_chr1to6_{arm}_molecules.tsv.gz"),
                                                 col="barcode")

    elif arm in ("kinnex_bulk", "kinnex_sc"):
        bam = o(f"{d}_chr1to6_{arm}_segmented.bam")
        n, w, on, off, unattr = filter_kinnex(
            a, keep, arm,
            f(f"{d}_full_{arm}_segmented.bam"),
            f(f"{d}_full_{arm}_arrays.tsv.gz"),
            bam, o(f"{d}_chr1to6_{arm}_arrays.tsv.gz"), arm)
        meta.update(segments_full=n, segments_written=w, fraction=round(w / n, 4),
                    molecules_on_target=on, molecules_off_target=off,
                    molecule_purity=round(on / (on + off), 4) if on + off else None,
                    segments_unattributed=unattr)
        meta["pbi"] = pbindex(a, bam)
        subset_transcripts(f(f"{d}_full_{arm}_transcripts.tsv"),
                           o(f"{d}_chr1to6_{arm}_transcripts.tsv"), keep)
        if arm == "kinnex_sc":
            for side in ("cells", "clonotypes"):
                p = f(f"{d}_full_{side}.tsv")
                if os.path.exists(p):
                    shutil.copyfile(p, o(f"{d}_chr1to6_{side}.tsv"))
    else:
        raise SystemExit(f"unknown arm {arm!r}; one of {ARMS}")

    meta["runtime_s"] = round(time.time() - t0)
    jpath = o(f"{d}_chr1to6_{arm}.json")
    with open(jpath, "w") as fh:
        json.dump(meta, fh, indent=2)
        fh.write("\n")
    log(f"{d} {arm}: {json.dumps({k: v for k, v in meta.items() if k not in ('chroms',)})}")
    return jpath


def subset_transcripts(src, dst, keep):
    if not os.path.exists(src):
        return
    n = 0
    with open(src) as fh, open(dst, "w") as out:
        out.write(fh.readline())
        for line in fh:
            if keep[rec_num(line.split("\t", 1)[0])]:
                out.write(line)
                n += 1
    log(f"    transcript manifest: {n:,} records -> {os.path.basename(dst)}")


def barcodes(path, col="barcode"):
    p, raw = pigz_in(path, threads=3)
    cols = raw.readline().decode().rstrip("\n").split("\t")
    i = cols.index(col)
    seen = set()
    for line in raw:
        seen.add(line.split(b"\t")[i])
    close_in(p)
    return len(seen)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    m = sub.add_parser("map", help="build the record -> chromosome map for a dataset")
    m.add_argument("--dataset", required=True)
    m.add_argument("--release", required=True)
    m.add_argument("--catalog-out", required=True)
    m.add_argument("--gtf", required=True)
    m.add_argument("--out-dir", required=True)

    g = sub.add_parser("filter", help="derive one arm's chr1to6 library")
    g.add_argument("--arm", required=True, choices=ARMS)
    g.add_argument("--dataset", required=True)
    g.add_argument("--release", required=True)
    g.add_argument("--map", required=True)
    g.add_argument("--out-dir", default=None,
                   help="default: alongside the full library, which is where the label distinguishes them")
    g.add_argument("--samtools-cmd", default="samtools {args}")
    g.add_argument("--pbindex-cmd", default="")

    a = ap.parse_args()
    if a.cmd == "map":
        build_map(a)
    else:
        run_filter(a)


if __name__ == "__main__":
    main()
