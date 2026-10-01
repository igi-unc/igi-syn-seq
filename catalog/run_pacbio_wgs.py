#!/usr/bin/env python3
"""PacBio HiFi whole-genome library for one dataset, one library, one chromosome.

    python3 run_pacbio_wgs.py --design design.yaml --paths paths.yaml --dataset IGI-SYN-SEQ-01 \
        --chrom chr8 --library tumor --depth 30 --out DIR --work DIR

This is the first read path that builds from `rearrange.py`. The exome path rebuilds each capture interval
in reference coordinates and represents a rearrangement as a junction contig, which is the right
compromise for 150 bp reads. A HiFi read is long enough to span a breakpoint and carry flanking sequence
on both sides, so the tumour genome has to exist as a derived chromosome: segments in their rearranged
order, with germline and somatic edits applied, and a coordinate map from derived position back to
reference so the truth stays checkable.

Depth follows the clone model exactly as the exome does, so an amplified region is deeper and a lost one
shallower; there is no capture step, so no GC or bait model applies.
"""
import argparse, json, os, subprocess, sys, time
from collections import defaultdict

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from igi_catalog.designer import build_env
from igi_catalog.genome_build import read_events
from igi_catalog.pipeline import depth_denominator, source_plan
from igi_catalog import rearrange


def sv_rows(catalog_dir, dataset, chrom):
    """Designed structural variants on this chromosome, as rearrange.plan_chromosome expects them."""
    import csv
    path = os.path.join(catalog_dir, f"{dataset}.svs.tsv")
    if not os.path.exists(path):
        return []
    out = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("chrom") != chrom or r.get("svtype") == "TRA":
                continue        # translocations need both partners; handled separately
            try:
                out.append({"event_id": r["event_id"], "svtype": r["svtype"],
                            "start": int(r["start"]), "end": int(r["end"]),
                            "clone": r.get("clone") or "T",
                            "haplotype": int(r.get("haplotype") or 0),
                            "inserted_sequence": r.get("inserted_sequence") or ""})
            except (ValueError, KeyError):
                continue
    return out


def mei_sequence_factory(rng):
    """Sequence for a mobile-element insertion. The catalog records a family and a length, not bases."""
    def make(event):
        return "".join(rng.choice("ACGT") for _ in range(300))
    return make


def split_fasta(path, chunk_bp, work, tag, overlap=60_000):
    """Split one derived chromosome into overlapping chunks so simulation parallelises past one task.

    The overlap is one maximum read length, and only the chunk that owns a region emits reads for it, so
    molecules spanning a boundary are not lost and none is duplicated.
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
    out = []
    pos = 0
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", required=True)
    ap.add_argument("--paths", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--library", choices=["tumor", "normal"], required=True)
    ap.add_argument("--depth", type=float, required=True)
    ap.add_argument("--catalog-dir", default="output")
    ap.add_argument("--background", default=None)
    ap.add_argument("--extra-events", nargs="*", default=[])
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    design = yaml.safe_load(open(a.design))
    paths = yaml.safe_load(open(a.paths))
    env = build_env(paths, design, a.dataset)
    dcfg = design["datasets"][a.dataset]
    purity = dcfg["purity"]
    import random
    rng = random.Random(f"{a.seed}:{a.dataset}:{a.chrom}:{a.library}")

    tables = [os.path.join(a.catalog_dir, f"{a.dataset}.snv_indel.tsv")]
    tables += [t for t in ([a.background] + list(a.extra_events)) if t and os.path.exists(t)]
    events = read_events(*tables)
    print(f"[{a.dataset} {a.chrom} {a.library}] somatic tables: "
          f"{', '.join(os.path.basename(t) for t in tables)}", flush=True)

    work = os.path.join(a.work, f"{a.dataset}_{a.chrom}_{a.library}")
    os.makedirs(work, exist_ok=True)
    os.makedirs(a.out, exist_ok=True)
    chrom_len = env.genome.lengths[a.chrom]
    denom = depth_denominator(env.clones, purity)
    tumor = a.library == "tumor"
    svs = sv_rows(a.catalog_dir, a.dataset, a.chrom) if tumor else []

    t0 = time.time()
    pieces, plan_rows = [], []
    for src, hap, kind, w in source_plan(env.clones, purity, a.chrom, chrom_len // 2, tumor=tumor):
        clone = "T" if (src == "NORMAL" or kind == "pre_cna") else src
        pre_only = kind == "pre_cna"
        mine = [s for s in svs
                if s["haplotype"] == hap and env.clones.is_descendant(clone, s["clone"])] if src != "NORMAL" else []
        # plan_chromosome returns (segments, applied, dropped). `dropped` is truth, not noise: an SV that
        # could not be placed -- overlapping another, or too near a chromosome end -- is absent from the
        # reads, and a truth table that still claims it would be wrong in the way F9 was.
        segs, applied, dropped = rearrange.plan_chromosome(a.chrom, chrom_len, mine)
        if dropped:
            print(f"    {src} hap{hap} {kind}: {len(dropped)} SV(s) not placeable: "
                  f"{[d[-1].get('event_id') if isinstance(d, tuple) else d for d in dropped][:4]}", flush=True)
        seq = rearrange.realise(env.genome, env.germline, events, env.clones, clone, hap, segs,
                                mei_sequence=mei_sequence_factory(rng), only_pre_cna=pre_only)
        if len(seq) < 10_000:
            continue
        tag = f"{src}_hap{hap}_{kind}"
        fa = rearrange.write_derived(os.path.join(work, f"{tag}.fa"), f"{a.chrom}_{tag}", seq)
        cmap = rearrange.write_coordinate_map(
            os.path.join(work, f"{tag}.coordmap.tsv"),
            [(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]) for r in rearrange.coordinate_map(segs)])
        cov = max(0.01, a.depth * w / denom)
        plan_rows.append({"source": src, "haplotype": hap, "copy_kind": kind, "weight": round(w, 5),
                          "coverage": round(cov, 4), "derived_length": len(seq),
                          "reference_length": chrom_len, "segments": len(segs),
                          "svs_offered": [s["event_id"] for s in mine],
                          "svs_applied": applied,
                          "svs_dropped": [(d[-1].get("event_id") if isinstance(d, tuple) else str(d))
                                          for d in dropped],
                          "fasta": fa, "coordinate_map": cmap})
        pieces.append((tag, fa, cov))
        print(f"    {tag}: {len(segs)} segments, derived {len(seq):,} bp "
              f"(reference {chrom_len:,}), {len(applied)}/{len(mine)} SVs applied, "
              f"coverage {cov:.3f}x", flush=True)
    print(f"  derived chromosomes built in {time.time() - t0:.0f}s", flush=True)

    # Badread rather than pbsim3 plus ccs. The ccs route is correct -- consensus over ten simulated
    # passes gave rq 0.99754 against the real 0.99823 -- but it took 3.61 h for 3x of chr21, which is
    # 3,247 h per library per dataset at the design's 30x. Badread reaches the same length distribution
    # 4.2x faster in one pass, so no consensus step is needed. Its identity is calibrated against the
    # real reads rather than taken from the flag, because the qscore model floors the error rate: asking
    # for 99.823 yielded 99.908, which is half the true error and would hand a caller cleaner HiFi than
    # exists.
    #
    # Each derived chromosome is split into chunks so the work parallelises past one task per chromosome.
    # Reads are independent, so the only cost is molecules that would have spanned a chunk boundary; the
    # chunks overlap by one maximum read length and the overlap region is simulated once, in the chunk
    # that owns it.
    t1 = time.time()
    br = paths["badread_cmd"]
    ident = paths.get("pacbio_identity", "99.4,99.8,0.5")
    lmean = paths.get("pbsim_length_mean", 16689)
    lsd = paths.get("pbsim_length_sd", 4593)
    chunk_bp = int(paths.get("pacbio_chunk_bp", 20_000_000))
    fq = []
    for i, (tag, fa, cov) in enumerate(pieces):
        for j, sub in enumerate(split_fasta(fa, chunk_bp, work, f"{tag}_c{j:03d}")):
            out = os.path.join(work, f"br_{tag}_{j:03d}.fq.gz")
            args = (f"simulate --reference {sub} --quantity {cov:.4f}x "
                    f"--length {lmean},{lsd} --identity {ident} "
                    f"--error_model {paths.get('pacbio_error_model', 'pacbio2021')} "
                    f"--qscore_model {paths.get('pacbio_qscore_model', 'pacbio2021')} "
                    f"--seed {a.seed + i * 100 + j} --start_adapter_seq '' --end_adapter_seq ''")
            with open(out, "wb") as fh:
                pr = subprocess.run(f"{br.format(args=args)} 2>/dev/null | gzip -1",
                                    shell=True, stdout=fh)
            if pr.returncode != 0 or os.path.getsize(out) < 100:
                raise RuntimeError(f"badread produced nothing for {tag} chunk {j}")
            fq.append(out)
    print(f"  badread done in {time.time() - t1:.0f}s, {len(fq)} chunk(s)", flush=True)

    # HiFi is delivered as an unaligned BAM, which is what a Revio run yields
    out_bam = os.path.join(a.out, f"{a.dataset}_{a.chrom}_{a.library}_hifi.bam")
    sam = paths["samtools_cmd"]
    cat = os.path.join(work, "all.fq.gz")
    with open(cat, "wb") as fh:
        subprocess.run(["cat"] + fq, stdout=fh, check=True)
    subprocess.run(sam.format(args=f"import -0 {cat} -o {out_bam}"), shell=True, check=True)
    n_reads = int(subprocess.run(sam.format(args=f"view -c {out_bam}"), shell=True,
                                 capture_output=True, text=True).stdout.strip() or 0)

    meta = {"dataset": a.dataset, "chrom": a.chrom, "library": a.library, "depth": a.depth,
            "platform": "PacBio Revio HiFi", "simulator": "badread pacbio2021",
            "reads": n_reads, "bam": out_bam, "sources": plan_rows,
            "identity": ident, "chunk_bp": chunk_bp, "chunks": len(fq),
            "runtime_s": round(time.time() - t0)}
    with open(os.path.join(a.out, f"{a.dataset}_{a.chrom}_{a.library}_hifi.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    print(f"  {n_reads:,} HiFi reads -> {out_bam} "
          f"({os.path.getsize(out_bam) / 1e6:.0f} MB)", flush=True)


if __name__ == "__main__":
    main()
