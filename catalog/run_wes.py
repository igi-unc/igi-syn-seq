#!/usr/bin/env python3
"""Build one exome library (tumour or normal) for one dataset and one chromosome.

    python3 run_wes.py --design design.yaml --paths paths.yaml --dataset IGI-SYN-SEQ-01 \
        --chrom chr1 --library tumor --depth 150 --out /path/out
"""
import argparse, json, os, random, subprocess, time
import yaml

from igi_catalog.designer import build_env
from igi_catalog.genome_build import read_events
from igi_catalog.junctions import fusion_junctions, sv_junctions, viral_junctions
from igi_catalog.pipeline import (GcBias, WesBuilder, merged_capture, off_target_bands,
                                 quality_model, write_record_map)
from igi_catalog.readnames import NAME_SPACE_STRIDE, shuffle_and_rename


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", required=True)
    ap.add_argument("--paths", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--library", choices=["tumor", "normal"], required=True)
    ap.add_argument("--depth", type=float, required=True)
    ap.add_argument("--catalog", required=True, help="the dataset's snv_indel.tsv")
    ap.add_argument("--background", default=None, help="the dataset's background.tsv (passenger mutations)")
    ap.add_argument("--extra-events", nargs="*", default=[],
                    help="further somatic tables to apply, e.g. splice_causal.tsv. Generating such a "
                         "table is not enough: it has to be handed to the builder or the variants it "
                         "describes never reach the sequence.")
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--max-intervals", type=int, default=None)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--no-off-target", action="store_true",
                    help="build only the bait set; the library then has no reads outside it")
    a = ap.parse_args()

    design = yaml.safe_load(open(a.design))
    paths = yaml.safe_load(open(a.paths))
    env = build_env(paths, design, a.dataset)
    purity = design["datasets"][a.dataset]["purity"]
    tables = [t for t in ([a.catalog, a.background] + list(a.extra_events)) if t and os.path.exists(t)]
    events = read_events(*tables)
    print(f"  somatic tables applied: {', '.join(os.path.basename(t) for t in tables)}", flush=True)
    capture = merged_capture(paths["exome_bed"], {a.chrom})
    work = os.path.join(a.work, f"{a.dataset}_{a.chrom}_{a.library}")
    os.makedirs(work, exist_ok=True)
    os.makedirs(a.out, exist_ok=True)

    t0 = time.time()
    gc = GcBias(paths.get("gc_bias_curve"))
    print(f"  GC bias model: {'measured, ' + str(len(gc.factors)) + ' bins' if gc.enabled else 'flat'}", flush=True)
    wb = WesBuilder(env, events, purity, work, gc_bias=gc)
    # the chromosome ordinal namespaces this run's record ids so the merged record map is unambiguous
    _ord = {f"chr{i}": i for i in range(1, 23)}
    _ord.update({"chrX": 23, "chrY": 24, "chrM": 25})
    wb.record_prefix = f"{_ord.get(a.chrom, 0):02d}"
    # structural variants, fusion breakpoints and viral integrations reach the DNA as junction contigs
    if a.library == "tumor":
        import csv as _csv

        def _rows(name):
            path = os.path.join(os.path.dirname(a.catalog), f"{a.dataset}.{name}.tsv")
            if not os.path.exists(path):
                return []
            with open(path) as fh:
                return list(_csv.DictReader(fh, delimiter="\t"))

        jn = (sv_junctions(env.genome, _rows("svs"))
              + fusion_junctions(env.genome, _rows("fusions"))
              + viral_junctions(env.genome, paths.get("viral_fasta"), _rows("viruses")))
        jn = [j for j in jn if j.get("chrom") == a.chrom or j.get("end_chrom") == a.chrom]
        wb.index_junctions(jn)
        print(f"  {len(jn)} junction contigs on {a.chrom}", flush=True)
    sources, n_iv = wb.write_source_fastas(capture, tumor=(a.library == "tumor"),
                                           max_intervals=a.max_intervals)
    print(f"[{a.dataset} {a.chrom} {a.library}] {n_iv} intervals, {len(sources)} sources, "
          f"{time.time() - t0:.0f}s", flush=True)
    for key, (path, w) in sorted(sources.items()):
        print(f"    {key[0]} hap{key[1]} {key[2]} cn{key[3]} gc{key[4]:02d}: weight {w:.4f}  "
              f"{os.path.getsize(path) / 1e6:.1f} MB", flush=True)

    t1 = time.time()
    qual = quality_model(paths)
    print(f"  quality model: {qual}", flush=True)
    pieces = wb.simulate(sources, paths["art_cmd"], a.depth,
                         f"{a.dataset}_{a.chrom}_{a.library}", seed=a.seed, qual=qual)

    # Off-target: coverage decaying away from each bait plus a thin genome-wide background. Without it
    # every base outside the bait set sits at exactly zero depth, which no real capture library shows.
    bands = []
    if not a.no_off_target:
        rng = random.Random(f"{a.seed}:{a.dataset}:{a.chrom}:offtarget")
        for name, rel, spans in off_target_bands(capture, a.chrom, env.genome.lengths[a.chrom], rng):
            bsrc, b_iv = wb.write_source_fastas({a.chrom: spans}, tumor=(a.library == "tumor"),
                                                max_intervals=a.max_intervals, tag=f"_{name}")
            if not bsrc:
                continue
            pieces += wb.simulate(bsrc, paths["art_cmd"], a.depth * rel,
                                  f"{a.dataset}_{a.chrom}_{a.library}_{name}",
                                  seed=a.seed + 1000 + len(bands), qual=qual)
            # count only the intervals actually written, which `--max-intervals` may have capped
            used = spans[:b_iv] if a.max_intervals else spans
            bases = sum(e - s + 1 for s, e in used)
            bands.append({"band": name, "relative_depth": rel, "intervals": b_iv, "bases": bases})
            print(f"  off-target {name}: {b_iv} intervals, {bases / 1e6:.1f} Mb, "
                  f"{a.depth * rel:.2f}x", flush=True)
    print(f"  ART done in {time.time() - t1:.0f}s", flush=True)

    # concatenate, then shuffle and give Illumina-style names in one disk-based pass, writing the map
    # from each read name back to its source record (owner decision D2)
    r1 = os.path.join(a.out, f"{a.dataset}_{a.chrom}_{a.library}_R1.fastq.gz")
    r2 = os.path.join(a.out, f"{a.dataset}_{a.chrom}_{a.library}_R2.fastq.gz")
    rmap = os.path.join(a.out, f"{a.dataset}_{a.chrom}_{a.library}_readmap.tsv.gz")
    recmap = os.path.join(a.out, f"{a.dataset}_{a.chrom}_{a.library}_records.tsv.gz")
    cat1 = os.path.join(work, "all_1.fq")
    cat2 = os.path.join(work, "all_2.fq")
    for idx, out in ((0, cat1), (1, cat2)):
        with open(out, "w") as fh:
            subprocess.run(["cat"] + [p[idx] for p in pieces if os.path.exists(p[idx])],
                           stdout=fh, check=True)
    # each chromosome takes its own slice of the Illumina name space, so the per-chromosome libraries
    # concatenate into one library without colliding names
    ordinals = {f"chr{i}": i for i in range(1, 23)}
    ordinals.update({"chrX": 23, "chrY": 24, "chrM": 25})
    offset = ordinals.get(a.chrom, 0) * NAME_SPACE_STRIDE
    n_reads = shuffle_and_rename(cat1, cat2, r1, r2, rmap, work, seed=a.seed, index_offset=offset)
    write_record_map(wb.record_map, recmap)
    if wb.junctions_used:
        jpath = os.path.join(a.out, f"{a.dataset}_{a.chrom}_{a.library}_junctions.tsv")
        with open(jpath, "w") as fh:
            fh.write("event_id\tchrom\tinterval_start\tinterval_end\tclone\thaplotype\n")
            for row in wb.junctions_used:
                fh.write("\t".join(str(x) for x in row) + "\n")
        print(f"  {len(wb.junctions_used)} junctions placed -> {os.path.basename(jpath)}", flush=True)
    for f in (cat1, cat2):
        if os.path.exists(f):
            os.remove(f)
    print(f"  {n_reads:,} pairs, Illumina names, maps in {os.path.basename(rmap)} and "
          f"{os.path.basename(recmap)}", flush=True)

    # The seed is recorded because the off-target windows are derived from it, so a checker that is told
    # a different seed tests membership in windows this library was never built from. Acceptance was run
    # with --seed 5000 while chr13 was built with 5013, and expected three events the builder could not
    # have placed; the check passed anyway because it compared with >=.
    meta = {"dataset": a.dataset, "chrom": a.chrom, "library": a.library, "depth": a.depth,
            "seed": a.seed,
            "intervals": n_iv, "off_target_bands": bands, "quality_model": qual,
            "gc_bias": "measured" if gc.enabled else "flat",
            "capture_bases": sum(e - s + 1 for s, e in capture.get(a.chrom, [])), "sources": {f"{k[0]}_hap{k[1]}_{k[2]}_cn{k[3]}": w for k, (_p, w) in sources.items()},
            "coverage_per_source": [{"source": s, "hap": h, "coverage": c} for _r1, _r2, s, h, c in pieces],
            "r1": r1, "r2": r2, "read_map": rmap, "record_map": recmap,
            "pairs_written": n_reads, "runtime_s": round(time.time() - t0)}
    with open(os.path.join(a.out, f"{a.dataset}_{a.chrom}_{a.library}.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    for p in pieces:
        for f in p[:2]:
            if os.path.exists(f):
                os.remove(f)
    print(f"  wrote {r1} ({os.path.getsize(r1) / 1e6:.0f} MB)", flush=True)


if __name__ == "__main__":
    main()
