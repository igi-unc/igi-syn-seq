#!/usr/bin/env python3
"""Generate the passenger background mutation table for one dataset.

    python3 run_background.py --design design.yaml --paths paths.yaml --dataset IGI-SYN-SEQ-01 \
        --catalog output/IGI-SYN-SEQ-01.snv_indel.tsv --out output/ [--chroms chr1,chr6]

Writes <dataset>.background.tsv next to the designed catalog. The designed catalog is read only to reserve
the spans it occupies, so a passenger never collides with an event under test.
"""
import argparse, json, os, random, sys, time
from collections import Counter

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from igi_catalog.designer import build_env, write_table
from igi_catalog.background import BackgroundDesigner, channel

RESOURCE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources",
                        "cosmic_v3.4_sbs_subset.tsv")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", required=True)
    ap.add_argument("--paths", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--catalog", nargs="*", default=[], help="designed tables whose spans to reserve")
    ap.add_argument("--out", required=True)
    ap.add_argument("--chroms", default=None, help="comma-separated; default autosomes plus sex chromosome(s)")
    ap.add_argument("--signatures", default=RESOURCE)
    ap.add_argument("--cache", default=None,
                    help="netMHCpan cache file; defaults to a background-specific one so a "
                         "concurrent catalog run cannot clobber it")
    a = ap.parse_args(argv)

    design = yaml.safe_load(open(a.design))
    paths = yaml.safe_load(open(a.paths))
    dcfg = design["datasets"][a.dataset]
    rng = random.Random(f"{design['seed']}:{a.dataset}:background")
    env = build_env(paths, design, a.dataset)
    # the designer may be running against the shared cache; a whole-file rewrite from either
    # side would discard the other's entries, so background scoring keeps its own
    env.netmhc.cache_path = a.cache or os.path.join(env.work, "netmhcpan_cache.background.json")
    env.netmhc.cache = {}
    if os.path.exists(env.netmhc.cache_path):
        env.netmhc.cache = json.load(open(env.netmhc.cache_path))

    if a.chroms:
        chroms = a.chroms.split(",")
    else:
        chroms = [f"chr{i}" for i in range(1, 23)] + (["chrX"] if dcfg["sex"] == "female" else ["chrX", "chrY"])
    chroms = [c for c in chroms if c in env.genome.lengths]

    bg = BackgroundDesigner(env, a.dataset, dcfg, design, rng, a.signatures)
    n_res = bg.reserve_designed(a.catalog)
    print(f"[background] {a.dataset}: {dcfg['background_mut_per_mb']} mut/Mb, signatures "
          f"{dcfg['signatures']}, {n_res} reserved positions from {len(a.catalog)} designed table(s)", flush=True)
    t0 = time.time()
    log = lambda m: print(m, flush=True)
    events = bg.design(chroms, list(dcfg["clones"]),
                       indel_fraction=float(design.get("background_indel_fraction", 0.08)),
                       truncal_fraction=float(design.get("background_truncal_fraction", 0.60)),
                       log=log)
    env.netmhc.flush()
    os.makedirs(a.out, exist_ok=True)
    path = write_table(events, os.path.join(a.out, f"{a.dataset}.background.tsv"))
    summ = {
        "dataset": a.dataset, "n": len(events),
        "runtime_s": round(time.time() - t0),
        "mut_per_mb": dcfg["background_mut_per_mb"],
        "signatures": dcfg["signatures"],
        "by_class": dict(Counter(e["class"] for e in events)),
        "by_clone": dict(Counter(e["clone"] for e in events)),
        "by_clonality_tier": dict(Counter(e["clonality_tier"] for e in events)),
        "by_consequence": dict(Counter(e["consequence"] for e in events)),
        "by_binding_tier": dict(Counter(e["binding_tier"] for e in events)),
        "coding": sum(1 for e in events if e["consequence"] != "non_coding"),
        "chr1to6": sum(1 for e in events if e["chr1to6"]),
        "spectrum": dict(Counter(channel(e["trinuc"], e["ref"], e["alt"])
                                 for e in events if e["class"] == "snv" and e["trinuc"])),
    }
    with open(os.path.join(a.out, f"{a.dataset}.background.summary.json"), "w") as fh:
        json.dump(summ, fh, indent=2)
    print(json.dumps({k: v for k, v in summ.items() if k != "spectrum"}, indent=2))
    print(f"[done] {len(events)} background mutations -> {path}")


if __name__ == "__main__":
    main()
