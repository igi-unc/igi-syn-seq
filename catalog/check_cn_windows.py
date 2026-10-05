#!/usr/bin/env python3
"""Assert that no WGS window straddles a copy-number boundary.

    python3 check_cn_windows.py --design design.yaml --paths paths.yaml [--chroms chr1,...]

`run_wgs.py` tiles a chromosome into windows and calls `source_plan` once per window, at the window's
midpoint. That makes one plan stand for a few megabases of reads, which is only sound if copy number is
constant across the window. With plain 5 Mb tiles it was not, and the consequence was not subtle: a 0.3 Mb
deletion has about a 6 % chance of containing the midpoint of the window it lands in, so focal events were
simply absent from the library. Acceptance measured PTEN_homdel at a depth ratio of 1.0 against an
expected 0.30, and RB1_homdel at 1.03 against 0.379.

The check is boundary containment, not sampled states. A three-point check on each window -- start,
midpoint, end -- reports zero straddled windows for exactly the case that is broken, because a focal event
sits in the interior and none of the three points touches it. That is how a first version of this check
passed the defect it was written to find.
"""
import argparse
import sys

import yaml

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from igi_catalog.designer import build_env      # noqa: E402
from run_wgs import whole_chromosome, WINDOW    # noqa: E402

ALL = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", default="design.yaml")
    ap.add_argument("--paths", required=True)
    ap.add_argument("--chroms", default=",".join(ALL))
    ap.add_argument("--window", type=int, default=WINDOW)
    a = ap.parse_args()
    paths = yaml.safe_load(open(a.paths))
    design = yaml.safe_load(open(a.design))

    n_bad = n_win = 0
    for ds in design["datasets"]:
        env = build_env(paths, design, ds)
        for chrom in a.chroms.split(","):
            if chrom not in env.genome.lengths:
                continue
            length = env.genome.lengths[chrom]
            cuts = [c for c in env.clones.cn_boundaries(chrom) if 1 < c <= length]
            wins = whole_chromosome(chrom, length, a.window, clones=env.clones)[chrom]
            n_win += len(wins)
            for s, e in wins:
                inside = [c for c in cuts if s < c <= e]
                if inside:
                    n_bad += 1
                    print(f"  FAIL {ds} {chrom} window {s:,}-{e:,} straddles "
                          f"{len(inside)} boundary/ies at {inside[:3]}")
    print(f"\n  {n_win:,} windows checked, {n_bad} straddle a copy-number boundary")
    return 1 if n_bad else 0


if __name__ == "__main__":
    sys.exit(main())
