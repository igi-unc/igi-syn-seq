"""Illumina-style read naming with a lookup table back to the source.

Simulators name reads after the sequence they came from, which would let a caller read the answer off a
read. Delivered FASTQs therefore carry ordinary Illumina names, and the truth bundle carries the map from
each name to the record it came from. Both are derived from the seed, so a rerun reproduces them.

Shuffling and renaming happen in one disk-based pass, so memory does not scale with library size.
"""
import gzip
import hashlib
import os
import subprocess

INSTRUMENT = "IGISYN"
FLOWCELL = "FC1"


def illumina_name(index, lane_count=4, run=1):
    """A well-formed Illumina name derived from a read's position in the shuffled library."""
    lane = index % lane_count + 1
    tile = 1101 + (index // lane_count) % 78
    x = 1000 + (index // (lane_count * 78)) % 20000
    y = 1000 + (index // (lane_count * 78 * 20000)) % 20000
    return f"{INSTRUMENT}:{run}:{FLOWCELL}:{lane}:{tile}:{x}:{y}"


def _key(name, seed):
    return hashlib.md5(f"{seed}:{name}".encode()).hexdigest()[:12]


def shuffle_and_rename(cat_r1, cat_r2, out_r1, out_r2, out_map, work, seed=1):
    """Shuffle a paired FASTQ and give it Illumina names, writing the name-to-source map.

    Records are keyed by a seeded digest and sorted on disk, so the peak memory is the sort buffer rather
    than the library.
    """
    keyed = os.path.join(work, "keyed.tsv")
    with open(cat_r1) as f1, open(cat_r2) as f2, open(keyed, "w") as out:
        while True:
            a = [f1.readline() for _ in range(4)]
            b = [f2.readline() for _ in range(4)]
            if not a[0]:
                break
            src = a[0][1:].split()[0].rsplit("-", 1)[0].rstrip("/1")
            out.write(_key(a[0].strip(), seed) + "\t" + src + "\t"
                      + "\t".join(x.rstrip("\n") for x in a[1:]) + "\t"
                      + "\t".join(x.rstrip("\n") for x in b[1:]) + "\n")
    sorted_path = os.path.join(work, "sorted.tsv")
    env = dict(os.environ, LC_ALL="C")
    subprocess.run(["sort", "-T", work, "-k1,1", "-o", sorted_path, keyed], check=True, env=env)
    n = 0
    with open(sorted_path) as fh, gzip.open(out_r1, "wt") as g1, gzip.open(out_r2, "wt") as g2, \
            gzip.open(out_map, "wt") as gm:
        gm.write("read_name\tsource_record\n")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            name = illumina_name(n)
            g1.write(f"@{name} 1:N:0:1\n{f[2]}\n+\n{f[4]}\n")
            g2.write(f"@{name} 2:N:0:1\n{f[5]}\n+\n{f[7]}\n")
            gm.write(f"{name}\t{f[1]}\n")
            n += 1
    for p in (keyed, sorted_path):
        if os.path.exists(p):
            os.remove(p)
    return n
