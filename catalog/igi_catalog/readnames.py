"""Illumina-style read naming with a lookup table back to the source.

Simulators name reads after the sequence they came from, which would let a caller read the answer off a
read. Delivered FASTQs therefore carry ordinary Illumina names, and the truth bundle carries the map from
each name to the record it came from. Both are derived from the seed, so a rerun reproduces them.

Shuffling and renaming happen in one disk-based pass, so memory does not scale with library size
-- but only because the sort buffer is bounded explicitly; see `sort_buffer_arg`.
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


# GNU sort sizes its default buffer from the machine's PHYSICAL memory and never reads the cgroup limit a
# SLURM step runs under, so a bare `sort` on a 503 GB node inside a 16 GB allocation is killed by the OOM
# killer no matter how the sort is otherwise set up. The earlier 320 GB request did not bound anything; it
# only happened to be larger than the buffer sort chose. Size the buffer from the allocation instead, so
# the step's footprint is a property of the job and not of whichever node it lands on.
SORT_BUFFER_FRACTION = 0.25
SORT_BUFFER_MIN_MB = 256
SORT_BUFFER_MAX_MB = 8192


def sort_buffer_mb(env=None):
    """Megabytes to allow GNU sort, derived from this step's SLURM memory allocation.

    A quarter of the allocation leaves room for the Python process, the gzip writers and sort's own
    per-thread overhead, which sit alongside the buffer. Outside SLURM the floor applies, which is slower
    but cannot be killed.
    """
    env = os.environ if env is None else env
    total = None
    per_node = env.get("SLURM_MEM_PER_NODE")
    if per_node and per_node.isdigit():
        total = int(per_node)
    else:
        per_cpu = env.get("SLURM_MEM_PER_CPU")
        cpus = env.get("SLURM_CPUS_ON_NODE") or env.get("SLURM_CPUS_PER_TASK")
        if per_cpu and per_cpu.isdigit() and cpus and cpus.isdigit():
            total = int(per_cpu) * int(cpus)
    if not total:
        return SORT_BUFFER_MIN_MB
    return max(SORT_BUFFER_MIN_MB, min(SORT_BUFFER_MAX_MB, int(total * SORT_BUFFER_FRACTION)))


def sort_buffer_arg(env=None):
    """`-S` argument for GNU sort, as a list ready to splice into a command."""
    return ["-S", f"{sort_buffer_mb(env)}M"]


# Names are a pure function of a read's index, so two runs starting from 0 produce the identical name
# sequence. A library built a chromosome at a time therefore had 100 % name collision between every pair
# of chromosomes, and merging them would have repeated each name 24 times and made the name-to-record map
# ambiguous -- destroying exactly the truth artefact owner decision D2 asks for. Each run is given a slice
# of the name space instead. The space is 4 lanes x 78 tiles x 20000 x 20000 = 1.25e11 names, so a stride
# of 4e9 holds 31 runs without overlap.
NAME_SPACE_STRIDE = 4_000_000_000


def shuffle_and_rename(cat_r1, cat_r2, out_r1, out_r2, out_map, work, seed=1, index_offset=0):
    """Shuffle a paired FASTQ and give it Illumina names, writing the name-to-source map.

    Records are keyed by a seeded digest and sorted on disk, so the peak memory is the sort buffer rather
    than the library, and that buffer is capped against the step's allocation rather than left to GNU
    sort's physical-memory heuristic. `index_offset` places this run's names in their own slice of the name space, so
    libraries built in pieces can be concatenated without collision.
    """
    keyed = os.path.join(work, "keyed.tsv")
    with open(cat_r1) as f1, open(cat_r2) as f2, open(keyed, "w") as out:
        while True:
            a = [f1.readline() for _ in range(4)]
            b = [f2.readline() for _ in range(4)]
            if not a[0]:
                break
            # `rstrip("/1")` would strip a *set* of characters, not a suffix, so every record id
            # ending in 1 lost its trailing 1s and its reads could not be traced back to a source.
            name = a[0][1:].split()[0]
            if name[-2:] in ("/1", "/2"):
                name = name[:-2]
            src = name.rsplit("-", 1)[0]
            out.write(_key(a[0].strip(), seed) + "\t" + src + "\t"
                      + "\t".join(x.rstrip("\n") for x in a[1:]) + "\t"
                      + "\t".join(x.rstrip("\n") for x in b[1:]) + "\n")
    sorted_path = os.path.join(work, "sorted.tsv")
    env = dict(os.environ, LC_ALL="C")
    subprocess.run(["sort", "-T", work, "-k1,1"] + sort_buffer_arg()
                   + ["-o", sorted_path, keyed], check=True, env=env)
    n = 0
    with open(sorted_path) as fh, gzip.open(out_r1, "wt") as g1, gzip.open(out_r2, "wt") as g2, \
            gzip.open(out_map, "wt") as gm:
        gm.write("read_name\tsource_record\n")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            name = illumina_name(index_offset + n)
            g1.write(f"@{name} 1:N:0:1\n{f[2]}\n+\n{f[4]}\n")
            g2.write(f"@{name} 2:N:0:1\n{f[5]}\n+\n{f[7]}\n")
            gm.write(f"{name}\t{f[1]}\n")
            n += 1
    for p in (keyed, sorted_path):
        if os.path.exists(p):
            os.remove(p)
    return n
