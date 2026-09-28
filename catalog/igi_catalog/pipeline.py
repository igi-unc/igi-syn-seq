"""Stages 2-5 driver: clone genomes -> reads -> a packaged release.

WES is built capture interval by capture interval. Each interval is rebuilt on every clone and haplotype
that retains it, which keeps all edits in reference coordinates and avoids having to map positions through
indels. Reads are drawn from each interval in proportion to the clone's share of sequenced molecules and
the haplotype's copy number there, so the realised allele fractions follow the truth table by construction.
"""
import os
import subprocess
from collections import defaultdict

from .genome_build import haplotype_sequence, read_events
from .simulate import clone_weights


def merged_capture(bed_path, chroms, pad=100, gap=200):
    """Capture intervals for the requested chromosomes, padded and merged."""
    by_chrom = defaultdict(list)
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            if f[0] not in chroms:
                continue
            by_chrom[f[0]].append((max(1, int(f[1]) + 1 - pad), int(f[2]) + pad))
    out = {}
    for c, ivs in by_chrom.items():
        ivs.sort()
        merged = [list(ivs[0])]
        for s, e in ivs[1:]:
            if s - merged[-1][1] <= gap:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out[c] = [tuple(m) for m in merged]
    return out


def source_plan(clones, purity, chrom, pos, tumor=True):
    """(source, haplotype, weight) triples for one locus.

    For the tumour library the weight of a clone is split across its haplotypes by copy number, so an
    amplified haplotype contributes proportionally more reads. For the normal library only the germline
    haplotypes contribute, equally.
    """
    if not tumor:
        return [("NORMAL", 0, 0.5), ("NORMAL", 1, 0.5)]
    w = clone_weights(clones, purity)
    plan = []
    for src, weight in w.items():
        if weight <= 0:
            continue
        if src == "NORMAL":
            plan += [("NORMAL", 0, weight * 0.5), ("NORMAL", 1, weight * 0.5)]
            continue
        a, b, _ = clones.cn(src, chrom, pos)
        total = a + b
        if total == 0:
            continue
        for hap, cn in enumerate((a, b)):
            if cn:
                plan.append((src, hap, weight * cn / total))
    return plan


class WesBuilder:
    """Write per-source interval FASTAs for an exome library, then simulate reads from them."""

    def __init__(self, env, events_by_chrom, purity, workdir):
        self.env = env
        self.events = events_by_chrom
        self.purity = purity
        self.workdir = workdir
        os.makedirs(workdir, exist_ok=True)

    def write_source_fastas(self, capture, tumor=True, max_intervals=None):
        """One FASTA per (source, haplotype); each record is one capture interval on that haplotype.

        Returns {(source, hap): (path, weight)} with weights averaged over intervals.
        """
        handles = {}
        weight_sum = defaultdict(float)
        n = 0
        for chrom, ivs in sorted(capture.items()):
            evs = self.events.get(chrom, [])
            for start1, end1 in ivs:
                mid = (start1 + end1) // 2
                for src, hap, w in source_plan(self.env.clones, self.purity, chrom, mid, tumor=tumor):
                    key = (src, hap)
                    if key not in handles:
                        path = os.path.join(self.workdir, f"{'tumor' if tumor else 'normal'}_{src}_hap{hap}.fa")
                        handles[key] = (path, open(path, "w"))
                    use = [] if src == "NORMAL" else evs
                    clone = "T" if src == "NORMAL" else src
                    seq, _stats = haplotype_sequence(self.env.genome, self.env.germline, use,
                                                     self.env.clones, chrom, hap, clone, start1, end1,
                                                     strict=False)
                    if len(seq) < 100:
                        continue
                    fh = handles[key][1]
                    fh.write(f">{chrom}_{start1}_{end1}_{src}_hap{hap}\n")
                    for i in range(0, len(seq), 60):
                        fh.write(seq[i:i + 60] + "\n")
                    weight_sum[key] += w
                n += 1
                if max_intervals and n >= max_intervals:
                    break
            if max_intervals and n >= max_intervals:
                break
        out = {}
        total = sum(weight_sum.values()) or 1.0
        for key, (path, fh) in handles.items():
            fh.close()
            out[key] = (path, weight_sum[key] / total)
        return out, n

    def simulate(self, sources, art_cmd, depth, out_prefix, read_len=150, seed=1):
        """Run ART per source at coverage proportional to its weight; returns the FASTQ pieces."""
        pieces = []
        for i, ((src, hap), (path, weight)) in enumerate(sorted(sources.items())):
            # weights sum to 1 across sources, so a source is simulated at its share of the total depth
            cov = max(1, round(depth * weight))
            pre = os.path.join(self.workdir, f"{os.path.basename(out_prefix)}_{src}_h{hap}_")
            cmd = art_cmd.format(args=(f"-ss HS25 -i {path} -p -l {read_len} -f {cov} -m 350 -s 60 "
                                       f"-rs {seed + i} -na -o {pre}"))
            subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
            pieces.append((f"{pre}1.fq", f"{pre}2.fq", src, hap, cov))
        return pieces


def concat_fastqs(pieces, out_r1, out_r2, gzip_output=True):
    """Concatenate per-source FASTQs into one library, optionally bgzip-compressed."""
    for idx, out in ((0, out_r1), (1, out_r2)):
        parts = [p[idx] for p in pieces if os.path.exists(p[idx])]
        if not parts:
            continue
        if gzip_output:
            with open(out, "wb") as fh:
                cat = subprocess.Popen(["cat"] + parts, stdout=subprocess.PIPE)
                gz = subprocess.Popen(["gzip", "-c"], stdin=cat.stdout, stdout=fh)
                cat.stdout.close()
                if gz.wait() or cat.wait():
                    raise RuntimeError(f"failed to write {out}")
        else:
            with open(out, "w") as fh:
                subprocess.run(["cat"] + parts, stdout=fh, check=True)
    return out_r1, out_r2
