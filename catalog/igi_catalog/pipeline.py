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
    """(source, haplotype, weight) triples for one locus, weights summing to 1.

    Read depth follows the number of copies present, not the number of cells: a clone contributes in
    proportion to its exclusive cell fraction times the absolute copy number of that haplotype at this
    locus. An amplified haplotype therefore contributes proportionally more reads and a lost one none,
    and the resulting allele fractions match the expected VAF computed by the clone model, which divides
    by the same mean local copy number.
    """
    if not tumor:
        return [("NORMAL", 0, "all", 0.5), ("NORMAL", 1, "all", 0.5)]
    raw = []
    for src, frac in clones.excl.items():
        if frac <= 0:
            continue
        a, b, _ = clones.cn(src, chrom, pos)
        for hap, cn in enumerate((a, b)):
            if not cn:
                continue
            # One copy of the haplotype carries everything this clone has: truncal variants from before
            # and after the copy-number changes, plus the clone's own. The remaining copies existed
            # before those later events, so they carry only the pre-CNA truncal variants. This is what
            # keeps multiplicity right under whole-genome doubling, where a haplotype has several copies.
            raw.append((src, hap, "all", purity * frac * 1))
            if cn > 1:
                raw.append((src, hap, "pre_cna", purity * frac * (cn - 1)))
    raw += [("NORMAL", 0, "all", (1.0 - purity)), ("NORMAL", 1, "all", (1.0 - purity))]
    total = sum(w for _s, _h, _k, w in raw) or 1.0
    return [(s, h, k, w / total) for s, h, k, w in raw]


def cn_profile(clones, chrom, pos):
    """A hashable copy-number state at a locus, used to group intervals that share a coverage plan."""
    return tuple(sorted((c, clones.cn(c, chrom, pos)[0], clones.cn(c, chrom, pos)[1]) for c in clones.clones))


class WesBuilder:
    """Write per-source interval FASTAs for an exome library, then simulate reads from them."""

    def __init__(self, env, events_by_chrom, purity, workdir):
        self.env = env
        self.events = events_by_chrom
        self.purity = purity
        self.workdir = workdir
        self._n_rec = 0
        self.record_map = []       # (record id, chrom, start, end, source, haplotype, copy kind)
        self.rejected = []         # designed somatic edits that could not be applied
        self.germline_rejected = 0
        os.makedirs(workdir, exist_ok=True)

    def write_source_fastas(self, capture, tumor=True, max_intervals=None):
        """One FASTA per (source, haplotype); each record is one capture interval on that haplotype.

        Returns {(source, hap): (path, weight)} with weights averaged over intervals.
        """
        handles = {}
        weights = {}
        n = 0
        for chrom, ivs in sorted(capture.items()):
            evs = self.events.get(chrom, [])
            for start1, end1 in ivs:
                mid = (start1 + end1) // 2
                prof = cn_profile(self.env.clones, chrom, mid) if tumor else ("normal",)
                pid = self._profile_id(prof)
                for src, hap, kind, w in source_plan(self.env.clones, self.purity, chrom, mid, tumor=tumor):
                    key = (src, hap, kind, pid)
                    if key not in handles:
                        path = os.path.join(self.workdir,
                                            f"{'tumor' if tumor else 'normal'}_{src}_hap{hap}_{kind}_cn{pid}.fa")
                        handles[key] = (path, open(path, "w"))
                        weights[key] = w
                    use = [] if src == "NORMAL" else evs
                    # a pre-CNA copy carries only truncal events that predate the copy-number changes
                    clone = "T" if (src == "NORMAL" or kind == "pre_cna") else src
                    pre_only = (kind == "pre_cna")
                    seq, stats = haplotype_sequence(self.env.genome, self.env.germline, use,
                                                    self.env.clones, chrom, hap, clone, start1, end1,
                                                    strict=False, only_pre_cna=pre_only)
                    # a designed somatic edit that fails to apply means the reads will not contain a
                    # change the truth table claims, so it is collected and raised rather than dropped
                    if stats.get("somatic_rejected"):
                        self.rejected.extend(stats["somatic_rejected"])
                    self.germline_rejected += stats.get("germline_rejected", 0)
                    if len(seq) < 100:
                        continue
                    fh = handles[key][1]
                    # A stable id, unique across every source file. The previous name omitted `kind`
                    # and the copy-number profile, so the "all" and "pre_cna" files produced identical
                    # ART read ids and 860k of 2.0M names collided.
                    self._n_rec += 1
                    fh.write(f">e{self._n_rec:09d}\n")
                    self.record_map.append((f"e{self._n_rec:09d}", chrom, start1, end1, src, hap, kind))
                    for i in range(0, len(seq), 60):
                        fh.write(seq[i:i + 60] + "\n")
                n += 1
                if max_intervals and n >= max_intervals:
                    break
            if max_intervals and n >= max_intervals:
                break
        out = {}
        for key, (path, fh) in handles.items():
            fh.close()
            out[key] = (path, weights[key])      # per-interval weight for this copy-number profile
        if self.rejected:
            uniq = sorted(set(self.rejected))
            raise RuntimeError(
                f"{len(uniq)} designed somatic edits could not be applied, so the reads would not "
                f"contain changes the truth table describes: {uniq[:8]}")
        return out, n

    def _profile_id(self, prof):
        if not hasattr(self, "_profiles"):
            self._profiles = {}
        if prof not in self._profiles:
            self._profiles[prof] = len(self._profiles)
        return self._profiles[prof]

    def simulate(self, sources, art_cmd, depth, out_prefix, read_len=150, seed=1):
        """Run ART per source at coverage proportional to its weight; returns the FASTQ pieces."""
        pieces = []
        for i, (key, (path, weight)) in enumerate(sorted(sources.items())):
            src, hap, kind, pid = key
            # weights sum to 1 within a copy-number profile, so every locus receives `depth` in total.
            # ART accepts fractional fold-coverage, and rounding to an integer distorts the smallest
            # clones most: a 2.33x target became 2x, which is a 14% shortfall on the deepest subclone.
            cov = max(0.01, depth * weight)
            pre = os.path.join(self.workdir, f"{os.path.basename(out_prefix)}_{src}_h{hap}_{kind}_cn{pid}_")
            cmd = art_cmd.format(args=(f"-ss HS25 -i {path} -p -l {read_len} -f {cov:.5f} -m 350 -s 60 "
                                       f"-rs {seed + i} -na -o {pre}"))
            subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
            pieces.append((f"{pre}1.fq", f"{pre}2.fq", src, hap, round(cov, 5)))
        return pieces


def write_record_map(record_map, path):
    """The map from each source record to where it came from, for the truth bundle."""
    import gzip
    with gzip.open(path, "wt") as fh:
        fh.write("record\tchrom\tstart\tend\tsource\thaplotype\tcopy_kind\n")
        for row in record_map:
            fh.write("\t".join(str(x) for x in row) + "\n")
    return path


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
