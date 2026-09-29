"""MHC class I binding via netMHCpan 4.1.

This build scores roughly 40 peptide-allele pairs per second per process and starts a fresh process for
each allele and length, so the designer keeps the work small and parallel:

* only peptides that span a changed residue are submitted (peptide-list input, `-p`), not whole windows;
* the requested lengths are all scored (8-11mers by default, the full class I range);
* every (allele, length) combination runs as its own concurrent process, since netMHCpan forks one
  process per combination anyway;
* every peptide-allele result is cached on disk, so re-runs and overlapping windows are free.
"""
import hashlib
import json
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor


def spanning_peptides(protein, first_changed, n_changed=1, lengths=(8, 9, 10, 11)):
    """Peptides of the given lengths that overlap [first_changed, first_changed + n_changed)."""
    out = set()
    if first_changed is None or not protein:
        return out
    last = first_changed + max(1, n_changed) - 1
    for L in lengths:
        # a peptide [start, start+L-1] overlaps [first_changed, last] when it starts at or before `last`
        # and ends at or after `first_changed`; requiring it to *contain* the whole changed region
        # instead yields nothing as soon as the region is longer than the peptide
        lo = max(0, first_changed - L + 1)
        hi = min(last, len(protein) - L)
        for start in range(lo, hi + 1):
            pep = protein[start:start + L]
            if len(pep) == L and "*" not in pep and "X" not in pep:
                out.add(pep)
    return out


class NetMHCpan:
    """Peptide-level binding predictions with an on-disk cache keyed by (allele, peptide)."""

    def __init__(self, cmd_template, workdir, cache_path=None, threads=6):
        self.cmd = cmd_template
        self.workdir = workdir
        os.makedirs(workdir, exist_ok=True)
        self.cache_path = cache_path
        self.threads = threads
        self.cache = {}
        if cache_path and os.path.exists(cache_path):
            try:
                with open(cache_path) as fh:
                    self.cache = json.load(fh)
            except (ValueError, OSError):
                self.cache = {}
        self._dirty = False

    @staticmethod
    def _fmt(allele):
        return allele.replace("*", "")

    @staticmethod
    def _key(allele, peptide):
        return f"{allele}|{peptide}"

    def predict(self, peptides, alleles):
        """Score every (peptide, allele) pair not already cached.

        Returns {peptide: {allele: (rank_el, rank_ba, aff_nM)}}.
        """
        peptides = sorted({p for p in peptides if p})
        # one job per (allele, length): netMHCpan forks a process per combination anyway, so running them
        # concurrently costs nothing extra and keeps wall time flat as more lengths are requested
        jobs = []
        for a in alleles:
            by_len = {}
            for p in peptides:
                if self._key(a, p) not in self.cache:
                    by_len.setdefault(len(p), []).append(p)
            jobs += [(a, L, ps) for L, ps in sorted(by_len.items())]
        if jobs:
            with ThreadPoolExecutor(max_workers=min(self.threads, len(jobs))) as ex:
                for (allele, _L, submitted), rows in zip(jobs, ex.map(lambda j: self._run_one(*j), jobs)):
                    for pep, vals in rows.items():
                        self.cache[self._key(allele, pep)] = vals
                    for pep in submitted:
                        self.cache.setdefault(self._key(allele, pep), None)
            self._dirty = True
            self.flush()
        out = {}
        for p in peptides:
            per = {}
            for a in alleles:
                v = self.cache.get(self._key(a, p))
                if v:
                    per[a] = tuple(v)
            if per:
                out[p] = per
        return out

    def _run_one(self, allele, length, peps):
        """One netMHCpan process for one allele and one peptide length."""
        rows = {}
        tag = hashlib.md5((allele + str(length) + "".join(peps)).encode()).hexdigest()[:16]
        path = os.path.join(self.workdir, f"pep_{tag}.txt")
        with open(path, "w") as fh:
            fh.write("\n".join(peps) + "\n")
        cmd = self.cmd.format(fasta=path, alleles=self._fmt(allele)) + f" -p -l {length}"
        try:
            out = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(f"netMHCpan failed for {allele} length {length}: {exc.stderr[-400:]}") from exc
        finally:
            if os.path.exists(path):
                os.remove(path)
        for line in out.splitlines():
            f = line.split()
            if len(f) < 16 or not f[0].isdigit() or not f[1].startswith("HLA-"):
                continue
            try:
                rows[f[2]] = (float(f[12]), float(f[14]), float(f[15]))
            except ValueError:
                continue
        return rows

    def flush(self):
        if self.cache_path and self._dirty:
            tmp = self.cache_path + ".tmp"
            with open(tmp, "w") as fh:
                json.dump(self.cache, fh)
            os.replace(tmp, self.cache_path)
            self._dirty = False

    @staticmethod
    def best(pep_alleles, among=None):
        """Lowest %rank_EL across alleles -> (rank, allele, affinity).

        `among` restricts the search to a subset of alleles, which is how the retained-allele tier is
        scored: a tumour that has lost an HLA haplotype cannot present through the alleles on it, so the
        best rank over the surviving alleles is a different and equally real number (owner decision D4).
        Returns (None, None, None) when no allele in `among` was scored for this peptide.
        """
        items = pep_alleles if among is None else {a: v for a, v in pep_alleles.items() if a in among}
        if not items:
            return None, None, None
        allele, vals = min(items.items(), key=lambda kv: kv[1][0])
        return vals[0], allele, vals[2]

    @staticmethod
    def best_over(res, peptides, among=None):
        """Best (rank, allele, affinity, peptide) over a set of peptides, optionally restricted by allele."""
        best = None
        for pep in peptides:
            al = res.get(pep)
            if not al:
                continue
            r, allele, aff = NetMHCpan.best(al, among=among)
            if r is None:
                continue
            if best is None or r < best[0]:
                best = (r, allele, aff, pep)
        return best or (None, None, None, None)

    @staticmethod
    def tier(rank_el):
        if rank_el is None:
            return "na"
        if rank_el <= 0.5:
            return "strong"
        if rank_el <= 2.0:
            return "weak"
        return "non"
