"""The per-cell model: one roster of cells shared by every single-cell assay.

Four libraries in the design see the same suspension: 10x 5' gene expression, 10x TCR, Kinnex
single-cell, and (if built) ONT single-cell. They must agree cell by cell. A barcode that is a tumour cell
of clone A1 in the gene-expression library cannot be a CD8 T cell in the TCR library, and a clonotype
observed in TCR has to belong to a barcode that the expression library also calls a T cell. That only
holds if one roster is built once and every assay reads it, which is why this is its own module rather
than something each builder decides for itself.

The roster is deterministic from the design seed, so a rerun reproduces it and the assays stay consistent
across separate jobs.
"""
import bisect
import csv
import random

# Cell types and their share of the suspension, from design section 8.
DEFAULT_COMPOSITION = [
    ("tumor", 0.35), ("CD8_T", 0.18), ("CD4_T", 0.09), ("Treg", 0.03), ("B", 0.05),
    ("plasma", 0.03), ("macrophage", 0.12), ("DC", 0.02), ("CAF", 0.08),
    ("endothelial", 0.04), ("NK", 0.01),
]
# Tumour cells are split by clone using the clone model's exclusive fractions, never a list written here.
# Design section 8 gives "T-only 35 %, A 40 %, A1 12 % of A, B 25 %", which sums to 1.12 and is ambiguous
# about whether A's share includes A1. The clone model resolves it: for dataset 01 the exclusive fractions
# are T 0.35, A 0.28, A1 0.12, B 0.25, and 0.28 + 0.12 = 0.40, so the prose's "A 40 %" is A inclusive of
# its subclone. Deriving from `excl` makes per-cell genotypes reproduce the CCFs by construction, which is
# what section 8 asks for, and it is right per dataset: dataset 02's A is 0.25, not 0.28. Hardcoding the
# numbers here made every clone low by a factor of 0.893 and was the same defect as the PIK3CA driver bug,
# two declarations of one quantity with only one of them wired up.
T_CELL_TYPES = ("CD8_T", "CD4_T", "Treg")
# Rough per-cell transcript yield by type: a plasma cell is a secretory factory, a resting T cell is not.
UMI_SCALE = {"tumor": 1.0, "CD8_T": 0.45, "CD4_T": 0.45, "Treg": 0.40, "B": 0.55,
             "plasma": 1.6, "macrophage": 1.1, "DC": 0.8, "CAF": 1.2,
             "endothelial": 0.9, "NK": 0.5}


def load_whitelist(path, limit=None):
    """10x barcode whitelist, optionally truncated for a smoke test."""
    out = []
    opener = open
    if path.endswith(".gz"):
        import gzip
        opener = gzip.open
    with opener(path, "rt") as fh:
        for line in fh:
            b = line.strip()
            if b:
                out.append(b)
            if limit and len(out) >= limit:
                break
    return out


def power_law_sizes(n_cells, n_clonotypes, rng, top_share=0.30, top_n=10):
    """Clonotype sizes following a power law, calibrated so the top `top_n` hold `top_share` of cells.

    The design asks for ~300 clonotypes over ~1,200 T cells with the top ten holding about 30 %. Rather
    than pick an exponent and hope, the exponent is solved for numerically against that constraint, so the
    realised share matches what was specified.
    """
    lo, hi = 0.1, 4.0
    for _ in range(60):
        mid = (lo + hi) / 2
        w = [(i + 1) ** -mid for i in range(n_clonotypes)]
        tot = sum(w)
        share = sum(w[:top_n]) / tot
        if share < top_share:
            lo = mid
        else:
            hi = mid
    w = [(i + 1) ** -((lo + hi) / 2) for i in range(n_clonotypes)]
    tot = sum(w)
    sizes = [max(1, int(round(n_cells * x / tot))) for x in w]
    # correct rounding drift so the sizes sum to n_cells exactly
    i = 0
    while sum(sizes) > n_cells and any(s > 1 for s in sizes):
        if sizes[i % n_clonotypes] > 1:
            sizes[i % n_clonotypes] -= 1
        i += 1
    while sum(sizes) < n_cells:
        sizes[i % n_clonotypes] += 1
        i += 1
    return sizes


class CellRoster:
    """Cells, their barcodes, types, clones and clonotypes; written once and shared by every assay."""

    def __init__(self, n_cells, whitelist, seed, clones, composition=None,
                 doublet_rate=0.06, ambient_rate=0.03, n_clonotypes=300, tcr_cells=1200,
                 alpha_dual_rate=0.10):
        self.rng = random.Random(f"{seed}:cells")
        self.n_cells = n_cells
        self.doublet_rate = doublet_rate
        self.ambient_rate = ambient_rate
        comp = composition or DEFAULT_COMPOSITION
        tumor_clones = sorted(clones.excl.items())
        tot = sum(v for _k, v in tumor_clones)
        if abs(tot - 1.0) > 1e-6:
            raise ValueError(f"clone exclusive fractions sum to {tot}, not 1; the clone tree is wrong")
        if len(whitelist) < n_cells * 2:
            raise ValueError(f"whitelist has {len(whitelist)} barcodes, too few for {n_cells} cells")
        bcs = self.rng.sample(whitelist, n_cells)

        # cell types, allocated by count rather than by independent draws so the realised composition
        # matches the design instead of merely having the right expectation
        counts = self._allocate(n_cells, comp)
        types = []
        for t, k in counts.items():
            types += [t] * k
        self.rng.shuffle(types)

        tumor_idx = [i for i, t in enumerate(types) if t == "tumor"]
        clone_counts = self._allocate(len(tumor_idx), tumor_clones)
        clone_seq = []
        for c, k in clone_counts.items():
            clone_seq += [c] * k
        self.rng.shuffle(clone_seq)

        self.cells = []
        ci = 0
        for i, (bc, t) in enumerate(zip(bcs, types)):
            clone = ""
            if t == "tumor":
                clone = clone_seq[ci]
                ci += 1
            self.cells.append({"barcode": bc, "cell_type": t, "clone": clone,
                               "umi_scale": UMI_SCALE.get(t, 1.0),
                               "clonotype": "", "trb": "", "tra": "", "tra2": "",
                               "is_doublet": False, "partner_barcode": ""})

        self._assign_clonotypes(n_clonotypes, tcr_cells, alpha_dual_rate)
        self._assign_doublets()

    @staticmethod
    def _allocate(n, shares):
        """Integer counts matching the requested shares, with the remainder given to the largest groups."""
        total = sum(s for _k, s in shares)
        raw = [(k, n * s / total) for k, s in shares]
        out = {k: int(v) for k, v in raw}
        rem = n - sum(out.values())
        for k, v in sorted(raw, key=lambda kv: -(kv[1] - int(kv[1])))[:rem]:
            out[k] += 1
        return out

    def _assign_clonotypes(self, n_clonotypes, tcr_cells, alpha_dual_rate):
        t_idx = [i for i, c in enumerate(self.cells) if c["cell_type"] in T_CELL_TYPES]
        self.rng.shuffle(t_idx)
        take = t_idx[:min(tcr_cells, len(t_idx))]
        n_clonotypes = max(1, min(n_clonotypes, len(take)))
        sizes = power_law_sizes(len(take), n_clonotypes, self.rng)
        self.clonotypes = []
        pos = 0
        for ci, size in enumerate(sizes, start=1):
            name = f"clonotype{ci}"
            trb = self._cdr3(self.rng, "TRB")
            tra = self._cdr3(self.rng, "TRA")
            members = take[pos:pos + size]
            pos += size
            for i in members:
                c = self.cells[i]
                c["clonotype"], c["trb"], c["tra"] = name, trb, tra
                if self.rng.random() < alpha_dual_rate:
                    c["tra2"] = self._cdr3(self.rng, "TRA")
            self.clonotypes.append({"clonotype": name, "size": len(members), "trb": trb, "tra": tra})

    @staticmethod
    def _cdr3(rng, chain):
        """A CDR3 amino-acid sequence of plausible length and composition.

        These are synthetic sequences, not real V(D)J recombinants: the design uses clonotypes to test
        clustering and expansion, not germline gene assignment, and inventing a real-looking V-J pairing
        would imply a fidelity this does not have. Documented rather than hidden.
        """
        n = rng.randint(10, 18)
        mid = "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(n - 3))
        return ("CAS" if chain == "TRB" else "CAV") + mid + "F"

    def _assign_doublets(self):
        n = int(round(self.n_cells * self.doublet_rate))
        idx = self.rng.sample(range(len(self.cells)), min(n * 2, len(self.cells)))
        for a, b in zip(idx[::2], idx[1::2]):
            self.cells[a]["is_doublet"] = True
            self.cells[a]["partner_barcode"] = self.cells[b]["barcode"]

    # ---------------------------------------------------------------- output
    def write(self, path):
        cols = ["barcode", "cell_type", "clone", "umi_scale", "clonotype", "trb", "tra", "tra2",
                "is_doublet", "partner_barcode"]
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(self.cells)
        return path

    def write_clonotypes(self, path):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["clonotype", "size", "trb", "tra"], delimiter="\t")
            w.writeheader()
            w.writerows(sorted(self.clonotypes, key=lambda c: -c["size"]))
        return path

    def summary(self):
        from collections import Counter
        t = Counter(c["cell_type"] for c in self.cells)
        cl = Counter(c["clone"] for c in self.cells if c["clone"])
        sizes = sorted((c["size"] for c in self.clonotypes), reverse=True)
        n_t = sum(1 for c in self.cells if c["clonotype"])
        return {
            "n_cells": len(self.cells),
            "by_cell_type": dict(t),
            "cell_type_fraction": {k: round(v / len(self.cells), 4) for k, v in t.items()},
            "tumor_by_clone": dict(cl),
            "tumor_clone_fraction": {k: round(v / max(1, sum(cl.values())), 4) for k, v in cl.items()},
            "n_clonotypes": len(self.clonotypes),
            "t_cells_with_clonotype": n_t,
            "top10_clonotype_share": round(sum(sizes[:10]) / max(1, n_t), 4),
            "largest_clonotype": sizes[0] if sizes else 0,
            "alpha_dual": sum(1 for c in self.cells if c["tra2"]),
            "doublets": sum(1 for c in self.cells if c["is_doublet"]),
            "doublet_fraction": round(sum(1 for c in self.cells if c["is_doublet"]) / len(self.cells), 4),
        }
