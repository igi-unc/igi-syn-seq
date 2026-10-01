"""Per-cell molecule sampling, shared by every single-cell assay.

The four single-cell libraries -- 10x 5' gene expression, 10x TCR, Kinnex single cell and ONT single cell
-- are different reads of the same suspension. They must agree cell by cell, and they must agree about
which *molecules* each cell contained: a UMI seen in gene expression for a given barcode and gene should
be the same molecule Kinnex sees, not an independently invented one. That only holds if the molecule pool
is drawn once and each assay samples from it, which is what this module provides.

A molecule is a (cell, transcript, UMI) triple. Expression per cell comes from the cell's type and, for a
tumour cell, its clone: a clone's copy number at a locus scales the transcript's abundance, so per-cell
genotypes and dosage reproduce the CCFs the clone model declares. Which allele a molecule came from is
recorded, because allele-specific expression is observable in long reads and is part of the truth.
"""
import bisect
import csv
import gzip
import random

UMI_LEN = 10
# 10x 5' v2: the cell barcode and UMI sit on read 1 with the TSO; cDNA is on read 2.
R1_STRUCTURE = "16bp barcode + 10bp UMI + TSO"


def umi(rng):
    return "".join(rng.choice("ACGT") for _ in range(UMI_LEN))


class MoleculePool:
    """Molecules per cell, drawn once so every single-cell assay samples the same pool."""

    def __init__(self, roster, transcripts, clones, seed, mean_molecules=8000, dispersion=0.35):
        """`transcripts` is [(id, gene, length, abundance, chrom, pos, hap, clone_or_empty)].

        `mean_molecules` is molecules captured per cell, not transcripts present: 10x captures a small
        and variable fraction, so the count is drawn per cell around the mean with the given dispersion
        and scaled by the cell type's yield.
        """
        self.rng = random.Random(f"{seed}:molecules")
        self.roster = roster
        self.clones = clones
        self.mean_molecules = mean_molecules
        self.dispersion = dispersion
        # one sampling table per cell type is wrong: a tumour cell's abundances depend on its clone's copy
        # number, so the table is built per clone and shared by cells of that clone
        self.tx = transcripts
        self._tables = {}

    def _table(self, key, scale_fn):
        """Cumulative abundance table for one cell class, built once and reused."""
        if key not in self._tables:
            w, cum, tot = [], [], 0.0
            for t in self.tx:
                a = max(0.0, scale_fn(t))
                tot += a
                cum.append(tot)
            self._tables[key] = (cum, tot)
        return self._tables[key]

    def _scale_for(self, cell):
        ct, clone = cell["cell_type"], cell["clone"]
        if ct != "tumor":
            # non-tumour cells are diploid and carry no designed somatic events
            return lambda t: t[3] * (1.0 if not t[7] else 0.0)

        def f(t):
            a = t[3]
            if t[7] and t[7] != clone and not self.clones.is_descendant(clone, t[7]):
                return 0.0          # a designed transcript of a clone this cell is not descended from
            if t[4] and t[5]:
                cn = self.clones.cn(clone, t[4], int(t[5]))[:2]
                hap = int(t[6]) if str(t[6]).isdigit() else 0
                base = sum(self.clones.base)
                a *= (cn[hap] * 2.0 / base) if base else 1.0
            return a
        return f

    def draw(self, cell):
        """Molecules for one cell: [(transcript_index, umi)]."""
        n = max(50, int(self.rng.gauss(self.mean_molecules * cell["umi_scale"],
                                       self.mean_molecules * cell["umi_scale"] * self.dispersion)))
        key = ("tumor", cell["clone"]) if cell["cell_type"] == "tumor" else ("other",)
        cum, tot = self._table(key, self._scale_for(cell))
        if tot <= 0:
            return []
        out = []
        for _ in range(n):
            i = bisect.bisect_left(cum, self.rng.random() * tot)
            if i < len(self.tx):
                out.append((i, umi(self.rng)))
        return out

    def write(self, path, cells=None, limit_per_cell=None):
        """The molecule truth table: which cell held which transcript under which UMI.

        This is what makes a single-cell truth set checkable. Without it a caller's barcode-UMI-gene
        matrix cannot be compared against anything, only inspected.
        """
        cells = cells if cells is not None else self.roster.cells
        n = 0
        with gzip.open(path, "wt") as fh:
            fh.write("barcode\tcell_type\tclone\ttranscript\tgene\tumi\thaplotype\n")
            for c in cells:
                mols = self.draw(c)
                if limit_per_cell:
                    mols = mols[:limit_per_cell]
                for i, u in mols:
                    t = self.tx[i]
                    fh.write(f"{c['barcode']}\t{c['cell_type']}\t{c['clone']}\t{t[0]}\t{t[1]}\t{u}\t{t[6]}\n")
                    n += 1
        return n
