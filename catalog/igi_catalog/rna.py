"""Bulk RNA construction: per-clone, per-haplotype transcripts at abundance-proportional depth.

Read depth per transcript has to follow its abundance. ART applies one fold-coverage to every sequence in
a FASTA, so transcripts are bucketed into log-spaced abundance bins and each bin is simulated at its own
coverage. Sixty bins over the roughly four-and-a-half logs of TPM in a tumour reproduce each transcript's
intended depth to within about 8 %, at the cost of sixty simulator invocations rather than one.

Allele-specific expression comes from the copy number of each haplotype in each clone, so a lost
haplotype contributes nothing and a gained one contributes proportionally more. A designed event can
override the split for its own gene (balanced, silenced or dominant).
"""
import math
import os
from collections import defaultdict

from .genome_build import EditSet, germline_edits, somatic_edits
from .transcriptome import exon_sequence

ASE_FRACTION = {"balanced": 0.5, "silenced": 0.1, "dominant": 0.9}


def haplotype_split(clones, clone, chrom, pos, ase=None, mutant_hap=None):
    """Fraction of a gene's expression coming from each haplotype in one clone.

    Copy number sets the default split, so loss of heterozygosity silences a haplotype and a gain raises
    it. An allelic-expression setting on a designed event overrides the split for that gene, with the
    stated fraction going to the haplotype carrying the mutation.
    """
    a, b, _ = clones.cn(clone, chrom, pos)
    total = a + b
    if total == 0:
        return (0.0, 0.0)
    frac = (a / total, b / total)
    if ase and ase in ASE_FRACTION and mutant_hap in (0, 1) and frac[mutant_hap] > 0:
        m = ASE_FRACTION[ase]
        frac = (m, 1 - m) if mutant_hap == 0 else (1 - m, m)
    return frac


def dosage(clones, clone, chrom, pos, baseline_ploidy=2):
    """Overall expression scaling from total copy number at a locus."""
    a, b, _ = clones.cn(clone, chrom, pos)
    return (a + b) / baseline_ploidy if baseline_ploidy else 1.0


class RnaBuilder:
    def __init__(self, env, ds_cfg, clone_weights, min_tpm=0.5):
        self.env = env
        self.cfg = ds_cfg
        self.weights = clone_weights
        self.min_tpm = min_tpm
        self.records = []          # (id, clone, hap, sequence, abundance)
        self.ase_by_gene = {}
        self.events_by_gene = defaultdict(list)

    def index_events(self, events_rows):
        """Remember each gene's allelic-expression setting and the haplotype its mutation sits on."""
        for r in events_rows:
            g = r.get("gene")
            if not g:
                continue
            self.events_by_gene[g].append(r)
            ase = r.get("allelic_expression")
            if ase and g not in self.ase_by_gene:
                self.ase_by_gene[g] = (ase, int(r.get("haplotype", 0)))

    def add_transcripts(self, transcripts, events_by_chrom):
        """Reference transcripts, rebuilt on both haplotypes of every clone."""
        for t in transcripts:
            base = self.env.expr.transcript(t.tid, t.gene_id)
            if base < self.min_tpm and t.gene_name not in self.events_by_gene:
                continue
            mid = (t.start + t.end) // 2
            ase, mut_hap = self.ase_by_gene.get(t.gene_name, (None, None))
            for src, w in self.weights.items():
                clone = "T" if src == "NORMAL" else src
                if src == "NORMAL":
                    split, dose = (0.5, 0.5), 1.0
                else:
                    split = haplotype_split(self.env.clones, clone, t.chrom, mid, ase, mut_hap)
                    dose = dosage(self.env.clones, clone, t.chrom, mid)
                for hap in (0, 1):
                    share = split[hap]
                    if share <= 0:
                        continue
                    es = EditSet(t.chrom)
                    es.edits += germline_edits(self.env.germline, t.chrom, hap, t.start, t.end).edits
                    if src != "NORMAL":
                        es.edits += somatic_edits(events_by_chrom.get(t.chrom, []), self.env.clones,
                                                  t.chrom, hap, clone, t.start, t.end).edits
                    seq = exon_sequence(self.env.genome, t, es)
                    if len(seq) < 150:
                        continue
                    self.records.append({
                        "id": f"{t.tid}|{src}|hap{hap}", "clone": src, "hap": hap, "source": "reference",
                        "gene": t.gene_name, "sequence": seq,
                        "abundance": base * dose * share * w,
                    })

    def add_designed(self, records):
        """Fusion, ERV, splice-isoform, CTA and viral transcripts.

        A transcript belonging to a clone is expressed by that clone and every clone descended from it,
        so a truncal fusion appears in all tumour cells rather than only in the truncal population.
        """
        for rec in records:
            origin = rec.get("clone", "T")
            tpm = float(rec.get("tpm") or 0)
            if tpm <= 0 or not rec.get("sequence"):
                continue
            for src, w in self.weights.items():
                if src == "NORMAL":
                    if not rec.get("in_normal"):
                        continue
                    share = float(rec.get("normal_tpm", 0)) / tpm if tpm else 0
                    if share <= 0:
                        continue
                elif not self.env.clones.is_descendant(src, origin):
                    continue
                else:
                    share = 1.0
                hap = int(rec.get("hap", 0) or 0)
                self.records.append({
                    "id": f"{rec['id']}|{src}", "clone": src, "hap": hap,
                    "source": rec.get("source", "designed"), "gene": rec.get("gene", ""),
                    "sequence": rec["sequence"], "abundance": tpm * share * w,
                })

    # ------------------------------------------------------------------ depth
    def coverage_plan(self, target_pairs, read_len, n_bins=60):
        """Per-record coverage proportional to abundance, grouped into log-spaced bins.

        Returns [(bin_index, coverage, [records])] and the realised pair count.
        """
        live = [r for r in self.records if r["abundance"] > 0 and r["sequence"]]
        denom = sum(r["abundance"] * len(r["sequence"]) for r in live)
        if denom <= 0:
            return [], 0
        k = target_pairs * 2 * read_len / denom
        for r in live:
            r["coverage"] = k * r["abundance"]
        covs = [r["coverage"] for r in live]
        lo, hi = max(min(covs), 1e-6), max(covs)
        if hi <= lo:
            return [(0, max(hi, 0.01), live)], target_pairs
        step = (math.log10(hi) - math.log10(lo)) / n_bins
        bins = defaultdict(list)
        for r in live:
            idx = min(n_bins - 1, int((math.log10(max(r["coverage"], lo)) - math.log10(lo)) / step))
            bins[idx].append(r)
        plan = []
        pairs = 0
        for idx in sorted(bins):
            recs = bins[idx]
            # the bin's coverage is the abundance-weighted mean, so total reads are preserved
            w = sum(x["coverage"] * len(x["sequence"]) for x in recs)
            bp = sum(len(x["sequence"]) for x in recs)
            cov = w / bp if bp else 0
            if cov <= 0:
                continue
            for r in recs:
                r["bin"] = idx
                r["bin_coverage"] = cov     # what the simulator was actually given
            plan.append((idx, cov, recs))
            pairs += cov * bp / (2 * read_len)
        return plan, int(pairs)

    def write_bin(self, recs, path):
        with open(path, "w") as fh:
            for i, r in enumerate(recs):
                # names are opaque: an aligner must not be able to read the truth off a read name
                fh.write(f">r{abs(hash((r['id'], i))) % (10 ** 12):012d}\n")
                s = r["sequence"]
                for j in range(0, len(s), 60):
                    fh.write(s[j:j + 60] + "\n")
        return path

    def write_manifest(self, path):
        cols = ["id", "source", "gene", "clone", "hap", "length", "abundance",
                "ideal_coverage", "bin", "bin_coverage"]
        with open(path, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in self.records:
                fh.write("\t".join(str(x) for x in [
                    r["id"], r["source"], r.get("gene", ""), r["clone"], r["hap"],
                    len(r["sequence"]), f"{r['abundance']:.8g}",
                    f"{r.get('coverage', 0):.8g}", r.get("bin", ""),
                    f"{r.get('bin_coverage', 0):.8g}"]) + "\n")
        return path
