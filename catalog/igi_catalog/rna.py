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
    # An allelic-expression setting redistributes between haplotypes that exist. It can never give
    # expression to a haplotype with no copies, and "balanced" means "no override" rather than "50/50",
    # so a gene in a one-haplotype region keeps the copy-number split.
    if ase in (None, "", "balanced") or mutant_hap not in (0, 1):
        return frac
    if frac[0] <= 0 or frac[1] <= 0:
        return frac
    m = ASE_FRACTION.get(ase)
    if m is None:
        return frac
    return (m, 1 - m) if mutant_hap == 0 else (1 - m, m)


def dosage(clones, clone, chrom, pos):
    """Expression scaling from local copy number, relative to this clone's own baseline ploidy.

    Total RNA per cell does not scale with whole-genome doubling, so the reference point is the clone's
    own ploidy rather than a fixed 2. Only departures from that baseline, gains and losses, change a
    gene's expression.
    """
    a, b, _ = clones.cn(clone, chrom, pos)
    base = sum(clones.base)
    return (a + b) / base if base else 1.0


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
        """Index allelic-expression settings by gene and by the clone that acquired the event.

        A gene can carry several events that disagree; the one with the highest expected VAF wins, since
        that is the allele most of the tumour's RNA would come from. The setting applies only in clones
        descended from the clone that acquired it.
        """
        best = {}
        for r in events_rows:
            g = r.get("gene")
            if not g:
                continue
            self.events_by_gene[g].append(r)
            ase = r.get("allelic_expression")
            if not ase:
                continue
            vaf = float(r.get("expected_vaf_tumor") or 0)
            key = g
            if key not in best or vaf > best[key][0]:
                best[key] = (vaf, ase, int(r.get("haplotype", 0) or 0), r.get("clone", "T"))
        self.ase_by_gene = {g: (v[1], v[2], v[3]) for g, v in best.items()}

    def add_transcripts(self, transcripts, events_by_chrom):
        """Reference transcripts, rebuilt on both haplotypes of every clone."""
        for t in transcripts:
            base = self.env.expr.transcript(t.tid, t.gene_id)
            if base < self.min_tpm and t.gene_name not in self.events_by_gene:
                continue
            mid = (t.start + t.end) // 2
            ase, mut_hap, ase_clone = self.ase_by_gene.get(t.gene_name, (None, None, None))
            for src, w in self.weights.items():
                clone = "T" if src == "NORMAL" else src
                if src == "NORMAL":
                    split, dose = (0.5, 0.5), 1.0
                else:
                    carries = ase_clone is not None and self.env.clones.is_descendant(clone, ase_clone)
                    split = haplotype_split(self.env.clones, clone, t.chrom, mid,
                                            ase if carries else None, mut_hap)
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
                        # the single-cell assays scale a molecule's abundance by the cell's clone copy
                        # number at the locus, so the locus has to travel with the record
                        "chrom": t.chrom, "pos": (t.start + t.end) // 2,
                    })

    def add_designed(self, records, events_by_chrom=None):
        """Fusion, ERV, splice-isoform, CTA and viral transcripts.

        A transcript belonging to a clone is expressed by that clone and every clone descended from it,
        so a truncal fusion appears in all tumour cells rather than only in the truncal population.

        A record may carry a `transcript` instead of a fixed `sequence`, in which case the sequence is
        built inside the clone and haplotype loop with that clone's germline and somatic edits applied.
        A cancer-testis antigen built from bare reference exons carries none of the variants designed
        into it, so its reads would contradict the truth table.
        """
        for rec in records:
            origin = rec.get("clone", "T")
            tpm = float(rec.get("tpm") or 0)
            t = rec.get("transcript")
            if tpm <= 0 or not (rec.get("sequence") or t is not None):
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
                # Allelic expression of a designed class. A cancer-testis antigen or an ERV locus is
                # de-repressed epigenetically, which acts on both alleles, so emitting it from one
                # haplotype only makes every germline heterozygous site inside it read as homozygous. A
                # splice isoform caused by a somatic splice-site variant is genuinely confined to that
                # variant's haplotype; an isoform switch with no genomic cause is not.
                haps = rec.get("haplotypes")
                if haps is None:
                    haps = [int(rec.get("hap", 0) or 0)]
                if len(haps) > 1 and src != "NORMAL":
                    split = haplotype_split(self.env.clones, src, rec["chrom"], int(rec["pos"]))
                elif len(haps) > 1:
                    split = (0.5, 0.5)
                else:
                    split = None
                for hap in haps:
                    hshare = share * (split[hap] if split else 1.0)
                    if hshare <= 0:
                        continue
                    self._emit(rec, t, src, hap, tpm, hshare, events_by_chrom)

    def _emit(self, rec, t, src, hap, tpm, share, events_by_chrom):
        """One record for this designed transcript on one clone and one haplotype."""
        seq = rec.get("sequence")
        if t is not None:
            es = EditSet(t.chrom)
            es.edits += germline_edits(self.env.germline, t.chrom, hap, t.start, t.end).edits
            if src != "NORMAL":
                es.edits += somatic_edits((events_by_chrom or {}).get(t.chrom, []), self.env.clones,
                                          t.chrom, hap, src, t.start, t.end).edits
            seq = exon_sequence(self.env.genome, t, es)
        if not seq or len(seq) < 150:
            return
        self.records.append({
            "id": f"{rec['id']}|{src}|hap{hap}", "clone": src, "hap": hap,
            "source": rec.get("source", "designed"), "gene": rec.get("gene", ""),
            "sequence": seq, "abundance": tpm * share,
            # as above: carried for the single-cell assays' copy-number scaling
            "chrom": rec.get("chrom") or (t.chrom if t is not None else ""),
            "pos": rec.get("pos") or ((t.start + t.end) // 2 if t is not None else 0),
            "reconcile": rec.get("reconcile", "add"),
            "reconcile_gene": rec.get("reconcile_gene") or rec.get("gene", ""),
            "normal_gene_tpm": rec.get("normal_gene_tpm"),
        })

    def reconcile_designed(self, log=print):
        """Make a designed transcript take its abundance from its gene instead of adding to it.

        Designed records were appended while the gene's reference transcripts stayed at full baseline, so
        a gene with a designed splice isoform expressed more in total than the same gene without one, and
        a cancer-testis antigen was emitted twice: once as the designed record at its tumour tier and
        once as the gene's own reference transcript at the cohort baseline, including in the normal
        fraction, which defeats the tumour specificity the class exists to test.

        Three behaviours, declared per record:

        `replace`     the designed record is the gene's tumour expression, so the reference transcripts
                      of that gene drop to zero in tumour clones and to the gene's normal-tissue level in
                      the normal fraction. Cancer-testis antigens.
        `redistribute` the designed record takes a share of the gene's existing budget, so the reference
                      transcripts scale down by that amount and the gene's total is unchanged. Splice
                      isoforms and the 5' partner of a fusion.
        `add`         a locus of its own rather than a gene's transcript, so nothing is taken away.
                      Endogenous retroviruses and viral transcripts.
        """
        ref, des = defaultdict(list), defaultdict(list)
        for r in self.records:
            if r["source"] == "reference":
                ref[(r.get("gene") or "", r["clone"])].append(r)
            else:
                # a fusion's display gene is "A-B", which matches no reference gene; the share comes off
                # the 5' partner, so grouping uses the gene the record declares it draws from
                des[(r.get("reconcile_gene") or r.get("gene") or "", r["clone"])].append(r)
        stats = defaultdict(int)
        for key, drecs in des.items():
            gene, clone = key
            if not gene:
                continue
            refs = ref.get(key)
            if not refs:
                continue
            have = sum(r["abundance"] for r in refs)
            if have <= 0:
                continue
            replace = [d for d in drecs if d.get("reconcile") == "replace"]
            redistribute = sum(d["abundance"] for d in drecs if d.get("reconcile") == "redistribute")
            if replace:
                if clone == "NORMAL":
                    # keep the gene at what normal tissue really shows, not at the cohort tumour median
                    normal = max((d.get("normal_gene_tpm") or 0) for d in replace)
                    factor = min(1.0, float(normal) / have) if normal else 0.0
                else:
                    factor = 0.0
                stats[f"replace/{clone}"] += 1
            elif redistribute > 0:
                factor = max(0.0, (have - redistribute) / have)
                stats[f"redistribute/{clone}"] += 1
            else:
                continue
            for r in refs:
                r["abundance"] *= factor
        if stats:
            log("  reconciled designed transcripts against their genes: "
                + ", ".join(f"{k} x{v}" for k, v in sorted(stats.items())))
        return dict(stats)

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

    def assign_record_ids(self):
        """Give every record a stable id, reproducible from the seed rather than from `hash()`."""
        for i, r in enumerate(sorted(self.records, key=lambda x: (x["id"], x["clone"], x["hap"]))):
            r["rec"] = f"r{i:09d}"

    def write_bin(self, recs, path):
        """FASTA for one abundance bin. Record names are stable ids, not descriptions, so nothing about
        the source can be read off a read; the manifest maps the id back to clone, haplotype and gene."""
        with open(path, "w") as fh:
            for r in recs:
                fh.write(f">{r['rec']}\n")
                s = r["sequence"]
                for j in range(0, len(s), 60):
                    fh.write(s[j:j + 60] + "\n")
        return path

    def attach_counts(self, read_map_path, read_len):
        """Record how many pairs each transcript actually produced, alongside how many were planned.

        Fragment ends cost reads on short molecules, so the realised count is below the plan and the gap
        widens as abundance falls. Consumers of the truth bundle should use the observed column.
        """
        import gzip
        from collections import Counter
        counts = Counter()
        with gzip.open(read_map_path, "rt") as fh:
            next(fh, None)
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) > 1:
                    counts[parts[1]] += 1
        for r in self.records:
            r["planned_pairs"] = r.get("bin_coverage", 0) * len(r["sequence"]) / (2 * read_len)
            r["observed_pairs"] = counts.get(r.get("rec", ""), 0)

    def write_manifest(self, path):
        cols = ["rec", "id", "source", "gene", "clone", "hap", "length", "abundance",
                "ideal_coverage", "bin", "bin_coverage", "planned_pairs", "observed_pairs"]
        with open(path, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in self.records:
                fh.write("\t".join(str(x) for x in [
                    r.get("rec", ""), r["id"], r["source"], r.get("gene", ""), r["clone"], r["hap"],
                    len(r["sequence"]), f"{r['abundance']:.8g}",
                    f"{r.get('coverage', 0):.8g}", r.get("bin", ""),
                    f"{r.get('bin_coverage', 0):.8g}",
                    f"{r.get('planned_pairs', 0):.1f}", r.get("observed_pairs", "")]) + "\n")
        return path
