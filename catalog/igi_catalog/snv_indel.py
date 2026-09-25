"""SNV and small-indel designer: builds a scored candidate pool from the annotation and fills the evidence grid."""
import random
from collections import defaultdict

from .annotation import CodingModel
from .binding import NetMHCpan, spanning_peptides
from .expression import Expression
from .context import ContextAnnotator
from . import diffcards

PRIMARY = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
CHR1TO6 = {f"chr{i}" for i in range(1, 7)}
AA3 = {}  # not needed; single-letter codes used throughout


class Candidate(dict):
    """A candidate somatic change with everything needed to place it in a tier cell."""


class SnvIndelDesigner:
    def __init__(self, env, ds_name, ds_cfg, cfg, rng):
        self.env = env            # namespace: genome, tx, rep, expr, germline, ctx, clones, netmhc, hla
        self.ds = ds_name
        self.dcfg = ds_cfg
        self.cfg = cfg
        self.rng = rng
        self.driver_genes = {d[0] for evs in ds_cfg.get("drivers", {}).values() for d in evs}
        self.excluded = self._excluded_regions()
        self.models = {}
        self.events = []
        self.cards = []
        self.n_id = defaultdict(int)

    # ------------------------------------------------------------------ helpers
    def _excluded_regions(self):
        """Regions where random events are not placed: chromothripsis arm, homozygous deletions."""
        cm = self.env.clones
        out = []
        ct = self.dcfg.get("chromothripsis")
        if ct:
            out.append(cm.region(ct["region"]))
        out += [(c, s, e) for c, s, e, _ in cm.regions_with(lambda a, b: a + b == 0)]
        return out

    def _excluded(self, chrom, pos1):
        return any(c == chrom and s < pos1 <= e for c, s, e in self.excluded)

    def model(self, t):
        if t.tid not in self.models:
            self.models[t.tid] = CodingModel(t, self.env.genome)
        return self.models[t.tid]

    def next_id(self, prefix):
        self.n_id[prefix] += 1
        return f"{self.ds}-{prefix}-{self.n_id[prefix]:04d}"

    # ------------------------------------------------------------------ candidate generation
    def sample_sites(self, n_genes, per_gene, kinds, region_filter=None, gene_filter=None):
        """Random CDS positions in representative transcripts; returns candidates classified by consequence."""
        genes = [t for t in self.env.rep.values() if t.chrom in PRIMARY and t.gene_name not in self.driver_genes
                 and (gene_filter is None or gene_filter(t))]
        self.rng.shuffle(genes)
        out = []
        for t in genes:
            if len(out) >= n_genes * per_gene:
                break
            cm = self.model(t)
            if len(cm.protein) < 60 or "*" in cm.protein:
                continue
            for _ in range(per_gene):
                i = self.rng.randrange(3, len(cm.map) - 3)
                gpos = cm.map[i]
                if self._excluded(t.chrom, gpos):
                    continue
                if region_filter and not region_filter(t.chrom, gpos):
                    continue
                cand = self._make_candidate(t, cm, gpos, kinds)
                if cand:
                    out.append(cand)
        return out

    def _make_candidate(self, t, cm, gpos, kinds):
        g = self.env.genome
        ref = g.seq(t.chrom, gpos - 1, gpos)
        alts = [b for b in "ACGT" if b != ref]
        self.rng.shuffle(alts)
        tries = []
        if "snv" in kinds:
            tries += [(ref, a) for a in alts]
        if "indel" in kinds:
            L = self.rng.choice([1, 1, 2, 3, 4, 5, 6, 9, 12, 15, 21, 30])
            anchor = g.seq(t.chrom, gpos - 1, gpos + L)          # deletion: ref = anchor + L bases
            tries.append((anchor, anchor[0]))
            ins = "".join(self.rng.choice("ACGT") for _ in range(L))
            tries.append((ref, ref + ins))                       # insertion after the anchor base
        for r, a in tries:
            # a somatic allele overlapping a germline variant cannot be written unambiguously in
            # reference coordinates, and the haplotype edit would not apply; skip such sites
            if self.env.germline.overlaps_variant(t.chrom, gpos, len(r)):
                continue
            mprot, k, cons = cm.mutate_protein(gpos, r, a)
            if mprot is None:
                continue
            want = {"snv": ("missense", "nonsense", "synonymous", "start_loss", "stop_loss"),
                    "indel": ("frameshift", "inframe_insertion", "inframe_deletion")}
            if cons not in sum((want[k] for k in kinds), ()):
                continue
            return self._annotate(t, cm, gpos, r, a, mprot, k, cons)
        return None

    def _annotate(self, t, cm, gpos, ref, alt, mprot, k, cons):
        f = self.env.ctx.features(t.chrom, gpos, cm)
        tpm = self.env.expr.gene(t.gene_id)
        c = Candidate(chrom=t.chrom, pos=gpos, ref=ref, alt=alt, gene=t.gene_name, gene_id=t.gene_id, transcript=t.tid,
                      consequence=cons, aa_index=k, wt_protein=cm.protein, mut_protein=mprot,
                      gene_tpm=round(tpm, 3), expression_tier=Expression.tier(tpm),
                      context_features=f, context=ContextAnnotator.stratum(f),
                      nmd="escape" if (cons in ("nonsense", "frameshift") and cm.is_last_coding_exon(gpos)) else
                          ("sensitive" if cons in ("nonsense", "frameshift") else "na"))
        c["aa_change"] = self._aa_change(cm.protein, mprot, k, cons)
        c["neo_window"] = self._neo_window(cm.protein, mprot, k, cons)
        return c

    @staticmethod
    def _aa_change(wt, mut, k, cons):
        if k is None:
            return ""
        if cons == "missense":
            return f"{wt[k]}{k + 1}{mut[k]}"
        if cons == "nonsense":
            return f"{wt[k]}{k + 1}*"
        if cons == "synonymous":
            return f"{wt[k] if k < len(wt) else '?'}{k + 1}="
        if cons == "frameshift":
            return f"{wt[k] if k < len(wt) else '?'}{k + 1}fs*{max(0, len(mut) - k)}"
        if cons in ("inframe_insertion", "inframe_deletion"):
            return f"{wt[k] if k < len(wt) else '?'}{k + 1}{'ins' if cons.endswith('insertion') else 'del'}"
        if cons == "stop_loss":
            return f"*{len(wt) + 1}ext*{len(mut) - len(wt)}"
        return cons

    @staticmethod
    def _neo_window(wt, mut, k, cons):
        """Mutant protein fragment containing every peptide (8-11mer) that spans a changed residue."""
        if k is None or cons in ("synonymous", "nonsense"):
            return ""
        if cons in ("missense", "inframe_insertion", "inframe_deletion", "stop_loss"):
            n_changed = 1 if cons == "missense" else max(1, len(mut) - len(wt)) if cons != "inframe_deletion" else 1
            return mut[max(0, k - 10):k + n_changed + 10]
        if cons == "frameshift":
            return mut[max(0, k - 10):]          # whole neo-ORF
        return ""

    # ------------------------------------------------------------------ binding
    def _n_changed(self, c):
        """How many mutant residues a peptide must span to be a neoepitope."""
        cons = c["consequence"]
        if cons == "missense":
            return 1
        if cons == "inframe_insertion":
            return max(1, (len(c["alt"]) - len(c["ref"])) // 3)
        if cons == "inframe_deletion":
            return 1                      # the junction residue
        if cons in ("frameshift", "stop_loss", "start_loss"):
            mut, k = c["mut_protein"], c["aa_index"]
            return max(1, len(mut) - k) if (mut and k is not None) else 1
        return 0

    def score(self, cands, lengths=(8, 9, 10, 11)):
        """Assign binding tiers: score every mutation-spanning peptide, discard peptides also present in
        the wild-type protein, and take the best remaining %rank_EL across the patient's alleles."""
        want = {}
        for c in cands:
            n = self._n_changed(c)
            if not n or c["aa_index"] is None or not c.get("mut_protein"):
                c["_mut_peps"], c["_wt_peps"] = set(), set()
                continue
            mp = spanning_peptides(c["mut_protein"], c["aa_index"], n, lengths)
            wp = spanning_peptides(c["wt_protein"], c["aa_index"], 1, lengths) if c.get("wt_protein") else set()
            c["_mut_peps"], c["_wt_peps"] = mp, wp
            for p in mp:
                want[p] = True
        if want:
            self.env.netmhc.predict(list(want), self.env.hla)
        for c in cands:
            neo = sorted(c.get("_mut_peps", set()) - c.get("_wt_peps", set()))
            res = self.env.netmhc.predict(neo, self.env.hla) if neo else {}
            best = None
            n_lt2 = 0
            for pep, al in res.items():
                r, allele, aff = NetMHCpan.best(al)
                if r <= 2.0:
                    n_lt2 += 1
                if best is None or r < best[0]:
                    best = (r, allele, aff, pep)
            if best:
                c["best_rank"], c["best_allele"], c["best_aff_nM"], c["best_peptide"] = round(best[0], 3), best[1], round(best[2], 1), best[3]
                c["binding_tier"] = NetMHCpan.tier(best[0])
            else:
                c["best_rank"] = c["best_allele"] = c["best_aff_nM"] = c["best_peptide"] = None
                c["binding_tier"] = "na"
            c["n_neopeptides_rank_lt2"] = n_lt2
            c.pop("_mut_peps", None)
            c.pop("_wt_peps", None)
        self._score_wt_counterparts(cands)
        return cands

    def _wt_counterpart(self, c):
        """The wild-type peptide occupying the same window as the best mutant peptide.

        Defined only where mutant and wild-type proteins stay in register (substitutions); an indel or
        frameshift shifts everything downstream, so no counterpart exists.
        """
        if c.get("consequence") not in ("missense", "start_loss", "stop_loss"):
            return None
        pep, mut, wt = c.get("best_peptide"), c.get("mut_protein"), c.get("wt_protein")
        if not pep or not mut or not wt:
            return None
        s = mut.find(pep)
        if s < 0 or s + len(pep) > len(wt):
            return None
        cand = wt[s:s + len(pep)]
        return cand if cand != pep and "*" not in cand and "X" not in cand else None

    def _score_wt_counterparts(self, cands):
        """Rank the wild-type version of each best mutant peptide, to separate a genuinely novel epitope
        from one whose wild-type counterpart binds just as well (low agretopicity)."""
        pairs = {}
        for c in cands:
            wtp = self._wt_counterpart(c)
            if wtp:
                pairs[c["best_peptide"]] = wtp
        if pairs:
            self.env.netmhc.predict(sorted(set(pairs.values())), self.env.hla)
        for c in cands:
            wtp = pairs.get(c.get("best_peptide"))
            c["wt_peptide"] = wtp
            c["wt_peptide_rank"] = c["agretopicity"] = None
            if not wtp:
                continue
            res = self.env.netmhc.predict([wtp], self.env.hla).get(wtp)
            if not res:
                continue
            allele = c.get("best_allele")
            vals = res.get(allele) or min(res.values(), key=lambda v: v[0])
            c["wt_peptide_rank"] = round(vals[0], 3)
            if c.get("best_aff_nM") and vals[2]:
                c["agretopicity"] = round(vals[2] / c["best_aff_nM"], 3)   # >1 means the mutant binds better

    # ------------------------------------------------------------------ placement
    def place(self, c, clone, subclass, flagpost=False, hap=None, pre_cna=None, ase="balanced", pair_with=None):
        cm = self.env.clones
        if hap is None:
            retained = cm.retained_haplotypes(clone, c["chrom"], c["pos"])
            hap = self.rng.choice(retained) if retained else 0
        if pre_cna is None:
            pre_cna = (clone == "T") and (self.rng.random() < 0.5 if cm.wgd else True)
        vaf, mults, tcn = cm.vaf(c["chrom"], c["pos"], hap, clone, pre_cna=pre_cna)
        a, b, cn_label = cm.cn(clone, c["chrom"], c["pos"])
        ev = {
            "event_id": self.next_id("SNV" if len(c["ref"]) == 1 and len(c["alt"]) == 1 else "IND"),
            "dataset": self.ds, "class": "snv" if len(c["ref"]) == 1 and len(c["alt"]) == 1 else "indel", "subclass": subclass,
            "chrom": c["chrom"], "pos": c["pos"], "ref": c["ref"], "alt": c["alt"], "gene": c["gene"], "transcript": c["transcript"],
            "consequence": c["consequence"], "aa_change": c["aa_change"], "nmd": c["nmd"],
            "haplotype": hap, "clone": clone, "ccf": cm.ccf(clone), "timing": "pre_cna" if pre_cna else "post_cna",
            "clonality_tier": cm.clonality_tier(clone, c["chrom"], c["pos"], hap),
            "expected_vaf_tumor": round(vaf, 4), "multiplicity": mults.get(clone, 0), "tumor_cn_at_locus": round(tcn, 2), "cn_segment": cn_label,
            "gene_tpm": c["gene_tpm"], "expression_tier": c["expression_tier"], "allelic_expression": ase,
            "binding_tier": c.get("binding_tier", "na"), "best_rank_el": c.get("best_rank"), "best_allele": c.get("best_allele"),
            "wt_peptide": c.get("wt_peptide"), "wt_peptide_rank_el": c.get("wt_peptide_rank"), "agretopicity": c.get("agretopicity"),
            "best_peptide": c.get("best_peptide"), "best_aff_nM": c.get("best_aff_nM"), "n_neopeptides_rank_lt2": c.get("n_neopeptides_rank_lt2", 0),
            "context": c["context"], "flagpost": flagpost, "chr1to6": c["chrom"] in CHR1TO6, "paired_event": pair_with or "",
        }
        ev.update({f"ctx_{k}": v for k, v in c["context_features"].items()})
        self._card(ev, c)
        self.events.append(ev)
        return ev

    def _card(self, ev, c):
        dna, wt_w, mut_w = diffcards.dna_card(self.env.genome, self.env.germline, ev["chrom"], ev["pos"], ev["ref"], ev["alt"], ev["haplotype"])
        prot, wt_p, mut_p = diffcards.protein_card(c["wt_protein"], c["mut_protein"], c["aa_index"], c["consequence"])
        ev["wt_hap_window"], ev["mut_hap_window"], ev["wt_protein_window"], ev["mut_protein_window"] = wt_w, mut_w, wt_p, mut_p
        title = f"{ev['gene']} {ev['aa_change']} {ev['consequence']}  clone={ev['clone']} tier={ev['clonality_tier']} VAF={ev['expected_vaf_tumor']} expr={ev['expression_tier']} bind={ev['binding_tier']} ctx={ev['context']}"
        extra = {"peptide": f"{ev['best_peptide']} {ev['best_allele']} rank {ev['best_rank_el']}" if ev["best_peptide"] else "none"}
        if ev.get("wt_peptide"):
            extra["wt pep"] = (f"{ev['wt_peptide']} rank {ev['wt_peptide_rank_el']}"
                               f"  (agretopicity {ev['agretopicity']})")
        self.cards.append(diffcards.render(ev["event_id"], title, dna, prot, extra))

    # ------------------------------------------------------------------ high-level design
    def design(self, log=print):
        cfg = self.cfg["counts"]
        tiers = self.cfg["tiers"]
        clones = list(self.env.clones.clones)
        # 1. hotspots (flagposts)
        for gene, change in self.cfg.get("hotspots", []):
            c = self.hotspot(gene, change)
            if c:
                self.place(c, "T", "hotspot", flagpost=True, pre_cna=True)
            else:
                log(f"  hotspot {gene} {change}: not resolvable in representative transcript, skipped")
        # 2. missense pool, scored
        ps = float(self.cfg.get("pool_scale", 1.0))
        P = lambda n: max(20, int(n * ps))
        log(f"  sampling missense pool (pool_scale={ps})")
        pool = [c for c in self.sample_sites(P(6000), 2, ("snv",)) if c["consequence"] == "missense"]
        # extra sampling inside LOH and amplified regions so those clonality tiers can be filled
        cm = self.env.clones
        loh = cm.regions_with(lambda a, b: (a == 0) != (b == 0))
        amp = cm.regions_with(lambda a, b: max(a, b) >= 4)
        def in_regions(regs):
            return lambda chrom, pos: any(c == chrom and s < pos <= e for c, s, e, _ in regs)
        pool += [c for c in self.sample_sites(P(1500), 3, ("snv",), region_filter=in_regions(loh)) if c["consequence"] == "missense"]
        # amplified regions hold only a few dozen genes, so sample many sites per gene there
        pool += [c for c in self.sample_sites(P(400), 40, ("snv",), region_filter=in_regions(amp)) if c["consequence"] == "missense"]
        log(f"  missense candidates: {len(pool)}; scoring binding")
        self.score(pool)
        used = set()
        def take(pred, n, prefer16=True):
            picks = []
            avail = [c for c in pool if id(c) not in used and pred(c)]
            self.rng.shuffle(avail)
            if prefer16:
                avail.sort(key=lambda c: 0 if c["chrom"] in CHR1TO6 else 1)
                half = [c for c in avail if c["chrom"] in CHR1TO6][: (n + 1) // 2]
                rest = [c for c in avail if id(c) not in {id(x) for x in half}][: n - len(half)]
                picks = half + rest
            else:
                picks = avail[:n]
            for c in picks:
                used.add(id(c))
            return picks
        # 3. core grid: clonality x expression x binding, clean context
        reps = cfg["snv_missense_grid_replicates"]
        for ct in tiers["clonality"]:
            for et in tiers["expression"]:
                for bt in tiers["binding"]:
                    def pred(c, et=et, bt=bt, ct=ct):
                        if c["context"] != "clean" or c["expression_tier"] != et or c["binding_tier"] != bt:
                            return False
                        need = {"T_LOH": lambda a, b: (a == 0) != (b == 0), "T_amp": lambda a, b: max(a, b) >= 4}.get(ct)
                        a, b, _ = cm.cn("T", c["chrom"], c["pos"])
                        if need:
                            return need(a, b)
                        return (a >= 1 and b >= 1) and max(a, b) < 4
                    picks = take(pred, reps)
                    clone = "T" if ct.startswith("T_") else ct
                    for c in picks:
                        hap = None
                        if ct == "T_LOH":
                            a, b, _ = cm.cn("T", c["chrom"], c["pos"]); hap = 0 if a >= 1 else 1
                        if ct == "T_amp":
                            a, b, _ = cm.cn("T", c["chrom"], c["pos"]); hap = 0 if a > b else 1
                        ase = self.rng.choice(["balanced", "balanced", "silenced", "dominant"]) if et != "T0" else "balanced"
                        self.place(c, clone, "grid_missense", hap=hap, pre_cna=True if ct.startswith("T_") else None, ase=ase)
                    if len(picks) < reps:
                        log(f"  grid cell {ct}/{et}/{bt}: only {len(picks)}/{reps} filled")
        # 4. context strata (missense, any expression, any binding)
        for stratum, n in cfg["snv_missense_context"].items():
            if stratum == "phased_pair":
                continue
            if stratum == "deep_intronic":
                picks = self.deep_intronic(n)
                for c in picks:
                    self.place(c, self.rng.choice(clones), "context_deep_intronic")
                continue
            picks = take(lambda c, s=stratum: c["context"] == s, n)
            if len(picks) < n and stratum in ("homopolymer", "segdup", "near_germline_het", "exon_edge"):
                extra = [c for c in self.sample_sites(P(3000), 3, ("snv",)) if c["consequence"] == "missense" and c["context"] == stratum]
                self.score(extra); pool.extend(extra); picks += take(lambda c, s=stratum: c["context"] == s, n - len(picks))
            for c in picks:
                self.place(c, self.rng.choice(clones), f"context_{stratum}")
            log(f"  context {stratum}: {len(picks)}/{n}")
        # 5. phased pairs: second SNV within 150 bp on the same haplotype, same clone
        n_pairs = cfg["snv_missense_context"].get("phased_pair", 0)
        made = 0
        for c in take(lambda c: c["context"] == "clean", n_pairs * 2):
            if made >= n_pairs:
                break
            partner = self.partner(c)
            if partner:
                clone = self.rng.choice(clones)
                e1 = self.place(c, clone, "phased_pair")
                self.place(partner, clone, "phased_pair", hap=e1["haplotype"], pre_cna=e1["timing"] == "pre_cna", pair_with=e1["event_id"])
                e1["paired_event"] = self.events[-1]["event_id"]
                made += 1
        log(f"  phased pairs: {made}/{n_pairs}")
        # 6. other SNV consequences
        other = self.sample_sites(P(6000), 2, ("snv",))
        for cons, key in (("nonsense", "snv_nonsense"), ("synonymous", "snv_synonymous")):
            cs = [c for c in other if c["consequence"] == cons]
            self.rng.shuffle(cs)
            for c in cs[:cfg[key]]:
                self.place(c, self.rng.choice(clones), cons)
        ssl = [c for c in other if c["consequence"] in ("start_loss", "stop_loss")]
        self.score([c for c in ssl if c["consequence"] == "stop_loss"])
        for c in ssl[:cfg["snv_start_stop_loss"]]:
            self.place(c, self.rng.choice(clones), c["consequence"])
        # 7. indels
        ind = self.sample_sites(P(8000), 2, ("indel",))
        fs_del = [c for c in ind if c["consequence"] == "frameshift" and len(c["ref"]) > len(c["alt"])]
        fs_ins = [c for c in ind if c["consequence"] == "frameshift" and len(c["ref"]) < len(c["alt"])]
        if_del = [c for c in ind if c["consequence"] == "inframe_deletion"]
        if_ins = [c for c in ind if c["consequence"] == "inframe_insertion"]
        for lst, key in ((fs_del, "indel_frameshift_del"), (fs_ins, "indel_frameshift_ins"), (if_del, "indel_inframe_del"), (if_ins, "indel_inframe_ins")):
            self.rng.shuffle(lst)
            sel = lst[:cfg[key]]
            self.score(sel)
            for c in sel:
                self.place(c, self.rng.choice(clones), key.replace("indel_", ""))
            log(f"  {key}: {len(sel)}/{cfg[key]}")
        hp = [c for c in ind if c["context"] == "homopolymer"]
        if len(hp) < cfg["indel_homopolymer"]:
            hp += [c for c in self.sample_sites(P(6000), 3, ("indel",)) if c["context"] == "homopolymer"]
        self.rng.shuffle(hp)
        sel = hp[:cfg["indel_homopolymer"]]
        self.score(sel)
        for c in sel:
            self.place(c, self.rng.choice(clones), "homopolymer_indel")
        log(f"  homopolymer indels: {len(sel)}/{cfg['indel_homopolymer']}")
        return self.events

    def hotspot(self, gene, change):
        """Resolve e.g. ('TP53','R248Q') to a genomic SNV in the representative transcript."""
        wt_aa, num, mut_aa = change[0], int(change[1:-1]), change[-1]
        # legacy hotspot numbering (e.g. BRAF V600E) may refer to a non-MANE isoform: try the representative first, then every complete CDS of the gene
        rep = next((t for t in self.env.rep.values() if t.gene_name == gene), None)
        others = [t for t in self.env.tx.values() if t.gene_name == gene and t.cds and t.gene_type == "protein_coding"
                  and "cds_start_NF" not in t.tags and "cds_end_NF" not in t.tags and (rep is None or t.tid != rep.tid)]
        t = None
        for cand_t in ([rep] if rep else []) + others:
            cmx = self.model(cand_t)
            if num - 1 < len(cmx.protein) and cmx.protein[num - 1] == wt_aa:
                t = cand_t
                break
        if t is None:
            return None
        cm = self.model(t)
        for off in range(3):
            i = (num - 1) * 3 + off
            gpos = cm.map[i]
            ref = self.env.genome.seq(t.chrom, gpos - 1, gpos)
            for alt in "ACGT":
                if alt == ref:
                    continue
                mprot, k, cons = cm.mutate_protein(gpos, ref, alt)
                if mprot and k == num - 1 and cons == "missense" and mprot[k] == mut_aa:
                    c = self._annotate(t, cm, gpos, ref, alt, mprot, k, cons)
                    self.score([c])
                    return c
        return None

    def partner(self, c):
        """A second missense SNV within 150 bp downstream in the same transcript."""
        t = self.env.tx[c["transcript"]]
        cm = self.model(t)
        i = cm.cds_index(c["pos"])
        if i is None:
            return None
        for _ in range(30):
            j = i + self.rng.randrange(1, 60)
            if j >= len(cm.map) - 3:
                continue
            gpos = cm.map[j]
            if abs(gpos - c["pos"]) > 150:
                continue
            ref = self.env.genome.seq(t.chrom, gpos - 1, gpos)
            for alt in self.rng.sample("ACGT", 4):
                if alt == ref:
                    continue
                mprot, k, cons = cm.mutate_protein(gpos, ref, alt)
                if mprot and cons == "missense":
                    p = self._annotate(t, cm, gpos, ref, alt, mprot, k, cons)
                    self.score([p])
                    return p
        return None

    def deep_intronic(self, n):
        """Non-coding SNVs > 50 bp inside introns of expressed genes and outside the capture BED (WGS-only sites)."""
        out = []
        genes = [t for t in self.env.rep.values() if t.chrom in PRIMARY and t.n_exons() > 2]
        self.rng.shuffle(genes)
        for t in genes:
            if len(out) >= n:
                break
            introns = [(t.exons[i][1] + 1, t.exons[i + 1][0] - 1) for i in range(len(t.exons) - 1) if t.exons[i + 1][0] - t.exons[i][1] > 400]
            if not introns:
                continue
            s, e = self.rng.choice(introns)
            gpos = self.rng.randrange(s + 150, e - 150)
            if self._excluded(t.chrom, gpos):
                continue
            f = self.env.ctx.features(t.chrom, gpos)
            if f["on_target"]:
                continue
            ref = self.env.genome.seq(t.chrom, gpos - 1, gpos)
            alt = self.rng.choice([b for b in "ACGT" if b != ref])
            tpm = self.env.expr.gene(t.gene_id)
            c = Candidate(chrom=t.chrom, pos=gpos, ref=ref, alt=alt, gene=t.gene_name, gene_id=t.gene_id, transcript=t.tid,
                          consequence="intronic", aa_index=None, wt_protein=None, mut_protein=None, aa_change="", neo_window="",
                          gene_tpm=round(tpm, 3), expression_tier=Expression.tier(tpm), context_features=f, context="deep_intronic", nmd="na")
            out.append(c)
        return out
