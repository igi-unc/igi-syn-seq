"""Gene-fusion designer: partner selection, exon junction choice, genomic breakpoints, and junction peptides."""
from .binding import NetMHCpan, spanning_peptides
from .expression import Expression
from . import fusion_core as fc
from . import diffcards

CHR1TO6 = {f"chr{i}" for i in range(1, 7)}

# Partner pairs with cancer relevance; the designer picks the exon junction and the mechanism.
KNOWN_PAIRS = [
    ("ETV6", "NTRK3"), ("MYB", "NFIB"), ("ESR1", "CCDC170"), ("SEC16A", "NOTCH1"), ("MAGI3", "AKT3"),
    ("FGFR3", "TACC3"), ("EWSR1", "FLI1"), ("BCR", "ABL1"), ("TMPRSS2", "ERG"), ("EML4", "ALK"),
    ("CD74", "ROS1"), ("KIF5B", "RET"), ("PAX3", "FOXO1"), ("NAB2", "STAT6"), ("FGFR2", "BICC1"),
]
# Adjacent same-strand gene pairs: transcribed as read-through, no DNA breakpoint (tumor-associated).
READTHROUGH_PAIRS = [("SLC45A3", "ELK4"), ("CTBS", "GNG5"), ("NFATC3", "PLA2G15"), ("TSNAX", "DISC1")]


class FusionDesigner:
    def __init__(self, env, ds_name, ds_cfg, cfg, rng, next_id):
        self.env = env
        self.ds = ds_name
        self.dcfg = ds_cfg
        self.cfg = cfg
        self.rng = rng
        self.next_id = next_id
        self.by_name = {}
        for t in env.rep.values():
            self.by_name.setdefault(t.gene_name, t)
        self.events = []
        self.cards = []

    # ---------------------------------------------------------------- junction search
    def find_junction(self, t5, t3, want_in_frame=True, min_prot=120):
        """Pick (n5, from3) giving a fusion with the requested frame and a usable protein."""
        segs5, segs3 = len(fc.cds_segments(t5)), len(fc.cds_segments(t3))
        opts = []
        for i in range(1, segs5):
            for j in range(1, segs3):
                f = fc.build_fusion(self.env.genome, t5, i, t3, j)
                if f["stop_before_junction"]:
                    continue
                if f["in_frame"] != want_in_frame:
                    continue
                if want_in_frame and f["protein_len"] < min_prot:
                    continue
                if not want_in_frame and f["protein_len"] < f["junction_codon"] + 5:
                    continue          # out-of-frame needs at least a short neo-ORF past the junction
                opts.append((i, j, f))
        if not opts:
            return None
        # prefer junctions near the middle of both partners
        opts.sort(key=lambda o: abs(o[0] / segs5 - 0.5) + abs(o[1] / segs3 - 0.5))
        return opts[0] if len(opts) < 4 else self.rng.choice(opts[:4])

    @staticmethod
    def mechanism(t5, t3):
        """DNA event implied by the partners' positions and orientations."""
        if t5.chrom != t3.chrom:
            return "TRA"
        if t5.strand != t3.strand:
            return "INV"
        upstream_first = (t5.start < t3.start) if t5.strand == "+" else (t5.start > t3.start)
        return "DEL" if upstream_first else "DUP"

    # ---------------------------------------------------------------- event construction
    def make(self, g5, g3, kind, clone, flagpost=False, near_exon=False):
        t5, t3 = self.by_name.get(g5), self.by_name.get(g3)
        if t5 is None or t3 is None or not t5.cds or not t3.cds:
            return None
        want_if = kind in ("in_frame", "in_frame_novel", "readthrough")
        hit = self.find_junction(t5, t3, want_in_frame=want_if)
        if hit is None:
            return None
        n5, from3, f = hit
        if kind == "readthrough":
            mech, bp5, bp3 = "none", None, None
        else:
            mech = self.mechanism(t5, t3)
            bp5 = fc.intron_breakpoint(t5, n5 - 1, "after", self.rng, near_exon=near_exon)
            bp3 = fc.intron_breakpoint(t3, from3, "before", self.rng, near_exon=near_exon)
            if bp5 is None or bp3 is None:
                return None
        cm = self.env.clones
        chrom = t5.chrom
        pos = bp5[1] if bp5 else t5.start
        hap = self.rng.choice(cm.retained_haplotypes(clone, chrom, pos) or [0])
        vaf, mults, tcn = cm.vaf(chrom, pos, hap, clone, pre_cna=(clone == "T"))
        # junction neopeptides: 8-11mers spanning the junction codon, minus any peptide present in either parent
        peps = spanning_peptides(f["protein"], f["junction_codon"], 1, (8, 9, 10, 11))
        parent5 = self._parent_peptides(t5)
        parent3 = self._parent_peptides(t3)
        neo = sorted(p for p in peps if p not in parent5 and p not in parent3)
        res = self.env.netmhc.predict(neo, self.env.hla) if neo else {}
        best = None
        for pep, al in res.items():
            r, allele, aff = NetMHCpan.best(al)
            if best is None or r < best[0]:
                best = (r, allele, aff, pep)
        tpm5 = self.env.expr.gene(t5.gene_id)
        tpm3 = self.env.expr.gene(t3.gene_id)
        ev = {
            "event_id": self.next_id("FUS"),
            "dataset": self.ds, "class": "fusion", "subclass": kind,
            "gene_5p": g5, "gene_3p": g3, "transcript_5p": t5.tid, "transcript_3p": t3.tid,
            "chrom_5p": t5.chrom, "chrom_3p": t3.chrom, "strand_5p": t5.strand, "strand_3p": t3.strand,
            "breakpoint_5p": bp5[1] if bp5 else "", "breakpoint_3p": bp3[1] if bp3 else "",
            "mechanism": mech, "exons_5p": n5, "exon_start_3p": from3 + 1,
            "in_frame": f["in_frame"], "fusion_protein_len": f["protein_len"], "junction_codon": f["junction_codon"],
            "junction_window": f["junction_window"],
            "haplotype": hap, "clone": clone, "ccf": cm.ccf(clone),
            "clonality_tier": cm.clonality_tier(clone, chrom, pos, hap),
            "expected_vaf_dna": round(vaf, 4) if mech != "none" else "",
            "tumor_cn_at_locus": round(tcn, 2),
            "gene_tpm_5p": round(tpm5, 3), "gene_tpm_3p": round(tpm3, 3),
            "expression_tier": Expression.tier(min(tpm5, tpm3)),
            "binding_tier": NetMHCpan.tier(best[0]) if best else "na",
            "best_rank_el": round(best[0], 3) if best else None,
            "best_allele": best[1] if best else None,
            "best_peptide": best[3] if best else None,
            "n_junction_neopeptides": len(neo),
            "wes_visible": self._captured(bp5) or self._captured(bp3),
            "flagpost": flagpost,
            "chr1to6": t5.chrom in CHR1TO6 and t3.chrom in CHR1TO6,
        }
        self._card(ev, f, t5, t3, n5, from3)
        self.events.append(ev)
        return ev

    def _captured(self, bp, flank=50):
        """Whether a breakpoint falls inside a capture region once the kit's flank is allowed for."""
        if not bp:
            return False
        chrom, pos = bp
        return self.env.ctx.exome.any(chrom, max(0, pos - 1 - flank), pos + flank)

    def _parent_peptides(self, t):
        from .annotation import CodingModel
        cm = CodingModel(t, self.env.genome)
        prot = cm.protein
        return {prot[i:i + L] for L in (8, 9, 10, 11) for i in range(max(0, len(prot) - L + 1))}

    def _card(self, ev, f, t5, t3, n5, from3):
        lines = [f"  fusion    {ev['gene_5p']}(ex1-{n5}) :: {ev['gene_3p']}(ex{from3 + 1}-end)  {ev['mechanism']}"
                 f"  {'in-frame' if ev['in_frame'] else 'out-of-frame'}"]
        if ev["breakpoint_5p"]:
            lines.append(f"  DNA       {ev['chrom_5p']}:{ev['breakpoint_5p']}({ev['strand_5p']}) -> {ev['chrom_3p']}:{ev['breakpoint_3p']}({ev['strand_3p']})")
        else:
            lines.append("  DNA       none (read-through transcript; tumor-associated)")
        lines += fc.junction_diff(self.env.genome, t5, n5, t3, from3)
        prot = [f"  protein   fusion length {f['protein_len']} aa, junction at aa {f['junction_codon'] + 1}",
                f"  junction  {f['junction_window']}",
                f"            {' ' * min(10, f['junction_codon'])}^"]
        extra = {"peptide": f"{ev['best_peptide']} {ev['best_allele']} rank {ev['best_rank_el']}" if ev["best_peptide"] else "none",
                 "VAF": ev["expected_vaf_dna"] or "n/a (no DNA breakpoint)",
                 "clone": f"{ev['clone']} ({ev['clonality_tier']})"}
        title = f"{ev['gene_5p']}-{ev['gene_3p']} {ev['subclass']}"
        self.cards.append(diffcards.render(ev["event_id"], title, lines, prot, extra))

    # ---------------------------------------------------------------- top level
    def design(self, n_total, clones, log=print):
        plan = []
        n_if = max(1, int(round(n_total * 0.50)))
        n_oof = max(1, int(round(n_total * 0.20)))
        n_rt = max(1, int(round(n_total * 0.10)))
        pairs = list(KNOWN_PAIRS)
        self.rng.shuffle(pairs)
        plan += [(p, "in_frame") for p in pairs[:n_if]]
        plan += [(p, "out_of_frame") for p in pairs[n_if:n_if + n_oof]]
        rt = list(READTHROUGH_PAIRS)
        self.rng.shuffle(rt)
        plan += [(p, "readthrough") for p in rt[:n_rt]]
        # remaining slots: novel partner pairs drawn from expressed genes
        expressed = [t for t in self.env.rep.values()
                     if self.env.expr.gene(t.gene_id) > 10 and len(fc.cds_segments(t)) >= 4]
        self.rng.shuffle(expressed)
        i = 0
        while len(plan) < n_total and i + 1 < len(expressed):
            plan.append(((expressed[i].gene_name, expressed[i + 1].gene_name), "in_frame_novel"))
            i += 2
        made = 0
        for (g5, g3), kind in plan:
            clone = "T" if made < max(1, n_total // 3) else self.rng.choice(clones)
            ev = self.make(g5, g3, kind, clone, flagpost=(made < 4 and kind == "in_frame"),
                           near_exon=(made % 5 < 2))
            if ev:
                made += 1
            else:
                log(f"  fusion {g5}-{g3} ({kind}): no usable junction, skipped")
            if made >= n_total:
                break
        log(f"  fusions: {made}/{n_total}")
        return self.events
