"""Structural variant and viral-integration designer.

Structural variants are specified in reference coordinates with a type, a size, the haplotype and clone
that carry them, and their effect on any gene they touch. A subset is placed deliberately inside coding
sequence, so the exon-level consequence is part of the truth rather than an accident of position.

Viral events are either integrations, which add a junction in the tumour genome and a host-virus fusion
transcript, or episomes, which add copies of the viral genome without a host junction.
"""
import random

CHR1TO6 = {f"chr{i}" for i in range(1, 7)}
SV_TYPES = ["DEL", "DUP", "INV", "TRA", "INS"]
MEI_FAMILIES = ["L1HS", "AluY", "SVA_E"]

# log-spaced size classes, in bp, for the non-fusion SV set
SIZE_CLASSES = [(50, 500), (500, 5_000), (5_000, 50_000), (50_000, 500_000), (500_000, 5_000_000)]


class StructuralDesigner:
    def __init__(self, env, ds_name, ds_cfg, cfg, rng, next_id):
        self.env = env
        self.ds = ds_name
        self.dcfg = ds_cfg
        self.cfg = cfg
        self.rng = rng
        self.next_id = next_id
        self.events = []
        self.viral = []
        self.cards = []
        self.genes = [t for t in env.rep.values()
                      if t.chrom in {f"chr{i}" for i in list(range(1, 23))} | {"chrX"} and len(t.cds) >= 3]

    # ------------------------------------------------------------------ helpers
    def _pick_gene(self, expressed=True):
        for _ in range(200):
            t = self.rng.choice(self.genes)
            tpm = self.env.expr.gene(t.gene_id)
            if (tpm > 5) == expressed:
                return t, tpm
        return None, 0.0

    def _exon_span(self, t, n_exons):
        """Genomic span covering n_exons consecutive coding exons, and the exons chosen."""
        segs = t.cds if t.strand == "+" else t.cds[::-1]
        if len(segs) <= n_exons:
            return None
        i = self.rng.randrange(0, len(segs) - n_exons)
        chosen = segs[i:i + n_exons]
        lo = min(s for s, _ in chosen)
        hi = max(e for _, e in chosen)
        return lo, hi, i + 1, i + n_exons

    def _coding_effect(self, t, start1, end1, svtype):
        """Effect of removing or duplicating [start1, end1] within a transcript's CDS."""
        removed = sum(max(0, min(e, end1) - max(s, start1) + 1) for s, e in t.cds)
        if removed == 0:
            return "intronic_or_utr", 0
        total = sum(e - s + 1 for s, e in t.cds)
        if svtype == "DEL" and removed >= total * 0.95:
            return "whole_gene_loss", removed
        frame = "in_frame" if removed % 3 == 0 else "out_of_frame"
        kind = {"DEL": "exon_loss", "DUP": "exon_duplication", "INV": "exon_inversion"}.get(svtype, "exon_disruption")
        return f"{kind}_{frame}", removed

    # ------------------------------------------------------------------ SVs
    def make_sv(self, svtype, clone, coding=False, size_class=None):
        cm = self.env.clones
        gene = None
        if coding:
            gene, _tpm = self._pick_gene(expressed=True)
            if gene is None:
                return None
            span = self._exon_span(gene, self.rng.choice([1, 1, 2, 3]))
            if span is None:
                return None
            start1, end1, ex_from, ex_to = span
            # widen into the flanking introns so the breakpoints are not exactly at exon boundaries
            start1 = max(1, start1 - self.rng.randrange(30, 400))
            end1 = end1 + self.rng.randrange(30, 400)
            chrom = gene.chrom
        else:
            lo, hi = size_class or self.rng.choice(SIZE_CLASSES)
            size = self.rng.randrange(lo, hi)
            chrom = self.rng.choice(sorted({t.chrom for t in self.genes}))
            limit = self.env.genome.lengths[chrom] - size - 1_000_000
            if limit <= 1_000_000:
                return None
            start1 = self.rng.randrange(1_000_000, limit)
            end1 = start1 + size
        if svtype == "TRA":
            other = self.rng.choice([c for c in sorted({t.chrom for t in self.genes}) if c != chrom])
            end_chrom, end_pos = other, self.rng.randrange(1_000_000, self.env.genome.lengths[other] - 1_000_000)
        else:
            end_chrom, end_pos = chrom, end1
        hap = self.rng.choice(cm.retained_haplotypes(clone, chrom, start1) or [0])
        vaf, _m, tcn = cm.vaf(chrom, start1, hap, clone, pre_cna=(clone == "T"))
        effect, removed = ("na", 0)
        if gene is not None and svtype != "TRA":
            effect, removed = self._coding_effect(gene, start1, end1, svtype)
        ins_seq = ""
        if svtype == "INS":
            fam = self.rng.choice(MEI_FAMILIES)
            ins_len = {"L1HS": 6000, "AluY": 300, "SVA_E": 2000}[fam] + self.rng.randrange(-100, 100)
            effect = f"mobile_element_{fam}"
            ins_seq = f"<{fam}:{ins_len}>"
        ev = {
            "event_id": self.next_id("SV"),
            "dataset": self.ds, "class": "sv", "subclass": ("coding" if coding else "non_coding"),
            "svtype": svtype, "chrom": chrom, "start": start1,
            "end_chrom": end_chrom, "end": end_pos,
            "size": (end_pos - start1 if svtype != "TRA" else ""),
            "gene": gene.gene_name if gene else "", "transcript": gene.tid if gene else "",
            "exons_affected": (f"{ex_from}-{ex_to}" if coding and gene else ""),
            "coding_effect": effect, "cds_bp_affected": removed,
            "inserted_sequence": ins_seq,
            "haplotype": hap, "clone": clone, "ccf": cm.ccf(clone),
            "clonality_tier": cm.clonality_tier(clone, chrom, start1, hap),
            "expected_vaf_dna": round(vaf, 4), "tumor_cn_at_locus": round(tcn, 2),
            "gene_tpm": round(self.env.expr.gene(gene.gene_id), 3) if gene else "",
            "wes_visible": self._captured(chrom, start1) or self._captured(end_chrom, end_pos),
            "chr1to6": chrom in CHR1TO6 and end_chrom in CHR1TO6,
        }
        self.events.append(ev)
        return ev

    # must match the padding used when the exome intervals are built for generation
    # (pipeline.merged_capture pad=100); gap-merging there can admit a few more breakpoints, so this
    # flag is a conservative predictor and the junctions table written by the run is authoritative
    CAPTURE_PAD = 100

    def _captured(self, chrom, pos, flank=CAPTURE_PAD):
        return self.env.ctx.exome.any(chrom, max(0, pos - 1 - flank), pos + flank)

    # ------------------------------------------------------------------ viruses
    def make_viruses(self, clones):
        cm = self.env.clones
        for spec in self.dcfg.get("viruses", []):
            mode = spec["mode"]
            clone = spec.get("clone", "T")
            ev = {
                "event_id": self.next_id("VIR"),
                "dataset": self.ds, "class": "virus", "subclass": mode,
                "virus": spec["name"], "accession": spec["accession"],
                "clone": clone, "ccf": cm.ccf(clone),
                "copies_per_cell": spec.get("copies", 1),
                "expressed": spec.get("expressed", False),
                "normal_trace_copies": spec.get("normal_trace_copies", 0),
            }
            if mode == "integrated":
                chrom, pos = spec["site"].split(":")
                pos = int(pos)
                hap = self.rng.choice(cm.retained_haplotypes(clone, chrom, pos) or [0])
                vaf, _m, tcn = cm.vaf(chrom, pos, hap, clone, pre_cna=(clone == "T"))
                ev.update({
                    "chrom": chrom, "integration_pos": pos, "haplotype": hap,
                    "expected_vaf_dna": round(vaf, 4), "tumor_cn_at_locus": round(tcn, 2),
                    "host_junction": f"{chrom}:{pos}|{spec['accession']}:1",
                    "wes_visible": self._captured(chrom, pos),
                    "chr1to6": chrom in CHR1TO6,
                })
            else:
                ev.update({"chrom": "", "integration_pos": "", "haplotype": "", "expected_vaf_dna": "",
                           "tumor_cn_at_locus": "", "host_junction": "", "wes_visible": False, "chr1to6": True})
            self.viral.append(ev)
        return self.viral

    # ------------------------------------------------------------------ top level
    def design(self, n_total, clones, log=print):
        plan = []
        n_coding = max(1, n_total // 2)
        for i in range(n_coding):
            plan.append((self.rng.choice(["DEL", "DEL", "DUP", "INV"]), True, None))
        per_type = {"DEL": 25, "DUP": 20, "INV": 15, "TRA": 15, "INS": 15}
        scale = (n_total - n_coding) / sum(per_type.values())
        for t, n in per_type.items():
            for _ in range(max(1, int(round(n * scale)))):
                plan.append((t, False, self.rng.choice(SIZE_CLASSES)))
        self.rng.shuffle(plan)
        made = 0
        for svtype, coding, size_class in plan:
            if made >= n_total:
                break
            clone = "T" if made < n_total // 2 else self.rng.choice(clones)
            if self.make_sv(svtype, clone, coding=coding, size_class=size_class):
                made += 1
        log(f"  structural variants: {made}/{n_total}")
        v = self.make_viruses(clones)
        log(f"  viral events: {len(v)}")
        return self.events
