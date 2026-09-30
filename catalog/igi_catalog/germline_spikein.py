"""Pathogenic germline alleles added to a baseline genotype.

The design gives each dataset a germline driver, and until now it existed only as a note in the driver
list. A germline allele is not a neoantigen source, but it has to be in the sequence: it is present in
both the tumour and the normal library, it is the reason the tumour's somatic LOH at that locus matters,
and a caller that reports it as somatic is making a mistake the benchmark should be able to see.

The allele is specified by coding position rather than by genomic coordinate, so the coordinate is derived
from the annotation and checked against the expected codon, wild-type residue and consequence before it is
written. A hard-coded position would be silently wrong on a different annotation release.
"""
from .annotation import CodingModel
from .genome import left_align


class GermlineSpikein:
    def __init__(self, env, ds_name, ds_cfg):
        self.env = env
        self.ds = ds_name
        self.dcfg = ds_cfg
        self.by_name = {}
        for t in env.rep.values():
            self.by_name.setdefault(t.gene_name, t)

    def _retained_haplotype(self, label, chrom):
        """The haplotype a labelled truncal CNA keeps, so the pathogenic copy is the surviving one."""
        for r in self.dcfg["cna"].get("T", []):
            if r.get("label") == label:
                kept = [h for h, cn in enumerate(r["cn"]) if cn > 0]
                if len(kept) != 1:
                    raise ValueError(f"{label}: expected exactly one retained haplotype, got cn={r['cn']}")
                return kept[0]
        raise KeyError(f"no truncal CNA labelled {label} in {self.ds}")

    def resolve(self, spec):
        """Turn a coding-level specification into a checked, left-aligned VCF record."""
        t = self.by_name.get(spec["gene"])
        if t is None or not t.cds:
            raise KeyError(f"{spec['gene']}: no representative coding transcript")
        cm = CodingModel(t, self.env.genome)
        c_from, c_to = spec["cds_delete"]          # 1-based inclusive coding positions
        idx = list(range(c_from - 1, c_to))
        if idx[-1] >= len(cm.map):
            raise ValueError(f"{spec['gene']}: c.{c_from}_{c_to} is past the coding sequence")
        deleted_cds = cm.cds[idx[0]:idx[-1] + 1]
        gpos = sorted(cm.map[i] for i in idx)
        if gpos[-1] - gpos[0] != len(idx) - 1:
            raise ValueError(f"{spec['gene']}: c.{c_from}_{c_to} spans an intron; not representable as one record")
        # the consequence comes from the coding-derived coordinates, which always start at the first
        # deleted base; the VCF record is left-aligned afterwards and may sit further left
        mprot, k, cons = cm.mutate_protein(gpos[0], self.env.genome.seq(t.chrom, gpos[0] - 1, gpos[-1]), "")
        if mprot is None:
            raise ValueError(f"{spec['id']}: {cons} applying c.{c_from}_{c_to} to {t.tid}")
        anchor = gpos[0] - 1
        ref = self.env.genome.seq(t.chrom, anchor - 1, gpos[-1])
        pos1, ref, alt = left_align(self.env.genome, t.chrom, anchor, ref, ref[0])
        exp = spec.get("expect", {})
        codon = idx[0] // 3 + 1
        problems = []
        if "codon" in exp and codon != exp["codon"]:
            problems.append(f"codon {codon} != expected {exp['codon']}")
        if "wt_aa" in exp and cm.protein[codon - 1] != exp["wt_aa"]:
            problems.append(f"wild-type residue {cm.protein[codon - 1]} != expected {exp['wt_aa']}")
        if "consequence" in exp and cons != exp["consequence"]:
            problems.append(f"consequence {cons} != expected {exp['consequence']}")
        if problems:
            raise ValueError(f"{spec['id']}: {'; '.join(problems)}")
        if self.env.germline.overlaps_variant(t.chrom, pos1, len(ref)):
            raise ValueError(f"{spec['id']}: the baseline genotype already has a variant over "
                             f"{t.chrom}:{pos1}-{pos1 + len(ref) - 1}")
        hap = (self._retained_haplotype(spec["cna_label"], t.chrom)
               if spec["haplotype"] == "retained_after" else int(spec["haplotype"]))
        return {
            "event_id": spec["id"], "dataset": self.ds, "class": "germline", "subclass": "pathogenic",
            "gene": spec["gene"], "transcript": t.tid, "chrom": t.chrom, "pos": pos1,
            "ref": ref, "alt": alt, "haplotype": hap, "genotype": ("1|0" if hap == 0 else "0|1"),
            "cds_change": f"c.{c_from}_{c_to}del{deleted_cds}",
            "codon": codon, "wt_aa": cm.protein[codon - 1], "consequence": cons,
            "protein_change": f"p.{cm.protein[codon - 1]}{codon}{mprot[k] if k < len(mprot) else ''}fs",
            "wt_protein_len": len(cm.protein), "mut_protein_len": len(mprot), "aa_index": k,
            "wt_protein": cm.protein, "mut_protein": mprot,
            "somatic_loh_label": spec.get("cna_label", ""),
            "note": spec.get("note", ""),
        }
