"""Which HLA alleles a tumour has lost, and the two binding tiers that follow (owner decision D4).

LENS ranks candidates against all of a patient's alleles today, but a tumour that has lost an HLA
haplotype cannot present through the alleles on it. Both numbers are therefore real and the truth bundle
carries both: the tier over all six alleles, which is what the evidence grid is filled on because it
matches how LENS ranks, and the tier over the retained alleles only, plus a flag when the best allele is
one the tumour lost.

The arrangement of alleles onto haplotypes is not something a SNV-level phased VCF settles, so it is
declared in the design and justified there by linkage disequilibrium rather than inferred here.
"""


class HlaLoss:
    """Allele-to-haplotype assignment, and the copy number of each allele in each clone."""

    def __init__(self, ds_cfg, clones):
        self.alleles = list(ds_cfg["hla"])
        self.clones = clones
        raw = ds_cfg.get("hla_haplotypes") or {}
        self.by_hap = {int(k): list(v) for k, v in raw.items()}
        self.hap_of = {}
        for hap, alist in self.by_hap.items():
            for a in alist:
                # a homozygous allele sits on both haplotypes; record every haplotype it is on
                self.hap_of.setdefault(a, set()).add(hap)
        self.region = self._region(ds_cfg)

    @staticmethod
    def _region(ds_cfg):
        """The chromosome and span of the truncal CNA that changes HLA copy number, if there is one."""
        for clone, regions in (ds_cfg.get("cna") or {}).items():
            for r in regions:
                if "HLA" not in (r.get("label") or ""):
                    continue
                reg = r["region"]
                chrom, rest = reg.split(":", 1)
                if "-" in rest:
                    lo, hi = (int(x) for x in rest.split("-"))
                else:                       # a whole arm; use its midpoint as the probe position
                    lo = hi = None
                return {"clone": clone, "chrom": chrom, "start": lo, "end": hi,
                        "cn": list(r["cn"]), "label": r["label"], "region": reg}
        return None

    def probe(self):
        """A position inside the HLA region, for asking the clone model about copy number."""
        if not self.region or self.region["start"] is None:
            return None
        return self.region["chrom"], (self.region["start"] + self.region["end"]) // 2

    def allele_copies(self, clone):
        """Copies of each allele in `clone`, from the haplotype copy number at the HLA locus."""
        p = self.probe()
        if p is None:
            return {a: 2 if len(self.hap_of.get(a, ())) == 2 else 1 for a in set(self.alleles)}
        chrom, pos = p
        out = {}
        for a in set(self.alleles):
            haps = self.hap_of.get(a)
            if not haps:                    # arrangement not declared; assume it survives
                out[a] = None
                continue
            cn = self.clones.cn(clone, chrom, pos)[:2]
            out[a] = sum(cn[h] for h in sorted(haps))
        return out

    def retained(self, clone=None):
        """Alleles with at least one copy in `clone` (the clone carrying the loss by default)."""
        clone = clone or (self.region or {}).get("clone") or "T"
        cop = self.allele_copies(clone)
        return sorted(a for a, n in cop.items() if n is None or n > 0)

    def lost(self, clone=None):
        clone = clone or (self.region or {}).get("clone") or "T"
        cop = self.allele_copies(clone)
        return sorted(a for a, n in cop.items() if n == 0)

    def rows(self):
        """One row per allele per clone, for `<dataset>.hla_loh.tsv`."""
        out = []
        for clone in self.clones.clones:
            cop = self.allele_copies(clone)
            for a in sorted(set(self.alleles)):
                out.append({
                    "allele": a, "gene": "HLA-" + a.split("*")[0].replace("HLA-", "")[0],
                    "haplotypes": ",".join(str(h) for h in sorted(self.hap_of.get(a, ()))) or "",
                    "clone": clone, "ccf": self.clones.ccf(clone),
                    "copies": "" if cop[a] is None else cop[a],
                    "lost": cop[a] == 0,
                    "loh_label": (self.region or {}).get("label", ""),
                    "loh_region": (self.region or {}).get("region", ""),
                })
        return out
