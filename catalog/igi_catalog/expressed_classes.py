"""Cancer-testis antigen, ERV, and splice-variant designers.

These three classes are defined by expression rather than by a change in the genome, so each event says
what the transcriptome builder must produce and what a caller should therefore be able to see.

Tumour-specific versus tumour-associated is the distinction that matters for ERVs and splice variants:
  tumour-specific  - the locus or junction is absent from the normal reference panel, and for splice
                     variants is caused by a designed somatic variant at the splice site;
  tumour-associated - present in normals at low level, raised in the tumour, with no somatic cause.
"""
import csv
import gzip

from .annotation import translate
from .binding import NetMHCpan, spanning_peptides
from .expression import Expression
from .genome import revcomp

CHR1TO6 = {f"chr{i}" for i in range(1, 7)}
SPLICE_MECHANISMS = ["exon_skip", "intron_retention", "cryptic_5p", "cryptic_3p", "novel_exon"]


def load_gene_list(path):
    with open(path) as fh:
        return [l.strip() for l in fh if l.strip() and not l.startswith("#")]


def load_pan_normal(quant_sf):
    """Name -> normal-tissue 95th-percentile TPM from the LENS pan-normal reference.

    The table carries transcripts and ERV loci alike (31,292 `Hsap38.*` entries), so one lookup serves
    cancer-testis antigens and ERVs both.
    """
    out = {}
    op = gzip.open if str(quant_sf).endswith(".gz") else open
    with op(quant_sf, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        try:
            i_name, i_tpm = hdr.index("Name"), hdr.index("TPM")
        except ValueError:
            return out
        for line in fh:
            f = line.rstrip("\n").split("\t")
            name = f[i_name].split("|")[0]
            tpm = float(f[i_tpm])
            out[name] = tpm
            # transcripts carry a version suffix that callers strip; ERV locus ids are dotted throughout
            # and must not be truncated, so only ENST/ENSG-style names get an alias
            if name.startswith(("ENST", "ENSG")) and "." in name:
                out.setdefault(name.split(".")[0], tpm)
    return out


class ExpressedClassDesigner:
    def __init__(self, env, ds_name, ds_cfg, cfg, rng, next_id, paths):
        self.env = env
        self.ds = ds_name
        self.dcfg = ds_cfg
        self.cfg = cfg
        self.rng = rng
        self.next_id = next_id
        self.paths = paths
        self.events = []
        self.by_name = {}
        for t in env.rep.values():
            self.by_name.setdefault(t.gene_name, t)
        self.pan_normal = load_pan_normal(paths["pan_normal_quant"]) if paths.get("pan_normal_quant") else {}

    # ------------------------------------------------------------------ shared
    def _score_window(self, protein, first_changed, n_changed):
        peps = spanning_peptides(protein, first_changed, n_changed, (8, 9, 10, 11))
        if not peps:
            return None, None, None, "na", 0
        res = self.env.netmhc.predict(sorted(peps), self.env.hla)
        best = None
        n_lt2 = 0
        for pep, al in res.items():
            r, allele, aff = NetMHCpan.best(al)
            if r <= 2.0:
                n_lt2 += 1
            if best is None or r < best[0]:
                best = (r, allele, aff, pep)
        if not best:
            return None, None, None, "na", 0
        return round(best[0], 3), best[1], best[3], NetMHCpan.tier(best[0]), n_lt2

    @staticmethod
    def longest_orf(seq, min_aa=8):
        """Longest methionine-started ORF over the three forward frames of `seq`."""
        best = ""
        for frame in range(3):
            prot = translate(seq[frame:])
            for piece in prot.split("*"):
                i = piece.find("M")
                if i >= 0 and len(piece) - i > len(best):
                    best = piece[i:]
        return best if len(best) >= min_aa else ""

    def _clone_fields(self, clone, chrom, pos):
        cm = self.env.clones
        hap = self.rng.choice(cm.retained_haplotypes(clone, chrom, pos) or [0])
        return {"clone": clone, "ccf": cm.ccf(clone), "haplotype": hap,
                "clonality_tier": cm.clonality_tier(clone, chrom, pos, hap)}

    # ------------------------------------------------------------------ CTAs
    def design_ctas(self, n, clones, log=print):
        names = load_gene_list(self.paths["cta_gene_list"])
        self.rng.shuffle(names)
        tiers = ["T0", "T1", "T10", "T100", "T1000"]
        made = 0
        for name in names:
            if made >= n:
                break
            t = self.by_name.get(name)
            if t is None or not t.cds:
                continue
            from .annotation import CodingModel
            cm = CodingModel(t, self.env.genome)
            if len(cm.protein) < 50 or "*" in cm.protein:
                continue
            target = tiers[made % len(tiers)]
            clone = "T" if made % 3 else self.rng.choice(clones)
            mid = len(cm.protein) // 2
            rank, allele, pep, tier, n_lt2 = self._score_window(cm.protein, mid, 1)
            ev = {
                "event_id": self.next_id("CTA"), "dataset": self.ds, "class": "cta",
                "subclass": "clone_restricted" if clone != "T" else "tumor_wide",
                "gene": name, "transcript": t.tid, "chrom": t.chrom, "start": t.start, "end": t.end,
                "target_expression_tier": target,
                "normal_tissue_tpm_p95": round(self.pan_normal.get(t.tid.split(".")[0], 0.0), 3),
                "baseline_tumor_tpm": round(self.env.expr.gene(t.gene_id), 3),
                "best_rank_el": rank, "best_allele": allele, "best_peptide": pep,
                "binding_tier": tier, "n_peptides_rank_lt2": n_lt2,
                "chr1to6": t.chrom in CHR1TO6,
            }
            ev.update(self._clone_fields(clone, t.chrom, t.start))
            self.events.append(ev)
            made += 1
        log(f"  CTAs: {made}/{n}")
        return made

    # ------------------------------------------------------------------ ERVs
    def design_ervs(self, n, clones, log=print):
        loci = []
        with open(self.paths["erv_loci"]) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("chr") and row.get("start"):
                    loci.append(row)
        extra = self.paths.get("erv_all_loci")
        if extra and len(loci) < n * 4:
            with open(extra) as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    if row.get("chr") and row.get("start"):
                        loci.append(row)
        self.rng.shuffle(loci)
        n_spec = n // 3
        n_assoc = n // 3
        n_spec_made = n_assoc_made = 0
        made = 0
        for row in loci:
            if made >= n:
                break
            chrom = row["chr"] if row["chr"].startswith("chr") else f"chr{row['chr']}"
            if chrom not in self.env.genome.lengths:
                continue
            start, end = int(row["start"]), int(row["end"])
            prot = row.get("protein") or ""
            # status follows the measured pan-normal level of this locus, not the loop index
            normal_tpm = self.pan_normal.get(row.get("ID", ""))
            if normal_tpm is None:
                continue                      # no measurement: cannot classify it honestly
            normal_tpm = round(normal_tpm, 3)
            if normal_tpm <= 0.05:
                status = "tumor_specific" if n_spec_made < n_spec else "unexpressed_negative"
                if status == "tumor_specific":
                    n_spec_made += 1
            elif normal_tpm <= 5.0:
                if n_assoc_made >= n_assoc:
                    continue
                status = "tumor_associated"
                n_assoc_made += 1
            else:
                continue                      # too highly expressed in normals to serve any tier
            rank = allele = pep = None
            tier, n_lt2 = "na", 0
            orf_source = "annotated"
            if not prot or len(prot) < 9:
                seq = self.env.genome.seq(chrom, start - 1, end)
                if row.get("strand") == "-":
                    seq = revcomp(seq)
                prot = self.longest_orf(seq)
                orf_source = "derived_longest_orf" if prot else "none"
            if prot and len(prot) >= 9:
                rank, allele, pep, tier, n_lt2 = self._score_window(prot, len(prot) // 2, 1)
            clone = "T" if made % 4 else self.rng.choice(clones)
            ev = {
                "event_id": self.next_id("ERV"), "dataset": self.ds, "class": "erv",
                "subclass": status, "locus_id": row.get("ID", ""), "chrom": chrom,
                "start": start, "end": end, "strand": row.get("strand", "+"),
                "repeat_family": row.get("Repbase", ""), "orf_aa_length": row.get("AA_length", ""),
                "has_annotated_orf": orf_source == "annotated", "orf_source": orf_source,
                "orf_aa": len(prot) if prot else 0,
                "normal_panel_tpm": normal_tpm, "normal_panel_source": "pan_normal_95th_percentile",
                "target_tumor_tpm": 0.0 if status == "unexpressed_negative" else round(self.rng.uniform(5, 120), 2),
                "best_rank_el": rank, "best_allele": allele, "best_peptide": pep,
                "binding_tier": tier, "n_peptides_rank_lt2": n_lt2,
                "chr1to6": chrom in CHR1TO6,
            }
            ev.update(self._clone_fields(clone, chrom, start))
            self.events.append(ev)
            made += 1
        log(f"  ERVs: {made}/{n}")
        return made

    # ------------------------------------------------------------------ splice
    def _splice_peptides(self, t, exons):
        """Protein of a novel isoform and the residue where it first differs from the reference.

        The isoform is translated from its longest ORF, which is what a caller working from assembled
        RNA would see; the comparison point is the first residue that differs from the reference protein.
        """
        from .annotation import CodingModel
        from .transcriptome import exon_sequence
        stub = type(t)(tid=t.tid + ".iso", gene_id=t.gene_id, gene_name=t.gene_name, gene_type=t.gene_type,
                       chrom=t.chrom, strand=t.strand, start=min(s for s, _ in exons),
                       end=max(e for _, e in exons), exons=exons, cds=[], tags=[], level=t.level)
        iso = self.longest_orf(exon_sequence(self.env.genome, stub))
        if not iso:
            return None, None
        ref = CodingModel(t, self.env.genome).protein
        k = 0
        while k < min(len(iso), len(ref)) and iso[k] == ref[k]:
            k += 1
        if k >= len(iso):
            return iso, None
        return iso, k

    def design_splice(self, n, clones, snv_designer, log=print):
        """Half the events are caused by a designed somatic splice-site variant (tumour-specific);
        half are isoform switches with no genomic cause (tumour-associated)."""
        n_specific = n // 2
        made = 0
        genes = [t for t in self.env.rep.values() if len(t.cds) >= 4 and self.env.expr.gene(t.gene_id) > 5]
        self.rng.shuffle(genes)
        for t in genes:
            if made >= n:
                break
            specific = made < n_specific
            mech = SPLICE_MECHANISMS[made % len(SPLICE_MECHANISMS)]
            from .transcriptome import skip_exon, retain_intron, cryptic_site, novel_exon
            feasible = {"exon_skip": lambda i: skip_exon(t, i), "intron_retention": lambda i: retain_intron(t, i),
                        "cryptic_5p": lambda i: cryptic_site(t, i, "5p"), "cryptic_3p": lambda i: cryptic_site(t, i, "3p"),
                        "novel_exon": lambda i: novel_exon(t, i)}[mech]
            # index in EXON space, not CDS space: the two differ whenever a transcript has untranslated
            # exons, and modifying a non-coding exon leaves the protein unchanged
            exons = t.exons if t.strand == "+" else t.exons[::-1]
            cds_lo = min(a for a, _b in t.cds)
            cds_hi = max(b for _a, b in t.cds)
            cand = [j for j in range(2, len(exons))
                    if exons[j - 1][1] >= cds_lo and exons[j - 1][0] <= cds_hi and feasible(j)]
            if not cand:
                continue
            i = self.rng.choice(cand) - 1
            exon_s, exon_e = exons[i]
            clone = "T" if made % 3 else self.rng.choice(clones)
            causal = None
            if specific:
                # a somatic base change at the canonical donor or acceptor dinucleotide of this exon
                donor = t.strand == "+"
                site = (exon_e + 1) if donor else (exon_s - 2)
                ref = self.env.genome.seq(t.chrom, site - 1, site + 1)
                if len(ref) < 2 or self.env.germline.overlaps_variant(t.chrom, site, 2):
                    continue
                # pick a base that leaves a dinucleotide no spliceosome recognises. "A" at a GT donor
                # would give AT, which is the U12 donor, so the change would not unambiguously break
                # the site. Donors on a minus-strand gene read as their reverse complement here.
                recognised = {"GT", "GC", "AT", "AC", "AG", "CT"}
                alt = next((b for b in "CATG" if b != ref[0] and (b + ref[1]) not in recognised), None)
                if alt is None:
                    continue
                causal = {"chrom": t.chrom, "pos": site, "ref": ref[0], "alt": alt}
            iso_prot, first_diff = self._splice_peptides(t, feasible(i + 1) or [])
            if iso_prot is None or first_diff is None:
                continue          # the isoform does not change the protein: not a usable antigen source
            rank = allele = pep = None
            tier, n_lt2 = "na", 0
            if iso_prot and first_diff is not None:
                n_new = max(1, len(iso_prot) - first_diff)
                rank, allele, pep, tier, n_lt2 = self._score_window(iso_prot, first_diff, min(n_new, 25))
            ev = {
                "event_id": self.next_id("SPL"), "dataset": self.ds, "class": "splice",
                "subclass": "tumor_specific" if specific else "tumor_associated",
                "mechanism": mech, "gene": t.gene_name, "transcript": t.tid, "chrom": t.chrom,
                "exon_index": i + 1, "exon_start": exon_s, "exon_end": exon_e,
                "causal_variant": (f"{causal['chrom']}:{causal['pos']}{causal['ref']}>{causal['alt']}" if causal else ""),
                # A junction created by a somatic splice-site change cannot be present in normals. An
                # isoform switch is an existing minor isoform, so its normal level is a small share of the
                # gene's own normal expression rather than an invented number.
                "normal_panel_junction_tpm": 0.0 if specific else round(
                    0.08 * self.pan_normal.get(t.tid.split(".")[0], 0.0), 3),
                "normal_panel_source": "none_by_construction" if specific else "0.08x gene pan-normal TPM",
                "target_tumor_junction_tpm": round(self.rng.uniform(5, 80), 2),
                "gene_tpm": round(self.env.expr.gene(t.gene_id), 3),
                "expression_tier": Expression.tier(self.env.expr.gene(t.gene_id)),
                "isoform_orf_aa": len(iso_prot) if iso_prot else 0,
                "first_changed_aa": (first_diff + 1) if first_diff is not None else "",
                "best_rank_el": rank, "best_allele": allele, "best_peptide": pep,
                "binding_tier": tier, "n_peptides_rank_lt2": n_lt2,
                "chr1to6": t.chrom in CHR1TO6,
            }
            ev.update(self._clone_fields(clone, t.chrom, exon_s))
            self.events.append(ev)
            made += 1
        log(f"  splice variants: {made}/{n} ({n_specific} tumour-specific)")
        return made
