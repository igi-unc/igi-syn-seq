"""Background (passenger) somatic mutations carrying the dataset's mutational signature.

The designed catalog supplies the events a neoantigen caller is meant to find. On its own it gives a
tumour genome whose only somatic differences are the ones under test, which is not what a caller sees and
makes any false-positive rate measured against it meaningless. This module adds the passenger background:
mutations drawn from a COSMIC signature mixture, placed across the callable genome, assigned to clones and
haplotypes by the same rules as the designed events.

Placement follows the signature's own convention. A COSMIC profile is a distribution over the 96
pyrimidine-centred trinucleotide channels of *observed* mutations, so a channel is drawn first and a
position with that context second; the resulting spectrum reproduces the profile without needing a
context-frequency correction.

Background mutations that land in coding sequence are annotated and scored like designed ones. A passenger
missense variant in a well expressed gene is a real neoepitope, and a benchmark that did not record it
would count a caller correct answer as a false positive.
"""
import bisect
import csv
import random

from .annotation import CodingModel
from .binding import NetMHCpan, spanning_peptides
from .expression import Expression

CHR1TO6 = {f"chr{i}" for i in range(1, 7)}
COMP = str.maketrans("ACGT", "TGCA")
# consequences whose mutant protein differs from the wild type, and so can yield a neoepitope
PROTEIN_CHANGING = {"missense", "nonsense", "stop_loss", "start_loss", "frameshift",
                    "inframe_insertion", "inframe_deletion"}


def revcomp(s):
    return s.translate(COMP)[::-1]


def channel(trinuc, ref, alt):
    """The COSMIC 96-channel label for a plus-strand substitution, folded onto the pyrimidine strand.

    Events are recorded on the plus strand, because that is what a VCF has to say; the signature
    convention is pyrimidine-centred, so the reported spectrum folds purine-centred calls back.
    """
    if ref in "AG":
        trinuc, ref, alt = revcomp(trinuc), revcomp(ref), revcomp(alt)
    return f"{trinuc[0]}[{ref}>{alt}]{trinuc[2]}"


class SignatureMix:
    """A weighted mixture of COSMIC SBS profiles, sampled as (context, ref, alt) triples."""

    def __init__(self, tsv, weights):
        rows = []
        with open(tsv) as fh:
            reader = csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")
            missing = [s for s in weights if s not in reader.fieldnames]
            if missing:
                raise KeyError(f"{tsv} has no column(s) {missing}; available: {reader.fieldnames}")
            for r in reader:
                rows.append(r)
        if len(rows) != 96:
            raise ValueError(f"{tsv}: expected 96 channels, found {len(rows)}")
        total_w = sum(weights.values())
        self.channels, probs = [], []
        for r in rows:
            ch = r["channel"]                     # e.g. A[C>A]A
            ctx = ch[0] + ch[2] + ch[6]
            self.channels.append((ctx, ch[2], ch[4]))
            probs.append(sum(w * float(r[s]) for s, w in weights.items()) / total_w)
        s = sum(probs)
        self.cum = []
        acc = 0.0
        for p in probs:
            acc += p / s
            self.cum.append(acc)
        self.cum[-1] = 1.0
        self.weights = dict(weights)

    def draw(self, rng):
        """(context, ref, alt) with the context written 5'->3' on the pyrimidine strand."""
        return self.channels[bisect.bisect_left(self.cum, rng.random())]


class CdsIndex:
    """Position lookup into representative coding transcripts, one sorted segment list per chromosome."""

    def __init__(self, rep):
        self.by_chrom = {}
        for t in rep.values():
            if not t.cds:
                continue
            self.by_chrom.setdefault(t.chrom, []).append((min(s for s, _ in t.cds),
                                                          max(e for _, e in t.cds), t))
        for chrom, items in self.by_chrom.items():
            items.sort(key=lambda x: x[0])
            self.by_chrom[chrom] = (items, [i[0] for i in items])

    def transcript_at(self, chrom, pos1):
        """A representative transcript whose CDS covers pos1, or None."""
        entry = self.by_chrom.get(chrom)
        if not entry:
            return None
        items, starts = entry
        i = bisect.bisect_right(starts, pos1)
        for lo, hi, t in reversed(items[max(0, i - 40):i]):
            if lo <= pos1 <= hi and any(s <= pos1 <= e for s, e in t.cds):
                return t
        return None


class BackgroundDesigner:
    def __init__(self, env, ds_name, ds_cfg, cfg, rng, sig_tsv):
        self.env = env
        self.ds = ds_name
        self.dcfg = ds_cfg
        self.cfg = cfg
        self.rng = rng
        self.sig = SignatureMix(sig_tsv, ds_cfg["signatures"])
        self.cds = CdsIndex(env.rep)
        self.events = []
        self.n = 0
        self._reserved = set()

    # ------------------------------------------------------------------ setup
    def reserve_designed(self, tsvs):
        """Reference spans already used by designed events, so a passenger cannot land on one."""
        for tsv in tsvs:
            try:
                fh = open(tsv)
            except FileNotFoundError:
                continue
            with fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    if row.get("class") not in ("snv", "indel"):
                        continue
                    pos, ref = int(row["pos"]), row["ref"]
                    for p in range(pos - 1, pos + len(ref) + 1):
                        self._reserved.add((row["chrom"], p))
        return len(self._reserved)

    def clone_weights(self, clones, truncal_fraction):
        """Share of the passenger burden per clone: truncal mutations plus a per-branch share."""
        sub = [c for c in clones if c != "T"]
        w = {"T": truncal_fraction}
        tot = sum(self.env.clones.ccf(c) for c in sub) or 1.0
        for c in sub:
            w[c] = (1.0 - truncal_fraction) * self.env.clones.ccf(c) / tot
        return w

    # ------------------------------------------------------------------ placement
    def callable_bases(self, chrom):
        """Non-N bases on a chromosome."""
        seq = self.env.genome.seq(chrom, 0, self.env.genome.lengths[chrom]).upper()
        return len(seq) - seq.count("N")

    def _site_for(self, chrom, ctx, ref, alt, tries=4000):
        """A position whose trinucleotide matches the channel, on either strand."""
        L = self.env.genome.lengths[chrom]
        rc_ctx, rc_ref, rc_alt = revcomp(ctx), revcomp(ref), revcomp(alt)
        for _ in range(tries):
            pos = self.rng.randrange(2, L)          # 1-based; needs one base either side
            tri = self.env.genome.seq(chrom, pos - 2, pos + 1).upper()
            if tri == ctx:
                r, a = ref, alt
            elif tri == rc_ctx:
                r, a = rc_ref, rc_alt
            else:
                continue
            if (chrom, pos) in self._reserved:
                continue
            if self.env.germline.overlaps_variant(chrom, pos, 1):
                continue
            return pos, r, a
        return None

    def _indel(self, chrom, tries=4000):
        """A short passenger indel: length 1-10, deletion or insertion, outside designed and germline spans."""
        L = self.env.genome.lengths[chrom]
        size = min(10, 1 + int(self.rng.expovariate(0.6)))
        is_del = self.rng.random() < 0.6
        for _ in range(tries):
            pos = self.rng.randrange(2, L - size - 2)
            span = self.env.genome.seq(chrom, pos - 1, pos + size + 1).upper()
            if "N" in span:
                continue
            if any((chrom, p) in self._reserved for p in range(pos - 1, pos + size + 2)):
                continue
            if self.env.germline.overlaps_variant(chrom, pos, size + 1):
                continue
            if is_del:
                return pos, span[:size + 1], span[0]
            ins = "".join(self.rng.choice("ACGT") for _ in range(size))
            return pos, span[0], span[0] + ins
        return None

    # ------------------------------------------------------------------ annotation
    def _coding(self, chrom, pos, ref, alt):
        """Consequence and, for protein-changing variants, the mutant protein and first changed residue."""
        t = self.cds.transcript_at(chrom, pos)
        if t is None:
            return None
        cm = CodingModel(t, self.env.genome)
        mprot, k, cons = cm.mutate_protein(pos, ref, alt)
        if cons in ("non_coding", "ref_mismatch", "frame_unresolved") or mprot is None:
            return None
        return {"gene": t.gene_name, "transcript": t.tid, "consequence": cons,
                "wt_protein": cm.protein, "mut_protein": mprot, "aa_index": k, "t": t}

    # ------------------------------------------------------------------ top level
    def design(self, chroms, clones, indel_fraction=0.08, truncal_fraction=0.60, log=print):
        rate = float(self.dcfg["background_mut_per_mb"])
        cw = self.clone_weights(clones, truncal_fraction)
        names, cums = list(cw), []
        acc = 0.0
        for c in names:
            acc += cw[c]
            cums.append(acc)
        cums[-1] = 1.0
        cm_model = self.env.clones
        coding_hits = []
        for chrom in chroms:
            mb = self.callable_bases(chrom) / 1e6
            want = int(round(rate * mb))
            made = failed = 0
            while made < want:
                if self.rng.random() < indel_fraction:
                    hit = self._indel(chrom)
                    kind = "indel"
                else:
                    ctx, ref, alt = self.sig.draw(self.rng)
                    hit = self._site_for(chrom, ctx, ref, alt)
                    kind = "snv"
                if hit is None:
                    failed += 1
                    if failed > 200:
                        break
                    continue
                pos, r, a = hit
                for p in range(pos - 1, pos + len(r) + 1):
                    self._reserved.add((chrom, p))
                clone = names[bisect.bisect_left(cums, self.rng.random())]
                haps = cm_model.retained_haplotypes(clone, chrom, pos) or [0]
                hap = self.rng.choice(haps)
                pre = clone == "T"
                vaf, _m, tcn = cm_model.vaf(chrom, pos, hap, clone, pre_cna=pre)
                self.n += 1
                ev = {
                    "event_id": f"BG-{self.n:06d}",
                    "dataset": self.ds, "class": kind, "subclass": "background",
                    "chrom": chrom, "pos": pos, "ref": r, "alt": a,
                    "trinuc": (self.env.genome.seq(chrom, pos - 2, pos + 1).upper() if kind == "snv" else ""),
                    "gene": "", "transcript": "", "consequence": "non_coding",
                    "haplotype": hap, "clone": clone, "ccf": cm_model.ccf(clone),
                    "timing": "pre_cna" if pre else "post_cna",
                    "clonality_tier": cm_model.clonality_tier(clone, chrom, pos, hap),
                    "expected_vaf_dna": round(vaf, 4), "tumor_cn_at_locus": round(tcn, 2),
                    "gene_tpm": "", "expression_tier": "", "binding_tier": "na",
                    "best_rank_el": "", "best_allele": "", "best_peptide": "",
                    "binding_tier_retained": "na", "best_rank_el_retained": "",
                    "best_allele_retained": "", "best_peptide_retained": "",
                    "best_allele_is_lost": False,
                    "flagpost": False, "chr1to6": chrom in CHR1TO6,
                }
                c = self._coding(chrom, pos, r, a)
                if c:
                    tpm = self.env.expr.gene(c["t"].gene_id)
                    ev.update({"gene": c["gene"], "transcript": c["transcript"],
                               "consequence": c["consequence"], "gene_tpm": round(tpm, 3),
                               "expression_tier": Expression.tier(tpm)})
                    if c["consequence"] in PROTEIN_CHANGING:
                        coding_hits.append((ev, c))
                self.events.append(ev)
                made += 1
            log(f"  {chrom}: {made} background mutations over {mb:.1f} callable Mb"
                + (f" ({failed} placement failures)" if failed else ""))
        self._score(coding_hits, log)
        return self.events

    def _score(self, hits, log):
        """Binding for passengers that change a protein, so real passenger neoepitopes are in the truth."""
        if not hits:
            return
        want = {}
        for ev, c in hits:
            if c["consequence"] in ("missense", "inframe_insertion", "inframe_deletion", "start_loss"):
                n = max(1, len(c["mut_protein"]) - c["aa_index"]) if c["consequence"] != "missense" else 1
            else:
                n = max(1, len(c["mut_protein"]) - c["aa_index"])   # frameshift / readthrough tail
            peps = spanning_peptides(c["mut_protein"], c["aa_index"], n, (8, 9, 10, 11))
            wt = spanning_peptides(c["wt_protein"], c["aa_index"], 1, (8, 9, 10, 11)) if c["wt_protein"] else set()
            want[ev["event_id"]] = sorted(peps - wt)
        allp = sorted({p for ps in want.values() for p in ps})
        log(f"  scoring {len(allp)} passenger neopeptides over {len(hits)} protein-changing background events")
        res = self.env.netmhc.predict(allp, self.env.hla)
        for ev, _c in hits:
            peps = want[ev["event_id"]]
            r, allele, _aff, pep = NetMHCpan.best_over(res, peps)
            if r is not None:
                ev.update({"binding_tier": NetMHCpan.tier(r), "best_rank_el": round(r, 3),
                           "best_allele": allele, "best_peptide": pep})
            # owner decision D4: the same peptides over the alleles the tumour retains
            rr, ral, _ra, rpep = NetMHCpan.best_over(res, peps, among=self.env.hla_retained)
            ev.update({"binding_tier_retained": NetMHCpan.tier(rr),
                       "best_rank_el_retained": round(rr, 3) if rr is not None else None,
                       "best_allele_retained": ral, "best_peptide_retained": rpep,
                       "best_allele_is_lost": bool(allele and allele in self.env.hla_lost)})
