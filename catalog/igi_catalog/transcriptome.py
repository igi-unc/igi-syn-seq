"""Stage 3: build the per-clone transcriptome that read simulation samples from.

The output is a set of transcript records, each with a sequence and a TPM per clone:

  * reference transcripts, rebuilt on each haplotype so germline and somatic small variants appear in RNA;
  * fusion transcripts, from the designed exon junctions;
  * novel isoforms for designed splice events (exon skip, intron retention, cryptic sites, novel exon);
  * ERV transcripts, whose loci are outside the GENCODE model;
  * viral transcripts for expressed integrations and episomes.

Expression starts from the cohort baseline, is overridden by any designed event on that gene, and is then
split between haplotypes according to the event's allelic-expression setting.
"""
import csv
from collections import defaultdict

from .genome import revcomp


class TranscriptRecord(dict):
    """id, source, gene, chrom, sequence, and tpm per clone."""


def exon_sequence(genome, t, edits=None):
    """Spliced transcript sequence in transcript orientation, with optional per-haplotype edits applied."""
    parts = []
    for s, e in t.exons:
        seq = genome.seq(t.chrom, s - 1, e)
        if edits is not None:
            seq, _n, _rej = edits.apply(seq, s, strict=False)
        parts.append(seq)
    full = "".join(parts)
    return full if t.strand == "+" else revcomp(full)


def skip_exon(t, index):
    """Exon list with one internal exon removed (1-based index in transcript order)."""
    ex = t.exons if t.strand == "+" else t.exons[::-1]
    if index <= 1 or index >= len(ex):
        return None
    kept = [tuple(e) for i, e in enumerate(ex, start=1) if i != index]
    return sorted(kept)


def retain_intron(t, index):
    """Exon list with the intron following exon `index` retained (exons merged)."""
    ex = t.exons if t.strand == "+" else t.exons[::-1]
    if index < 1 or index >= len(ex):
        return None
    a, b = sorted([tuple(ex[index - 1]), tuple(ex[index])])
    merged = (a[0], b[1])
    kept = [tuple(e) for i, e in enumerate(ex, start=1) if i not in (index, index + 1)] + [merged]
    return sorted(kept)


def cryptic_site(t, index, side, shift=30, min_exon=30):
    """Exon list with one exon boundary moved, modelling a cryptic donor or acceptor.

    The shift is reduced to fit short exons; an exon too short to move is refused rather than inverted.
    """
    ex = [tuple(e) for e in sorted(t.exons)]
    order = list(range(len(ex))) if t.strand == "+" else list(reversed(range(len(ex))))
    if index < 1 or index > len(ex):
        return None
    i = order[index - 1]
    s, e = ex[i]
    if e - s + 1 < min_exon:
        return None
    shift = max(3, min(shift, (e - s + 1) // 2 - 5))
    if side == "5p":
        e = e - shift if t.strand == "+" else e
        s = s + shift if t.strand == "-" else s
    else:
        s = s + shift if t.strand == "+" else s
        e = e - shift if t.strand == "-" else e
    if e - s < 20:
        return None
    out = list(ex)
    out[i] = (s, e)
    return out


def novel_exon(t, index, length=90):
    """Exon list with a new exon inserted into an intron, searching outward from `index` for one that
    is long enough to hold the new exon plus flanking intron."""
    ex = [tuple(e) for e in sorted(t.exons)]
    if len(ex) < 2:
        return None
    need = length + 200
    order = sorted(range(1, len(ex)), key=lambda i: abs(i - index))
    for i in order:
        gap_s, gap_e = ex[i - 1][1] + 1, ex[i][0] - 1
        if gap_e - gap_s >= need:
            mid = (gap_s + gap_e) // 2
            return sorted(ex + [(mid, mid + length - 1)])
    return None


class TranscriptomeBuilder:
    def __init__(self, env, ds_name, clones, expression):
        self.env = env
        self.ds = ds_name
        self.clones = clones
        self.expr = expression
        self.records = []

    # ------------------------------------------------------------------ expression
    def _clone_tpm(self, base_tpm, event=None):
        """TPM per clone: the baseline everywhere, replaced by a target in clones carrying the event."""
        out = {c: base_tpm for c in self.clones.clones}
        if event is None:
            return out
        target = event.get("target_tpm")
        if target is None:
            return out
        origin = event.get("clone", "T")
        for c in self.clones.clones:
            if self.clones.is_descendant(c, origin):
                out[c] = float(target)
        return out

    # ------------------------------------------------------------------ builders
    def add_reference_transcripts(self, chrom, hap, edits_by_clone, limit=None):
        """Reference transcripts on one chromosome and haplotype, one sequence per clone edit set."""
        n = 0
        for t in self.env.rep.values():
            if t.chrom != chrom:
                continue
            base = self.expr.transcript(t.tid, t.gene_id)
            for clone, edits in edits_by_clone.items():
                seq = exon_sequence(self.env.genome, t, edits)
                self.records.append(TranscriptRecord(
                    id=f"{t.tid}|hap{hap}|{clone}", source="reference", gene=t.gene_name,
                    chrom=chrom, hap=hap, clone=clone, sequence=seq, tpm=base))
            n += 1
            if limit and n >= limit:
                break
        return n

    def add_splice_isoform(self, event):
        """One novel isoform for a designed splice event."""
        t = self.env.tx.get(event["transcript"])
        if t is None:
            return None
        idx = int(event["exon_index"])
        mech = event["mechanism"]
        exons = {
            "exon_skip": lambda: skip_exon(t, idx),
            "intron_retention": lambda: retain_intron(t, idx),
            "cryptic_5p": lambda: cryptic_site(t, idx, "5p", 30),
            "cryptic_3p": lambda: cryptic_site(t, idx, "3p", 30),
            "novel_exon": lambda: novel_exon(t, idx),
        }[mech]()
        if not exons:
            return None
        stub = type(t)(tid=f"{t.tid}.{mech}", gene_id=t.gene_id, gene_name=t.gene_name,
                       gene_type=t.gene_type, chrom=t.chrom, strand=t.strand,
                       start=min(s for s, _ in exons), end=max(e for _, e in exons),
                       exons=exons, cds=[], tags=[], level=t.level)
        seq = exon_sequence(self.env.genome, stub)
        rec = TranscriptRecord(
            id=f"{event['event_id']}|{t.tid}|{mech}", source="splice_isoform", gene=t.gene_name,
            chrom=t.chrom, hap=int(event.get("haplotype", 0)), clone=event.get("clone", "T"),
            sequence=seq, tpm=float(event.get("target_tumor_junction_tpm", 0) or 0),
            # an isoform switch moves expression within the gene rather than adding to it
            reconcile="redistribute", reconcile_gene=t.gene_name)
        self.records.append(rec)
        return rec

    def add_fusion_transcript(self, event, fusion_core):
        """Transcript for a designed fusion, from its exon junction."""
        t5 = self.env.tx.get(event["transcript_5p"])
        t3 = self.env.tx.get(event["transcript_3p"])
        if t5 is None or t3 is None:
            return None
        n5, from3 = int(event["exons_5p"]), int(event["exon_start_3p"]) - 1
        seq = (fusion_core.five_prime_cds(self.env.genome, t5, n5)
               + fusion_core.three_prime_cds(self.env.genome, t3, from3))
        tpm = min(float(event.get("gene_tpm_5p", 0) or 0), float(event.get("gene_tpm_3p", 0) or 0))
        rec = TranscriptRecord(
            id=f"{event['event_id']}|{event['gene_5p']}-{event['gene_3p']}", source="fusion",
            gene=f"{event['gene_5p']}-{event['gene_3p']}", chrom=event["chrom_5p"],
            hap=int(event.get("haplotype", 0)), clone=event.get("clone", "T"),
            sequence=seq, tpm=max(1.0, tpm * 0.3),
            # the fusion transcript is transcribed from the 5' partner's promoter, so its share comes
            # off that gene rather than being added on top of full wild-type expression
            reconcile="redistribute", reconcile_gene=event["gene_5p"])
        self.records.append(rec)
        return rec

    def add_erv_transcript(self, event):
        """Transcript for a designed ERV locus, which is outside the GENCODE model."""
        chrom, s, e = event["chrom"], int(event["start"]), int(event["end"])
        if chrom not in self.env.genome.lengths:
            return None
        seq = self.env.genome.seq(chrom, s - 1, e)
        if event.get("strand") == "-":
            seq = revcomp(seq)
        rec = TranscriptRecord(
            id=f"{event['event_id']}|{event.get('locus_id') or chrom + ':' + str(s)}", source="erv",
            gene=event.get("locus_id", ""), chrom=chrom, hap=int(event.get("haplotype", 0)),
            clone=event.get("clone", "T"), sequence=seq,
            tpm=float(event.get("target_tumor_tpm", 0) or 0),
            reconcile="add")      # a locus in its own right, not one of a gene's transcripts
        self.records.append(rec)
        return rec

    # ------------------------------------------------------------------ output
    def write_fasta(self, path):
        with open(path, "w") as fh:
            for r in self.records:
                if not r["sequence"]:
                    continue
                fh.write(f">{r['id']} source={r['source']} gene={r['gene']} clone={r['clone']} hap={r['hap']} tpm={r['tpm']}\n")
                seq = r["sequence"]
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i + 60] + "\n")
        return path

    def write_table(self, path):
        cols = ["id", "source", "gene", "chrom", "hap", "clone", "length", "tpm"]
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for r in self.records:
                row = dict(r)
                row["length"] = len(r["sequence"])
                w.writerow(row)
        return path


def read_table(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))
