"""GENCODE GTF parsing to a compact transcript model, with CDS/protein extraction and coordinate mapping."""
import gzip, json, os, re
from dataclasses import dataclass, field, asdict
from .genome import revcomp

CODON = {a + b + c: aa for (a, b, c), aa in zip(
    [(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
def translate(cds):
    return "".join(CODON.get(cds[i:i + 3], "X") for i in range(0, len(cds) - len(cds) % 3, 3))

@dataclass
class Transcript:
    tid: str; gene_id: str; gene_name: str; gene_type: str; chrom: str; strand: str
    start: int; end: int  # 1-based inclusive transcript span
    exons: list = field(default_factory=list)   # [(start,end)] 1-based inclusive, genomic order
    cds: list = field(default_factory=list)     # [(start,end)] 1-based inclusive, genomic order (no stop codon in GENCODE CDS)
    tags: list = field(default_factory=list)
    level: int = 0
    def cds_len(self): return sum(e - s + 1 for s, e in self.cds)
    def n_exons(self): return len(self.exons)

_ATTR = re.compile(r'(\S+) "([^"]*)"')
def parse_gtf(path, keep_types=("protein_coding", "lncRNA")):
    """Return {tid: Transcript} for genes of keep_types (None = all)."""
    op = gzip.open if path.endswith(".gz") else open
    tx = {}
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if f[2] not in ("transcript", "exon", "CDS"): continue
            a = dict(_ATTR.findall(f[8]))
            gt = a.get("gene_type", "")
            if keep_types and gt not in keep_types: continue
            tid = a["transcript_id"]
            if f[2] == "transcript":
                tags = re.findall(r'tag "([^"]+)"', f[8])
                tx[tid] = Transcript(tid, a["gene_id"], a.get("gene_name", ""), gt, f[0], f[6], int(f[3]), int(f[4]), tags=tags, level=int(a.get("level", 0)))
            elif tid in tx:
                (tx[tid].exons if f[2] == "exon" else tx[tid].cds).append((int(f[3]), int(f[4])))
    for t in tx.values():
        t.exons.sort(); t.cds.sort()
    return tx

def load_or_build(gtf, cache):
    if cache and os.path.exists(cache):
        with gzip.open(cache, "rt") as fh:
            return {k: Transcript(**v) for k, v in json.load(fh).items()}
    tx = parse_gtf(gtf)
    if cache:
        with gzip.open(cache, "wt") as fh:
            json.dump({k: asdict(v) for k, v in tx.items()}, fh)
    return tx

def representative_transcripts(tx):
    """One coding transcript per protein-coding gene: MANE_Select, else Ensembl_canonical, else longest CDS among 'basic'."""
    by_gene = {}
    for t in tx.values():
        if t.gene_type != "protein_coding" or not t.cds or "cds_start_NF" in t.tags or "cds_end_NF" in t.tags: continue
        by_gene.setdefault(t.gene_id, []).append(t)
    rep = {}
    for g, ts in by_gene.items():
        def key(t): return (("MANE_Select" in t.tags), ("Ensembl_canonical" in t.tags), ("basic" in t.tags), t.cds_len())
        rep[g] = max(ts, key=key)
    return rep

class CodingModel:
    """CDS sequence, protein, and genomic<->CDS coordinate mapping for one transcript."""
    def __init__(self, t: Transcript, genome):
        self.t = t; self.g = genome
        segs = t.cds if t.strand == "+" else t.cds[::-1]
        self.map = []  # genomic 1-based positions in CDS order
        for s, e in segs:
            rng = range(s, e + 1) if t.strand == "+" else range(e, s - 1, -1)
            self.map.extend(rng)
        seq = "".join(genome.seq(t.chrom, s - 1, e) for s, e in t.cds)
        coding = seq if t.strand == "+" else revcomp(seq)
        # GENCODE CDS excludes the stop codon, so the model had no way to know where the stop was:
        # `stop_loss` could not be sampled and was instead inferred from protein length, which fires on
        # frame artefacts. The three downstream bases are appended so the stop is addressable like any
        # other codon, and their genomic positions join the coordinate map.
        stop_positions = self._stop_codon_positions(t, genome)
        self.stop_index = len(coding) // 3 if stop_positions else None
        if stop_positions:
            self.map.extend(stop_positions)
            coding += "".join(genome.seq(t.chrom, p - 1, p) if t.strand == "+"
                              else revcomp(genome.seq(t.chrom, p - 1, p)) for p in stop_positions)
        self.cds = coding
        full = translate(self.cds)
        stop_at = full.find("*")
        self.protein = full[:stop_at] if stop_at >= 0 else full
        self.has_stop = stop_at >= 0
        self.pos_index = {p: i for i, p in enumerate(self.map)}

    @staticmethod
    def _stop_codon_positions(t, genome):
        """Genomic positions of the three bases immediately after the CDS, in transcript order."""
        if t.strand == "+":
            last = t.cds[-1][1]
            pos = [last + 1, last + 2, last + 3]
            return pos if pos[-1] <= genome.lengths.get(t.chrom, 0) else []
        first = t.cds[0][0]
        pos = [first - 1, first - 2, first - 3]
        return pos if pos[-1] >= 1 else []
    def cds_index(self, gpos): return self.pos_index.get(gpos)
    def codon_number(self, gpos):
        i = self.cds_index(gpos); return None if i is None else i // 3
    def exon_number(self, gpos):
        for n, (s, e) in enumerate(self.t.exons if self.t.strand == "+" else self.t.exons[::-1], start=1):
            if s <= gpos <= e: return n
        return None
    def is_last_coding_exon(self, gpos):
        segs = self.t.cds if self.t.strand == "+" else self.t.cds[::-1]
        s, e = segs[-1]; return s <= gpos <= e
    def mutate_protein(self, gpos, ref, alt):
        """Apply a simple substitution/indel at genomic position (1-based, ref/alt on + strand) and return
        (mutant_protein, first_changed_aa_index, consequence). Handles SNV, insertion, deletion within CDS."""
        i = self.cds_index(gpos)
        if i is None: return None, None, "non_coding"
        cds = self.cds
        if self.t.strand == "-":
            ref_t, alt_t = revcomp(ref), revcomp(alt)
            # on minus strand the anchor base is the last base of the reversed allele
            i = i - (len(ref_t) - 1)
        else:
            ref_t, alt_t = ref, alt
        if cds[i:i + len(ref_t)] != ref_t: return None, None, "ref_mismatch"
        mut = cds[:i] + alt_t + cds[i + len(ref_t):]
        # extend with downstream genomic sequence so frameshifts can run into the 3'UTR
        tail = self._downstream(3000)
        mprot = translate(mut + tail)
        wprot = self.protein
        stop = mprot.find("*"); mprot = mprot[:stop] if stop >= 0 else mprot
        k = 0
        while k < min(len(mprot), len(wprot)) and mprot[k] == wprot[k]: k += 1
        d = len(alt_t) - len(ref_t)
        if d == 0 and len(ref_t) == 1:
            codon = i // 3
            # consequence is decided by which codon was hit, not by how the protein lengths compare
            if self.stop_index is not None and codon == self.stop_index:
                cons = "stop_loss" if len(mprot) > len(wprot) else "synonymous"
            elif codon == 0:
                cons = "start_loss" if mprot[:1] != "M" else "synonymous"
            elif k >= len(wprot) and len(mprot) == len(wprot): cons = "synonymous"
            elif len(mprot) < len(wprot): cons = "nonsense"
            elif len(mprot) > len(wprot):
                # a substitution outside the stop codon cannot lengthen the protein; this means the
                # transcript's reading frame does not close cleanly, so the candidate is unusable rather
                # than mislabelled as missense
                return None, None, "frame_unresolved"
            else: cons = "missense"
        elif d % 3 == 0: cons = "inframe_insertion" if d > 0 else "inframe_deletion"
        else: cons = "frameshift"
        return mprot, k, cons
    def _downstream(self, n):
        """Sequence after the coding region, starting past the stop codon.

        The stop codon is part of `self.cds`, so the tail must begin after it; otherwise a readthrough
        immediately meets a duplicated copy of the stop it just destroyed and appears to extend by one
        residue.
        """
        t = self.t
        skip = 3 if self.stop_index is not None else 0
        if t.strand == "+":
            last = t.cds[-1][1] + skip
            return self.g.seq(t.chrom, last, last + n)
        first = t.cds[0][0] - skip
        return revcomp(self.g.seq(t.chrom, first - 1 - n, first - 1))
