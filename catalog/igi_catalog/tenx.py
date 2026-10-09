"""10x Chromium 5' v2 short-read library construction.

Cell Ranger expects two reads per fragment, and they carry different things:

    R1  26 bp   16 bp cell barcode + 10 bp UMI
    R2  90 bp   cDNA

So R1 is not sequenced from the template at all -- it is the barcode and UMI, which this module knows
exactly because it drew them. Only R2 needs a read simulator, and it gets the same fitted NovaSeq X quality
profile the WES, WGS and bulk RNA arms use, via ART, so the four Illumina assays share one error model
rather than each inventing its own.

Two properties of the real assay that the naive construction misses:

**PCR duplicates.** The design asks for ~40,000 reads per cell over ~8,000 captured molecules, so a
molecule is sequenced about five times on average. Emitting one read per molecule would produce a library
with no duplicates at all, and UMI deduplication -- the thing the UMIs exist for -- would have nothing to
collapse. Reads per molecule are drawn around the mean so the duplicate-rate distribution is realistic
rather than uniform.

**5' bias.** 5' chemistry captures the transcript's 5' end, so R2 clusters there instead of spanning the
whole cDNA. Only the 5'-proximal window is handed to the simulator, which both reproduces that bias and
keeps the intermediate FASTA to a third of its full-length size.

The barcode and UMI come from the shared `MoleculePool`, so a barcode-UMI-transcript triple seen here is
the same molecule Kinnex single cell and ONT single cell see. That only holds because `MoleculePool.draw`
is a pure function of the cell; it used to depend on how many cells had been drawn before it.
"""
from .genome import revcomp
import bisect
import json
import os

R1_LEN = 26          # 16 bp barcode + 10 bp UMI, the 5' v2 read-1 length
R2_LEN = 90          # cDNA
FIVE_PRIME_WINDOW = 450   # bp from the transcript 5' end that R2 can reach


def r1_sequence(barcode, umi_seq):
    """Read 1: barcode then UMI, exactly as Cell Ranger parses it."""
    return (barcode + umi_seq)[:R1_LEN]


_R1Q_PATH = os.path.join(os.path.dirname(__file__), "..", "resources", "tenx_r1_quality.json")
_R1Q = {}


def r1_profile(kind="gex"):
    """Per-cycle (characters, cumulative weights) for R1, measured from the real 10x libraries."""
    if not _R1Q:
        with open(os.path.abspath(_R1Q_PATH)) as fh:
            _R1Q.update(json.load(fh))
    key = f"_cum_{kind}"
    if key not in _R1Q:
        cyc = []
        for cycle in _R1Q[kind]["cycles"]:
            chars, cum, tot = [], [], 0.0
            for ch, w in cycle:
                tot += w
                chars.append(ch)
                cum.append(tot)
            cyc.append((chars, cum, tot))
        _R1Q[key] = cyc
    return _R1Q[key]


def r1_quality(rng, length=R1_LEN, kind="gex"):
    """An R1 quality string drawn per cycle from the real per-cycle distribution.

    R1 is not simulated from a template -- inventing read errors in the barcode would silently move cells
    between barcodes, which is a different experiment -- so only the QUALITIES are modelled, never the
    bases. They were previously a flat two-level string, the top bin 98 % of the time and one lower bin
    otherwise, which gave 663 distinct quality strings in 20,000 reads where R2 had 19,461. Real NovaSeq X
    R1 is four-level binned (`#*9I`, so Q2/Q9/Q24/Q40) with a per-cycle shape: the first cycle is the
    worst at 83 % top bin, and later cycles run above 93 %. Drawing per cycle from the measurement
    reproduces both the alphabet and the shape.
    """
    cyc = r1_profile(kind)
    out = []
    for i in range(length):
        chars, cum, tot = cyc[i % len(cyc)]
        u = rng.random() * tot
        out.append(chars[bisect.bisect_left(cum, u) if u <= cum[-1] else len(chars) - 1])
    return "".join(out)


def five_prime_window(seq, window=FIVE_PRIME_WINDOW):
    """The 5'-proximal slice of a cDNA that 5' chemistry actually sequences."""
    return seq[:window] if len(seq) > window else seq


def art_read_strands(aln_path):
    """read_id -> '+' or '-' from ART's ALN file: which strand of the window ART drew.

    ART samples either strand with equal probability, and `-na` suppresses exactly the file that records
    which. The GEX and TCR builders passed `-na` and wrote whatever ART emitted straight out as R2, so
    both libraries came out UNSTRANDED -- Cell Ranger measured 45,140 sense against 44,751 antisense on
    ds-02 GEX and refused the run, because it infers chemistry from R2 strand bias and a coin flip matches
    no chemistry:

        Unable to distinguish between [SC5P-R2, SC3Pv2] chemistries based on the R2 read mapping

    Recovering the strand from the ALN rather than by testing whether the read is a substring of the
    window: a substring test fails on any read carrying a sequencing error, which at HS25 is most of them.
    The ALN is exact, and it costs nothing -- reads are BIT-IDENTICAL with and without `-na` at the same
    `-rs`, verified before this was written.

    ALN record format, three lines per read:

        >ref_id  read_id  aln_start  strand
        <reference segment, gapped>
        <read, gapped>
    """
    strands = {}
    if not os.path.exists(aln_path):
        return strands
    with open(aln_path) as fh:
        for line in fh:
            if line.startswith(">"):
                f = line[1:].split()
                if len(f) >= 4:
                    strands[f[1]] = f[-1]
    return strands


def art_chunk_index(read_id):
    """The chunk index encoded in an ART read id, e.g. `m12-1` -> 12.

    This is not cosmetic. Both builders write the ART input FASTA as `>m{i}` with `i` the index within
    the chunk, but SKIP any record whose window is shorter than the read length -- and then paired the
    k-th read in the FASTQ with `chunk[k]`. One skipped short transcript therefore shifted every
    subsequent read onto the wrong barcode, UMI and truth-map record. Keying on the id ART carries
    through removes the assumption that the FASTQ is dense.
    """
    head = read_id.split("-")[0]
    return int(head[1:]) if head[1:].isdigit() else None


def orient_antisense(seq, qual, strand):
    """Orient one R2 antisense to its transcript window, as 10x 5' chemistry produces.

    In 5' chemistry the barcode and UMI are incorporated at the transcript's 5' end by template
    switching, so R1 reads the 5' end in the sense direction and R2 reads back toward it -- antisense.
    Cell Ranger's own chemistry definitions encode this: SC5P-R2 expects antisense R2, SC3Pv2 expects
    sense, and its detector decides between them on the observed bias.

    An ART read drawn from `+` is the sense window, so it is reverse complemented here; one drawn from
    `-` is already the antisense strand and is left alone. The quality string is reversed with the
    sequence, never complemented.
    """
    if strand == "+":
        return revcomp(seq), qual[::-1]
    return seq, qual


def reads_for_molecule(rng, mean_reads):
    """How many times one molecule is sequenced. At least one, so no molecule silently vanishes."""
    n = int(rng.gauss(mean_reads, max(1.0, mean_reads * 0.6)))
    return max(1, n)


def cellranger_names(out_dir, sample, lane=1):
    """The filenames Cell Ranger requires: it parses the sample, lane and read from the name."""
    base = os.path.join(out_dir, f"{sample}_S1_L{lane:03d}")
    return f"{base}_R1_001.fastq.gz", f"{base}_R2_001.fastq.gz"
