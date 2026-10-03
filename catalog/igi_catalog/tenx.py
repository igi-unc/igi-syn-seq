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

R1_LEN = 26          # 16 bp barcode + 10 bp UMI, the 5' v2 read-1 length
R2_LEN = 90          # cDNA
FIVE_PRIME_WINDOW = 450   # bp from the transcript 5' end that R2 can reach


def r1_sequence(barcode, umi_seq):
    """Read 1: barcode then UMI, exactly as Cell Ranger parses it."""
    return (barcode + umi_seq)[:R1_LEN]


def r1_quality(rng, length=R1_LEN, hi="I", lo="F"):
    """A plausible R1 quality string.

    R1 is 26 cycles on a fitted-profile instrument, so its qualities are high and nearly flat; the
    four-level binning of current chemistry means almost every base is the top bin. A fitted per-cycle
    profile is not used here because R1 is not simulated from a template -- inventing read errors in the
    barcode would silently move cells between barcodes, which is a different experiment.
    """
    return "".join(hi if rng.random() > 0.02 else lo for _ in range(length))


def five_prime_window(seq, window=FIVE_PRIME_WINDOW):
    """The 5'-proximal slice of a cDNA that 5' chemistry actually sequences."""
    return seq[:window] if len(seq) > window else seq


def reads_for_molecule(rng, mean_reads):
    """How many times one molecule is sequenced. At least one, so no molecule silently vanishes."""
    n = int(rng.gauss(mean_reads, max(1.0, mean_reads * 0.6)))
    return max(1, n)


def cellranger_names(out_dir, sample, lane=1):
    """The filenames Cell Ranger requires: it parses the sample, lane and read from the name."""
    import os
    base = os.path.join(out_dir, f"{sample}_S1_L{lane:03d}")
    return f"{base}_R1_001.fastq.gz", f"{base}_R2_001.fastq.gz"
