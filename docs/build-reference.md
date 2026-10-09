# IGI-SYN-SEQ build reference

How each deliverable is actually produced, so that one aspect can be corrected without rebuilding
everything. For *why* a choice was made see `IGI-SYN-SEQ-design.md` (specification) and
`catalog-design-notes.md` (design rationale and measurements). This document is the operational map:
which code runs, what it reads, what it writes, what it costs, and what to edit to change a given
behaviour.

Two things to know before using it:

- **The build is sbatch arrays over Python drivers, not Nextflow.** Section 13 of the design
  specification describes a Nextflow DSL2 pipeline. That is the intended end state; it does not exist.
  What produced the current release is `catalog/run_*.py` launched by the job scripts in `jobs/`, whose
  `README.md` gives the run order. Where the two disagree, this document is correct.
- **Stages are independently re-runnable and that is the point.** Every stage writes a file and the
  next stage reads it. Changing a read simulator does not re-run netMHCpan; changing the catalog does
  not re-measure a quality profile. The dependency table in §7 says what a given edit forces.

## 1. The fifteen sample types

Per dataset (`IGI-SYN-SEQ-01` from HG002, `IGI-SYN-SEQ-02` from IPISRC044). "Built" means release-scale
reads exist and passed acceptance; "tested" means a single-chromosome run completed and was checked.

| # | Sample | Builder | Status |
|---|---|---|---|
| 1 | tumor WES 150x | `run_wes.py` | **built** |
| 2 | normal WES 100x | `run_wes.py` | **built** |
| 3 | tumor bulk RNA 80M pairs | `run_rna.py` | **built** |
| 4 | tumor bulk WGS 30x Illumina | `run_wgs.py` | tested (chr21) |
| 5 | normal bulk WGS 30x Illumina | `run_wgs.py` | tested (chr21) |
| 6 | tumor WGS 30x PacBio HiFi | `run_pacbio_wgs.py` | tested (chr21) |
| 7 | normal WGS 30x PacBio HiFi | `run_pacbio_wgs.py` | tested (chr21) |
| 8 | tumor Kinnex bulk RNA | `run_kinnex_bulk.py` | **validated on chr21** |
| 9 | tumor Kinnex scRNA (MAS 16-mer, 10x 5' v2) | `run_kinnex_sc.py` | written, adapter layout unresolved |
| 10 | tumor 10x 5' GEX | not written | — |
| 11 | tumor 10x 5' TCR | not written | — |
| 12 | tumor ONT scRNA | not written | — |
| 13 | tumor ONT bulk RNA | not written | — |
| 14 | tumor ONT WGS | not written | — |
| 15 | normal ONT WGS | not written | — |

Earlier counts of "thirteen" and "fourteen" in the channel log were wrong; 15 is the count after the
2026-10-01 decisions to add Illumina bulk WGS and to mirror the four PacBio sample types in ONT. The
design specification's §3 table lists 12 and predates the ONT bulk RNA and ONT WGS decisions.

Samples 9–12 share one per-cell roster, so a barcode means the same cell in all four. That coupling
lives in `cells.py`/`molecules.py` and is described in §5.

## 2. Configuration: three files, three different jobs

| File | Tracked | Holds | Change it to |
|---|---|---|---|
| `catalog/design.yaml` | yes | what the datasets *are*: clone trees, purity, CNAs, event counts per tier, hotspot list, HLA haplotypes, seed | change the biology or the catalog size |
| `config/sites/unc-lccc.paths.yaml` | yes | where inputs live *and* every measured simulator constant | move to another site, or recalibrate a simulator |
| `catalog/paths.example.yaml` | yes | the same keys with placeholder values, as the documented contract | add a new input key |

Site files contain absolute paths, which is why they are per-site and why no other tracked file may
contain one. A new site copies `paths.example.yaml`; `validate_catalog.py` fails on a missing key
rather than silently defaulting.

Keys read by each driver, extracted from the code rather than from memory:

| Driver | `design.yaml` keys | `paths.yaml` keys |
|---|---|---|
| `igi_catalog.designer` | `baseline` `clones` `hla` `sex` `seed` `tiers` `counts` `hotspots` | `gtf` `annotation_cache` `reference_fasta` `expression_tsv` `germline_vcf` `exome_bed` `rmsk_bed` `segdups_bed` `netmhcpan_cmd` `netmhcpan_threads` `erv_loci` `erv_all_loci` `cta_gene_list` `pan_normal_quant` `arms_bed` |
| `run_background.py` | `background_mut_per_mb` `signatures` `clones` `sex` | as designer, plus the COSMIC subset |
| `run_germline_spikein.py` | `baseline` `sex` `cna` | `germline_vcf` `germline_vcf_baseline` |
| `run_splice_causal.py` | — | `reference_fasta` |
| `run_wes.py` | `purity` `clones` | `art_cmd` `exome_bed` `gc_bias_curve` `art_profile_r1` `art_profile_r2` `viral_fasta` `reference_fasta` |
| `run_wgs.py` | `purity` `clones` | `art_cmd` `art_profile_r1` `art_profile_r2` `viral_fasta` `reference_fasta` |
| `run_rna.py` | `purity` `clones` | `art_cmd` `viral_fasta` `reference_fasta` `expression_tsv` |
| `run_pacbio_wgs.py` | `purity` `clones` | `badread_cmd` `samtools_cmd` `pacbio_length_mean` `pacbio_length_sd` `pacbio_identity` `pacbio_error_model` `pacbio_qscore_model` `pacbio_chunk_bp` `pacbio_chunk_timeout_s` |

## 3. Stage 1 — the catalog

`python3 -m igi_catalog.designer --design design.yaml --paths <site> --dataset <D> --out output/`
(`jobs/01_catalog.sbatch`, array 1-2, 24 cores / 64 GB / ~1 h per dataset, dominated by netMHCpan.)

Writes to `catalog/output/`, all tracked, all deterministic given `seed`:

| File | Content |
|---|---|
| `<D>.snv_indel.tsv` | designed SNVs and indels with clonality / expression / binding / context tiers |
| `<D>.snv_indel.summary.json` | the grid fill report: which (clonality x expression x binding) cells were satisfiable |
| `<D>.fusions.tsv`, `.svs.tsv`, `.viruses.tsv`, `.expressed.tsv`, `.hla_loh.tsv` | the other antigen classes |
| `<D>.splice_causal.tsv` | splice-altering variants; a *genomic* edit, so read builders take it via `--extra-events` |
| `<D>.germline_spikein.tsv` | designed germline alleles, from `run_germline_spikein.py` |
| `<D>.background.tsv` | passenger mutations, from `run_background.py` |
| `*.diffcards.txt` | per-event debugging cards |

The three auxiliary drivers exist separately because they are expensive and orthogonal: the background
model is a whole-genome signature draw (~2 min, 2,926 mutations at 1.0 mut/Mb), the germline spike-in
rewrites a VCF, and the splice designer needs only the genome. Re-running the designer does not re-run
them, so **a designer change requires re-running `run_background.py` only if the reserved-position set
changes** — the background model reserves designed positions so the two never collide.

Binding is netMHCpan `%rank_EL` over the patient's own alleles, and HLA loss is applied twice: a
`by_binding_tier` before loss and `by_binding_tier_retained` after. For `IGI-SYN-SEQ-02`, 36 events
have their best allele on the lost haplotype and 22 change tier as a result. Anything that reads a
binding tier must say which of the two it means.

## 4. Stages 2-4 — short-read builders

`run_wes.py` and `run_wgs.py` share `pipeline.py`; `run_rna.py` uses `rna.py`. All three follow the
same shape and it is worth stating once:

1. Build the clone/haplotype set from `clones.py` (`excl` exclusive fractions, `cn()` per region,
   pre- vs post-CNA event timing).
2. For each target region, rebuild the sequence per (clone, haplotype) that retains it, in **reference
   coordinates** — `genome_build.haplotype_sequence()` applies germline, designed, background and extra
   events to a fresh window rather than mapping positions through indels. This is why allele fractions
   come out right by construction instead of being corrected afterwards.
3. Group identical requests into *sources*, keyed `(src, hap, kind, pid, gc_bin)`, write one FASTA per
   source, and give ART a per-source depth equal to the clone's molecule share times haplotype copy
   number times GC factor.
4. Shuffle and rename (`readnames.py`) into instrument-style names, emitting a read map that ties every
   read back to its source record.

**Regions.** WES uses capture intervals from `exome_bed` plus off-target bands
(`pipeline.off_target_bands()`, design §10.1). WGS uses `whole_chromosome()` in 5 Mb windows with a
flat `GcBias(None)` — the GC curve is *capture efficiency* and does not apply without a capture.

**Read names must not collide across the 96 independent array tasks.** Two mechanisms, both needed:
`NAME_SPACE_STRIDE = 4_000_000_000` offsets the read index by chromosome ordinal, and the record id
carries a per-chromosome `record_prefix`. The second was added after the first release: without it
record ids restarted at 1 on every chromosome and 667,090 of 1,000,832 were duplicated — and the
integrity check passed anyway, because it was checking that names resolve, not that they are unique.
`jobs/11_merge.sbatch` now proves uniqueness by sorting, and fails the job if it does not hold.

**Cost.** WES at release scale: 96 tasks, 8 cores / 48 GB, under 16 h wall; 215M read pairs total.
RNA: 2 tasks, 16 cores / **320 GB**, 6h31m and 6h45m, 76.5M pairs each. The RNA memory is high because
`--min-tpm 0.01` keeps roughly three times the records that the 0.5 default would.

**`--min-tpm` is a benchmark-correctness parameter, not a performance knob.** LENS filters expression at
`-p 50 --exclude-zeros`, the median of the sample's own non-zero TPM, so truncating the low tail raises
that bar and shifts every downstream expression threshold. A 0.5 floor puts p50 at 0.895 against a real
0.115 (7.8x inflated); 0.01 puts it at 0.168 (1.5x). The floor is set so truncation comes from
sequencing depth, not from the generator.

## 5. Single-cell coupling (unbuilt, designed)

`cells.py` builds the roster: `power_law_sizes()` solves the exponent numerically for a target top-10
clone share, and the tumour split is derived from `clones.excl` with a hard check that the exclusive
fractions sum to 1. It does **not** read clone fractions from prose — an earlier version hardcoded them
and silently disagreed with `design.yaml`.

`molecules.py` builds one pool of (cell, transcript, UMI) molecules, abundance scaled by the clone's
copy number at that locus; `_scale_for()` returns 0 for transcripts of non-descendant clones and for
zero-copy haplotypes. Kinnex sc, 10x GEX, 10x TCR and ONT scRNA each sample that one pool
independently, which is what makes a barcode mean the same cell in all four.

## 6. Long-read builders

### 6.1 Kinnex bulk RNA (MAS-seq 8-mer)

A Kinnex library ligates several full-length cDNAs into one long molecule separated by known adapters and
sequences that as a single HiFi read; `skera split` cuts it back into segments. So the builder constructs
*arrays*, not segments, and hands the splitting to the real skera -- segmenting it ourselves would make the
pre-skera BAM pointless, since the reason to ship it is that skera's own output is reproduced.

    transcript records (shared with bulk RNA via rna_assembly.assemble)
      -> molecules sampled by abundance x size selection
      -> arrays A0 S0 A1 S1 ... S7 A8, forward frame only
      -> badread, one full-length read per array
      -> pre-skera BAM with zm/np/rq  ->  real `skera split`  ->  segmented BAM

Each link was tested rather than assumed:

| Assumption | Test | Result |
|---|---|---|
| badread emits a whole short reference as one read | 30 synthetic arrays | all 32 reads started at 0 and ran to the end |
| skera accepts a synthetic BAM | `samtools import` output | **fails**: "Bam record missing the read quality tag". `pacbio_bam.py` writes `zm`/`np`/`rq` |
| skera recovers what was built | 22 arrays | 22/22 segment counts exact, 153/159 lengths within 2 % |
| the adapter set is in one frame | build all-forward, ask skera | coherent ascending/descending runs per ZMW |

**The adapter set was recovered from real data**, not transcribed from a kit document, which also guarantees
the simulated arrays are segmentable by the same skera that produced the reference set. skera writes a
MessagePack `ds` tag on every segment carrying each adapter's label and observed sequence;
`jobs/measure/fit_mas_adapters.py` decodes it. The tags also give the grammar, and one of them does not mean
what it looks like: `dl` is the **left adapter index**, not a segment count, and `di` is the segment index.

Two modelling decisions matter more than they look:

1. **Molecules are sampled by abundance alone, not abundance x length.** Kinnex yields one read per molecule
   and TPM is already length-normalised. On the release manifest the two weights give mean cDNA lengths of
   1,788 bp and 4,197 bp, a 2.3x error in which transcripts dominate the library.
2. **Size selection is fitted and it is doing heavy lifting.** Even with the right weight the molecules run
   short against real segments (p10 498 bp against 1,355), so
   `jobs/measure/fit_kinnex_size_selection.py` fits an efficiency per length bin, normalised to an
   abundance-weighted mean of 1. The corrections reach 50x suppression at 500-750 bp, and they are not pure
   chemistry: the real segments are HG002 and the molecules come from a TCGA-BRCA basal baseline, so the
   ratio also carries a transcriptome difference and the low-abundance short tail a 0.01 TPM floor keeps.
   Checking that the corrected density matches the real one is an identity, not a validation.

Orientation is supplied by badread choosing a strand per read, so arrays are built in the forward frame
only. Randomising it in both places was a real bug: arrays built with the array-level flip disabled still
came back from skera 60 % reverse.

**Validated on chr21** (40,000 molecules, 5,791 arrays, 7m53s): the whole chain ran and the output matches
real Kinnex on every statistic that was not an input.

| | simulated | real |
|---|---|---|
| segment length mean / sd | 2,263 / 900 | 2,223 / 935 |
| median / p90 | 1,854 / 3,283 | 1,974 / 3,374 |
| reverse-oriented segments | 50.8 % | 52.5 % |
| complete 8-segment arrays | 74.6 % | 77.1 % |
| mean segments per array | 6.882 | 6.908 |

The length agreement is the validation the size-selection curve needed, since the identity check could not
provide one: without the curve the molecules would average 1,788 bp against a real 2,223.

Two residual notes. The p10 comes out at 1,637 bp against a real 1,355, so the curve slightly over-corrects
the short end. And of 38,903 segments, **none of the 33,728 in complete arrays violate the adapter grammar**;
530 of the 5,175 in incomplete arrays do, which is expected rather than a defect -- an array that lost a
middle segment does not follow `dl = k - di` any more than it follows `dl = di`.

The fill passes behaved as predicted to within a few tenths of a percent: 36.94 %, 13.73 % and 4.96 % of
arrays still missing after each pass against a predicted 36.8 %, 13.5 % and 5.0 %, leaving 2.09 % never
sequenced, whose molecules the array map records as such.

### 6.3 Kinnex single-cell RNA (MAS 16-mer, 10x 5' v2)

Same array machinery as bulk, with the molecules being 10x 5' v2 cDNAs drawn from the shared per-cell pool
so a barcode and UMI mean the same molecule here as in 10x gene expression, 10x TCR and ONT single cell.

**The 10x segment architecture was recovered from the real arrays**, not from a kit document, by orienting
segments against the barcode whitelist and taking a column-wise consensus. 96.1 % of real segments carry a
whitelist barcode at offset 22 and the layout that follows is exact:

| offset | content | evidence |
|---|---|---|
| 0-22 | `CTACACGACGCTCTTCCGATCT` TruSeq R1 primer | 100 % column agreement |
| 22-38 | 16 bp cell barcode | in the whitelist |
| 38-48 | 10 bp UMI | 27 % agreement, i.e. random |
| 48-61 | `TTTCTTATATGGG` 10x 5' v2 TSO | found at offset 48 = 22+16+10 in 96.8 % |
| 61.. | cDNA, then polyA (median 29 bp), then `GTACTCTGCGTTGATACCACTGCTT` | modal 25 bp tail |

**Sixteen cDNAs per array, and how I got it wrong first.** A real array read carries a median of 16 TSOs
(mean 15.1, p90 16) over 17,092 bp, against a single-cDNA median of 952 bp. That is direct and settles it.

I first concluded the array was an 8-mer and edited the design specification to say so. The evidence for
that was real but insufficient: all nine bulk adapters appear about once per read, and 83 % of the segments
skera returns carry exactly one TSO, one R1 primer and one barcode. What I missed is that 12 % of skera's
output is whole un-segmented arrays averaging **8.64 TSOs** -- the other eight cDNAs were in there. A base
check would have caught it earlier: skera's segments cover 96 % of the read bases, so 7.6 segments per read
at a 952 bp median cannot be the whole story. My first attempt at that check used a 2,000-read file against
a 3,093-read skera run and so reported an impossible 148 %, which I took as a measurement artifact instead
of pursuing.

**The adapter layout is not resolved, and that is this assay's open item.** Nine adapters are identified and
no tenth is detectable: inside the un-segmented blocks every high-multiplicity motif is per-cDNA 10x
structure -- the R1 primer, the TSO, the polyA-SMART junction at about six copies per block -- rather than a
distinct array adapter. Nine adapters cannot bracket sixteen segments, so some cDNA boundaries must be
adapter-free, joined where one cDNA's SMART primer meets the next one's R1 primer. `MasArrays.build_sc`
therefore brackets the first eight boundaries with the known adapters and joins the rest directly, which
**reproduces what skera demonstrably does to the real data** -- about six clean single-cDNA segments per read
plus one large block -- rather than asserting a chemistry. Built arrays come out at 16,924-18,122 bp against
a real 17,092. Substituting the true MAS-16 adapter list, from PacBio or from any reference segmented BAM
for the single-cell kit, would remove the assumption.

**Single-cell size selection is fitted separately** (`kinnex_sc_size_selection.json`), because the 10x cDNA
distribution differs sharply from bulk: median 952 bp against 1,974, mean 1,072 against 2,223, max 4,876.
The curve is fitted only on segments carrying exactly one TSO and one R1 primer. Fitting it on skera's raw
output instead put 12 % of density at 7-15 kb and demanded a 6x boost there, which would have reproduced the
un-segmented-array artifact as though it were biology.

**Validated against the real single-cell arrays**, on the metrics that were not inputs:

| | synthetic | real |
|---|---|---|
| segments per read (skera) | 7.93 | 7.607 |
| mean segment length | 2,276 | 2,155 |
| full arrays (skera) | 96.3 % | 92.5 % |
| single-cDNA fraction | 76.7 % | 81.6 % |
| single-cDNA median / mean length | 1,036 / 1,140 | 952 / 1,072 |
| whitelist barcode at offset 22 | 93.3 % | 96.1 % |

That last row was 0.0 % on the first attempt, and it is the one worth remembering. A 1-2 bp random spacer
inserted at every adapter boundary shifted each segment's contents by 1-2 bp, so the barcode sat at 23 or 24
and a tool indexing it at a fixed offset would have found nothing -- while 89 % of the barcodes were present
and correct, so every other check passed. skera reports that spacer, but at array ENDS; it does not belong at
interior boundaries. The residual 23-25 offsets after the fix are HiFi indels shifting the frame, 21 segments
of 971.

**Captured molecule sets nest across assays rather than overlapping partially.** `mean_molecules` is
molecules *captured*, so it is legitimately assay-specific -- 10x gene expression uses 8,000 per cell and
Kinnex single cell 10,600, the latter set so that a 12 M post-skera segment target is reachable. Because
both draw from one RNG stream keyed on the barcode, the smaller capture is an exact prefix of the larger:
verified at 100 % overlap, with the 8,000-molecule set a strict subset of the 10,600-molecule set.

That is a simplification, and in the generous direction. Two independent library preps of the same cell
would overlap partially, not nest; perfect nesting maximises cross-assay agreement, which is what a truth
set wants but is not what two real captures do. A consumer comparing UMI counts per cell between the two
assays will see a systematic difference that is capture depth rather than biology.

**A defect fixed in the shared machinery.** `MoleculePool.draw` drew from one pool-wide random stream, so the
molecules a cell received depended on how many cells had been drawn before it. Gene expression and Kinnex
asking for the same cell got different molecules, and drawing the same cell twice from one pool did not even
agree with itself -- the exact opposite of the cross-assay consistency the module exists to provide. It is
now seeded per barcode and is a pure function of the cell, verified across call orders.

### 6.2 PacBio HiFi WGS

`run_pacbio_wgs.py`. Shape:

1. `rearrange.py` produces the derived chromosome per clone/haplotype, including SVs.
2. `split_fasta()` chunks it at `pacbio_chunk_bp` (20 Mb).
3. A `ThreadPoolExecutor` runs Badread over the chunks, 8 concurrent, each with
   `pacbio_chunk_timeout_s` and up to 3 retries.
4. Chunks are combined newline-safely with a framing assertion, then converted to a `hifi_reads` BAM.

Four lessons are baked into that list and each cost real time:

- A blocking `subprocess.run` loop managed 5 chunks/hour. The thread pool does ~98/hour.
- Piping Badread to `gzip` deadlocked, leaving six shells at 0% CPU holding executor slots. The gzip
  was removed (the output is converted to BAM anyway) and the timeout added so a hang cannot be silent.
- `cat` over chunks is not newline-safe: 6 of 54 chunk FASTQs lacked a trailing newline and
  concatenation merged records. Hence the explicit framing assertion.
- Running this on the login node alongside other work dropped throughput 18x. It runs under sbatch.

**Badread's flags do not mean what they appear to.** Realized identity is offset from the requested
value by **+0.43 points for PacBio and −0.22 for ONT cDNA** — not a constant, so each model is
calibrated separately against real reads and the calibrated value is what `paths.yaml` stores.
`ERRHMM-ONT-HQ` measured marginally *worse* than plain `nanopore2023`. pbsim3 was rejected: its
`--accuracy-mean` is ignored, both methods returning 0.964 regardless.

## 7. Measured constants: provenance and re-measurement

Everything below is fitted from real data. Each row says what to re-run to change it.

| Constant(s) | Value | Measured from | Re-measure with |
|---|---|---|---|
| `art_profile_r1/r2` | 4-level Q3/Q13/Q25/Q41 | 4M pairs, IPISRC044 blood-normal WGS, NovaSeq X LH00499, 150 cycles | `jobs/measure/art_profile.sbatch` |
| `gc_bias_curve` | 0.40 at 30-35% GC to 1.25 at 60-65%, mean-1 normalised | 161,618 capture intervals, IPISRC044 WES normal | `jobs/measure/gc_bias.sbatch`, which calls `fit_gc_bias.py` |
| `pacbio_length_mean/sd` | 16689 / 4593 | 200k reads, HG002 Revio HiFi | `jobs/measure/longread_measure.sbatch` |
| *realised* `sd` | ~4,870 | measured across 29 release chromosomes | see the note below |
| `pacbio_identity` | `99.4,99.8,0.5` (targets 0.99823) | same | same job, then apply the +0.43 offset |
| `kinnex_segment_length_mean/sd` | 2224 / 933 | 200k post-skera segmented reads | same job |
| `ont_identity` | request `98.41` for a 0.9822 target | IPISRC044 R10.4.1 SUP scRNA, 2025-07-17 | `measure_longread.py`; the ONT cDNA model runs ~0.19 points low |
| `ont_length_mean/sd` | 910 / 461 | same | **applies to ONT scRNA only; see the caveat below** |
| `ont_wgs_length_mean/sd`, `..._identity_target` | 19117 / 15530, 0.98676 | ONT open-data `giab_2025.01/basecalling/sup/HG002/PAW70337` | `measure_longread.py`; the BAM is remote, see `ont_wgs_bam` |
| expression baseline | per-transcript median TPM | 191 TCGA-BRCA basal-like tumours, UCSC Xena Toil recompute | not yet tracked; the original job is in scratch |
| `background_mut_per_mb`, `signatures` | 1.0 mut/Mb; SBS3 0.6 / SBS1 0.15 / SBS5 0.15 / SBS13 0.1 | COSMIC v3.4, TNBC-typical | edit `design.yaml` |

Two notes on these:

- The ART profile comes from WGS, not WES, because IPISRC044 has no 150 bp exome (BostonGene exomes are
  101 bp normal / 140 bp tumour). A per-cycle quality profile is a property of the run, not the capture,
  and capture is modelled separately by the GC curve.
- GIAB's own HG002 ONT is all R9.4-era (2D reads from 2016, Guppy v2/v3) and unusable as an R10.4.1
  target, which is why the ONT WGS parameters come from ONT's open-data bucket instead.
- `ont_wgs_length_sd` ≈ `mean` with a long tail. **Verified**: Badread reproduces it. At a requested
  19,117/15,530 the realised distribution is mean 18,554 and sd 15,378 -- within 2.9 % and 1.0 % -- with a
  median of 14,134, a p90 of 38,554 and a maximum of 129,812, which is the shape real ONT gives. The gamma
  length model therefore handles sd ≈ mean, which had been the open risk for this assay.
- **ONT RNA read length is set by the molecule, not by a parameter, and only the single-cell target is
  measured.** For the ONT RNA assays one molecule is one read, so the read-length distribution comes from
  the transcript and the size-selection curve rather than from `ont_length_mean`. Against the real
  IPISRC044 ONT scRNA median of 910 bp, the simulated single-cell reads come out at 1,129 bp -- 24 % long,
  because the size curve is fitted on Kinnex single cell rather than on ONT. The bulk ONT arm has **no real
  reference at all**: the only IPISRC044 ONT data is single cell, so bulk uses the Kinnex bulk size curve as
  the best available proxy and comes out at 2,216 bp against nothing to compare it with. Neither figure
  should be quoted as validated; obtaining real ONT bulk cDNA would settle it.
- ONT genomic identity is calibrated like PacBio's and in the opposite direction: Badread's realised value
  runs about **0.49 points below** the request (98.409 for a requested 98.9), against PacBio's +0.43 above.
  `ont_wgs_identity` is therefore 99.17 to land on the measured 0.98676. Three models, three different
  offsets, none of them the flag's face value.
- **Badread's realised length spread runs about 6 % wide.** Across 29 release chromosomes the mean read
  length is 16,768 bp against a target 16,689 (+0.5 %) and the accuracy 0.99830 against 0.99823, both well
  within tolerance -- **but that accuracy is derived from the quality strings, not from an alignment, and
  the two differ by 5.4x.** Measured against the reference on reads simulated from it, the true accuracy is
  0.99089 (error 9.11e-03/bp) while the quality strings claim 0.99833 (1.67e-03/bp). Real HiFi is
  1.77e-03/bp, so the delivered reads carry about 5.1x the error rate of real HiFi. `pacbio_identity` is a
  mean request of 99.4 % and the reads deliver on it; the request was tuned until the quality-derived number
  matched real, which drove the wrong quantity to agreement. See CHANNEL.md [23].
  inside tolerance, but the standard deviation is consistently 4,870 against 4,593. It is a systematic
  property of Badread's gamma length model rather than sampling noise -- every one of the 29 chromosomes is
  high, across both datasets and both libraries -- so a caller that is sensitive to the read-length tail
  will see a slightly broader distribution than a real Revio run gives. Correcting it would mean narrowing
  the requested `pacbio_length_sd` below the measured value, which trades a documented bias for an
  undocumented one; it is left as a recorded characteristic.
- `fit_gc_bias.py` reproduces the committed `ipisrc044_wes_gc_bias.json` exactly from the saved
  `bedcov.txt`, all 15 bins and both the factors and the captured-base shares. Before it was written the
  curve was a one-off computation with no script, so it could not be re-derived or corrected.

## 8. What to change, and what it forces

| To change | Edit | Then re-run | Cost |
|---|---|---|---|
| clone tree, purity, CNAs | `design.yaml` `clones`/`purity`/`cna` | designer, background, **all read builders** | full rebuild |
| event counts or tier grid | `design.yaml` `counts`/`tiers` | designer, then all read builders | full rebuild |
| a hotspot / flagpost | `design.yaml` `hotspots` | designer, then all read builders | full rebuild |
| passenger rate or signature mix | `design.yaml` `background_*`/`signatures` | `run_background.py`, then read builders | ~2 min + rebuild |
| HLA alleles or which haplotype is lost | `design.yaml` `hla` | designer (re-scores binding), then read builders | full rebuild |
| Illumina quality model | `art_profile_r1/r2` | affected read builders only | per-assay |
| capture GC behaviour | `gc_bias_curve` | `run_wes.py` only | WES only |
| exome target set | `exome_bed` | designer (context annotation) and `run_wes.py` | WES + catalog |
| WES or WGS depth | `--depth` in the job script | that builder only | per-assay |
| RNA depth or expression floor | `--pairs` / `--min-tpm` | `run_rna.py` only | RNA only |
| long-read length or identity | `pacbio_*` / `ont_*` in site paths | that long-read builder only | per-assay |
| long-read chunking or timeout | `pacbio_chunk_bp` / `_timeout_s` | that builder; output is unchanged | per-assay |
| read-name style | `readnames.py` | affected builders **and** `release_merge` uniqueness check | per-assay |

The column that matters is the third. Only four kinds of edit force a full rebuild, and all four are
changes to what the datasets *are*. Everything about how they are *sequenced* is per-assay.

## 9. Validation

`acceptance.py` runs per chromosome against an aligned BAM and the truth bundle; `class_presence.py`
checks that every designed antigen class produced reads.

Current state: **44/44 on both datasets, zero failures** — 40 per-chromosome checks (ten chromosomes x
four) plus four read-level checks. Re-established after fixing defects 12 and 13, and worth stating plainly
that the earlier reported figure was not fully earned: the junction checks in it were comparing a merged
library's 87 events against a per-chromosome expectation of 3, under the wrong seed. With correct per-chromosome
seeds all eighteen junction checks pass on exact equality, and all four of the failures the fix exposed
were the checker's seed rather than the generator -- in each case the expected set under the build's own
seed matches what was placed exactly.

The read-level checks: no duplicate names in 400,000 sampled, no name carrying truth information, every one
of 68.1 M and 68.2 M reads tracing to a source record, and 78.2 % on target, matching an independent 78.0 %
prediction. Alongside: 4/4 alignments at 99.99 % mapped, RNA truth-vs-realised correlation 0.994, and every
designed record of all six classes represented (fusion 39/39, erv 102/102, splice 132/132, cta 69/69,
virus 9/9). The three headline LOH events verified at 0.897/0.980 (TP53 17p), 0.859/0.925 (BRCA1 17q) and
0.700/0.725 (RB1 13q).

Fourteen defects have been found in the acceptance suite itself, four of them visible only at release
scale. They are listed in `catalog-design-notes.md` §14. The pattern worth carrying forward: **six of the
fourteen were checks that passed on data that was wrong.** An interval-count floor of 20 skipped an
11-interval amplicon; a depth baseline required the same chromosome, which chr13p and chr17 cannot
satisfy; a baseline neutrality probe tested only a window midpoint and crossed into LOH; the read-map
integrity check verified resolution rather than uniqueness; and the junction check compared a merged
library's 87 events against 3 expected on one chromosome with `>=`, which hid a seed mismatch that had
made the expected and placed sets completely disjoint. When adding a check, state what it would fail on,
then confirm that it does.

The two assays validated only on chr21 were each checked against their real-data target rather than
against themselves:

| Assay | Measured | Target | Source of target |
|---|---|---|---|
| PacBio HiFi | length 16,788 sd 4,916, accuracy 0.99830 (Q27.7) | 16,689 sd 4,593, 0.99823 (Q27.5) | HG002 Revio HiFi |
| Illumina WGS | 3,703,101 pairs, 0.974 of callable-base expectation | — | chr21 is 18.6 % N |

The WGS figure is the one to be careful with: measured against raw chromosome length it looks like a 21 %
shortfall, and it is not one. 8,671,409 of chr21's 46,709,983 bases are N, and yield follows callable
bases. The remaining 2.6 % is ART declining reads near N runs and window edges. Judge any WGS yield
against callable bases.

**Short-read WGS has an acceptance arm** (`acceptance.py --assay wgs`, run by `jobs/22_accept_wgs.sbatch`).
It measures depth over the reference's non-N runs instead of capture intervals, does not apply the GC curve,
and requires every designed breakpoint on the chromosome to be placed, which is stricter than the exome's
test. It has not been run at release scale because there is no release-scale WGS yet. **The long-read assays
still have no acceptance arm**, which is now the largest open gap.

## 10. Open items


- 10x 5' GEX and TCR builders; ONT scRNA, bulk RNA and WGS builders.
- ONT WGS identity calibration, and verification that Badread's length model handles sd ~ mean with a
  551 kb tail (§7).
- Release-scale runs of Illumina WGS and PacBio HiFi; both are validated only on chr21.
- Whether tumour WES should rise from 150x. At 150x the lowest-CCF subclones are marginal to invisible:
  IGI-SYN-SEQ-02 clone A1 gives 3.4 expected alt reads and clones A and B about 10-12, which is why
  acceptance excludes those tiers from the allele-fraction check as too thin. Tumour WGS has been raised
  from 30x to 100x for the same reason; the exome has not.
- An acceptance arm for the long-read assays. The short-read WGS arm exists (`--assay wgs`, §9) but has
  not yet been run at release scale, because there is no release-scale WGS to run it against.
- Resolve HG002 haplotype parentage against HG003/HG004.
- The input manifest's checksums (`catalog-freeze.txt` covers the catalog inputs; the 28-row
  `input-manifest.tsv` does not yet carry them).
- Push the local commits; request GitHub GC for the rewritten history; close the forks and PRs.

---

## 11. Review findings of 2026-10-05 and what each one changed

Eight defects were found by review of the delivered data. Four were in the packaging and were repaired in
place; three needed the reads simulated again; one was a measurement that turned out to be right once the
correct target was used. The distinction is worth keeping: a packaging defect costs one streaming pass
over the library, a read defect costs the whole arm.

| # | Defect | Root cause | Fix | Cost |
|---|---|---|---|---|
| B1 | Every ONT read header carried the source reference, strand and position, or the source molecule id | Badread's default description, never stripped; `readnames` fixed this for the Illumina arms only | `longread.combine_fastq(map_path=...)` strips it into `*_read_map.tsv.gz`; `repair_longread.py --what ont` for delivered files | 1 pass over 396 GB |
| B2 | HiFi BAMs had 11 fields, no tags, UUID names, no `@RG`, no `.pbi` | `run_pacbio_wgs.py` built them with `samtools import` while `pacbio_bam.write_hifi_bam` sat unused next to it | `write_hifi_bam` in the builder; `repair_longread.py --what hifi` for the 4 delivered BAMs | 1 pass over 290 GB |
| B3 | Kinnex segment names were `/1/ccs/17_2371` | `@RG` had no `PU`, and skera takes the movie from `PU` | `pacbio_bam.header` sets `PU`/`PM`; repair then re-runs `skera split` | 1 pass + 4 skera runs |
| B4 | `np:i:8` constant across a library while `rq` varied per read | `np_passes` pinned to the profile's mean | `draw_np` from `resources/pacbio_np_model.json` | same pass as B2/B3 |
| B5 | ONT bulk RNA reads averaged 910 bp in a library the design calls "cDNA, full length" | `run_ont.py` applied the single-cell model, fitted from 10x 5' sc cDNA, to `bulk_rna` as well | `ont_bulk_rna_length_mean/_sd` = 2091/843, measured on real HG002 Kinnex FLNC | **re-simulate 2 libraries** |
| B6 | 10x read names were `@IGI-SYN-SEQ-01:1`; R1 qualities had 663 distinct strings in 20,000 reads | the two 10x builders never used `readnames.illumina_name`; `r1_quality` was a flat two-level draw | `illumina_name` with a name-space slice; `r1_quality` draws per cycle from `resources/tenx_r1_quality.json` | **re-simulate 4 libraries** |
| B7 | CDR3s were random peptides with internal cysteines, unrelated to the V/J the truth table named | `cells.CellRoster._cdr3` drew uniformly over 20 amino acids | `vdj.recombine` from `resources/vdj_germline.json`; clonotype table gains V/J/C, `cdr3_nt` and `flagpost_antigen` | **re-simulate 2 libraries** |

### 11.1 The R1 quality target was R1, not R2

The review compared 10x R1's 663 distinct quality strings against R2's 19,461 and read the gap as a
defect in R1. Half right. Measuring the real IPISRC044 library gives **2,819** distinct strings in 20,000
reads over the first 26 cycles, because real R1 is only 26 four-level-binned cycles (`#*9I`, so
Q2/Q9/Q24/Q40) and genuinely has low variety. R2's 19,461 was never the right target. The per-cycle model
now gives 2,947 against that real 2,819, and the old 663 was indeed wrong.

### 11.2 Why the roster did not have to be rebuilt for B7

`CellRoster` has one shared random stream, and after the clonotype loop it uses it to assign doublets and
partner barcodes. Those draws are already inside four delivered libraries -- Kinnex single cell, ONT
single cell and 10x GEX, both datasets -- which share the roster. Swapping the junction generator without
replacing its consumption of that stream would have shifted every later draw and silently changed which
cells are doublets.

`CellRoster._burn_legacy_cdr3` therefore still takes the old generator's draws and throws the value away.
Verified by rebuilding the roster and diffing it against the delivered `IGI-SYN-SEQ-01_full_cells.tsv`:
`barcode`, `cell_type`, `clone`, `umi_scale`, `clonotype`, `is_doublet` and `partner_barcode` differ in
**0 of 4,000** cells, while `trb`, `tra` and `tra2` differ in all of them. Only the 10x TCR arm needed
rebuilding.

### 11.3 Fitted resources added

| File | Measured from | Holds |
|---|---|---|
| `resources/pacbio_np_model.json` | 300k reads of `HG002_PacBio-Revio_m84039_230928_213653_s3.hifi_reads.bam`; 200k of the Kinnex `segmented.bam` | empirical `np` CDF and `ec/np` ratio, separately for HiFi WGS and Kinnex |
| `resources/tenx_r1_quality.json` | 500k R1 reads each of `DSCOLAB_IPISRC044_T3_SCG1` and `..._T3_TCR1` | per-cycle quality distributions for GEX and TCR |
| `resources/vdj_germline.json` | GENCODE v37 TR segments over GRCh38 | per V the CDR3 nucleotides from the conserved Cys, per J those to the conserved Phe, plus framework |

### 11.4 Verification

The repaired HiFi BAM was put through the PacBio toolchain rather than inspected by eye: `pbindex` indexes
it (949 of 949 reads), `pbindexdump` reads the index, `extracthifi` keeps 949 of 949, `zmwfilter` parses
the hole numbers and `bam2fastq` round-trips it. `skera split` on a repaired Kinnex BAM produces
`m84000_260101_000000_s0/1/ccs/17_1856` with the segment coordinates unchanged from the original run, so
the arrays still line up with their truth rows.

---

## 12. Focal copy-number events, and why they were missing

Found by acceptance, not by review: `PTEN_homdel` came back with an observed depth ratio of **1.0** against
an expected 0.30, and `RB1_homdel` **1.03** against 0.379. The expectation was right -- 0.30 is 1 minus
purity, which is exactly what a homozygous deletion should leave -- so the libraries were wrong.

Three builders had the same defect in three different sizes.

| Arm | Copy number evaluated | Consequence |
|---|---|---|
| Illumina WGS | once per 5 Mb window, at the midpoint | a 0.3 Mb event has a 6 % chance of containing the midpoint, so focal events were absent |
| PacBio HiFi WGS | **once per chromosome**, at `chrom_len // 2` | no CNA represented at all unless it spans the chromosome midpoint |
| ONT WGS | **once per chromosome**, at `chrom_len // 2` | same |

The windows and the whole-chromosome plan were both shortcuts that assumed copy number is constant over
the span they stand for. That assumption is what was never checked.

### 12.1 The fix

`Clones.cn_boundaries(chrom)` returns every 1-based coordinate where copy number can change, across all
clones and both the CNA and designed-SV interval indexes.

- **Illumina**: `run_wgs.whole_chromosome` splits at those boundaries before tiling. Every window is now
  copy-number-uniform, which is the property the midpoint shortcut always assumed. chr10 goes from 27
  windows to 29, and the new 0.40 Mb window at 87,700,001 is `PTEN_homdel` itself.
- **Long read**: `derived.cn_slices` cuts the *derived* chromosome into copy-number-uniform pieces and
  `build_derived` computes a coverage per piece. A piece whose weight is zero is not simulated, which is
  what a deletion is. Mapping reference boundaries onto derived coordinates goes through
  `rearrange.coordinate_map`, so it survives inversions and insertions.

`run_pacbio_wgs.py` also stopped carrying its own copy of the derived-chromosome loop. The copies had
drifted: this fix landed in `derived.py` and would silently have missed PacBio.

### 12.2 Verification

Summed plan weight at the event over the same at a control locus on the same chromosome:

| Event | at event | control | ratio | expected | old (midpoint, whole chromosome) |
|---|---|---|---|---|---|
| PTEN_homdel | 0.3000 | 1.0000 | 0.300 | 0.30 | 1.0000 |
| RB1_homdel | 0.3793 | 1.0000 | 0.379 | 0.379 | 1.0000 |
| MYC_amp | 4.5000 | 0.6500 | 6.923 | — | 1.0000 |

End to end on chr10, PacBio: `T_hap0_all_s000` spans 1-87,700,000 and stops exactly at the boundary;
slice `s001` is skipped for T, A, A1 and B on both haplotypes; `NORMAL` keeps every slice. 46 slices
across the ten derived copies, 38 simulated, 8 with no copies.

`check_cn_windows.py` asserts the Illumina property permanently. It tests boundary **containment**
(`start < c <= end`), not sampled copy-number states: a first version sampled each window's start,
midpoint and end and reported zero straddled windows on chr10 and chr13 -- the two chromosomes whose
focal deletions were entirely missing -- because the event sits in the window's interior and none of the
three points touches it. On the old tiling the containment check fails 14 of 264 windows over five
chromosomes; on the new one, 0 of 275.

### 12.3 What this cost

Rebuilding all three WGS arms is 504 task-hours, measured from the previous runs rather than estimated:
48 Illumina tumour tasks at 4.23 h, 48 normal at 1.30 h, 96 PacBio at 1.41 h, 96 ONT at 1.08 h. Wall
clock is set by how much of the 105-job association limit is free, not by the work: 29.7 h at 17
concurrent slots, 7.2 h at 70.

`50_merge_longread.sbatch` had to be fixed before any of it could land. All three merges were guarded with
`[ -s "$out" ] ||`, so with the previous merged libraries in place the merge would have skipped all twelve
tasks and reported success, leaving the deliverables identical to the ones the rebuild existed to replace
-- and acceptance would then have passed against the old data. `stale()` rebuilds when the output is
missing, empty, or older than its newest input.


### 12.4 Memory requests, and why they are a throughput decision

The QOS caps total memory per user (4,500 G here), not just the job count, so an over-request does not buy
safety -- it buys fewer concurrent tasks. Measured peak RSS over every prior task of each arm:

| Arm | tasks measured | mean RSS | peak RSS | requested | now |
|---|---|---|---|---|---|
| Illumina WGS | 362 | 10.3 G | 45.8 G | 48 G | 48 G, unchanged -- the tail is already close |
| PacBio HiFi WGS | 282 | 6.0 G | 11.6 G | 64 G | **24 G** |
| ONT WGS | 314 | 7.2 G | 39.6 G | 64 G | **48 G** |
| ONT bulk RNA | 3 | 25.7 G | 27.7 G | 320 G | **48 G** |
| ONT single-cell RNA | 3 | 23.1 G | 23.2 G | 320 G | **48 G** |
| bulk RNA | 3 | 2.7 G | 4.5 G | 320 G | **16 G** |
| 10x GEX | 4 | 7.8 G | 7.9 G | 320 G | **24 G** |
| Kinnex bulk | 2 | 8.0 G | 8.2 G | 320 G | **24 G** |
| Kinnex single-cell | 4 | 33.6 G | 46.6 G | 320 G | **64 G** |

The ONT RNA rows were written here before the scripts were changed, and for a while the table recorded an
intention rather than a fact: `41_ont_bulk_rna.sbatch` and `42_ont_sc_rna.sbatch` still asked for 320 G
after this section said they asked for 48 G. They now match. A table that describes the build is only
useful if it is checked against the build, so the numbers above were regenerated from `sacct` accounting
rather than copied forward.

The cost of getting this wrong was concrete twice in one day. The ONT bulk RNA job sat at
`QOSMaxMemoryPerUser` indefinitely asking for 320 G against a 29 G peak, and an ETA was reported for a job
that had never started. Then the ONT WGS array went two and a half hours without being given a single slot
while PacBio ran beside it, both asking 64 G, because memory was saturated at 4,496 of 4,500 G.

Size the request from the measured peak plus headroom for the tail, not from the mean and not from a round
number.


## 13. The chr1to6 subset

The subset exists so an end-to-end LENS run takes hours rather than days. Design §12 asked for it to be
derived from the full release "by alignment, not re-simulated"; it is derived from the full release by
**truth map**, and §12 now records that departure and the reasons. The short version: the per-chromosome
arms' subset is already defined by source locus, so defining the RNA arms' subset by where an aligner put
the reads would have left the subset release using two different meanings of "on chr1-6"; the truth map is
exact where an aligner is not; and the reads already carry the sidecars that make it free.

Two jobs and one script:

| Step | Runs | Produces |
|---|---|---|
| `jobs/90_subset_map.sbatch` | 2 tasks, ~1 min | `$IGI_WORK/subset/<ds>.chr1to6_records.tsv.gz` |
| `jobs/91_merge_chr1to6.sbatch` | 16 tasks | `<arm>_chr1to6/<ds>/<ds>_chr1to6_*` for WES, WGS, ONT WGS, PacBio |
| `jobs/90_subset_chr1to6.sbatch` | 12 tasks | `<arm>/<ds>/<ds>_chr1to6_<arm>*` for the six RNA arms |

### 13.1 The record map

`catalog/subset_chr1to6.py map` writes one row per transcript record: `rec`, `source`, the chromosomes it
resolves to, and the keep flag. Resolution is per source class -- GENCODE transcript for `reference` and
`splice_isoform`, the designed event's locus for `erv`/`cta`/`fusion`, and unconditional keep for `virus`.
A record that resolves to **nothing aborts the map**, because an unresolved source class would silently
shrink the subset; that check is why the map is trustworthy rather than merely produced.

Cross-checked against the truth tables' own `chr1to6` flags: 75 of 75 expressed events agree in both
datasets. Result: 442,380 of 1,226,495 records on chr1-6 for ds-01 (36.1 %) and 463,011 of 1,259,951 for
ds-02 (36.7 %).

### 13.2 Streaming, and the lockstep assertion

Every arm streams its deliverable and its truth sidecar **in lockstep and asserts the read names agree**,
read by read, rather than trusting that they do. The sidecars were written in emission order, but an arm
that was ever resumed or re-merged could have broken that, and a silent misalignment would mislabel every
read downstream. All reads and writes go through `pigz`; the gzip module is the bottleneck on a 49 GB
FASTQ.

### 13.3 Kinnex: which molecule is a segment?

An array is built `adapter_0 m_0 adapter_1 m_1 ... adapter_k`, so `skera split` returns one segment per
bracketed molecule and the segment's molecule index is **`min(dl, dr)`** of its two adapter tags. Taking
the minimum rather than `dl` is what makes it orientation-agnostic: about half the arrays are sequenced
reverse and come back with `dl > dr`. Verified on both libraries, both orientations and on short arrays.

Two things this surfaced, neither of which a naive implementation would have noticed:

- **4.2 % of ZMWs have fewer segments than molecules.** `skera` moved the rest to `non_passing.bam`.
  Harmless -- `min(dl, dr)` still identifies the molecule for the segments that remain, which is why
  counting segments in order would have been wrong.
- **0.019 % of segments carry an adapter index past the end of their array.** HiFi error on a 17 bp
  adapter occasionally makes adapter *j* look like *j+1*. Such a segment cannot be attributed, so it is
  **kept** and counted (`segments_unattributed`): keeping it adds a little background, while dropping it
  would remove reads nothing can account for, and an unexplained hole in a benchmark library is worse than
  an unexplained extra read.

### 13.4 Verification before launch

Each of the four code paths was run against truncated copies of the real deliverables and checked by an
independent reimplementation, not by eye:

| Arm | Checked | Result |
|---|---|---|
| bulk RNA | R1, R2 and readmap against the keep set, record by record | 38,074 of 100,000 pairs, identical order and content |
| Kinnex bulk | output BAM names against a separate recomputation | 76,996 of 199,995 segments, identical; `@RG`/`PU` preserved |
| Kinnex sc | molecule purity, barcodes | 85,342 of 199,995 segments, purity 0.494, 4,000 barcodes |
| ONT sc RNA | reads, molecules, barcodes | 69,170 of 200,000 reads, 4,000 barcodes |

At full scale the kept fractions are 34-38 % across the six arms. One number from the truncated runs did
not survive: 10x GEX retains **3,972 / 3,977 of 4,000** barcodes rather than all of them, because a few
cells have no chr1-6 molecule at all. §12 now says 99.3 % instead of "all".

### 13.5 Naming is a correctness requirement

`preflight_resolve_inputs` walks every search directory recursively and matches on the basename, so
`inputs/fastqs/<ds>/full` and `.../chr1to6` are one namespace. Every chr1to6 name therefore carries a
`chr1to6` label the full name lacks, and `make_lens_manifest.py` **fails** if any `File_Prefix` is a
prefix of another, within or across releases. 10x TCR is the single row whose two releases name the same
file, and it is linked under `full/` only.


## 14. The PacBio error floor, and the quality strings

F24 found that the delivered HiFi reads carried 5.1x the error their quality strings claimed. The follow-up
established two separate facts, and only one of them is fixable here.

### 14.1 Every long-read identity setting had been calibrated against the wrong thing

`jobs/measure/measure_longread.py` derives accuracy from the base-quality string, and from `rq` for a real
HiFi BAM. For a REAL read that is a fair measurement: the quality string is the instrument's own estimate.
For a SIMULATED read it is not a measurement at all. Badread builds sequence from an error model and
qualities from a separate qscore model, and nothing couples them.

`jobs/measure/alignment_identity.py` measures identity from NM over the aligned length, on reads simulated
from an unmodified reference window so the alignment is ground truth.
`jobs/81_calibrate_identity.sbatch` sweeps every setting.

ONT came out fine, and errs in the harmless direction:

| setting | request | true identity | target | verdict |
|---|---|---|---|---|
| `ont_wgs_identity` | 99.17 | ~0.9883 | 0.98676 | ~12 % too clean |
| `ont_identity` (cDNA) | 98.41 | ~0.9812 | 0.98220 | ~6 % too noisy |

For both, the quality strings UNDERSTATE accuracy (ratio 0.54-0.93), which does not touch the sequence.
No ONT change was made.

### 14.2 PacBio has an error floor, not a calibration offset

Across six requests the EXCESS over the nominal error rate is constant rather than proportional:

| request | nominal error | true error | excess | indel share |
|---|---|---|---|---|
| 99.4 | 6.0e-03 | 9.12e-03 | 3.12e-03 | 91 % |
| 99.7 | 3.0e-03 | 5.99e-03 | 2.99e-03 | 90 % |
| 99.85 | 1.5e-03 | 4.41e-03 | 2.91e-03 | 88 % |
| 99.95 | 0.5e-03 | 3.53e-03 | 3.03e-03 | 86 % |
| 99.99 | 0.1e-03 | 3.17e-03 | 3.07e-03 | 85 % |
| **100** | **0** | **3.05e-03** | **3.05e-03** | 84 % |

Asking for perfect reads still returns 3.05e-03/bp (Q25.2). Badread's `pacbio2021` model therefore cannot
produce real HiFi, whose error is 1.77e-03: the best it can do is **1.72x real**, and `--identity` has no
effect at the limit. `pacbio_identity` is set to the floor, `"100,100,0"`, which is 1.72x real instead of
the 5.15x that `99.4` delivered.

The floor is a per-base property, not a per-read one, which had to be checked rather than assumed: the
three arms driven by `pacbio_identity` deliver very different read lengths, and the Kinnex arms were
resubmitted at the new setting before this check returned. A 2 kb cDNA read floors in the same place as a
16.7 kb genomic one:

| case | read length | true error | true Q | true/claimed |
|---|---|---|---|---|
| genomic (`pacbio`) | 16.7 kb | 3.05e-03 | 25.15 | 5.152 |
| cDNA (`kinnex`) | 2.0 kb | **2.92e-03** | 25.34 | 4.989 |

4 % apart on a QA band of +/-0.0020, so one `pacbio_qual_error_scale` covers all three arms and no arm
needed re-queueing. Had the cDNA floor differed materially, the Kinnex libraries -- the longest-wall-time
arms in the build, at 20-44 h each -- would have had to be rebuilt a third time.

**Reaching real HiFi needs a different simulator, and that decision is not made here.** Note that the
pbsim3 + ccs route's figure of `rq 0.99754` came from the quality estimate too, by the same method §14.1
invalidates, so it is not established either; whichever route is tried next must have ITS true error
measured by alignment first. pbsim3's error models were also removed from this site as non-durable.

### 14.3 The quality strings were ours to fix, and the first fix did not work

The true/claimed error ratio is **5.15 at every request**, because the qscore model is keyed to the
REQUESTED identity rather than the realised one. So raising `--identity` improves the sequence and leaves
the claim exactly as optimistic as before -- reads at the floor would assert Q32.3 while carrying Q25.2.

That assertion is the F24 defect itself, and unlike the floor it is entirely in our hands, because we write
the BAM. `pacbio_bam.rescale_quality` multiplies each base's error probability by
`pacbio_qual_error_scale` before the read is written, which preserves the per-position structure of the
qscore model; `rq` then follows from the rescaled array rather than needing a correction of its own.

Setting that scale to the measured 5.15 was the obvious move and it was **not enough**. Measured on the
delivered ds-01 HiFi BAM after the first rebuild, the reads claimed 0.997802 -- with `rq` matching the
quality string exactly, so the wiring was right -- against a true 0.996948. Still 1.39x optimistic, not
1.00x.

**The correction saturates.** The error mass is concentrated on a handful of very bad bases:

| Phred | bases | share of claimed error mass |
|---|---|---|
| Q3 | 523 | 8.8 % |
| Q4 | 1,124 | 15.0 % |
| Q5 | 1,162 | 12.3 % |
| Q6 | 1,155 | 9.7 % |

**45.8 % of the claimed error mass sits on the 0.08 % of bases at Q3-Q7**, and those are exactly the bases
whose error probability cannot be multiplied by five without exceeding 1. Over 5.0 M bases:

| | mean implied error | identity |
|---|---|---|
| raw badread claim | 5.95e-04 | 0.999405 |
| x5.15, uncapped arithmetic | 3.07e-03 | 0.996935 |
| x5.15 as actually applied | **2.20e-03** | 0.997801 |

An effective 3.70x against the 5.15x asked for.

The scale that lands on the measured true error is **10.0**, solved numerically. Two independent quality
distributions -- 16.7 kb genomic and 2.0 kb cDNA -- solve to the same value, which is why one constant
still covers all three arms.

A second measured property decided the exact number. BAM Phred is integer and every input here is an
integer, so a constant scale shifts every base by the same whole number of decibels: the achievable
aggregate moves in ~1 dB steps and **no scale lands exactly on the target**. The two candidate plateaus:

| scale | genomic claimed/true | cDNA claimed/true |
|---|---|---|
| 7.5 - 8.91 | 0.93x (optimistic) | 0.95x |
| **9.0 - 11.0** | **1.05x (pessimistic)** | **1.07x** |

10.0 is taken from the pessimistic plateau, because a read understating its own accuracy is harmless where
overstating it is the defect being fixed -- the same rule already applied to the ONT arms in §14.1. Within
that plateau 10.0 is chosen because it is exactly 10 dB and so needs no rounding at all: the boundary
between the two plateaus sits at scale 8.9125, where a config float a few thousandths either way flips
every base by a whole Phred. Sitting a constant on that edge would have made the delivered qualities
depend on float rounding.

`rescale_quality` therefore applies the correction as the integer decibel shift it actually is, computed
once from the scale, rather than per base where it would sit on that edge. The rewrite was checked to be
bit-identical to the per-base form across the whole Phred range at every scale tried, so the only thing
that changed the output was the constant.

**Delivered result, measured on a rebuilt BAM rather than predicted:**

| | claimed error/bp | claimed identity | `rq` | vs true 3.05e-03 |
|---|---|---|---|---|
| scale 5.15 (first attempt) | 2.198e-03 | 0.997802 | 0.997802 | **0.72x -- optimistic** |
| scale 10.0, predicted | 3.210e-03 | 0.996790 | — | 1.05x |
| scale 10.0, **delivered** | **3.202e-03** | **0.996798** | 0.996801 | **1.049x -- pessimistic** |

Prediction and delivery agree to 0.3 %, and `rq` tracks the quality string to six decimals, so it is still
being derived from the rescaled array rather than carrying its own figure. The reads now assert slightly
more error than they carry, which is the direction that cannot mislead.

All three PacBio arms read `pacbio_identity` and all three are wired to the scale
(`run_pacbio_wgs.py`, `run_kinnex_bulk.py`, `run_kinnex_sc.py`), so all three were rebuilt again. They were
cancelled 1 h 20 m into their run to do it: the alternative was letting 20-44 h of Kinnex wall time
complete and then discarding it.

### 14.4 The QA accuracy target had to move, and why that is not moving the goalposts

`qa_release.py` derives accuracy from the quality strings, so with those calibrated the PacBio targets
become the simulator's own true accuracy, 0.99695, not real HiFi's 0.99823. The check now asks "do the
reads describe themselves correctly", which is answerable, instead of "did Badread reach real HiFi", which
it cannot. The gap is a simulator limit and belongs in this document, not in a check that fails on every
release.

The target stays at the measured true accuracy 0.99695 rather than moving to the 0.99679 the reads now
claim, because the point of the check is that the two agree. The ~1 dB quantization residue from §14.3 is
absorbed by the tolerance deliberately: delivered 0.996790 against target 0.99695 is 0.00016, well inside
the +/-0.0015 allowed for `pacbio` and the +/-0.0020 for the two Kinnex arms. If a future change puts a
delivered claim outside that band, the right response is to re-solve the scale, not to widen the band.

## 15. A VCF ALT column is not restricted to sequence

Reported from a LENS run: Cell Ranger refused the ds-01 10x GEX reads because they contained a literal
`*`. The report's diagnosis was exactly right, and the stray character turned out to be the mild half of
the defect.

`germline_edits` appended whatever allele the haplotype genotype selected straight into the sequence. The
HG002 Q100 VCF uses the spanning-deletion allele `*` at **157,318 of 5,945,526 records (2.65 %)**, always
as the second ALT of a multiallelic site. Keeping multiallelic records -- which was right, and which put
178,123 HG002 genotypes back into the sequence -- is what let these reach the builder. IPISRC044 has none,
because its VCF went through SHAPEIT5 normalisation, which is why ds-02 was clean and ds-01 was not.

### 15.1 Why the character was the smaller problem

`*` in VCF means "this allele is removed by a spanning deletion recorded elsewhere", so the correct action
is to apply **nothing**: the deletion's own edit already removes those bases. Instead `EditSet.apply`
replaced `len(REF)` reference bases with the single character `*`, so each affected site did two wrong
things at once:

```
 74,017  sites where a haplotype selected `*`
960,223  reference bases deleted across the two haplotypes, mean 13 bp, longest REF 63,066
  2.65%  of VCF records carry a `*` ALT; ds-02: zero
```

Delivered ds-01 reads carried the character at roughly 692 per million in bulk RNA, 7 in 10x GEX and 6 in
WGS, and 18 per 2 M records in the Kinnex bulk BAM. But nearly a megabase of spurious deletion per dataset
would read to a germline caller as **74,017 false-positive deletions absent from the truth VCF**. The
character blocked one tool loudly; the deletions would have quietly distorted every variant-calling result
on ds-01, which is the failure mode a truth set exists to prevent.

### 15.2 The fix, at one chokepoint

- `germline_edits` treats `*` as a no-op and counts them in `spanning_deletion_alleles`, so the skip is
  reported rather than silent.
- `EditSet.add` **rejects any allele outside `ACGTN`**. That is the single chokepoint every edit passes
  through, so no future source can reintroduce this. `rna.py` had spliced `.edits` in directly in four
  places, walking past the check; those now go through a validating `extend()`, which also carries the
  counter across.
- `qa_release.py` gained `check_alphabet`: every delivered base must be `ACGTN`, over both FASTQ and BAM.

The last point is the lesson worth keeping. The release already had checks on read length, accuracy, read
names, truth-map resolution and format conformance, and **not one of them asked whether the delivered
bases were bases**. A length check cannot see this; an accuracy check cannot either. The check was verified
both ways round -- it fires on the old ds-01 and passes on ds-02.

`EditSet.add` is also the model for §16.2: decide the constraint once, where everything must pass, rather
than at each call site.

### 15.3 Verification and blast radius

Verified on `chr1:750000-1000000`, the window carrying the first affected records: both haplotypes build
with zero non-base characters, with 2 and 4 spanning-deletion alleles correctly skipped, and the new guard
fires on a deliberately malformed allele.

All of ds-01 was regenerated except **10x TCR**, which is built from GENCODE V/J germline anchors rather
than the individual's genome; it was checked empirically for `*` across its whole library and had none.

**The reported blocker is resolved, and verified by complete count.** 10x GEX was the arm Cell Ranger
refused, so it is the one that had to be checked exhaustively rather than sampled -- Cell Ranger fails the
whole run on one character, so a sampled zero would have proved nothing:

```
IGI-SYN-SEQ-01-GEX_S1_L001_R1_001.fastq.gz   846 M   non-ACGTN over the FULL library: NONE
IGI-SYN-SEQ-01-GEX_S1_L001_R2_001.fastq.gz   4.8 G   non-ACGTN over the FULL library: NONE
```

Every read of all 140,980,551 pairs, both the barcode/UMI read and the cDNA read. The rebuild also wrote
**140,980,551** pairs against the pre-rebuild library's 140,980,551, so the seed derivation still
reproduces the same library and the rebuild changed the bases and nothing else. 4,000 cells, 7 h 10,
peak RSS 8.25 G against the 24 G this arm was reduced to.

**Confirmed on the worst-affected arm too.** ONT bulk RNA carried the highest rate of `*` of any deliverable,
4,147 non-ACGTN characters per 300,000 reads, which is about 13,800 per million. The rebuilt arm
(2026-10-08 18:15 to 23:42, 5 h 27, 49.8 GB) measures **zero across the full library** -- every read of
all 47 GB, not a sample.

Being exact about what that does and does not establish, since the two fixes landed hours apart: this build
contains the `*` fix (committed 16:40) but **not** the §17 ambiguity-code normalisation (22:19), because
Python loads its modules when the process starts. So the result validates the `*` fix and says nothing
about the other one. The arm's freedom from ambiguity codes has the same cause as 10x GEX's -- its reads
come from expressed transcripts, so one of the reference's 94 positions would have to fall inside an
expressed exon.

Its peak RSS was 20.6 G against the 320 G it was still requesting, which independently supports the 48 G
this arm was reduced to in §12.4.
ds-01's chr1to6 subsets and merges had to be redone as well, because they were derived from the corrupt
reads. `70_stage_inputs.sh` now holds back **any** chr1to6 deliverable older than what it was derived
from, which was not a hypothetical precaution: this fix invalidated every ds-01 library while that
dataset's chr1to6 subsets were still being written from the old reads.


## 16. A memory request does not bound GNU sort

§12.4 right-sized every arm from its measured peak. One of those cuts then failed: `igi_rna` task 1 for
ds-01 was OOM-killed at 16 G on 2026-10-08, and not inside the Python builder but inside GNU `sort`:

```
subprocess.CalledProcessError: Command '['sort', '-T', ..., '-k1,1', '-o', sorted.tsv, keyed.tsv]'
    died with <Signals.SIGKILL: 9>
slurmstepd: error: Detected 1 oom_kill event in StepId=11577347.batch
```

### 16.1 Why the request was never the bound

`readnames.shuffle_and_rename` keys every read by a seeded digest, sorts on disk, and renames in the sorted
order. Its docstring said "the peak memory is the sort buffer rather than the library". That was true as
written and false in effect, because nothing bounded the sort buffer.

GNU sort chooses its default buffer from the machine's **physical** memory. It does not read the cgroup
limit that a SLURM step runs under. On this cluster's 503 G nodes it will therefore try to hold the whole
input, whatever `--mem` said:

| | input | peak RSS | wall |
|---|---|---|---|
| `sort` (default buffer) | 0.93 G | **1003 M** -- the entire input | 6.15 s |
| `sort -S 128M` | 0.93 G | **134 M** | 7.44 s |

So the buffer tracked library size, and the 320 G request had not bounded anything either -- it merely
happened to exceed whatever sort chose. The accounting shows the same thing at release scale: the
full-library RNA arm peaked at **59.8 G**, against a keyed intermediate of about 51 G for 80 M pairs
(~638 B per record). The arm was not using 59.8 G to simulate reads; it was using it to sort them.

This also reframes §12.4. Four arms route through this sort -- bulk RNA, exome, Illumina WGS and all three
ONT arms via `longread.py` -- so their measured peaks were substantially the sort's appetite for the
library rather than the cost of simulation. Illumina WGS peaking at 45.8 G against a 48 G request was the
clearest symptom: that is a per-chromosome keyed file, not a working set.

### 16.2 The fix

Derive the buffer from the step's own allocation, so a step's footprint is a property of the job and not of
whichever node it lands on. A quarter of the allocation, floored at 256 M and capped at 8 G, leaves room
for the Python process, the gzip writers and sort's per-thread overhead, which all sit alongside the
buffer. Outside SLURM the floor applies -- slower, but it cannot be killed.

The same rule is written twice, because both a Python and a shell caller need it, and each is marked as
the other's counterpart:

- `sort_buffer_mb` / `sort_buffer_arg` in `catalog/igi_catalog/readnames.py`, spliced into the `sort` call.
- `sort_buffer_mb` / `sort_buffer_arg` in `jobs/lib.sh`, a new tracked, site-independent companion to the
  untracked `env.sh`, sourced by the three merge scripts whose read-name uniqueness proof also ran an
  unbounded `sort -u` over a whole library (`11_merge`, `50_merge_longread`, `91_merge_chr1to6`). Those had
  not failed yet; they were the same latent bug one node size away.

`EditSet.add` in §15.2 was a single chokepoint that every edit had to pass. This is the
same shape of fix: the buffer is decided in one place per language rather than at each call site.

### 16.3 Confirmed on the workload that failed

The rerun of ds-01 bulk RNA sorted the real file inside the allocation that had killed it:

```
keyed.tsv    48.4 GB     written 23:36   (the intermediate the unbounded sort wanted in RAM)
sorted.tsv   48.4 GB     written 23:44   (complete -- the sort took about 8 minutes)
peak RSS      8.84 G     against a 16 G allocation
```

So a 48.4 GB sort now runs with a step footprint of 8.84 G -- the 4 G buffer plus the Python process and
the gzip writers alongside it -- where before the buffer alone grew toward the size of the input and the
step was OOM-killed at 1 h 44. The 48.4 GB measured here is close to the ~51 GB predicted from 80 M pairs
at ~638 B per record, and it is the number that explains the old 59.8 G release-scale peak in §16.1.

The external merge is not the bottleneck anyone might fear: eight minutes for 48.4 GB across roughly a
dozen runs, against a 21 % penalty measured on the small case.

And at this scale the penalty is immaterial, because the sort was never the expensive part. The ds-02 arm
built with the *unbounded* sort took **6.89 h** end to end (24,812 s, 76.5 M pairs, 5.14 GB R1), and ds-01
is tracking the same shape: 27 min of record assembly, about 65 min of ART, **8 min of sorting**, and the
rest writing three gzip streams. The sort is roughly 2 % of the runtime, so a 21 % penalty on it is a
rounding error -- which is worth stating, because "bound the buffer and it will get slower" is the obvious
objection to §16.2 and the measurement says it does not matter here.

The write-out is what dominates, and that is a separate, untouched opportunity: `shuffle_and_rename`
writes R1, R2 and the read map through Python's `gzip` at default compression, single-threaded. Nothing
about the sort fix changed it, and it has not been optimised.

### 16.4 The fix changes no output, so nothing needs rebuilding

Buffer size could in principle change the order of records whose keys tie, which would change which read
receives which Illumina name. It does not. GNU sort's last-resort comparison falls back to the whole line
when the keys compare equal, and whole lines here are unique, so the order is total. Checked against a
400,000-record file deliberately collapsed onto 2,000 distinct keys -- about 200 ties per key:

```
-S 1M  -S 8M  -S 64M  -S 1024M  unbounded   ->  all md5 1ff0aeb126149fa19ae1e8a35c766514
```

Byte-identical at every buffer size. Libraries already built with the unbounded sort are therefore valid,
and only the OOM-killed ds-01 bulk RNA arm had to be rerun. Correctness of the rename itself was
re-verified independently: names unique, R1 and R2 in step, the read-name map aligned to the FASTQ order,
and every source record recovered.

### 16.5 What is still open

`run_rna.py` has no resume path -- it re-assembles records and re-simulates every bin on each invocation.
Rerunning the killed task therefore cost the full ~1 h 45 m rather than just the sort. That is a deliberate
trade for now: the alternative is caching 50 G intermediates, and this build has twice shipped or graded
stale data when a step trusted an output that merely existed (see the note at the top of `11_merge.sbatch`).

The post-fix peaks have not been re-measured. Every arm's request above was sized against a sort that was
free to take the library, so several are now larger than they need to be -- Illumina WGS most of all. They
should be re-measured once the current rebuild drains, and this table revised against that, not against
the figures here.

## 17. IUPAC ambiguity codes in delivered reads

Found while verifying that §15's alphabet check passes on the rebuilt ds-01 arms. It does -- rebuilt WES
and ONT WGS carry no `*` -- but the check surfaced a second, unrelated population of non-ACGTN characters
that had been in the delivered data all along.

Delivered ONT WGS reads carry IUPAC ambiguity codes, roughly 90 characters per 20,000 reads, in both
tumour and normal and across several chromosomes:

```
ds-01 chr10 normal   24 R  22 W  17 Y  13 M  12 K
ds-01 chr10 tumour   17 R  14 W  13 Y   6 M   6 K
ds-01 chr13 tumour   11 Y   5 R   5 K   3 M
ds-01 chr17 normal    8 R   7 Y
ds-02 chr1  tumour    6 Y   6 K   3 R   1 M
```

These are not ours. They are in the references:

| source | ambiguity codes |
|---|---|
| `reference_fasta` (GRCh38, what the builders read), chr1-chrY | **94** |
| of which chr10 | 36 -- R 13, Y 8, W 6, K 4, M 3, S 1, B 1 |
| chr17 12 · chr2 9 · chr3 7 · chr22 5 · chrX 5 · chr7 4 | chr9/12/13/21 3 · chr1 2 · chr6/chr16 1 |
| chr4, chr5, chr8, chr11, chr14, chr15, chr18-chr20, chrY | none |
| `virus_unmasked.02ec8.fa` | **497** -- Y 190, R 162, W 42, M 35, S 34, K 34 |

Note that this is `reference_fasta`, the full assembly the builders read via `Genome`, not the `noalt`
build used for alignment. The first draft of this section counted the noalt file; the primary-chromosome
totals happen to agree at 94, but they are different files and the one that matters is the one the
generator opens.

That both the tumour and the normal library carry them in proportion to coverage identifies the reference,
rather than a designed event, as the source, and the counts line up: chr10 holds 36 of the 94 positions and
a 20,000-read sample covers that chromosome about 2.8 times, which is the 88 characters observed.

**All of it is the reference, and none of it is the viral set.** This needed checking rather than
assuming, because the characters in the reads do not match the reference's at first glance -- chr13's reads
carry `R` and `M` where chr13's three positions are 2 Y and 1 K. The explanation is strand:
complement(Y) = R and complement(K) = M, and Badread emits minus-strand reads. Grouping both the reference
and the reads into complement classes resolves it exactly:

| | reference classes | observed in reads | classes with no reference source |
|---|---|---|---|
| chr3 | RY 4, W 2, BV 1 | RY 14, W 16, BV 17 | none |
| chr10 | RY 21, W 6, KM 7, S 1, BV 1 | RY 41, W 22, KM 25 | none |

So the viral reference's 497 codes are a *latent* source rather than a contributing one: nothing observed
in delivered reads requires them. They still matter, because `viral_records` feeds viral sequence into
expressed transcripts, and they are normalised along with everything else.

They show up in ONT WGS and not in the short-read arms for a simple reason: a 19 kb read has a few hundred
times the chance of spanning one of 94 isolated positions than a 150 bp read does, and the exome arms only
cover captured exons.

### 17.1 Why this had to be fixed even though nothing had failed

Appearing in only one assay is a property of this seed, not a safety margin:

- **No instrument emits them.** A real sequencer reading an ambiguous locus reports a definite base or N;
  the ambiguity lives in the reference's knowledge, not in the molecule. A synthetic read carrying `Y` is
  not realistic however faithful it is to the FASTA.
- **Strict consumers refuse them**, and that is exactly how §15 surfaced: Cell Ranger rejects any base
  outside ACGTN, and it fails the whole run on one character rather than dropping the read. A full scan of
  the ONT WGS and exome arms found that short reads are *not* immune: `IGI-SYN-SEQ-01_chr21_tumor_R1`
  carries one `M` in an 80,000-read sample, from one of chr21's three positions falling in a captured
  exon. A sampled zero therefore proves nothing about a library -- 300,000 reads of 10x GEX measured zero,
  but the delivered library is ~140 M reads, and the relevant question is whether it contains *one*
  character, not what its rate is. `rna_assembly.viral_records` also puts raw viral sequence into
  *expressed* transcript records, and the viral reference holds 497 codes.
- **`revcomp` does not complement them.** `genome.COMP` maps only `ACGTNacgtn`, so an ambiguity code
  survives reverse-complementing **unchanged** -- `revcomp("ACGTRY")` returns `"YRACGT"`, where `Y` should
  become `R` and `R` become `Y`. This is a defect in the helper rather than an observed corruption, and the
  distinction is worth keeping straight. The paths that would reach it are minus-strand transcripts
  (`transcriptome.py:33`), inverted SV segments (`rearrange.py:148`) and fusion partners, all of which feed
  delivered reads. But the one strand-flipped case actually found in the data runs the other way: chr3's
  single `B` appears as `V` on minus-strand reads, which is the *correct* complement, produced by Badread's
  own strand handling. So the latent defect is real and the delivered evidence for it is not; normalising
  at `Genome.seq` removes the class either way, because `revcomp` then only ever sees ACGTN.

  Tracing that `V` is also what caught an error in this section's first draft. `V` appears in no reference
  -- not GRCh38, not the viral set -- so it looked like a third source. It is chr3's `B` reverse
  complemented, and the two reads carrying it sit at one locus with identical flanks.

### 17.2 N, not a definite base

They become `N`. Resolving each code to one of its constituent bases would be more realistic, but it would
put a definite base at a position no truth table records, and a caller would report it as a variant absent
from the truth VCF. That is the §15 mistake in miniature -- writing something into the sequence that the
truth set does not account for. Callers skip N.

`genome.normalize_bases` does the translation, applied at `Genome.seq`, which is the single chokepoint for
reference access, plus the two viral fetches that do not go through it (`junctions.viral_junctions` and
`rna_assembly.viral_records`). Verified against the real reference: chr1 and chr3 both go from the counts
above to clean. This is the §15.2 shape of fix again -- one place that everything passes.

### 17.3 Why QA fails on any non-ACGTN character, in every arm

**This section previously argued for a WARN tier, and that argument was wrong. It is kept here, as the
reasoning it was, because the way it failed is the useful part.**

The original rule made the verdict depend on who reads the file:

| characters | arm | verdict (ORIGINAL, SUPERSEDED) |
|---|---|---|
| `*`, `<DEL>`, breakend notation | any | **FAIL** -- a symbolic allele; a generator defect |
| `R Y S W K M B D H V` | 10x GEX, 10x TCR | **FAIL** -- Cell Ranger refuses the run on one character |
| `R Y S W K M B D H V` | everything else | **WARN** -- reference ambiguity; clears on next rebuild |

The justification was that a FAIL nobody acts on teaches people to ignore failures, and that for the DNA
arms a FAIL would be exactly that, since "minimap2, pbmm2 and bwa all tolerate ambiguity codes, the rate
is ~2e-7 per base, and the arms clear whenever each is next built".

Every clause of that is true. The conclusion still did not hold, because the premise was **an enumeration
of the consumers I happened to know**, and a reference dataset does not get to decide what reads it. It
was falsified within a day, by §17.6: OptiType's `razers3` aborts while loading its first chunk --

```
seqan::ParseError: Unexpected character 'Y' found.
```

-- because SeqAn's reader accepts only A, C, G, T and N. It killed HLA typing on the ds-01 normal exome
through all seven retries, in the arm the table above had classified as tolerant, while every aligner in
the same run read the file without complaint.

The rule is therefore unconditional:

| characters | arm | verdict |
|---|---|---|
| anything outside `ACGTN` | **any** | **FAIL**, naming which of the two causes it is |

There is no WARN tier. The cost of the absolute guarantee is low -- `genome.normalize_bases` removes the
codes at source, so the only arms that can fail are those built before it, which need rebuilding anyway --
and the cost of the conditional one was a user's pipeline failing seven times on data we had checked and
passed. Tolerating a defect because the tools *we* tested survive it is how both §15 and §17 reached
someone else's run instead of being caught here.

### 17.4 What still carries them

Measured, not inferred. A scan of every ONT WGS and exome deliverable in both datasets (80,000-read sample
per file) found:

| arm | files carrying codes | note |
|---|---|---|
| ds-01 ONT WGS | 17 | chr2, 3, 6, 7, 9, 10, 12, 13, 16, 17, 21, 22, X |
| ds-02 ONT WGS | 19 | same mechanism, different seed |
| ds-01 exome | 1 | `chr21_tumor_R1`, a single `M` in 80,000 reads |
| ds-02 exome | 0 | at this sample size |

### 17.5 The normalisation works, and ds-01's arms are a per-task mixture

An accident of timing produced a clean natural experiment. The normalisation was committed at 22:19:57,
every job of the 2026-10-08 rebuild was submitted before it, and Python imports its modules when a process
starts -- so within a 48-task array, tasks scheduled before that minute ran the old code and tasks
scheduled after ran the new one. For the ds-01 PacBio array, **28 tasks started before and 20 after**.

chr10 appears twice in that array, once per library, on either side of the line. Same chromosome, same 36
reference positions, same reference file; the only difference is which code the process loaded:

```
chr10 tumour  (task 10, started 21:57:43, old code)   73 Y  50 R  36 W  31 M  27 K  6 S  2 V  2 B
chr10 normal  (task 34, started 22:26:46, new code)   none
```

227 characters against zero. That is the normalisation verified on delivered data rather than on a unit
test, and the character counts are themselves consistent with §17's complement-class argument -- Y/R at
73/50, M/K at 31/27, W and S self-complementing, and the B/V pair at 2/2.

It also means **no arm of ds-01 is uniformly normalised**, and which tasks are depends on when SLURM
happened to schedule them. When this was written, §17.3 made that a WARN rather than a defect outside the
10x arms, and the conclusion drawn here was that uniform normalisation "is not worth doing on its own
account". **§17.6 overturns that.** The mixture is a defect in every arm that carries a code, the dirty
units have to be rebuilt, and §17.6.3 records which they are and how they were identified -- by complete
measurement, not by this inference.

Illumina WGS and HiFi WGS were not scanned exhaustively; by the mechanism they carry codes at a rate
between the exome's and ONT's, since they are whole-genome but short or shorter.

The exome hit is the one that matters for the decision in §17.3, because it shows short reads are not
immune and therefore that a sampled zero on a 10x arm is not evidence of a clean library.

So the 10x arms were checked properly rather than sampled. **Both ds-02 GEX R2 libraries are clean over
their entire length** -- every read of the ~140 M-read full release and the ~50 M-read chr1to6 subset, zero
non-ACGTN characters:

```
IGI-SYN-SEQ-02-GEX_S1_L001_R2_001.fastq.gz           non-ACGTN: NONE   (full library, 5.0 GB)
IGI-SYN-SEQ-02-chr1to6-GEX_S1_L001_R2_001.fastq.gz   non-ACGTN: NONE   (full library, 1.9 GB)
```

Cell Ranger will therefore run on ds-02 as delivered, and §17.3's FAIL rule for the 10x arms is a guard
against a future build rather than a description of a current failure. The reason the exome can carry a
code while GEX does not is that a GEX read has to come from an *expressed* transcript, so one of the 94
positions must fall inside an expressed exon, where exome capture only requires it to fall inside a
captured one.

### 17.6 razers3, and three defects in the check that was supposed to catch this

Reported from downstream on 2026-10-09: `optitype_razers3` failed on `nd-wes-normal`, all seven attempts,
exit 1, aborting while loading its first 10M-read chunk with `seqan::ParseError: Unexpected character 'Y'
found`. The diagnosis that came with it was correct in every particular -- the codes come from the
reference rather than the VCF, `Genome.seq` returned the fetched sequence unchanged, the `*` fix's
`EditSet` guard validates VCF alleles and not reference bases, and the proposed fix was to translate
non-ACGTN to N at the `fetch(...).upper()` return in `genome.py`.

That is the fix that had already been committed, at that exact line, in `6ad1796` -- the report was
written against the pre-fix file. The useful content of the report was therefore not the fix but the
**evidence that two of my judgement calls were wrong**, and a LENS-side `tr` to N was correctly declined
as masking a data problem.

#### 17.6.1 The verdict was wrong

See §17.3. razers3 is the consumer the WARN tier assumed did not exist.

#### 17.6.2 The check sampled, and never read R2

`check_alphabet` had been reading 400,000 lines -- 100,000 reads -- off the front of each file. The report
put the exome rate at "one Y and one M in the first 2M reads", which is the whole problem: only 94
reference positions are ambiguous, a targeted assay covers few of them, and the resulting rate is about
**one offending read per million**. Measured on the library that broke razers3:

```
IGI-SYN-SEQ-01_chr10_normal_R1.fastq.gz   1,860,880 reads
  complete scan : R x1
  100k sample   : none          <-- what the check reported
```

A sample can witness contamination. It cannot establish absence at a rate below its own resolution, and
absence is precisely what every consumer of these files depends on. The check now counts completely.

Two further gaps found while fixing that one:

- **R2 was never examined.** The globs ended at `*_R1.fastq.gz` and `*_R1_001.fastq.gz`, so half of every
  paired library went unread -- including the 10x R2 that carries the cDNA and is the only mate Cell
  Ranger aligns. The earlier full-library 10x verifications in §17.5 happened to be done by hand on R2;
  the automated check could not have reproduced them.
- **A read failure was indistinguishable from a clean library.** The complete scan runs without `head` in
  the pipeline so `pipefail` can stay on, and a decompression or `samtools` error is now reported as
  `FAIL (could not read the file)` rather than as zero offences.

Completeness had to be made affordable or it would simply be skipped, which is how a check stops being a
check. The scan's inner loop tests `/[^ACGTNacgtn]/` before calling `gsub`, because `gsub` rebuilds the
string on every read it touches:

| inner loop | one exome library (1.86 M reads) | result |
|---|---|---|
| `gsub` on every read | 105 s | `R x1` |
| regex pre-test first | **11 s** | `R x1` |

I found this the honest way: my first pass over the 3.1 GB reference used `gsub` unconditionally and ran
for over ten minutes before I killed it; with the pre-test it is 32 s.

#### 17.6.3 Which units actually need rebuilding, by measurement

The 94 ambiguous positions are not spread evenly, so the repair can be targeted. Counted fresh from the
build reference (`Homo_sapiens.assembly38.fa`, the FASTA the generators read -- not the `noalt` build):

| chrom | n | breakdown | | chrom | n | breakdown |
|---|---|---|---|---|---|---|
| chr10 | 36 | R13 Y8 W6 K4 M3 B1 S1 | | chr12 | 3 | Y2 M1 |
| chr17 | 12 | Y5 R3 K2 S1 W1 | | chr13 | 3 | Y2 K1 |
| chr2 | 9 | Y4 W2 K1 M1 R1 | | chr21 | 3 | M2 R1 |
| chr3 | 7 | Y3 W2 B1 R1 | | chr9 | 3 | Y2 R1 |
| chr22 | 5 | R2 Y2 W1 | | chr1 | 2 | M1 R1 |
| chrX | 5 | Y2 R1 S1 W1 | | chr16 | 1 | R1 |
| chr7 | 4 | Y2 R1 S1 | | chr6 | 1 | Y1 |

**No ambiguous position exists on chr4, chr5, chr8, chr11, chr14, chr15, chr18, chr19, chr20 or chrY.**
The DNA arms are generated per chromosome, and each task reads only its own chromosome, so those ten
chromosomes cannot carry a code and their tasks do not need re-running -- 10 of 24 chromosomes, which is
42 % of the per-chromosome tasks in every DNA arm.

That is an inference, though, and §17.6.2 is a lesson about trusting inference over measurement. So the
rebuild set is decided by `jobs/92_scan_alphabet.sbatch`, which counts every base of all 612 delivered
files (~1.3 TB) 96-wide in about 40 minutes, against ~33 h serially. It scans the **per-chromosome units
rather than the merged libraries**, because those units are exactly the array tasks that build them: a
dirty file names the task to re-run. It also catches anything the chromosome inference would miss, such as
the viral reference's 497 codes reaching an RNA arm, where they do not arrive by chromosome at all.

## 18. Regression tests for the guards

`catalog/tests/test_guards.py` pins down the guards added on 2026-10-08 and 2026-10-09. It needs no
cluster and no data, takes 0.2 s, and runs either as `python3 catalog/tests/test_guards.py` or under
pytest. **10 tests.**

Each guard is tested because its entire value is doing what it claims, and this build has three times now
shipped or graded data behind a check that did not:

| test | what would otherwise go unnoticed |
|---|---|
| `check_alphabet` verdicts | a `*` passing, or an ambiguity code being tolerated in any arm -- the WARN tier razers3 falsified (§17.3) |
| the scan is complete, not sampled | the check reverting to a sample, which reported `none` on the library that aborted razers3 (§17.6.2) |
| both mates are examined | R2 going unread, as it did -- the only mate Cell Ranger aligns |
| `sort_buffer_mb` | the buffer drifting back to something the allocation does not bound (§16) |
| `normalize_bases` | an ambiguity code dropped from the translation table (§17) |
| `revcomp` after normalisation | the uncomplemented-ambiguity path reopening |
| `rescale_quality` is an integer dB shift | the per-base form creeping back, which sits on the knife edge at scale 8.9125 (§14.3) |
| `rescale_quality` saturates | the reasoning for 10.0 over the measured 5.15 being lost |

The suite is checked by mutation rather than only by running it green, because a test that cannot fail is
worth nothing. Each line below is a real edit to the source, the suite run, and the edit reverted:

```
baseline                              10/10
the scan samples by default again      9/10   <- the §17.6 defect, reintroduced
the globs stop at R1 again             9/10
the WARN tier restored for DNA arms    7/10
sort-buffer cap removed                9/10
normalisation table drops V and B      9/10
Phred shift forced to 0                8/10
```

The WARN mutation breaking three tests rather than one is the useful signal: the verdict is asserted per
arm, in aggregate, and in the wording that names razers3, so restoring the tolerant rule cannot pass by
satisfying one of them.

`revcomp("ACGTRY") == "YRACGT"` is asserted deliberately, documenting the untreated behaviour rather than
the desired one, so that anyone who fixes `COMP` properly is told that this test encodes the old contract.
