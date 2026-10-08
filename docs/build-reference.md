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

## 9. Review findings of 2026-10-05 and what each one changed

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

### 9.1 The R1 quality target was R1, not R2

The review compared 10x R1's 663 distinct quality strings against R2's 19,461 and read the gap as a
defect in R1. Half right. Measuring the real IPISRC044 library gives **2,819** distinct strings in 20,000
reads over the first 26 cycles, because real R1 is only 26 four-level-binned cycles (`#*9I`, so
Q2/Q9/Q24/Q40) and genuinely has low variety. R2's 19,461 was never the right target. The per-cycle model
now gives 2,947 against that real 2,819, and the old 663 was indeed wrong.

### 9.2 Why the roster did not have to be rebuilt for B7

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

### 9.3 Fitted resources added

| File | Measured from | Holds |
|---|---|---|
| `resources/pacbio_np_model.json` | 300k reads of `HG002_PacBio-Revio_m84039_230928_213653_s3.hifi_reads.bam`; 200k of the Kinnex `segmented.bam` | empirical `np` CDF and `ec/np` ratio, separately for HiFi WGS and Kinnex |
| `resources/tenx_r1_quality.json` | 500k R1 reads each of `DSCOLAB_IPISRC044_T3_SCG1` and `..._T3_TCR1` | per-cycle quality distributions for GEX and TCR |
| `resources/vdj_germline.json` | GENCODE v37 TR segments over GRCh38 | per V the CDR3 nucleotides from the conserved Cys, per J those to the conserved Phe, plus framework |

### 9.4 Verification

The repaired HiFi BAM was put through the PacBio toolchain rather than inspected by eye: `pbindex` indexes
it (949 of 949 reads), `pbindexdump` reads the index, `extracthifi` keeps 949 of 949, `zmwfilter` parses
the hole numbers and `bam2fastq` round-trips it. `skera split` on a repaired Kinnex BAM produces
`m84000_260101_000000_s0/1/ccs/17_1856` with the segment coordinates unchanged from the original run, so
the arrays still line up with their truth rows.

---

## 10. Focal copy-number events, and why they were missing

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

### 10.1 The fix

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

### 10.2 Verification

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

### 10.3 What this cost

Rebuilding all three WGS arms is 504 task-hours, measured from the previous runs rather than estimated:
48 Illumina tumour tasks at 4.23 h, 48 normal at 1.30 h, 96 PacBio at 1.41 h, 96 ONT at 1.08 h. Wall
clock is set by how much of the 105-job association limit is free, not by the work: 29.7 h at 17
concurrent slots, 7.2 h at 70.

`50_merge_longread.sbatch` had to be fixed before any of it could land. All three merges were guarded with
`[ -s "$out" ] ||`, so with the previous merged libraries in place the merge would have skipped all twelve
tasks and reported success, leaving the deliverables identical to the ones the rebuild existed to replace
-- and acceptance would then have passed against the old data. `stale()` rebuilds when the output is
missing, empty, or older than its newest input.


### 10.4 Memory requests, and why they are a throughput decision

The QOS caps total memory per user (4,500 G here), not just the job count, so an over-request does not buy
safety -- it buys fewer concurrent tasks. Measured peak RSS over every prior task of each arm:

| Arm | tasks measured | mean RSS | peak RSS | requested | now |
|---|---|---|---|---|---|
| Illumina WGS | 244 | 9.1 G | 45.5 G | 48 G | 48 G, unchanged -- the tail is already close |
| PacBio HiFi WGS | 120 | 5.4 G | 11.6 G | 64 G | **24 G** |
| ONT WGS | 158 | 6.0 G | 39.6 G | 64 G | **48 G** |
| ONT bulk RNA | 2 | — | 29.0 G | 320 G | **48 G** |

The cost of getting this wrong was concrete twice in one day. The ONT bulk RNA job sat at
`QOSMaxMemoryPerUser` indefinitely asking for 320 G against a 29 G peak, and an ETA was reported for a job
that had never started. Then the ONT WGS array went two and a half hours without being given a single slot
while PacBio ran beside it, both asking 64 G, because memory was saturated at 4,496 of 4,500 G.

Size the request from the measured peak plus headroom for the tail, not from the mean and not from a round
number.


## 11. The chr1to6 subset

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

### 11.1 The record map

`catalog/subset_chr1to6.py map` writes one row per transcript record: `rec`, `source`, the chromosomes it
resolves to, and the keep flag. Resolution is per source class -- GENCODE transcript for `reference` and
`splice_isoform`, the designed event's locus for `erv`/`cta`/`fusion`, and unconditional keep for `virus`.
A record that resolves to **nothing aborts the map**, because an unresolved source class would silently
shrink the subset; that check is why the map is trustworthy rather than merely produced.

Cross-checked against the truth tables' own `chr1to6` flags: 75 of 75 expressed events agree in both
datasets. Result: 442,380 of 1,226,495 records on chr1-6 for ds-01 (36.1 %) and 463,011 of 1,259,951 for
ds-02 (36.7 %).

### 11.2 Streaming, and the lockstep assertion

Every arm streams its deliverable and its truth sidecar **in lockstep and asserts the read names agree**,
read by read, rather than trusting that they do. The sidecars were written in emission order, but an arm
that was ever resumed or re-merged could have broken that, and a silent misalignment would mislabel every
read downstream. All reads and writes go through `pigz`; the gzip module is the bottleneck on a 49 GB
FASTQ.

### 11.3 Kinnex: which molecule is a segment?

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

### 11.4 Verification before launch

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

### 11.5 Naming is a correctness requirement

`preflight_resolve_inputs` walks every search directory recursively and matches on the basename, so
`inputs/fastqs/<ds>/full` and `.../chr1to6` are one namespace. Every chr1to6 name therefore carries a
`chr1to6` label the full name lacks, and `make_lens_manifest.py` **fails** if any `File_Prefix` is a
prefix of another, within or across releases. 10x TCR is the single row whose two releases name the same
file, and it is linked under `full/` only.


## 12. The PacBio error floor, and the quality strings

F24 found that the delivered HiFi reads carried 5.1x the error their quality strings claimed. The follow-up
established two separate facts, and only one of them is fixable here.

### 12.1 Every long-read identity setting had been calibrated against the wrong thing

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

### 12.2 PacBio has an error floor, not a calibration offset

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

**Reaching real HiFi needs a different simulator, and that decision is not made here.** Note that the
pbsim3 + ccs route's figure of `rq 0.99754` came from the quality estimate too, by the same method §12.1
invalidates, so it is not established either; whichever route is tried next must have ITS true error
measured by alignment first. pbsim3's error models were also removed from this site as non-durable.

### 12.3 The quality strings were ours to fix, and are fixed

The true/claimed error ratio is **5.15 at every request**, because the qscore model is keyed to the
REQUESTED identity rather than the realised one. So raising `--identity` improves the sequence and leaves
the claim exactly as optimistic as before -- reads at the floor would assert Q32.3 while carrying Q25.2.

That assertion is the F24 defect itself, and unlike the floor it is entirely in our hands, because we write
the BAM. `pacbio_bam.rescale_quality` multiplies each base's error probability by
`pacbio_qual_error_scale` (5.15, the measured ratio) before the read is written, which preserves the
per-position structure of the qscore model while making the aggregate honest; `rq` then follows from the
rescaled array rather than needing a correction of its own. Phred is clamped to [1, 93] -- 0 means "no
quality available" in SAM. Verified: a badread-shaped quality string claiming 0.99939 becomes 0.99695,
which is the measured true accuracy at the floor.

All three PacBio arms read `pacbio_identity`, so all three are affected and all three were wired to the
scale: `run_pacbio_wgs.py`, `run_kinnex_bulk.py`, `run_kinnex_sc.py`.

### 12.4 The QA accuracy target had to move, and why that is not moving the goalposts

`qa_release.py` derives accuracy from the quality strings, so with those calibrated the PacBio targets
become the simulator's own true accuracy, 0.99695, not real HiFi's 0.99823. The check now asks "do the
reads describe themselves correctly", which is answerable, instead of "did Badread reach real HiFi", which
it cannot. The gap is a simulator limit and belongs in this document, not in a check that fails on every
release.
