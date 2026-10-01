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
| 8 | tumor Kinnex bulk RNA | not written | — |
| 9 | tumor Kinnex scRNA (MAS 16-mer) | not written | — |
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

`run_pacbio_wgs.py` is the only long-read driver written. Shape:

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
| `pacbio_identity` | `99.4,99.8,0.5` (targets 0.99823) | same | same job, then apply the +0.43 offset |
| `kinnex_segment_length_mean/sd` | 2224 / 933 | 200k post-skera segmented reads | same job |
| `ont_length_mean/sd`, `ont_identity` | 910 / 461, `98.22,99.6,1.5` | IPISRC044 R10.4.1 SUP scRNA, 2025-07-17 | `measure_longread.py`; no ONT arm wired into the job yet |
| `ont_wgs_length_mean/sd`, `..._identity_target` | 19117 / 15530, 0.98676 | ONT open-data `giab_2025.01/basecalling/sup/HG002/PAW70337` | `measure_longread.py`; the BAM is remote, see `ont_wgs_bam` |
| expression baseline | per-transcript median TPM | 191 TCGA-BRCA basal-like tumours, UCSC Xena Toil recompute | not yet tracked; the original job is in scratch |
| `background_mut_per_mb`, `signatures` | 1.0 mut/Mb; SBS3 0.6 / SBS1 0.15 / SBS5 0.15 / SBS13 0.1 | COSMIC v3.4, TNBC-typical | edit `design.yaml` |

Two notes on these:

- The ART profile comes from WGS, not WES, because IPISRC044 has no 150 bp exome (BostonGene exomes are
  101 bp normal / 140 bp tumour). A per-cycle quality profile is a property of the run, not the capture,
  and capture is modelled separately by the GC curve.
- GIAB's own HG002 ONT is all R9.4-era (2D reads from 2016, Guppy v2/v3) and unusable as an R10.4.1
  target, which is why the ONT WGS parameters come from ONT's open-data bucket instead.
- `ont_wgs_length_sd` ≈ `mean` with a 551 kb tail. Badread's length model has not been verified to
  handle that; this is an open item, not a settled parameter.
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
checks that every designed antigen class produced reads. Current state: **40/40 + 5/5 on both
datasets**, 4/4 alignments at 99.99% mapped, 78.2% on-target against an independent 78.0% prediction,
RNA truth-vs-realised correlation 0.994, and every designed record of all six classes represented
(fusion 39/39, erv 102/102, splice 132/132, cta 69/69, virus 9/9). The three headline LOH events
verified at 0.897/0.980 (TP53 17p), 0.859/0.925 (BRCA1 17q), 0.700/0.725 (RB1 13q).

Eleven defects were found in the acceptance suite itself, four of them visible only at release scale.
They are listed in `catalog-design-notes.md` §13. The pattern worth carrying forward: **four of the
eleven were checks that passed on data that was wrong.** An interval-count floor of 20 skipped an
11-interval amplicon; a depth baseline required the same chromosome, which chr13p and chr17 cannot
satisfy; a baseline neutrality probe tested only a window midpoint and crossed into LOH; and the
read-map integrity check verified resolution rather than uniqueness. When adding a check, state what it
would fail on.

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

**The WGS and long-read assays have no acceptance arm yet.** The depth and allele-fraction checks
measure over capture intervals and there are none. This is the largest open gap.

## 10. Open items

- Acceptance arm for WGS and long reads (§9).
- Kinnex bulk (MAS 8-mer) and Kinnex sc (**16-mer**, not 8-mer) builders.
- 10x 5' GEX and TCR builders; ONT scRNA, bulk RNA and WGS builders.
- ONT WGS identity calibration and length-model verification (§7).
- Release-scale runs of Illumina WGS and PacBio HiFi; both are validated only on chr21.
- Catalog freeze: tag the catalog and both germline VCFs with checksums.
- Resolve HG002 haplotype parentage against HG003/HG004.
- Resolvable source URLs and a netMHCpan version pin in `data-sources.md`.
