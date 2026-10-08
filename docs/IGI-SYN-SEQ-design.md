# IGI-SYN-SEQ synthetic multi-assay tumor datasets: design specification

Status: draft v0.2 (2026-09-24); defaults marked **[default]** were accepted by the owner on 2026-09-24 (see Decisions log). Owner: Steven Vensko. Generator repo: `~/dev/igi-syn-seq`.
Decisions marked **[default]** were delegated and can be overridden before generation starts.

## 0. Decisions log

| Date | Decision |
|---|---|
| 2026-09-24 | All **[default]** choices below accepted (tumor architectures, catalog counts, 4,000-cell composition, depths, ~700 GB per full release pair). |
| 2026-09-24 | Tumor expression baseline: TCGA-BRCA basal-like samples (PanCanAtlas `Subtype_mRNA = Basal`, 193 samples, 191 present in the UCSC Xena Toil kallisto transcript-TPM table) rather than a new salmon run. |
| 2026-09-24 | ERV tumor-specific vs tumor-associated labels use the pan-normal reference shipped with lens-v2.0.0-dev-alt (`pan_normal_and_mtec_exp.95th_perc.homo_sapiens.quant.sf` and `erv_pep_exp_norm_and_mtec.homo_sapiens.tsv`) for now. |
| 2026-09-24 | IGI-SYN-SEQ-02 phasing: WhatsHap read-backed blocks scaffolded by SHAPEIT5 statistical phasing (1000 Genomes GRCh38 panel). |
| 2026-09-24 | Germline SVs are part of both baselines: Q100 `stvar` calls >= 50 bp for HG002 (~46,500, TRF-annotated by GIAB); for IPISRC044 a short-read caller on the normal WGS (Manta proposed, pending), since the IPISRC044 "ONT" data turned out to be single-cell cDNA, not genomic; both annotated for gene overlap and predicted consequence. |

## 1. Goals

- A fully synthetic, fully known-truth human tumor/normal dataset that exercises every LENS
  antigen source (SNV, InDel, fusion, SV, CTA/self, viral, ERV, splice) across every assay LENS
  consumes, including long-read and single-cell.
- Variants span a **spectrum of evidence** so the dataset serves both as positive control
  ("flagpost" events) and for setting detection thresholds (VAF, depth, expression,
  allelic expression, binding, sequence context, clonality).
- Every event is reflected consistently in every assay in which it is physically observable.
- Tumor architecture modeled on **triple-negative breast cancer (TNBC)**.
- A full-genome release plus a **chr1to6 subset** that mirrors the full release but runs fast.

## 2. Dataset naming and layout

Two independent datasets share one catalog design (the IGI-SYN-SEQ family) but differ in
germline baseline and tumor architecture. Each is its own RAFT dataset with its own manifest:

| Dataset | Germline baseline | Sex | Architecture |
|---|---|---|---|
| `IGI-SYN-SEQ-01` | HG002 (GIAB Q100 phased assembly calls) | male | diploid, purity 0.70, no WGD, HRD (SBS3) |
| `IGI-SYN-SEQ-02` | IPISRC044 (public; called + phased here) | male | WGD, purity 0.45, APOBEC-dominated |

Both baselines are male. TNBC in a male host is biologically unusual but has no simulation
cost; the benefit is that HG002 is the PacBio reference sample (public Revio HiFi WGS, Kinnex
bulk, Kinnex sc + matched 10x 5' Illumina) which is what makes accurate long-read error models
possible. If a female baseline is later wanted, HG004 (HG002's mother, GIAB trio-phased) is a
drop-in swap of the baseline VCF only. **[default: keep HG002 and IPISRC044]**

Each dataset has two releases: `full` and `chr1to6`. Manifest values, using `IGI-SYN-SEQ-01`
as the example: `Dataset` = `IGI-SYN-SEQ-01` or `IGI-SYN-SEQ-01-chr1to6`; `Patient_Name` =
`IGI-SYN-SEQ-01` (one synthetic patient per dataset).

RAFT workspace placement (relative to the RAFT workspace root, `$RAFT`), per dataset:

```
inputs/fastqs/IGI-SYN-SEQ-01/<full|chr1to6>/   Illumina FASTQs (WES T/N, bulk RNA, 10x GEX, 10x TCR)
inputs/bams/IGI-SYN-SEQ-01/<full|chr1to6>/     PacBio unaligned BAMs (+ .pbi) and derived FASTQs
inputs/metadata/IGI-SYN-SEQ-01/                manifests, truth bundle, design tables, README
```

## 3. Samples per dataset

| Sample | Assay | Platform | Depth / size [default] | Delivered as |
|---|---|---|---|---|
| tumor WES | WES, hg38_exome.bed | Illumina PE150 | 150x on-target | FASTQ R1/R2 |
| normal WES | WES | Illumina PE150 | 100x | FASTQ |
| tumor bulk RNA | RNA-Seq, polyA, stranded | Illumina PE150 | 80 M pairs | FASTQ |
| tumor bulk WGS | WGS, short read | Illumina PE150 | **100x** | FASTQ R1/R2 |
| normal bulk WGS | WGS, short read | Illumina PE150 | 30x | FASTQ |
| tumor WGS | HiFi WGS | PacBio Revio | 30x | `*.hifi_reads.bam` + `.pbi` + FASTQ |
| normal WGS | HiFi WGS | PacBio Revio | 30x | same |
| tumor Kinnex bulk RNA | Kinnex full-length (MAS 8-mer array) | PacBio Revio | 15 M segmented reads | pre-skera `hifi_reads.bam` **and** post-skera `segmented.bam` |
| tumor Kinnex scRNA | Kinnex single-cell (MAS 16-mer), 10x 5' v2 cDNA | PacBio Revio | 4,000 cells, ~12 M segmented reads | pre-skera + post-skera BAMs |
| tumor 10x 5' GEX | scRNA-Seq 5' v2 | Illumina | 4,000 cells, ~40k reads/cell | FASTQ (Cell Ranger naming) |
| tumor 10x 5' TCR | scTCR-Seq V(D)J | Illumina | ~5k reads/cell | FASTQ |
| tumor ONT scRNA | single-cell cDNA, 10x 5' structure | Oxford Nanopore R10.4.1 | 4,000 cells, ~10 M reads | FASTQ |
| tumor ONT bulk RNA | cDNA, full length | Oxford Nanopore R10.4.1 | ~20 M reads | FASTQ |
| tumor ONT WGS | WGS, long read | Oxford Nanopore R10.4.1 | 30x | FASTQ |
| normal ONT WGS | WGS, long read | Oxford Nanopore R10.4.1 | 30x | FASTQ |

Fifteen sample types per dataset. Earlier notes in this file and in the channel log said ten, eleven,
thirteen and fourteen; those predate one or more of the 2026-10-01 decisions to add Illumina bulk WGS
(samples 4-5) and to mirror all four PacBio sample types in ONT (samples 12-15). Build status per sample
type is tracked in `build-reference.md` §1, which is the one place to look.

**Confirmed 2026-10-02, after a wrong turn worth recording.** The single-cell array really does carry 16
cDNAs: a real array read contains a median of 16 10x TSOs (mean 15.1, p90 16) over 2,000 reads averaging
17,092 bp, against a single-cDNA median of 952 bp. On the way to that I briefly concluded it was an 8-mer,
from weaker evidence -- nine adapters each appearing once per read, and 83 % of skera's segments carrying a
single TSO -- and edited this table to say so. Both observations were real; the inference was not. 12 % of
skera's output is whole un-segmented arrays averaging 8.6 TSOs, which is where the other eight cDNAs were
hiding. The adapter layout remains unresolved and is the one open item for this assay; see
`build-reference.md` §6.3.

The Kinnex sc, 10x GEX and 10x TCR libraries derive from the same 4,000 cells and share
barcodes and UMIs (independent molecule sampling from a common per-cell molecule pool).
ONT scRNA samples that same pool, so a barcode means the same cell in all four single-cell libraries.
No normal RNA. LENS v2.0.0-dev-alt enters Kinnex sc at the segmented-BAM stage
(`lens.nf` line ~58), so the segmented BAM is the LENS input and the pre-skera BAM is for
validating skera itself.

**Tumour WGS depth, revised 2026-10-02.** It was 30x and is now 100x. 30x is a germline-grade depth and
could not see this dataset's own subclones: at 30x the truncal clone yields 10.5 expected alt reads in
IGI-SYN-SEQ-01 and 6.8 in IGI-SYN-SEQ-02, while clone A at CCF 0.40 yields 4.2 and clone A1 at CCF 0.10
yields 1.3 and 0.7. A dataset built to test subclonal neoantigen calling cannot ship a WGS arm that only
exercises truncal calling. At 100x the subclones give 4-14 alt reads and the truncal clone 22-35. Normal
stays at 30x. Tumour WES at 150x remains marginal for the lowest-CCF subclones -- IGI-SYN-SEQ-02 clone A1
sits at 3.4 expected alt reads -- and raising it to 250-300x is an open option, not yet taken.

Storage estimate (full, both datasets): **at least 680 GB** for the sample types costed below, dominated
by the twelve 30x WGS libraries; the long-read and single-cell RNA assays are on top of that and not yet
estimated. Derived from measured bytes per base rather than guessed: 0.761 for a HiFi BAM and 0.441 for
gzipped Illumina FASTQ, both taken from the chr21 validation runs.

| Component | Size |
|---|---|
| 4 PacBio HiFi WGS BAMs | 283 GB |
| 4 ONT WGS FASTQ | ~186 GB (at ~0.50 bytes/base, not yet measured) |
| 4 Illumina WGS FASTQ | 164 GB |
| 4 exomes, merged | 29 GB (measured) |
| 2 bulk RNA | 21 GB (measured) |
| Kinnex, 10x, ONT RNA | not yet estimated |

The previous figure of ~700 GB predates the six Illumina and ONT WGS libraries and counted only the HiFi
BAMs; it was coincidentally close for the wrong reason.

## 4. Germline baselines

- **IGI-SYN-SEQ-01 / HG002**: hg38-projected phased diploid VCF from the HPRC/GIAB HG002 assembly
  (SNV + indel + SV), plus HG002 HLA alleles from clinical typing, trio-phased and concordant with
  the diploid assembly (Chin et al. 2020, Nat Commun 11:4794, Supplementary Table 4):
  maternal A*01:01:01G, B*35:08:01G, C*04:01:01G, DRB1*10:01:01G, DQA1*01:01:01G, DQB1*05:01:01G;
  paternal A*26:01:01G, B*38:01:01G, C*12:03:01G, DRB1*04:02:01, DQA1*03:01:01G, DQB1*03:02:01G.
  Manifest `Alleles` (class I, two-field): HLA-A*01:01,HLA-A*26:01,HLA-B*35:08,HLA-B*38:01,HLA-C*04:01,HLA-C*12:03.
  The HLA LOH event in section 5.1 removes the paternal haplotype (A*26:01/B*38:01/C*12:03).
- **IGI-SYN-SEQ-02 / IPISRC044**: LENS's germline VCFs are exome-restricted, so a genome-wide set is
  built here: DeepVariant (WGS model) on the blood-normal Illumina WGS, then SHAPEIT5 statistical phasing
  (1000G GRCh38 panel) scaffolded by WhatsHap read-backed blocks. IPISRC044 has **no long-read DNA**: the
  "ONT" files are ONT sequencing of the T1 10x 5' single-cell cDNA library (polyA, 10x adapters, median
  read 291 bp), so read-backed phasing comes from short reads only and SHAPEIT5 carries the genome-wide
  phase. Germline SVs must come from a short-read caller on the normal WGS (Manta proposed), annotated
  for gene overlap. Procedure: `baseline-references.md`. HLA alleles A*01:01 homozygous, B*08:01 /
  B*27:05, C*01:02 / C*07:01: stated by the manifest without provenance, and independently confirmed by
  typing the patient's own recovered MHC reads (OptiType 1.3.5, all six calls identical). The MHC was
  initially empty because alt-aware alignment had diverted the reads onto the reference's 525 HLA
  contigs; realigning them to a primary-only reference restored 7,532 records in chr6:28-34 Mb against
  1,822 before, of which 46% are heterozygous.
- A pathogenic **BRCA1** germline frameshift is added to both baselines (TNBC/HRD
  realism; wild-type allele lost somatically, section 5).
- Germline variants are also deliberately placed in a subset of CTA and ERV ORFs
  (section 6) to exercise germline-aware peptide generation.

### 4.1 Designed germline alleles

Each dataset carries a pathogenic germline allele on top of its baseline genotype. It is not a
neoantigen source, but it has to be in the sequence: it is present in both the tumour and the normal
library, it is why the somatic loss of heterozygosity at that locus matters, and a caller that reports
it as somatic is making a mistake the benchmark should be able to see.

Both datasets carry BRCA1 c.68_69delAG (the 185delAG founder allele), p.Glu23ValfsTer17, phased onto the
haplotype that survives the truncal 17q loss so the tumour retains only the pathogenic copy. The allele
is specified by coding position and resolved against the annotation at build time, then checked against
the expected codon, wild-type residue and consequence before it is written; a hardcoded genomic
coordinate would be silently wrong on a different annotation release. The resulting VCF is the dataset's
own `germline.vcf.gz`, and the baseline it was built from is recorded alongside it.

## 5. Tumor architecture (TNBC)

Shared driver logic, dataset-specific parameters.

### 5.1 IGI-SYN-SEQ-01 / HG002 (purity 0.70, ploidy ~2.3, no WGD)

Clone tree: `T` (truncal, CCF 1.0) -> `A` (CCF 0.40) -> `A1` (CCF 0.12, nested in A);
`B` (CCF 0.25, sibling of A). Sum of sibling CCFs 0.65 <= 1.

| Clone | Drivers / arm events |
|---|---|
| T | TP53 R248Q + 17p LOH; RB1 frameshift + 13q LOH; PTEN homozygous focal deletion (SV, 10q23); BRCA1 wild-type allele lost (17q LOH); MYC 8q24 focal amp CN 12; 1q gain CN 3; 5q loss; 8p loss; 16q loss; HLA LOH on 6p (loses one haplotype: A*26:01/B*38:01/C*12:03); chromothripsis on 5p (~60 breakpoints); HRD tandem-duplication phenotype (~100 TDs, 1-10 kb, genome-wide) |
| A | PIK3CA H1047R; EGFR focal amp CN 8 |
| A1 | NF1 frameshift |
| B | KMT2C nonsense; one private in-frame fusion |

Signatures: SBS3 0.60, SBS1 0.15, SBS5 0.15, SBS13 0.10; indels ID6-dominated
(deletions with microhomology). Background genome-wide burden ~1 mut/Mb
(~3,000 SNV, ~300 indel), almost all non-coding.

Expected VAF tiers at purity 0.70: truncal LOH region ~0.54, truncal amplified (MYC
region) up to ~0.8, truncal het ~0.35, clone A ~0.14, clone B ~0.09, clone A1 ~0.04.

### 5.2 IGI-SYN-SEQ-02 / IPISRC044 (purity 0.45, WGD, ploidy ~3.6)

Same clone topology with four clones (T, A, A1, B) at CCF 1.0 / 0.35 / 0.10 / 0.30.
Drivers: TP53 frameshift + LOH (pre-WGD, so 0 wild-type copies out of 4 total at the
locus after WGD), RB1 whole-gene deletion, PTEN nonsense + LOH, BRCA1 germline + LOH,
MYC amp CN 20, EGFR amp in A, PIK3CA E545K in T, chromothripsis on chr3p, HLA LOH losing the B*27:05 / C*01:02 haplotype (hap1), matching
`design.yaml`; the retained haplotype carries B*08:01 / C*07:01. Signatures: SBS2+13 0.50, SBS3 0.20, SBS1+5 0.30; ~4 mut/Mb
background (~12,000 SNV, ~1,200 indel). Post-WGD private mutations have 1 of ~4 copies,
which pushes many events into the low-VAF regime (truncal het post-WGD ~0.12).

Tumor-in-normal contamination: 0 % **[default off; parameter available]**.

## 6. Designed variant catalog

### 6.1 Evidence-spectrum axes

Every designed event carries a tier label on each axis in the truth bundle.

| Axis | Tiers |
|---|---|
| clonality | T-LOH, T-amp, T-het, A, B, A1 |
| expression (gene TPM) | 0, ~1, ~10, ~100, ~1000 |
| allelic expression of mutant allele | balanced 0.5, silenced 0.1, dominant 0.9 (subset only) |
| MHC-I binding (patient HLA, netMHCpan %rank) | strong <=0.5, weak 0.5-2, non >2 |
| sequence context | clean-unique; homopolymer >=6; low-mappability/segdup; within 30 bp of germline het; somatic pair within one codon or <=150 bp; exon edge (<=10 bp from capture boundary); deep intronic/UTR (WGS-only) |
| capture | WES on-target vs off-target |

**Flagpost** = truncal, LOH/amplified (VAF >= 0.5), TPM >= 100, balanced, strong binder,
clean context, exon center. ~40 flagposts across all classes.

### 6.2 Counts per dataset

| Class | Designed count | Composition |
|---|---|---|
| SNV, missense core grid | 300 | 6 clonality x 5 expression x 3 binding cells, ~3-4 replicates each, clean context |
| SNV, missense context strata | 150 | 6 context strata x 25 |
| SNV, nonsense | 60 | incl. NMD-escape (last exon) vs NMD-sensitive |
| SNV, synonymous | 40 | negatives; a few create cryptic splice sites (count under splice) |
| SNV, start-loss / stop-loss | 20 | |
| SNV, hotspot flagposts | 12 | TP53 R248Q/R273H, PIK3CA H1047R/E545K, KRAS G12D, BRAF V600E, IDH1 R132H, etc. |
| InDel, frameshift del / ins | 80 / 60 | lengths 1-10, long tail 11-50; neo-ORF lengths span 1-60 aa |
| InDel, in-frame del / ins | 40 / 30 | |
| InDel, homopolymer context | 40 | PacBio-hard indels |
| Gene fusions | 20 | 4 flagposts; in-frame, out-of-frame, and read-through with no DNA breakpoint (tumor-associated); mechanism follows from the partners' positions and orientations. Exome visibility is not a quota: a breakpoint is visible when it falls in the padded, merged capture intervals used for generation, and each run's junctions table is authoritative. Achieved: -01 17 fusions with 8 visible, -02 20 with 9 |
| SVs (non-fusion designed) | 100 (achieved 98 per dataset; exome-visible 22 in -01 and 26 in -02, by the same interval rule as fusions) | DEL 25 (50 bp-5 Mb log-spaced), DUP 20, INV 15, TRA 15, INS 15 (10 L1/Alu/SVA MEIs, 5 novel sequence), complex 10; 50 land in coding sequence: exon-deleting in-frame 15, out-of-frame 15, whole-gene loss 10, intragenic exon dup 10 |
| SVs, structured background | ~160 | chromothripsis cluster (~60) + HRD tandem dups (~100) |
| Viruses | 3 | HPV16 integrated at 8q24 inside the MYC amplicon, E6/E7 expressed, host-virus fusion transcript, junction visible in WGS/RNA but not WES, the site being outside the bait set (D8); EBV episomal ~5 copies/cell, low expression, plus trace 0.05 copies/cell in the blood normal (realistic negative); HPV18 episomal in clone B only, unexpressed |
| CTAs | 15 | 5 expression tiers x 3; 4 carry germline coding variants, 3 carry somatic missense, 2 restricted to clone A (single-cell heterogeneity) |
| ERVs | 30 | split by the measured pan-normal reference rather than a fixed quota: tumor-specific at or below 0.05 TPM in normals, tumor-associated 0.05-5 TPM, the rest unexpressed negatives. Tier sizes are whatever the reference yields for the loci sampled, so these are targets not guarantees and each run records its achieved counts. Observed: 10/4/16 in -01, 10/5/15 in -02. A locus with no measurement is skipped rather than assigned |
| Splice variants | 30 | tumor-specific 15 = somatic splice-site/cryptic variants causing exon skip 5, intron retention 3, cryptic 5' 3, cryptic 3' 2, novel exon 2; tumor-associated 15 = annotated isoform switch 10 + novel junction with no DNA cause 5 |
| Negative controls | ~110 | 30 low-depth germline hets that mimic somatic; 30 A>I RNA-editing sites (RNA-only); 20 processed-pseudogene parent mismatches; 10 germline variants in CTA/ERV; 20 somatic SNVs in unexpressed genes |

Totals per dataset: ~620 designed coding SNVs, ~250 designed indels, 20 fusions,
~260 SVs, 3 viruses, 15 CTAs, 30 ERVs, 30 splice events, plus signature-driven background.
Coding burden is ~10x a real TNBC by design.

### 6.3 Placement rules

- >= 50 % of every class and every tier cell placed on chr1to6 so the subset preserves the
  spectrum. HLA (6p) and IGK (2p) are inside the subset; TRA/TRB/TRG (chr14/7) are not.
- Fusion partners chosen from LENS's fusion reference where possible plus novel pairs.
- CTA and ERV loci drawn from `references/homo_sapiens/cta_self` and the HERV annotation in
  `gencode.v37.annotation.with.hervs.gtf`.
- Peptide binding evaluated against each patient's HLA with netMHCpan (and mhcflurry as
  a second opinion) at design time; the truth bundle records both.

### 6.4 HLA loss and the two binding tiers (owner decision D4)

LENS ranks candidates against all of a patient's alleles, but a tumour that has lost an HLA haplotype
cannot present through the alleles on it. Both numbers are real, so the truth bundle carries both.

Per event: `binding_tier` / `best_rank_el` / `best_allele` / `best_peptide` over all six alleles;
`binding_tier_retained` / `best_rank_el_retained` / `best_allele_retained` / `best_peptide_retained`
over the alleles the tumour retains; and `best_allele_is_lost`, true when the best all-allele result
belongs to a lost allele. `<dataset>.hla_loh.tsv` lists every allele in every clone with its haplotype
assignment, copy number, lost flag and the CNA label responsible.

The evidence grid is filled on the **all-allele** tier, because that matches how LENS ranks today. The
retained tier is an additional column. Scoring a caller follows the same split: ranking against the
all-allele tier, and the warning a caller should emit about unusable alleles against the lost-allele flag.

Both datasets lose haplotype 1 over `chr6:28,500,000-33,500,000`:

| Dataset | Lost | Retained |
|---|---|---|
| IGI-SYN-SEQ-01 | A\*26:01, B\*38:01, C\*12:03 | A\*01:01, B\*35:08, C\*04:01 |
| IGI-SYN-SEQ-02 | B\*27:05, C\*01:02 | A\*01:01, B\*08:01, C\*07:01 |

Dataset 01 loses a full haplotype and therefore one allele of each class I gene. Dataset 02 keeps
A\*01:01 because it is homozygous, so only B and C are affected.

Which allele sits on which haplotype is not something a SNV-level phased VCF settles, so the arrangement
is declared in `design.yaml` with its basis recorded. For IGI-SYN-SEQ-02 it follows linkage
disequilibrium: A\*01:01-B\*08:01-C\*07:01 is the 8.1 ancestral haplotype and B\*27:05 travels with
C\*01:02. OptiType on the recovered normal returned exactly the six manifest alleles, which confirms the
genotype but not the phase. For IGI-SYN-SEQ-01 the B-C linkages are likewise established but which
haplotype is maternal is still unverified.

Because the lost alleles account for a large share of best hits -- about 50 % of events in dataset 01 and
36 % in dataset 02 -- the two tiers differ for a substantial fraction of the catalog, which is the point
of carrying both.

## 7. Cross-assay reflection matrix

| Event class | WES T | WES N | bulk RNA | HiFi WGS T | HiFi WGS N | Kinnex bulk | Kinnex sc | 10x GEX | 10x TCR |
|---|---|---|---|---|---|---|---|---|---|
| somatic SNV/indel (exonic) | yes | no | yes if expressed, ASE-weighted | yes, phased | no | yes if expressed | per cell, clone-restricted | per cell (5' bias) | no |
| somatic SNV/indel (deep intronic) | off-target only | no | intron retention only | yes | no | no | no | no | no |
| germline variants | yes | yes | yes | yes | yes | yes | yes | yes | no |
| CNA / LOH / HLA LOH | depth + BAF | baseline | ASE, dosage | depth + BAF, phased | baseline | dosage | dosage per cell (scevan) | dosage per cell | no |
| fusion (with DNA breakpoint) | if breakpoint captured | no | chimeric reads | split/spanning reads | no | full-length chimeric | per cell | 5'-end chimeric | no |
| read-through fusion | no | no | yes | no | no | yes | yes | yes | no |
| SV | if breakpoint captured | no | if transcribed | yes | no | if transcribed | if transcribed | rarely | no |
| virus integrated | no (see below) | no | host-virus + viral transcripts | junction + viral reads | no | yes | yes | yes | no |
| virus episomal | off-target | trace EBV | yes | yes | trace EBV | yes | yes | yes | no |
| CTA / ERV expression | germline vars only | germline vars only | yes | germline vars only | germline vars only | yes | per cell | per cell | no |
| splice (somatic-caused) | causal variant | no | junction | causal variant | no | full-length isoform | per cell | 5' only | no |
| splice (associated) | no | no | junction | no | no | full-length isoform | per cell | 5' only | no |
| TCR clonotypes | no | no | TRA/TRB reads at bulk level | germline loci | germline loci | some TCR transcripts | some TCR transcripts | 5' TCR reads | yes |

**Cancer-testis antigens are expressed regardless of their cohort baseline (owner decision D6).** A CTA
is silent in somatic tissue and de-repressed in tumour, and discovering that de-repression is the point of
including the class. Its tumour expression therefore comes from the designed `target_expression_tier`
(T1 1.5, T10 15, T100 120, T1000 600 TPM), never from the expression baseline, which reads essentially
zero for a testis-restricted gene in a breast cohort. A zero baseline suppresses only the gene's ordinary
reference transcript; it never suppresses the CTA record. Twelve of the fifteen CTAs per dataset are
expressed at their designed tier.

The remaining three per dataset sit in the T0 tier and are deliberately **not** expressed. They are
negative controls: a discovery tool that calls them is wrong, and they are there to catch that. They are
not a gap in coverage.

CTAs are expressed from both haplotypes, because de-repression is epigenetic and acts on both alleles.
Emitting from one would make every germline heterozygous site inside a CTA read as homozygous, and these
are variant-dense genes.

**Integrated virus and the exome (owner decision D8).** The HPV16 integration at chr8:127,400,000 sits
outside the bait set and its flanking bands, so it produces no exome reads at all and the row above reads
no rather than off-target. The site was chosen to sit inside the MYC amplicon, which is where such an
integration biologically belongs, and the owner declined to move it to make it capturable. Viral
integration is therefore detectable from WGS, RNA and the single-cell assays, and not from WES. Where a
designed convenience and realism conflict, this dataset prefers realism and says so.

## 8. Single-cell design (4,000 cells)

Composition **[default]**: tumor 35 % (split by clone: T-only 35 %, A 40 %, A1 12 % of A,
B 25 % of tumor cells, so per-cell genotypes reproduce the CCFs), CD8 T 18 % (exhausted,
effector, memory), CD4 T 9 %, Treg 3 %, B 5 %, plasma 3 %, macrophage/monocyte 12 %,
DC 2 %, CAF 8 %, endothelial 4 %, NK 1 %. Doublets 6 %, ambient RNA 3 %.
10x 5' v2 chemistry (16 bp barcode from `737K-august-2016.txt`, 10 bp UMI, TSO),
matching `references/single_cell/10x_5prime_v2_pacbio_kinnex`.

TCR: ~1,200 T cells with productive paired alpha/beta (10 % alpha-dual), ~300 clonotypes,
power-law expansion (top 10 clonotypes hold ~30 % of T cells); 5 flagpost clonotypes
mapped in the truth bundle to specific flagpost neoantigens.

Junctions are recombined from the real germline rather than invented. A CDR3 runs from the conserved
cysteine at IMGT 104, contributed by the V gene, to the phenylalanine of the J gene's `FGXG` motif at
IMGT 118; `catalog/resources/vdj_germline.json` holds exactly those nucleotides, plus the framework either
side, for the 58 functional V and 63 functional J genes GENCODE v37 carries over GRCh38. `vdj.recombine`
trims each end by whole codons, inserts a GC-rich N region free of stops and cysteines, and keeps the
result in frame. Measured over 3,000 draws: beta CDR3 length mean 13.9 (real ~14), alpha 13.5 (real
~13.5), 100 % anchored C..F, 0 % internal cysteine, 0 stops.

This replaced random peptides -- the literal string `CAS` or `CAV`, uniform draws over all twenty amino
acids, then `F`. That produced junctions like `CASCAWCQESIIYPATQVF`, which no selected repertoire
contains, and worse, left the junction unrelated to the V and J genes the truth table named, so a caller
that correctly recovered the V gene from the framework would have been scored wrong. The N-region
nucleotides are still drawn rather than taken from a real repertoire, so junction length and composition
are realistic but the clonotypes are not any individual's.

Tumor cells' expression reflects CNA dosage (for scevan) and clone-restricted CTAs.

### 7.1 Illumina short-read WGS (owner decision, 2026-10-01)

Samples 4 and 5 of fifteen. The design previously had WGS only as PacBio HiFi, so the only short-read DNA was
the exome pair. Most somatic callers and copy-number tools are tuned on short-read WGS and an exome cannot
substitute for one off-target, so both tumour and normal are added at 30x.

It is the exome path without the capture step. `WesBuilder` already rebuilds each interval on every clone
and haplotype that retains it, applies germline and somatic edits in reference coordinates, weights depth
by absolute copy number, substitutes a junction contig where a rearrangement sits, and names reads uniquely
across chromosomes; none of that is capture-specific. Three differences:

- Intervals are 5 Mb windows tiling the chromosome rather than bait regions, so there is no on-target or
  off-target distinction and no flanking bands. Windows keep a chromosome from being held in memory twice.
- **The GC model is flat.** The measured curve is *capture efficiency*, fitted from how much sequence the
  baits pulled down at each GC, and applying it without baits would impose a bias with no cause. Real WGS
  carries a milder PCR-driven GC bias which this does not model: a known limitation, recorded rather than
  approximated with the wrong curve.
- Depth is uniform, so copy number is the only thing modulating it.

Read names occupy a different slice of the flowcell coordinate space from the exome's, so an exome and a
WGS library from the same dataset can be held together without colliding.

Acceptance needs its own arm for this assay: the existing depth and allele-fraction checks measure over
capture intervals, and for WGS they have to work on genome windows instead.

### 8.1 ONT single-cell RNA (owner decision, 2026-10-01)

Sample 12 of fifteen. Kinnex single cell already covers long-read single cell, so what ONT adds is
**platform diversity**: its indel-heavy homopolymer errors break variant calling and isoform assignment
differently from HiFi's, and a pipeline tuned on one platform is exactly what a benchmark should expose.

The error model is fitted from IPISRC044's own ONT single-cell RNA, which is the best provenance in the
dataset: the same individual as IGI-SYN-SEQ-02 and the same assay being simulated, sequenced 2025-07-17 on
a PromethION with `dna_r10.4.1_e8.2_400bps_sup@v4.3.0`. Measured on 100,000 reads: mean length 910 (sd
461, median 843), accuracy 0.9822 (Q17.5).

pbsim3 cannot reach modern ONT at all. Both its ONT models produce about 84 % accuracy, a nine-fold higher
error rate than the real reads, and unlike HiFi there is no consensus step to recover it; the model named
"HQ" is marginally worse than the plain one. Badread's `nanopore2023` model reproduces the real reads
closely: length 917 (460), accuracy 0.9800 (Q17.0).

The real reads are used as the target and as parameters only, never as material to edit. Editing them
would give a perfect error profile but would carry that patient's own somatic variants into the
background, make the truth set incomplete, and redistribute real patient reads — a larger version of what
D12 required removing. Nothing of the individual's sequence reaches the release.

This assay shares the per-cell roster with Kinnex single cell, 10x 5' GEX and 10x TCR, so a barcode has
one cell type, one clone and one clonotype across all four libraries.

## 9. Expression baseline

Tumor-cell expression baseline: the median transcript-level TPM across TCGA-BRCA basal-like
tumors (PanCanAtlas `Subtype_mRNA = Basal`; 191 of 193 samples are in the UCSC Xena Toil
`tcga_Kallisto_tpm` table, GENCODE v23 transcript IDs, log2(TPM+0.001)). Transcript IDs are
mapped to the LENS annotation (GENCODE v37) by stable ENST ID, dropping version suffixes;
transcripts absent from v23 inherit their gene's mean. A single sample is not used because
the goal is a plausible basal-like profile, not a specific patient. Immune/stromal cell-type
profiles come from the IPISRC044 UCSF SCG1 data (real 10x 5' TME). Designed events override
the baseline for their genes (expression tier).

## 10. Read simulation and error models

| Assay | Simulator | Error/quality model source |
|---|---|---|
| Illumina WES/WGS-style | ART (art_illumina) with custom profiles + fragment model; GC bias and capture efficiency from real IPISRC044 WES | `art_profiler_illumina` on IPISRC044 WES reads |
| Illumina bulk RNA | molecule sampler (per-haplotype, per-clone transcriptome) -> fragmentation -> ART | IPISRC044 RNA reads; positional 3' bias fitted from real data |
| 10x 5' GEX / TCR | custom molecule sampler writing R1 (BC+UMI+TSO) and R2 with ART qualities | IPISRC044 SCG1/TCR reads |
| PacBio HiFi WGS | pbsim3 (quality-score model, multi-pass -> ccs) | HG002 Revio HiFi public data |
| Kinnex bulk / sc | Badread over MAS arrays built from full-length molecules; 8-mer array for both, 10x 5' v2 cDNA structure and polyA for single cell | HG002 Kinnex bulk `segmented.bam` (adapters and array structure decoded from skera's tags) and HG002 Kinnex sc `0-CCS` arrays (10x structure and segment lengths) |

PacBio outputs are written as unaligned BAM with Revio-style read names, `RG`/`PU`,
and HiFi tags (`np`, `rq`, `ec`), indexed with `pbindex`; skera is run to produce the
segmented BAM so the delivered pair is exactly what a Revio run yields.

Three constraints of that format are not optional and were each found by a tool refusing the output
rather than by reading a specification:

- The `@RG` `ID` must be **eight hexadecimal digits** (real ones look like `b0776b05`). pbbam parses it as
  a number, so `ID:synthetic` aborts `pbindex` with `ERROR: stoul` and leaves a 65-byte index reporting
  zero reads.
- The `@RG` must carry **`PU`**, the movie name. skera names each segment `<movie>/<zmw>/ccs/<qs>_<qe>`
  and takes the movie from `PU`, not from the read name it was handed; without it every segment name
  begins with a slash.
- Every field of the movie name after the instrument must be **decimal**, because pbbam parses the date
  and time fields as numbers.

**The HiFi error rate is about 5x real, and the quality scores do not say so.** Simulated from the plain
reference, where the alignment is ground truth, these reads measure 0.99089 accuracy (9.11e-03 error/bp)
while their quality strings claim 0.99833 (1.67e-03/bp). Real HiFi is 1.77e-03/bp. Badread's qscore model
floors the relationship between requested identity and reported quality, so the two cannot both be matched;
the identity request was tuned against the quality-derived figure, which is the reads' claim about
themselves rather than their content. The homopolymer indel enrichment is 3.24x against a real 4.34x, so
errors do concentrate in homopolymers but less than in real HiFi. Both are open; see CHANNEL.md [23].

**`np` is fabricated, and nothing should trust it as a pass count.** Badread applies an error model per
read with no pass structure, so there are no passes behind these reads at all. The `np` written here is
drawn from the real HG002 distribution so that a tool filtering on pass count sees a realistic
distribution rather than a single constant, but it is a plausible number, not a measurement: it does not
describe how the read was produced, and the read's accuracy was not obtained by consensus over that many
passes. `ec` inherits the same caveat, being derived from `np`. Any analysis that treats `np` as evidence
about consensus depth is reading something that was invented.

`np` is drawn per read from the empirical distribution measured off the real HG002 BAMs
(`catalog/resources/pacbio_np_model.json`: median 7, mean 7.60, p10 4, p90 13) and `ec` is `np` times the
measured `ec/np` ratio of 1.130. A single constant `np` leaves a caller filtering on pass count with
nothing to filter on. `np` is not derived from each read's own `rq`, even though that would be more
self-consistent, because Badread's per-read accuracy spread is narrower than Revio's -- sd 0.00089 against
a p10-p90 span 2.6x wider in the real data -- so deriving it would give a far too tight `np` distribution.
The read's `rq` remains computed from its own simulated qualities.

**Long-read FASTQ headers carry no truth.** Badread writes the source reference, strand and coordinates
into every genomic read's description and the source molecule id into every RNA read's, plus
`error-free_length` and `read_identity`, which do not exist in real data. The description is stripped and
written to `*_read_map.tsv.gz` instead, which is the same thing `readnames.shuffle_and_rename` does for
the Illumina arms.

All simulators run from containers (none are installed on the cluster; only wgsim,
samtools and bcftools are on PATH). Workload runs under SLURM.

### 10.1 Off-target capture coverage

A hybrid-capture library is not confined to its bait set. Probes pull down the sequence flanking each
bait, so coverage decays outward over a kilobase or two rather than stopping at the bait edge, and a
thin background covers the rest of the genome. Building only on-target reads would leave every base
outside the bait set at exactly zero depth, which no real exome shows and which makes off-target
copy-number signal, mapping artefacts near baits and off-target germline calls impossible to exercise.

The exome builder therefore emits three extra interval sets alongside the bait set, each at a depth
relative to the on-target depth:

| Set | Distance from a bait | Relative depth |
|---|---|---|
| proximal | 1-500 bp | 0.15 |
| mid | 501-2,000 bp | 0.02 |
| distal | random windows clear of every bait | 0.012 |

The depths are means over each band, not the peak at the bait edge, since real flank coverage decays
within a couple of hundred bases.

Measured on the full release (68,112,282 pairs per tumour library): 78.2 % of reads fall inside the
capture intervals the builder uses, 15.3 % proximal, 4.9 % mid and 1.6 % distal, identical to within a
tenth of a point across both datasets and matching the 78.0 % that the band sizes and depths predict
independently. An uncapped chr6 library gave 77.6 %, so the figure is stable from one chromosome to the
whole genome. That figure depends
entirely on which interval list it is measured against, and both numbers matter:

- Against the builder's capture intervals, which are the bait list padded by 100 bp and merged across
  200 bp gaps, the library is **77.6 % on target**. The padding is 2.19 Mb of the 6.42 Mb total, so a
  third of what this definition calls on target is really bait flank.
- Against the raw bait list (4.23 Mb over 11,091 baits on chr6), the same library is about **51 % on
  target**, which is where a real exome measured by Picard `HsMetrics` normally sits. This figure is
  calculated from the interval sizes rather than from an alignment, so treat it as an estimate until a
  BAM confirms it.

The distal set samples the genome-wide background rather than tiling it: any window it covers looks like
real off-target sequence, but the untiled remainder contributes nothing, so the library's total
off-target read count is lower than a real one. Local depth realism is the more useful property for a
benchmark, so it is the one preserved.

Every set runs through the same clone, haplotype and copy-number machinery as the bait set, so loss of
heterozygosity and amplification are visible off-target as well. A rearrangement junction is handed out
once across all four sets, because the bands overlap the bait flanks and the same junction would
otherwise be counted twice. The record map records which set each source record came from.

### 10.2 Passenger background

Designed events are the ones a caller is meant to find. On their own they give a tumour genome whose
only somatic differences are the ones under test, which is not what a caller sees and makes any
false-positive rate measured against it meaningless. Each dataset therefore carries a passenger
background at its configured `background_mut_per_mb`, drawn from a mixture of COSMIC v3.4 SBS profiles:
SBS3-dominant for IGI-SYN-SEQ-01 (HRD), APOBEC-dominant for IGI-SYN-SEQ-02.

A 96-channel is drawn first and a matching trinucleotide position second. This is the convention the
COSMIC profiles are written in, so the realised spectrum reproduces the profile without a
context-frequency correction. Passengers are assigned to clones by branch, avoid both germline variants
and the spans of designed events, and those landing in coding sequence are annotated and scored through
netMHCpan. A passenger missense in a well expressed gene is a genuine neoepitope; a truth set that
omitted it would score a caller's correct answer as a false positive.

## 11. Truth bundle (`inputs/metadata/IGI-SYN-SEQ-01/truth/`, likewise for -02)

- `germline.phased.vcf.gz` (SNV/indel/SV) with haplotype IDs.
- `somatic.vcf.gz` with INFO: clone, CCF, expected VAF per assay, tier labels, flagpost flag,
  chr1to6 flag, neoantigen peptides and binding ranks.
- `cna.seg` (per clone, per haplotype copy number), `hla_loh.tsv`.
- `sv.bedpe`, `fusions.tsv`, `viral_integrations.tsv`, `viral_copy_number.tsv`.
- `novel_isoforms.gtf` + `splice_events.tsv` (specific vs associated, mechanism).
- `erv_expression.tsv`, `cta_expression.tsv` (with tier and germline/somatic variant flags).
- `expression.tumor_bulk.tsv` (expected TPM per transcript per clone and mixed).
- `cells.tsv` (barcode, cell type, clone, doublet flag), `clonotypes.tsv` (barcode, TRA/TRB CDR3,
  V/J, flagpost antigen link), `negatives.tsv`.
- `manifests/` : LENS manifests for full and chr1to6, with HLA alleles filled.
- `design.yaml`, `seed.txt`, generator version.

## 12. chr1to6 subset

Derived from the full release, not re-simulated. Built in two halves, because the release is built
in two halves.

**Per-chromosome arms** (WES, Illumina WGS, PacBio HiFi WGS, ONT WGS) are simulated one
chromosome at a time, so their subset is the chr1-chr6 shards **merged into one library per
sample** (`jobs/91_merge_chr1to6.sbatch`). The merge is not cosmetic: lens-v2.0.0-dev's
`preflight_resolve_inputs` gathers every file whose name starts with `File_Prefix` and then fails
hard on more than one BAM, so PacBio's six shards cannot be addressed by any prefix; and with the
release's chromosome-first names (`<ds>_chr1_tumor_wgs_R1`) no prefix selects one assay's six
chromosomes without also selecting the other assays'.

**Whole-library arms** (bulk RNA, 10x GEX, Kinnex bulk, Kinnex sc, ONT bulk RNA, ONT sc RNA) are
filtered read by read by `catalog/subset_chr1to6.py`
(`jobs/90_subset_map.sbatch`, then `jobs/90_subset_chr1to6.sbatch`).

**Selection is by truth map, not by alignment**, which departs from this section's original
wording. Each library ships a read-level truth sidecar naming the transcript record every read
came from, and a read is kept when that record's locus is on chr1-6. Three reasons:

1. *It is the rule the other ten arms use.* A per-chromosome shard holds the reads whose SOURCE
   locus is on that chromosome. Selecting the RNA arms by where an aligner puts them would make
   the subset release use two different definitions of "on chr1-6" depending on the assay.
2. *It is exact* -- no MAPQ, multimapping or soft-clip judgement.
3. *It is cheap* -- aligning these twelve libraries (276 GB of reads) would cost more than
   simulating them did.

What alignment would have added is the read whose source is on chr7 but which maps into chr1-6
anyway. Those reads are excluded, and that is stated rather than discovered.

Per source class: reference and splice transcripts resolve through the GENCODE transcript, ERV and
CTA through their event's locus, a **fusion is kept if either partner** is on chr1-6 (it is one
molecule and half of it is in the subset), and **viral transcripts are kept unconditionally** --
there is no human chromosome to test, and this section keeps viral reads.

- 10x GEX and Kinnex sc keep all barcodes but only chr1to6 molecules. Measured: all 4,000
  barcodes retained in both; molecules drop to ~35 % of full, which is the UMI depth Cell Ranger
  cell calling has to work at.
- **10x TCR is kept complete** because TRA/TRB/TRG lie outside chr1to6. It is linked under
  `full/` only -- a second link under the same name made `<ds>-TCR_S1_L001` match four FASTQs --
  and the chr1to6 manifest points at that same file.
- Truth bundle filtered to chr1to6 events; manifests regenerated as
  `lens.<ds>.chr1to6.{with,no}-alleles.manifest` with `Dataset` = `<ds>-chr1to6`.

**Caveat, Kinnex single-cell.** A single-cell array carries sixteen cDNAs but the kit has only
nine adapters, so `skera` returns the last bracketed segment as molecules 7-15 concatenated (§10
records this as an open adapter-layout caveat). A segment is kept if **any** molecule in it is on
chr1-6, so that one segment almost always survives: measured molecule purity is **0.494**, against
1.000 for Kinnex bulk where one segment is one molecule. Half the molecules in the Kinnex sc
subset are therefore off-target. That is the adapter-layout caveat surfacing, not a new defect,
and each arm's `*_chr1to6_<arm>.json` reports the purity so the size of it stays visible.

**Caveat, event coverage.** §4 asks for >= 50 % of every class and tier cell on chr1to6. Measured
against the truth tables' own `chr1to6` flags, only SNVs (77 %/75 %) and viral integrations
(67 %) clear that; indels are 34 %/32 %, SVs 23 %/34 %, fusions 26 %/40 %, splice 23 %/37 %, CTA
20 %/27 %, ERV 40 %/40 %. The subset is a fast end-to-end exercise of the pipeline, not a
proportional miniature of the full truth set, and raising the non-SNV classes would mean
re-running the designer and invalidating every library built from it.

## 13. Generator pipeline

Nextflow DSL2 in `~/dev/igi-syn-seq`, containerized, seeded, parameterized by `design.yaml`
(so a different catalog or baseline can be regenerated). Stages:

0. `baseline`: fetch/construct phased germline per dataset.
1. `catalog`: seeded catalog design (Python): pick loci per tier, evaluate peptides with
   netMHCpan/mhcflurry, emit truth tables.
2. `genomes`: per-clone, per-haplotype genome FASTAs with CNA multiplicities and SVs;
   viral contigs and integrations.
3. `transcriptomes`: per-clone haplotype-aware transcript sets with designed isoforms,
   fusions, ERVs, CTAs, and expression vectors.
4. `simulate`: per-assay read simulation with clone/purity mixing.
5. `package`: BAM/FASTQ naming, pbindex, skera, manifests, README.
6. `subset`: chr1to6 derivation.
7. `validate`: align back, check VAF/depth/expression against truth, run Cell Ranger and
   skera/lima/isoseq dry-runs, run LENS chr1to6 end to end.

## 14. Open items

- Obtain a Kinnex single-cell read set (UCSF raw T1 BAM preferred; ENA LongBench as fallback).
- Cell Ranger for validating the 10x outputs: run from the community-built nf-core image
  (`quay.io/nf-core/cellranger:9.0.1`; 10.0.0 also available) under singularity, with the 10x GRCh38 2024-A
  GEX and 7.1.0 V(D)J references. 10x's EULA governs use; the image is not redistributed by us.
