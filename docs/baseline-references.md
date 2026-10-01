# Building the germline baselines for IGI-SYN-SEQ-01 and IGI-SYN-SEQ-02

The two IGI-SYN-SEQ datasets each have one synthetic patient, who needs a **haplotype-resolved diploid germline**
(phased SNVs and indels, ideally SVs too) on GRCh38, because every simulated assay, including
long reads and single cells, is generated from personalized haplotype sequences. This page
records how each baseline is produced so it can be reproduced on any cluster. Site-specific
values (paths, partitions, image locations) are deliberately absent.

Placeholders used below: `$REF` (GRCh38 FASTA with `.fai`), `$NORMAL_BAM` (short-read normal
WGS aligned to `$REF`), `$LR_FASTQ/` (long-read FASTQs), `$OUT/` (output directory), `$N` (CPUs).

## Tool versions

| Tool | Version | Image |
|---|---|---|
| DeepVariant | 1.8.0 | `google/deepvariant:1.8.0` |
| WhatsHap | 2.4 | `quay.io/biocontainers/whatshap:2.4--py310h184ae93_0` |
| minimap2 | 2.22 (2.28 recommended) | `quay.io/biocontainers/minimap2:2.28--h577a1d6_4` |
| samtools / bcftools | 1.21 | `quay.io/biocontainers/samtools:1.21--h50ea8bc_0`, `quay.io/biocontainers/bcftools:1.21--h8b25389_0` |

Run every command through `singularity exec <image>` (or apptainer); nothing is installed on hosts.

---

## IGI-SYN-SEQ-01: HG002 (GIAB Ashkenazi son, NA24385, male)

No variant calling is needed. NIST publishes an assembly-derived, haplotype-resolved call set
for HG002 projected onto GRCh38, plus raw Revio HiFi reads for error-model training.

### 1. Download
Base URL: `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/`

| Content | Path under base URL | Size |
|---|---|---|
| Small variants, phased (GT uses `\|`) | `analysis/NIST_HG002_DraftBenchmark_defrabbV0.020-20250117/GRCh38_HG2-T2TQ100-V1.1_smvar.vcf.gz` (+ `.tbi`) | 40 MB |
| Small-variant benchmark regions | `.../GRCh38_HG2-T2TQ100-V1.1_smvar.benchmark.bed` | 0.7 MB |
| Structural variants, phased | `.../GRCh38_HG2-T2TQ100-V1.1_stvar.vcf.gz` (+ `.tbi`) | 46 MB |
| SV benchmark regions | `.../GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed` | 0.2 MB |
| Revio HiFi WGS, raw unaligned movies | `HG002_NA24385_son/PacBio_HiFi-Revio_20231031/HG002_PacBio-Revio_m84039_230928_213653_s3.hifi_reads.bam`, `..._m84039_231005_222902_s1.hifi_reads.bam` | 37 + 41 GB (~48x) |
| README and md5 | `.../PacBio_HiFi-Revio_20231031/README_HG002-PacBio-Revio.md`, `checksums.md5` | |

```bash
wget -c "$BASE/analysis/NIST_HG002_DraftBenchmark_defrabbV0.020-20250117/GRCh38_HG2-T2TQ100-V1.1_smvar.vcf.gz"{,.tbi}
# ... remaining files as in the table
grep hifi_reads checksums.md5 | md5sum -c -
```
PacBio's own Kinnex HG002 releases (single-cell 10x 5' and full-length bulk RNA) are hosted only on
`downloads.pacbcloud.com`; see `data-sources.md` if that host is unreachable from your cluster.

### 2. Notes for use as a baseline
- The call set is assembly-based, so phasing is genome-wide (one phase set per chromosome), which is
  exactly what haplotype construction needs.
- Outside the benchmark BEDs the calls are lower confidence. Decision for IGI-SYN-SEQ: use all calls,
  but tag designed somatic events by whether they fall inside the benchmark regions (a tier axis).
- Sex: male. chrX outside the PARs and chrY are haploid; the VCF encodes this.
- HLA alleles (clinical, trio-phased; Chin et al. 2020 Nat Commun, Supplementary Table 4):
  maternal A*01:01:01G / B*35:08:01G / C*04:01:01G / DRB1*10:01:01G / DQA1*01:01:01G / DQB1*05:01:01G,
  paternal A*26:01:01G / B*38:01:01G / C*12:03:01G / DRB1*04:02:01 / DQA1*03:01:01G / DQB1*03:02:01G.
- Chromosome naming: `chr`-prefixed, matching `Homo_sapiens.assembly38`.

---

## IGI-SYN-SEQ-02: IPISRC044 (public sarcoma case, male)

Only exome-restricted, short-read-phased germline calls exist in the LENS outputs (156k records,
phase-block N50 138 bp), so a genome-wide call set is built from the blood-normal Illumina WGS.
**There is no long-read DNA for IPISRC044.** The files under `assays/ONT/` are ONT reads of the 10x 5'
single-cell cDNA libraries (median 291 bp, polyA and 10x adapters, exon depth in the 10^5 range vs ~7x
intergenic). They add read-backed phase only inside expressed genes; genome-wide phase therefore comes
from SHAPEIT5 statistical phasing, with WhatsHap short-read blocks as a scaffold.

### 1. Inputs
- `$NORMAL_BAM`: blood-normal Illumina WGS, bwa-mem aligned and sorted against `$REF` (the LENS
  alignment is reused; any equivalent alignment works). Index required.
- `$LR_FASTQ/` (optional): long-read **DNA** if available. For IPISRC044 there is none; the ONT cDNA reads
  were aligned once (32 of 97 files, 19x nominal but concentrated in exons) and contribute little.

### 2. Genome-wide small-variant calling
```bash
run_deepvariant --model_type WGS --ref $REF --reads $NORMAL_BAM \
  --output_vcf $OUT/SAMPLE.deepvariant.vcf.gz --output_gvcf $OUT/SAMPLE.deepvariant.g.vcf.gz \
  --num_shards $N --intermediate_results_dir $SCRATCH/dv
```
No `--regions`: the LENS run used the WES model restricted to the exome BED, which is what made it
unusable here. Expect roughly 4-5 M records for a 30x genome. Wall time is several hours at 64 CPUs.

### 3. Long-read alignment (optional; only useful with genomic long reads; no ONT data is simulated)
```bash
ls $LR_FASTQ/*.fastq.gz | sort -V | head -n 32 > files.txt
cat $(cat files.txt) \
 | minimap2 -t $N -ax map-ont --secondary=no -R '@RG\tID:ONT_T1\tSM:SAMPLE\tPL:ONT' $REF - \
 | samtools sort -@4 -m 2G -o $OUT/SAMPLE.ont.sorted.bam -
samtools index $OUT/SAMPLE.ont.sorted.bam
samtools coverage $OUT/SAMPLE.ont.sorted.bam     # confirm ~15-25x on autosomes
```
Use `lr:hq` instead of `map-ont` with minimap2 >= 2.27 for R10 super-accuracy reads.

### 4. Phasing
Per chromosome (WhatsHap is single-threaded, so run chromosomes as an array), then concatenate:
```bash
bcftools view -f PASS -r $CHR -Oz -o SAMPLE.$CHR.pass.vcf.gz $OUT/SAMPLE.deepvariant.vcf.gz && bcftools index -t SAMPLE.$CHR.pass.vcf.gz
whatshap phase --ignore-read-groups --indels --chromosome $CHR --reference $REF \
  -o SAMPLE.$CHR.phased.vcf SAMPLE.$CHR.pass.vcf.gz $NORMAL_BAM $OUT/SAMPLE.ont.sorted.bam
bgzip SAMPLE.$CHR.phased.vcf && tabix -p vcf SAMPLE.$CHR.phased.vcf.gz
whatshap stats --chromosome $CHR SAMPLE.$CHR.phased.vcf.gz
# after all chromosomes:
bcftools concat -Oz -o $OUT/SAMPLE.deepvariant.pass.phased.vcf.gz SAMPLE.chr{1..22}.phased.vcf.gz SAMPLE.chrX.phased.vcf.gz SAMPLE.chrY.phased.vcf.gz
bcftools index -t $OUT/SAMPLE.deepvariant.pass.phased.vcf.gz
whatshap stats --tsv $OUT/phasing.stats.tsv $OUT/SAMPLE.deepvariant.pass.phased.vcf.gz
```
`--ignore-read-groups` lets the tumor ONT BAM (different sample name) contribute reads.
WhatsHap requires each BAM's index to sit next to the BAM as `<name>.bam.bai`; if an index lives
elsewhere (Nextflow work dirs do this), stage BAM and index as side-by-side symlinks first.

### 5. Acceptance criteria
- After SHAPEIT5: every heterozygous PASS site phased, one phase set per chromosome.
- Read-backed WhatsHap stage alone (short reads): expect ~60 % of hets phased in blocks of a few hundred bp
  (observed chr8: 62 %, median block 126 bp, longest 36 kb); this is the scaffold, not the product.
- chrX (male) has essentially no heterozygous calls outside the PARs; chrY none.

### 6. Additional steps (decided 2026-09-24)
- **Unphased hets and block boundaries.** WhatsHap yields many blocks per chromosome. For genome
  construction, either (a) statistically phase across blocks with SHAPEIT5 against the 1000 Genomes
  GRCh38 panel and use WhatsHap's read-backed phase as a scaffold, or (b) assign block orientations
  at random and record the choice in the truth bundle. **Decision: (a), SHAPEIT5 5.1.1** `phase_common`
  per autosome against the 1000 Genomes 30x GRCh38 panel (3,202 samples, NYGC 2022-04-22 release from
  the EBI FTP; download the 23 per-chromosome VCFs in parallel, one stream is throttled to ~0.5 MB/s).
  WhatsHap output is **not** used as `--scaffold`: SHAPEIT5 treats a scaffold as chromosome-wide phase,
  which many short read-backed blocks are not. Practical requirements found on the first runs:
  - the target needs `INFO/AC` and `INFO/AN` (DeepVariant writes neither): `bcftools +fill-tags -- -t AC,AN`;
  - biallelic records only: `bcftools norm -m -any -f $REF | bcftools view -e 'ALT=="*"'`;
  - genetic maps in SHAPEIT format (`pos chr cM`): the `shapeit5` GitHub repository was unreachable
    ("repository access blocked"), the identical b38 maps ship in the `shapeit4` repository as
    `maps/genetic_maps.b38.tar.gz`; an HTML page saved as `chrN.b38.gmap.gz` makes SHAPEIT5 segfault
    after `GMAP parsing [n=0]`, so check the file type;
  - convert each panel chromosome to BCF once (`bcftools view -Ob` + index) for speed.
  ```bash
  bcftools view -f PASS -r $CHR calls.vcf.gz | bcftools norm -m -any -f $REF | bcftools view -e 'ALT=="*"' \
    | bcftools +fill-tags -Ob -o target.$CHR.bcf -- -t AC,AN && bcftools index target.$CHR.bcf
  SHAPEIT5_phase_common --input target.$CHR.bcf --reference 1kGP.$CHR.bcf --map $CHR.b38.gmap.gz \
    --region $CHR --thread 8 --output target.$CHR.shapeit5.bcf
  ```
  SHAPEIT5 emits only target sites present in the panel (chr22: 54,906 of 65,778, i.e. 83 %). The rest are
  placed by a combination step: for each remaining heterozygous site, use the WhatsHap block it belongs to,
  oriented to agree with SHAPEIT5 at the block's panel sites; otherwise a seeded random orientation. Every
  record gets `INFO/PHASE_SOURCE` in {shapeit5, whatshap_block, random, hom, haploid}; chrX outside the PARs
  and chrY are emitted haploid for a male. One phase set per chromosome results.
  Outcome for IPISRC044 (2026-09-24): 4.58 M PASS records; hets phased by SHAPEIT5 2,541,146, via WhatsHap
  blocks 73,017, random 130,424 (4.7 % of hets), homozygous 1,733,596, haploid chrX/chrY 109,641; zero
  unphased hets; `whatshap stats` reports one block per chromosome (chr1: 100 % phased, 248.9 Mb block).
- **Germline SVs.** With genomic long reads use sniffles2 (`sniffles --input lr.bam --vcf sv.vcf.gz
  --reference $REF --phase`). IPISRC044 has none, so SVs come from Manta 1.6.0 in germline mode on the
  short-read normal WGS, restricted to the primary chromosomes (lower insertion sensitivity accepted).
  **Decision: include**, then annotate with SnpEff (GENCODE v37 database) so the truth bundle documents
  each germline SV.
  ```bash
  awk -v OFS='\t' '$1 ~ /^chr([0-9]+|X|Y)$/ {print $1,0,$2}' $REF.fai | bgzip > primary.bed.gz && tabix -p bed primary.bed.gz
  configManta.py --bam $NORMAL_BAM --referenceFasta $REF --callRegions primary.bed.gz --runDir manta
  manta/runWorkflow.py -m local -j $N -g 60            # results/variants/diploidSV.vcf.gz
  bcftools view -f PASS diploidSV.vcf.gz | snpEff -dataDir $SNPEFF_DATA GRCh38.GENCODEv37 - | bgzip > germline_sv.snpeff.vcf.gz
  ```
  Manta needs the BAM index beside the BAM (stage symlinks as for WhatsHap). It ran the whole genome in a
  single candidate-generation bin despite 32 CPUs, so wall time was 3 h 08 min on one thread.
  Outcome for IPISRC044: 8,756 PASS SVs (4,510 DEL, 2,299 INS, 1,384 BND, 563 DUP); SnpEff finds 201
  gene-fusion-plus-frameshift and 80 gene-fusion consequences among them. Files: `germline_wgs/manta/`.
  For HG002 the Q100 `stvar` file contains all 6.27 M variants; keep the 46,504 records with
  `abs(strlen(ALT)-strlen(REF)) >= 50` (GIAB tandem-repeat annotations retained) and annotate the same
  way: 46 % intergenic, 25 % intronic, 77 transcript ablations, 60 exon losses. Files: `inputs/vcfs/HG002/germline/`.
- **HLA alleles.** From the manifest: A*01:01 homozygous, B*08:01/B*27:05, C*01:02/C*07:01.

---

## Provenance of the current baselines
- IGI-SYN-SEQ-01 (HG002) files were fetched from the GIAB mirror on 2026-09-24 and md5-verified.
- IGI-SYN-SEQ-02 (IPISRC044) calls were produced on 2026-09-24 with the versions above (DeepVariant 64 shards,
  2 h 45 min wall; minimap2 2.22 `map-ont` on 32 ONT files giving 19.2x autosomal depth, 12x chrX, 3 h;
  WhatsHap 2.4 per chromosome). Results, statistics and the site-specific
  paths are recorded in the untracked site config and the run logs, not in this repository.

## MHC recovery for an alt-aware normal (owner decision D3)

IPISRC044's blood normal was aligned alt-aware, which puts HLA reads on the 525 `HLA-*` contigs and leaves
the primary chr6 MHC at 7-11x. That is too shallow to genotype the class I loci, and the class I loci are
exactly where the designed HLA LOH has to be visible. If a baseline is built from such a BAM without this
step, the MHC is silently under-genotyped and the HLA LOH event has no heterozygous sites to be measured
against.

`jobs/03_mhc_recover.sbatch` does the recovery in five steps:

1. Build a primary-only reference: keep the primary chromosomes and the unplaced/random/decoy contigs, drop
   every `*_alt` and `HLA-*` contig so their reads have nowhere else to go. For GRCh38 this keeps 2,580 of
   3,366 contigs.
2. Extract the reads that could belong in the MHC — those parked on the HLA and alt contigs, those already on
   chr6:28-34 Mb, and the unplaced unmapped pool.
3. Realign them to the primary-only reference.
4. Recall chr6:28-34 Mb with DeepVariant.
5. Verify, and **assert** rather than print.

Measured on IPISRC044 (job 11552481, 1h26m):

| | before | after |
|---|---|---|
| HLA-A depth | 7-11x | 40.8x |
| HLA-B depth | 7-11x | 48.9x |
| PASS records in chr6:28-34 Mb | — | 7,532 |
| heterozygous | — | 3,991 |

Per locus: HLA-A 32 PASS / 0 het, HLA-B 93/64, HLA-C 105/100.

**HLA-A has no heterozygous calls because IPISRC044 is homozygous there**, `A*01:01/A*01:01`, which
`design.yaml` records. That is the correct result, not a failure of the recovery, and it has a consequence for
the benchmark: HLA-A LOH cannot be detected by allelic imbalance in IGI-SYN-SEQ-02, so a caller that reports
it there is wrong and a caller that stays silent is right. IGI-SYN-SEQ-01 is heterozygous at all three class I
loci and is the dataset that exercises HLA-A LOH detection. `hla_loss.py` already handles this correctly: for
dataset 02 it reports only `HLA-B*27:05` and `HLA-C*01:02` as lost, and does not claim HLA-A is.

One trap, because it nearly buried the whole result: the first version of step 5 ran its counts inside the
bwa/samtools image, which has no `bcftools`. Every count printed `command not found` followed by `0`, so the
job log said zero records had been recovered when 7,532 had. The counts now run in the bcftools image and are
asserted -- a recovery that genotypes under 1,000 PASS or under 500 heterozygous sites fails the job rather
than completing quietly.
