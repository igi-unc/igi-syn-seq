# External data sources for IGI-SYN-SEQ

Path placeholders: `$RAFT` = RAFT workspace root, `$DATASETS` = shared lab datasets directory.

Where a source has a stable public URL it is given. Where it does not -- PacBio's dataset portal reorganises
and several of these came from an institutional mirror -- the portal root and the exact dataset identifier
are given instead, which is what makes it findable. Nothing here is a reconstructed guess at a URL.

| Purpose | Source | Location on disk | Status (2026-10-01) |
|---|---|---|---|
| HG002 phased germline baseline (IGI-SYN-SEQ-01) | GIAB HG2-T2TQ100-V1.1 GRCh38 smvar + stvar (defrabb V0.020), from the GIAB FTP at `ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/` | `$DATASETS/HG002/GIAB_Q100_benchmark_defrabbV0.020/` | downloaded 2026-09-24 |
| HiFi WGS error model | GIAB HG002 Revio HiFi raw movies (m84039_230928_213653_s3, m84039_231005_222902_s1) | `$DATASETS/HG002/PacBio_HiFi-Revio_20231031/` | downloaded 2026-09-24 |
| Kinnex sc error / MAS-seq model | PacBio DATA-Revio-Kinnex-HG002-10x5p (pre-skera BAM + matched 10x Illumina) | `$RAFT/inputs/bams/PacBio_Kinnex_scRNA_HG002_10x5p_2024/` | present; obtained via the same institutional mirror after downloads.pacbcloud.com proved unreachable from the cluster |
| Kinnex bulk error / MAS-seq array model | PacBio Kinnex-full-length-RNA DATA-Revio-HG002-1 (post-skera `segmented.bam`, 38.7 M reads, with skera `ds`/`di`/`dl`/`zm` tags that let the 8-segment arrays be reconstructed per ZMW; `flnc.bam`, 37.2 M reads) | `$DATASETS/HG002/Kinnex_fullLengthRNA_HG002_Revio/` | copied 2026-09-24 from a local institutional mirror of downloads.pacbcloud.com |
| Extra bulk long-read RNA (HG002) | GIAB RNA effort MAS-seq, GM24385 refined FLNC FASTQ + lima/refine reports | `$DATASETS/HG002/GIAB_MASseq_RNA_GM24385/` | copied 2026-09-24, same mirror |
| IPISRC044 germline baseline (IGI-SYN-SEQ-02) | DeepVariant WGS on blood-normal BAM + WhatsHap (short reads) + SHAPEIT5 (1000G panel) | `$RAFT/inputs/vcfs/IPISRC044/germline_wgs/` | complete: `IPISRC044.germline.phased.final.vcf.gz`, 95 MB, checksummed in `catalog-freeze.txt` |
| IPISRC044 short-read error models | real IPISRC044 WES / RNA / 10x FASTQs | `$RAFT/inputs/fastqs/IPISRC044/` | present |
| ONT long-read single-cell RNA (IPISRC044 T1) | ONT PromethION R10.4.1 sup reads of the T1 10x 5' cDNA library (~58 % carry the 10x Read 1 adapter, 53 % the 5' TSO, 82 % polyA; median 291 bp) — a real long-read/short-read pair on one 10x library; not genomic, so unusable for phasing or SVs | `$RAFT/inputs/fastqs/IPISRC044/assays/ONT/` | present (T1, T2, T3) |
| ONT genomic error / length model (R10.4.1 SUP) | ONT open data, `giab_2025.01/basecalling/sup/HG002/PAW70337/calls.sorted.bam` at `https://ont-open-data.s3.eu-west-1.amazonaws.com/` | read remotely, not mirrored | present; see note below |
| Immune/stromal single-cell profiles | IPISRC044 UCSF SCG1 (10x 5') | `$RAFT/inputs/fastqs/IPISRC044/assays/ucsf/T*/FASTQ/` | present |
| Tumor expression baseline | UCSC Xena Toil `tcga_Kallisto_tpm` (transcript log2 TPM, all TCGA; `toil.xenahubs.net`) + PanCanAtlas subtypes `tcga_subtypes_syn8402849.csv` (BRCA PAM50 in `Subtype_mRNA`) | `$DATASETS/ucsc_xena/`, `$DATASETS/TCGA/tcga_subtypes_syn8402849/` | present |
| ERV / CTA pan-normal expression reference | lens-v2.0.0-dev-alt `pan_normal_and_mtec_exp.95th_perc.homo_sapiens.quant.sf`, `erv_pep_exp_norm_and_mtec.homo_sapiens.tsv` | `$RAFT/references/homo_sapiens/{expression,erv}/` | present (to be regenerated later per issue 31) |
| Statistical phasing panel | 1000 Genomes 30x GRCh38 phased panel (for SHAPEIT5), `20220422_3202_phased_SNV_INDEL_SV` under `ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/` | `$DATASETS/1000G_GRCh38_phased/` | present, 53 GB |

GIAB's own HG002 ONT data is all R9.4-era -- 2D reads from 2016 basecalled with Guppy v2/v3 -- so it cannot
calibrate an R10.4.1 model. The ONT open-data bucket is used instead, and only the length and accuracy
statistics are taken from it (19,117 bp mean, sd 15,530, accuracy 0.98676); the reads themselves are not
copied.

Note: LENS's own germline VCFs for IPISRC044 (lens-v2.0.0-dev, IPISRC044-WGS project) are
exome-restricted (`--regions hg38_exome.50bpflank.bed`, 156k records) and phased from
short reads only (block N50 138 bp), so they are not usable as a genome-wide phased baseline.

## Containers (all pulled into `$IMAGE_DIR` with `singularity pull docker://...`)

| Tool | Image |
|---|---|
| DeepVariant 1.8.0 | `google/deepvariant:1.8.0` |
| WhatsHap 2.4 | `quay.io/biocontainers/whatshap:2.4--py310h184ae93_0` |
| SHAPEIT5 5.1.1 | `quay.io/biocontainers/shapeit5:5.1.1--h34261f4_2` |
| sniffles 2.8.0 | `quay.io/biocontainers/sniffles:2.8.0--pyhdfd78af_1` |
| minimap2 2.28 / samtools 1.21 / bcftools 1.21 | `quay.io/biocontainers/{minimap2:2.28--h577a1d6_4,samtools:1.21--h50ea8bc_0,bcftools:1.21--h8b25389_0}` |
| Badread 0.4.1 | `badread_0.4.1.sif` — **the long-read simulator actually used**, for both PacBio HiFi and ONT |
| pbsim3 3.0.5 | `quay.io/biocontainers/pbsim3:3.0.5--h9948957_2` — **rejected**: `--accuracy-mean` is ignored, both methods returning 0.964 regardless of the requested value |
| ART 2016.06.05 | `quay.io/biocontainers/art:2016.06.05--h0704011_13` |
| skera 1.4.0 / pbtk 3.5.0 / lima 2.13 / isoseq 4.3 / pbmm2 1.17 / pigeon 1.4 | `quay.io/biocontainers/{pbskera:1.4.0--hdfd78af_0,pbtk:3.5.0--h9ee0642_0,lima:2.13.0--h9ee0642_0,isoseq:4.3.0--h9ee0642_0,pbmm2:1.17.0--h9ee0642_0,pbpigeon:1.4.0--h9948957_0}` |
| Cell Ranger 9.0.1 | `quay.io/nf-core/cellranger:9.0.1` (community build of the 10x tarball; 10x EULA applies) |
| netMHCpan 4.1b / mhcflurry 2.1.1 | LENS images `spvensko/netmhcpan:4.1b`, `spvensko/mhcflurry:2.1.1` |
