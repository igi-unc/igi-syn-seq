#!/bin/bash
#SBATCH -J igi_preflight
#SBATCH -c 8 --mem=48G -t 8:00:00
#SBATCH --array=1-14
# Push the first ~5 % of each delivered library through a tool LENS will actually use, and fail if the
# tool chokes on it.
#
# WHY THIS EXISTS. All three defects that reached a consumer were found by a consumer: the `*` allele by
# Cell Ranger, the IUPAC bases by OptiType's razers3, the unstranded R2 by Cell Ranger's chemistry
# detector. Every one had passed everything this release asserted at the time. The checks in
# qa_release.py each encode a property somebody thought to assert; a real tool's parser encodes every
# property its authors thought to assert, which is a much larger set and not one I have to enumerate.
#
# WHY 5 % AND NOT THE WHOLE FILE. The question is whether the consumer ACCEPTS the format, which the
# first bad record settles -- razers3 aborted while loading and never reached read two million. What 5 %
# canNOT do is prove the absence of a rare offender: the `Y` that stopped razers3 was one base in
# 4,858,789,950, and a 5 % sample would have missed it 95 % of the time. Absence is the complete
# alphabet scan's job (92_scan_alphabet, §17.6.2). The two are complementary:
#
#     92_scan_alphabet   every base, one cheap property        -- is it ACGTN
#     94_consumer_preflight   5 % of the bases, the real tool's entire opinion
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"

OUT=$IGI_WORK/release/preflight
mkdir -p "$OUT"
R=$IGI_WORK/release
I=$IGI_REF/imgs
# Consumer reference bundles come from the environment, not from a path baked in here.
P=${IGI_LENS_REFS:?set IGI_LENS_REFS in jobs/env.sh to the consumer reference directory}

# 48G is set by STAR alone: the CTAT index is 26.5 GiB of Genome + SA before any reads, and STAR was
# OOM-killed loading it under 16G. Every other tool here needs a fraction of that, but SLURM allocates
# per task, so the heaviest sets the request -- recorded here so it is not trimmed back by someone
# looking only at razers3.
SM="singularity exec -B $IGI_BIND $IGI_SAMTOOLS_IMG samtools {args}"
RAZERS="singularity exec -B $IGI_BIND $I/quay.io-biocontainers-optitype-1.3.5--0.img razers3 {args}"
SALMON="singularity exec -B $IGI_BIND $I/quay.io-biocontainers-salmon-1.11.4--h7f96273_0.img salmon {args}"
PBMM2="singularity exec -B $IGI_BIND $I/quay.io-biocontainers-pbmm2-1.17.0--h9ee0642_0.img pbmm2 {args}"
LIMA="singularity exec -B $IGI_BIND $I/quay.io-biocontainers-lima-2.13.0--h9ee0642_0.img lima {args}"
CRVDJ="singularity exec -B $IGI_BIND $I/nfcore-cellranger-7.1.0.img cellranger {args}"
STAR="singularity exec -B $IGI_BIND $I/quay.io-biocontainers-star-2.7.10a--h9ee0642_0.img STAR {args}"
# The STAR index that matches the GENCODE v37 annotation LENS uses, inside the CTAT bundle.
# versionGenome 2.7.4a, which STAR 2.7.10a accepts.
STARIDX=${IGI_STAR_INDEX:?set IGI_STAR_INDEX in jobs/env.sh to a STAR genome directory}

# label : tool : reads : reference
SPECS=(
  "ds01_wes_normal:razers3:$R/merged/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01_normal_R1.fastq.gz:"
  "ds01_wes_tumor:razers3:$R/merged/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01_tumor_R1.fastq.gz:"
  "ds02_wes_normal:razers3:$R/merged/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02_normal_R1.fastq.gz:"
  "ds02_wes_tumor:razers3:$R/merged/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02_tumor_R1.fastq.gz:"
  "ds01_rna_star:star:$R/rna/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01_full_rna_R1.fastq.gz,$R/rna/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01_full_rna_R2.fastq.gz:$STARIDX"
  "ds02_rna_star:star:$R/rna/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02_full_rna_R1.fastq.gz,$R/rna/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02_full_rna_R2.fastq.gz:$STARIDX"
  "ds01_rna:salmon:$R/rna/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01_full_rna_R1.fastq.gz,$R/rna/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01_full_rna_R2.fastq.gz:$P/GRCh38.GENCODEv37"
  "ds02_rna:salmon:$R/rna/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02_full_rna_R1.fastq.gz,$R/rna/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02_full_rna_R2.fastq.gz:$P/GRCh38.GENCODEv37"
  "ds01_tcr:cellranger_vdj:$R/tenx_tcr/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01-TCR_S1_L001_R1_001.fastq.gz,$R/tenx_tcr/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01-TCR_S1_L001_R2_001.fastq.gz:$P/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0"
  "ds02_tcr:cellranger_vdj:$R/tenx_tcr/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02-TCR_S1_L001_R1_001.fastq.gz,$R/tenx_tcr/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02-TCR_S1_L001_R2_001.fastq.gz:$P/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0"
  "ds01_hifi:pbmm2:$R/pacbio/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01_chr10_tumor_hifi.bam:$IGI_FASTA_NOALT"
  "ds02_hifi:pbmm2:$R/pacbio/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02_chr10_tumor_hifi.bam:$IGI_FASTA_NOALT"
  "ds01_kinnex:lima:$R/kinnex_bulk/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01_full_kinnex_bulk_segmented.bam:$P/10x_5kit_primers.fasta"
  "ds02_kinnex:lima:$R/kinnex_bulk/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02_full_kinnex_bulk_segmented.bam:$P/10x_5kit_primers.fasta"
)

i=$((SLURM_ARRAY_TASK_ID - 1))
[ "$i" -lt "${#SPECS[@]}" ] || { echo "no spec for task $SLURM_ARRAY_TASK_ID"; exit 0; }
spec=${SPECS[$i]}
label=${spec%%:*}; r1=${spec#*:}; tool=${r1%%:*}; r2=${r1#*:}; reads=${r2%%:*}; ref=${r2#*:}

case $tool in
  razers3)        CMD=$RAZERS ;;
  star)           CMD=$STAR ;;
  salmon)         CMD=$SALMON ;;
  pbmm2)          CMD=$PBMM2 ;;
  lima)           CMD=$LIMA ;;
  cellranger_vdj) CMD=$CRVDJ ;;
  *) echo "unknown tool $tool"; exit 1 ;;
esac

echo "=== $label via $tool (first 5 %) ==="
python3 "$(dirname "$IGI_ENV")/measure/consumer_preflight.py" \
  --tool "$tool" --label "$label" --reads "$reads" --ref "$ref" \
  --cmd "$CMD" --samtools-cmd "$SM" --fraction 0.05 --threads 8 \
  --out "$OUT/$label.$tool.json"
echo "PREFLIGHT $label DONE"
