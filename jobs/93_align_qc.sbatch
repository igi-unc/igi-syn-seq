#!/bin/bash
#SBATCH -J igi_align_qc
#SBATCH -c 8 --mem=24G -t 8:00:00
#SBATCH --array=1-18
# Align a sample of every delivered library and record what a real aligner makes of it.
#
# WHY THIS IS DIFFERENT FROM THE OTHER CHECKS. Every other check in this release asserts a property
# somebody thought to assert, and each was written AFTER a consumer hit the defect it catches: the `*`
# allele (§15), the IUPAC bases (§17), the unstranded R2 (§19). That pattern cannot catch the defect
# nobody has thought of yet.
#
# Alignment does not test a property. It hands the reads to the same class of tool the consumer will use
# and asks whether they behave like sequencing reads of this genome. A read reverse complemented when it
# should not be, built from the wrong reference, shifted off its coordinates, corrupted mid-record or
# paired with the wrong mate does not announce which invariant it broke -- it simply fails to align, or
# aligns with the wrong orientation, insert size or clipping.
#
# It would have caught the R2 defect outright: the rebuilt 10x TCR measures plus_strand_fraction 0.000,
# and the library Cell Ranger refused measured ~0.5.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"

OUT=$IGI_WORK/release/align_qc
mkdir -p "$OUT"
R=$IGI_WORK/release
BWA="singularity exec -B $IGI_BIND $IGI_REF/imgs/michaelfranklin-bwasamtools-0.7.17-1.10.img bwa {args}"
MM2="singularity exec -B $IGI_BIND $IGI_REF/imgs/quay.io_biocontainers_minimap2_2.28--h577a1d6_4.sif minimap2 {args}"
SM="singularity exec -B $IGI_BIND $IGI_SAMTOOLS_IMG samtools {args}"

# label:kind:reads  -- one task each. Paired libraries give R1,R2; a BAM is sampled through samtools.
# The per-chromosome arms are sampled on chr10, which carries the most IUPAC positions (36) and so is
# the chromosome most likely to show a reference-handling defect.
SPECS=(
  "ds01_wes_tumor:wes:$R/merged/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01_tumor_R1.fastq.gz,$R/merged/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01_tumor_R2.fastq.gz"
  "ds01_wes_normal:wes:$R/merged/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01_normal_R1.fastq.gz,$R/merged/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01_normal_R2.fastq.gz"
  "ds02_wes_tumor:wes:$R/merged/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02_tumor_R1.fastq.gz,$R/merged/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02_tumor_R2.fastq.gz"
  "ds02_wes_normal:wes:$R/merged/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02_normal_R1.fastq.gz,$R/merged/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02_normal_R2.fastq.gz"
  "ds01_wgs_tumor:wgs:$R/wgs/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01_chr10_tumor_wgs_R1.fastq.gz,$R/wgs/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01_chr10_tumor_wgs_R2.fastq.gz"
  "ds01_wgs_normal:wgs:$R/wgs/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01_chr10_normal_wgs_R1.fastq.gz,$R/wgs/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01_chr10_normal_wgs_R2.fastq.gz"
  "ds02_wgs_tumor:wgs:$R/wgs/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02_chr10_tumor_wgs_R1.fastq.gz,$R/wgs/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02_chr10_tumor_wgs_R2.fastq.gz"
  "ds02_wgs_normal:wgs:$R/wgs/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02_chr10_normal_wgs_R1.fastq.gz,$R/wgs/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02_chr10_normal_wgs_R2.fastq.gz"
  "ds01_ont_tumor:ont_wgs:$R/ont_wgs/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01_chr10_tumor_ont_wgs.fastq.gz"
  "ds02_ont_tumor:ont_wgs:$R/ont_wgs/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02_chr10_tumor_ont_wgs.fastq.gz"
  "ds01_hifi_tumor:pacbio:$R/pacbio/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01_chr10_tumor_hifi.bam"
  "ds02_hifi_tumor:pacbio:$R/pacbio/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02_chr10_tumor_hifi.bam"
  "ds01_gex:tenx_gex:$R/tenx_gex/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01-GEX_S1_L001_R2_001.fastq.gz"
  "ds02_gex:tenx_gex:$R/tenx_gex/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02-GEX_S1_L001_R2_001.fastq.gz"
  "ds01_tcr:tenx_tcr:$R/tenx_tcr/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01-TCR_S1_L001_R2_001.fastq.gz"
  "ds02_tcr:tenx_tcr:$R/tenx_tcr/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02-TCR_S1_L001_R2_001.fastq.gz"
  "ds01_ont_bulk_rna:ont_rna:$R/ont_bulk_rna/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01_full_ont_bulk_rna.fastq.gz"
  "ds02_ont_bulk_rna:ont_rna:$R/ont_bulk_rna/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02_full_ont_bulk_rna.fastq.gz"
)

i=$((SLURM_ARRAY_TASK_ID - 1))
[ "$i" -lt "${#SPECS[@]}" ] || { echo "no spec for task $SLURM_ARRAY_TASK_ID"; exit 0; }
spec=${SPECS[$i]}
label=${spec%%:*}; rest=${spec#*:}; kind=${rest%%:*}; reads=${rest#*:}

# Long reads are few and slow to align; short reads are many and fast. Sample accordingly.
case $kind in
  ont_wgs|pacbio|ont_rna) MAX=4000 ;;
  *)                      MAX=50000 ;;
esac
# minimap2 for anything long or spliced, bwa for short genomic DNA.
case $kind in
  ont_wgs|pacbio|ont_rna) ALN=$MM2 ;;
  *)                      ALN=$BWA ;;
esac

echo "=== $label ($kind), $MAX reads ==="
python3 "$(dirname "$IGI_ENV")/measure/align_qc.py" \
  --reads "$reads" --kind "$kind" --ref "$IGI_FASTA_NOALT" --label "$label" \
  --aligner-cmd "$ALN" --samtools-cmd "$SM" --max-reads "$MAX" --threads 8 \
  --out "$OUT/$label.json"
echo "ALIGN QC $label DONE"
