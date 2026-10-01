#!/bin/bash
#SBATCH -J igi_align
#SBATCH -c 24 --mem=64G -t 48:00:00
#SBATCH --array=1-4
# Align the merged exomes so acceptance runs on the release itself rather than a slice (D7 condition 1).
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
M=$IGI_WORK/release/merged; A=$IGI_WORK/release/bams
LB=(tumor normal)
i=$((SLURM_ARRAY_TASK_ID-1)); d=${DATASETS[$((i/2))]}; l=${LB[$((i%2))]}
mkdir -p "$A/$d"
for e in amb ann bwt pac sa; do
  [ -f "$IGI_FASTA_NOALT.$e" ] || { echo "missing bwa index $IGI_FASTA_NOALT.$e"; exit 1; }
done
B=$A/$d/${d}_${l}.bam
if [ ! -s "$B" ]; then
  # bwa's stderr is kept: suppressing it once turned a missing index into an unreadable
  # "failed to read header".
  singularity exec -B "$IGI_BIND" "$IGI_BWA_IMG" bash -c "
    bwa mem -t 24 -R '@RG\tID:${d}_${l}\tSM:${d}_${l}\tPL:ILLUMINA\tLB:${d}_${l}' \
      $IGI_FASTA_NOALT $M/$d/${d}_${l}_R1.fastq.gz $M/$d/${d}_${l}_R2.fastq.gz 2> $A/$d/${l}.bwa.err \
      | samtools sort -@6 -m 3G -o $B -" || { tail -20 "$A/$d/${l}.bwa.err"; exit 1; }
  singularity exec -B "$IGI_BIND" "$IGI_BWA_IMG" samtools index "$B"
fi
singularity exec -B "$IGI_BIND" "$IGI_BWA_IMG" samtools flagstat "$B" | head -6
echo "DONE $d $l $(date)"
