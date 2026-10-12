#!/bin/bash
#SBATCH -J igi_align_wgs_n
#SBATCH -c 24 --mem=64G -t 24:00:00
#SBATCH --array=1-2
# Align the merged NORMAL short-read WGS library, which nothing in the pipeline did.
#
# 22_accept_wgs aligns the tumour itself and 13_align covers the exomes, so wgs_bams held one BAM per
# dataset and acceptance only ever opened that one. That is how every IGI-SYN-SEQ-02 normal library
# shipped at 69 % of its designed depth (build-reference section 24): the arm the defect was in had never
# been aligned, so no check could have measured it even in principle. check_normal_depth now needs this
# BAM, and reports a FAILURE rather than a skip when it is missing, which is why the alignment has to
# exist as a stage rather than as something somebody remembers to run.
#
# SEPARATE FROM 22, not folded into it, because wall clock is the risk there. The tumour alignment plus
# ten chromosomes of acceptance took 21:09 and 21:42 against a 40 h limit; adding ~7 h of normal
# alignment to the same task would spend a third of the remaining margin, and a timeout costs the whole
# 28 h and is retried at the same limit. As its own array task it runs beside the tumour's and restarts
# on its own.
#
# 30x against the tumour's 100x, so this is roughly a third of 22's alignment work.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
M=$IGI_RELEASE/wgs_merged; A=$IGI_RELEASE/wgs_bams
d=${DATASETS[$((SLURM_ARRAY_TASK_ID-1))]}
mkdir -p "$A/$d"
for e in amb ann bwt pac sa; do
  [ -f "$IGI_FASTA_NOALT.$e" ] || { echo "missing bwa index $IGI_FASTA_NOALT.$e"; exit 1; }
done
R1=$M/$d/${d}_normal_wgs_R1.fastq.gz
R2=$M/$d/${d}_normal_wgs_R2.fastq.gz
for f in "$R1" "$R2"; do
  [ -s "$f" ] || { echo "missing merged library $f; run 50_merge_longread.sbatch first"; exit 1; }
done
# The read group's SM is ${d}_normal, matching 13_align and the tumour's ${d}_tumor, so the two BAMs a
# somatic caller is handed name the same sample pair the manifest does.
B=$A/$d/${d}_normal.bam
if [ ! -s "$B" ] || [ "$R1" -nt "$B" ]; then
  rm -f "$B" "$B.bai"
  singularity exec -B "$IGI_BIND" "$IGI_BWA_IMG" bash -c "
    bwa mem -t 24 -R '@RG\tID:${d}_normal_wgs\tSM:${d}_normal\tPL:ILLUMINA\tLB:${d}_normal_wgs' \
      $IGI_FASTA_NOALT $R1 $R2 2> $A/$d/normal.bwa.err \
      | samtools sort -@6 -m 3G -o $B -" || { tail -20 "$A/$d/normal.bwa.err"; exit 1; }
  singularity exec -B "$IGI_BIND" "$IGI_BWA_IMG" samtools index "$B"
else
  echo "  $B is newer than its library; keeping it"
fi
singularity exec -B "$IGI_BIND" "$IGI_BWA_IMG" samtools flagstat "$B" | head -6
echo "DONE $d normal wgs $(date)"
