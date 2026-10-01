#!/bin/bash
#SBATCH -J igi_rna
#SBATCH -c 16 --mem=320G -t 72:00:00
#SBATCH --array=1-2
# Bulk RNA, full genome: ERVs, splice isoforms and CTAs reach an RNA library here for the first time,
# which no chr6 slice could exercise.
#
# --min-tpm 0.01, not the 0.5 default. LENS filters on `-p 50 --exclude-zeros`, the median of the
# sample's own non-zero TPM, so truncating the low tail raises that bar and moves every downstream
# expression threshold. Measured on the transcript baseline: a 0.5 floor keeps 26.9 % of non-zero
# transcripts and puts p50 at 0.895 against a real 0.115 (7.8x); 0.01 keeps 88.5 % and p50 is 0.168
# (1.5x). The floor is set so truncation comes from sequencing depth -- at 80M pairs a transcript near
# 0.01 TPM yields under one pair -- rather than from the generator. Memory is high because this keeps
# roughly three times the records.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
O=$IGI_WORK/release/rna; W=$IGI_WORK/release/rna_work
d=${DATASETS[$((SLURM_ARRAY_TASK_ID-1))]}
CH=$(IFS=,; echo "${CHROMS_ALL[*]}")
mkdir -p "$O/$d"
cd "$IGI_REPO/catalog"
python3 run_rna.py --design design.yaml --paths "$IGI_PATHS" --dataset "$d" --chroms "$CH" \
  --label full --pairs 80000000 --min-tpm 0.01 \
  --out "$O/$d" --work "$W" --seed $((6000 + SLURM_ARRAY_TASK_ID))
echo "DONE $d $(date)"
