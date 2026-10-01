#!/bin/bash
#SBATCH -J igi_kinnex_bulk
#SBATCH -c 16 --mem=320G -t 72:00:00
#SBATCH --array=1-2
# Kinnex bulk RNA (PacBio MAS-seq, 8-mer array), one task per dataset. Produces both the pre-skera
# hifi_reads BAM, where one read is one whole array, and the segmented BAM from the real `skera split`.
#
# Memory matches 12_rna.sbatch: the same transcript record set is assembled, and --min-tpm 0.01 keeps about
# three times the records a 0.5 floor would.
#
# 15 M segmented reads against the real set's 38.7 M. One molecule yields one segment, so --segments is also
# the molecule count; at a mean 6.91 segments per array that is about 2.2 M arrays.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
O=$IGI_WORK/release/kinnex_bulk; W=$IGI_WORK/release/kinnex_bulk_work
d=${DATASETS[$((SLURM_ARRAY_TASK_ID-1))]}
CH=$(IFS=,; echo "${CHROMS_ALL[*]}")
mkdir -p "$O/$d"
cd "$IGI_REPO/catalog"
python3 run_kinnex_bulk.py --design design.yaml --paths "$IGI_PATHS" --dataset "$d" --chroms "$CH" \
  --segments 15000000 --min-tpm 0.01 --jobs 8 --label full \
  --out "$O/$d" --work "$W" --seed $((9000 + SLURM_ARRAY_TASK_ID))
echo "DONE $d $(date)"
