#!/bin/bash
#SBATCH -J igi_pacbio
#SBATCH -c 16 --mem=48G -t 48:00:00
#SBATCH --array=1-96
# PacBio HiFi WGS at 30x via Badread, tumour and normal. Same array layout as 10_wes.sbatch.
#
# Why this is a batch job and not a loop on the login node: running it beside two other jobs dropped
# throughput from 18 chunks/hour to 1. With dedicated cores a thread pool of 8 does ~98/hour.
#
# Badread's realized identity is offset from the requested value by about +0.43 points for PacBio (and
# -0.22 for ONT cDNA, so the offset is not a constant). pacbio_identity in the site paths file is already
# calibrated for that; it is not the target accuracy. Validated on chr21: length 16,788 sd 4,916 and
# accuracy 0.99830, against real HG002 Revio at 16,689 sd 4,593 and 0.99823.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
O=$IGI_WORK/release/pacbio; W=$IGI_WORK/release/pacbio_work
LB=(tumor normal)
i=$((SLURM_ARRAY_TASK_ID-1))
c=${CHROMS_ALL[$((i % 24))]}; j=$((i / 24)); l=${LB[$((j % 2))]}; d=${DATASETS[$((j / 2))]}
mkdir -p "$O/$d"
cd "$IGI_REPO/catalog"
python3 run_pacbio_wgs.py --design design.yaml --paths "$IGI_PATHS" --dataset "$d" --chrom "$c" \
  --library "$l" --depth 30 --jobs 8 \
  --catalog-dir output \
  --background output/$d.background.tsv \
  --extra-events output/$d.splice_causal.tsv \
  --out "$O/$d" --work "$W" --seed $((8000 + SLURM_ARRAY_TASK_ID))
echo "DONE $d $c $l $(date)"
