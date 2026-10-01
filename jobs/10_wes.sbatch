#!/bin/bash
#SBATCH -J igi_wes
#SBATCH -c 8 --mem=48G -t 16:00:00
#SBATCH --array=1-96
# Owner decision D7: full-genome release, no slice. 96 tasks = 24 chromosomes x 2 libraries x 2 datasets.
# The exome arm is kept separate from RNA so an RNA failure costs an RNA regeneration rather than both
# (D7 condition 3).
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
O=$IGI_WORK/release/wes; W=$IGI_WORK/release/wes_work
LB=(tumor normal)
i=$((SLURM_ARRAY_TASK_ID-1))
c=${CHROMS_ALL[$((i % 24))]}; j=$((i / 24)); l=${LB[$((j % 2))]}; d=${DATASETS[$((j / 2))]}
mkdir -p "$O/$d"
cd "$IGI_REPO/catalog"
DEPTH=150; [ "$l" = normal ] && DEPTH=100
python3 run_wes.py --design design.yaml --paths "$IGI_PATHS" --dataset "$d" --chrom "$c" --library "$l" \
  --depth $DEPTH \
  --catalog output/$d.snv_indel.tsv \
  --background output/$d.background.tsv \
  --extra-events output/$d.splice_causal.tsv \
  --out "$O/$d" --work "$W" --seed $((5000 + SLURM_ARRAY_TASK_ID))
echo "DONE $d $c $l $(date)"
