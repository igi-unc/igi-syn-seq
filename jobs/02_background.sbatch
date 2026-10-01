#!/bin/bash
#SBATCH -J igi_background
#SBATCH -c 24 --mem=32G -t 6:00:00
#SBATCH --array=1-2
# Passenger background mutations, genome wide. Runs AFTER 01_catalog so it can reserve the spans the
# designed catalog occupies and never collide with a designed event. ~2 min per dataset.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
d=${DATASETS[$((SLURM_ARRAY_TASK_ID-1))]}
cd "$IGI_REPO/catalog"
python3 run_background.py --design design.yaml --paths "$IGI_PATHS" --dataset "$d" \
  --catalog output/"$d".snv_indel.tsv output/"$d".expressed.tsv --out output/
echo "DONE $d $(date)"
