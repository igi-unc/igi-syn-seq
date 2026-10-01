#!/bin/bash
#SBATCH -J igi_catalog
#SBATCH -c 24 --mem=64G -t 12:00:00
#SBATCH --array=1-2
# Stage 1: the designed catalog, one array task per dataset. Deterministic given design.yaml's seed.
# Dominated by netMHCpan; ~1 h per dataset.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
d=${DATASETS[$((SLURM_ARRAY_TASK_ID-1))]}
cd "$IGI_REPO/catalog"
python3 -m igi_catalog.designer --design design.yaml --paths "$IGI_PATHS" --dataset "$d" --out output/
echo "DONE $d $(date)"
