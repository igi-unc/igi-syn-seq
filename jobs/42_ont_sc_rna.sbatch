#!/bin/bash
#SBATCH -J igi_ont_sc
#SBATCH -c 48 --mem=48G -t 48:00:00
#SBATCH --array=1-2
# ONT R10.4.1 single-cell cDNA, 10x 5' v2 structure. Seed MUST be 9100: this is the fourth assay sharing
# the cell roster and molecule pool with Kinnex sc, 10x GEX and 10x TCR, and a different seed gives a
# different roster.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
O=$IGI_WORK/release/ont_sc_rna; W=$IGI_WORK/release/ont_sc_rna_work
d=${DATASETS[$((SLURM_ARRAY_TASK_ID-1))]}
CH=$(IFS=,; echo "${CHROMS_ALL[*]}")
mkdir -p "$O/$d"
cd "$IGI_REPO/catalog"
python3 -u run_ont.py --design design.yaml --paths "$IGI_PATHS" --dataset "$d" --assay sc_rna \
  --chroms "$CH" --cells 4000 --reads 10000000 --mean-molecules 8000 --min-tpm 0.01 \
  --jobs 48 --label full --out "$O/$d" --work "$W" --seed 9100
echo "DONE $d $(date)"
