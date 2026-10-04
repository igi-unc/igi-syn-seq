#!/bin/bash
#SBATCH -J igi_ont_bulk
#SBATCH -c 48 --mem=320G -t 48:00:00
#SBATCH --array=1-2
# ONT R10.4.1 bulk cDNA. No MAS array: one molecule is one read, which is the whole difference from the
# Kinnex bulk arm. cDNA model constants are the IPISRC044 R10.4.1 SUP measurements (910/461), distinct from
# the genomic ones.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
O=$IGI_WORK/release/ont_bulk_rna; W=$IGI_WORK/release/ont_bulk_rna_work
d=${DATASETS[$((SLURM_ARRAY_TASK_ID-1))]}
CH=$(IFS=,; echo "${CHROMS_ALL[*]}")
mkdir -p "$O/$d"
cd "$IGI_REPO/catalog"
python3 -u run_ont.py --design design.yaml --paths "$IGI_PATHS" --dataset "$d" --assay bulk_rna \
  --chroms "$CH" --reads 20000000 --min-tpm 0.01 --jobs 48 --label full \
  --out "$O/$d" --work "$W" --seed $((9600 + SLURM_ARRAY_TASK_ID))
echo "DONE $d $(date)"
