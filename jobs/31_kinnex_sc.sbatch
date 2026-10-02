#!/bin/bash
#SBATCH -J igi_kinnex_sc
#SBATCH -c 16 --mem=320G -t 72:00:00
#SBATCH --array=1-2
# Kinnex single-cell RNA (PacBio MAS-seq, 10x 5' v2 cDNA), one task per dataset. Produces the pre-skera
# hifi_reads BAM, the segmented BAM from the real `skera split`, the cell roster and the clonotype table.
#
# The roster and the molecule pool are built from --seed, and every single-cell assay must use the SAME seed
# or a barcode stops meaning one cell across assays. 10x GEX, 10x TCR and ONT single cell, when written,
# take this same seed.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
O=$IGI_WORK/release/kinnex_sc; W=$IGI_WORK/release/kinnex_sc_work
d=${DATASETS[$((SLURM_ARRAY_TASK_ID-1))]}
CH=$(IFS=,; echo "${CHROMS_ALL[*]}")
mkdir -p "$O/$d"
cd "$IGI_REPO/catalog"
python3 -u run_kinnex_sc.py --design design.yaml --paths "$IGI_PATHS" --dataset "$d" --chroms "$CH" \
  --cells 4000 --segments 12000000 --min-tpm 0.01 --jobs 8 --label full \
  --out "$O/$d" --work "$W" --seed 9100
echo "DONE $d $(date)"
