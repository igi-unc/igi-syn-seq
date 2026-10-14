#!/bin/bash
#SBATCH -J igi_tenx_gex
#SBATCH -c 8 --mem=24G -t 48:00:00
#SBATCH --array=1-2
# 10x Chromium 5' v2 gene expression, one task per dataset.
#
# The seed MUST be 9100, the same as every other single-cell assay. The cell roster and the molecule pool
# are derived from it, and a different seed gives a different roster -- at which point a barcode no longer
# means one cell across assays, which is the property the whole single-cell design exists to provide.
# Verified on a 40-cell run: identical barcode sets, 99.5% of Kinnex molecules present in GEX, and no
# (barcode, UMI) pair mapping to more than one transcript.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
O=$IGI_RELEASE/tenx_gex; W=$IGI_RELEASE/tenx_gex_work
d=${DATASETS[$((SLURM_ARRAY_TASK_ID-1))]}
CH=$(IFS=,; echo "${CHROMS_ALL[*]}")
mkdir -p "$O/$d"
# A candidate release root is populated with symlinks into the live release for the arms that are
# not being rebuilt, and a write follows a link whichever direction you reasoned about it from.
# See assert_write_path in lib.sh for what this cost on 2026-10-14.
source "$(dirname "${IGI_ENV}")/lib.sh"
assert_write_path "$O/$d"
cd "$IGI_REPO/catalog"
python3 -u run_tenx_gex.py --design design.yaml --paths "$IGI_PATHS" --dataset "$d" --chroms "$CH" \
  --cells 4000 --reads-per-cell 40000 --mean-molecules 8000 --min-tpm 0.01 --label full \
  --out "$O/$d" --work "$W" --seed 9100
echo "DONE $d $(date)"
