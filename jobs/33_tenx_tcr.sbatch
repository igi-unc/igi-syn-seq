#!/bin/bash
#SBATCH -J igi_tenx_tcr
#SBATCH -c 16 --mem=96G -t 24:00:00
#SBATCH --array=1-2
# 10x Chromium 5' V(D)J (TCR), one task per dataset. The fifteenth and last sample type.
#
# Seed MUST be 9100: the roster is shared with Kinnex sc, 10x GEX and ONT sc. Verified at the shipping
# scale against the GEX run already on disk -- identical barcode sets and identical T-cell sets, 1,200 T
# cells over 300 clonotypes with a largest clonotype of 99.
#
# Transcripts are assembled V + CDR3 + J + C from real GRCh38 segments read straight from the GTF, because
# the annotation cache keeps 1 of 202 TR segments. A V(D)J caller infers the junction by aligning to
# germline segments, so invented framework would leave every read unassignable.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
O=$IGI_WORK/release/tenx_tcr; W=$IGI_WORK/release/tenx_tcr_work
d=${DATASETS[$((SLURM_ARRAY_TASK_ID-1))]}
mkdir -p "$O/$d"
cd "$IGI_REPO/catalog"
python3 -u run_tenx_tcr.py --design design.yaml --paths "$IGI_PATHS" --dataset "$d" \
  --cells 4000 --reads-per-cell 5000 --molecules-per-cell 120 --label full \
  --out "$O/$d" --work "$W" --seed 9100
echo "DONE $d $(date)"
