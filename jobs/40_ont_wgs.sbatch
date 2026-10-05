#!/bin/bash
#SBATCH -J igi_ont_wgs
#SBATCH -c 48 --mem=48G -t 48:00:00
#SBATCH --array=1-96
# ONT R10.4.1 WGS at 30x, tumour and normal. Same array layout as 10_wes.sbatch: 24 chromosomes x 2
# libraries x 2 datasets. Uses the same derived-chromosome construction as the PacBio arm, shared through
# igi_catalog/derived.py rather than reimplemented.
#
# Identity is calibrated: Badread's realised value for the ONT genomic model runs ~0.49 points BELOW the
# request, the opposite direction from PacBio's +0.43, so ont_wgs_identity is 99.17 to land on the measured
# 0.98676. The length model was the open risk for this assay and is verified: at a requested 19,117/15,530
# it returns mean 18,554 and sd 15,378, so sd ~ mean is reproduced.
# Memory is 48G against a measured peak of 39.6 G over 158 prior tasks; the mean is 6.0 G but the tail
# is what has to fit. See the note in 21_wgs_pacbio.sbatch for why the old 64G mattered.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
O=$IGI_WORK/release/ont_wgs; W=$IGI_WORK/release/ont_wgs_work
LB=(tumor normal)
i=$((SLURM_ARRAY_TASK_ID-1))
c=${CHROMS_ALL[$((i % 24))]}; j=$((i / 24)); l=${LB[$((j % 2))]}; d=${DATASETS[$((j / 2))]}
mkdir -p "$O/$d"
cd "$IGI_REPO/catalog"
# Chunks are 20 Mb of derived chromosome, so chr1 becomes ~13 per copy and all 48 slots are used. The
# first attempt left each copy whole: 14 chunks of 249 Mb, about 4 hours each against a 1 hour timeout,
# which lost 61 of 96 tasks. The timeout is raised as a second line of defence, not as the fix.
python3 -u run_ont.py --design design.yaml --paths "$IGI_PATHS" --dataset "$d" --assay wgs \
  --chrom "$c" --library "$l" --depth 30 --jobs 48 \
  --out "$O/$d" --work "$W" --seed $((9500 + SLURM_ARRAY_TASK_ID))
echo "DONE $d $c $l $(date)"
