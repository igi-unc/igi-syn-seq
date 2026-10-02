#!/bin/bash
#SBATCH -J igi_wgs
#SBATCH -c 8 --mem=48G -t 24:00:00
#SBATCH --array=1-96
# Illumina short-read WGS at 30x, tumour and normal (owner decision 2026-10-01). Same array layout as
# 10_wes.sbatch: 24 chromosomes x 2 libraries x 2 datasets.
#
# No GC curve here. The measured curve is CAPTURE efficiency, so applying it without a capture would invent
# a bias the assay does not have; run_wgs.py uses a flat GcBias(None) deliberately.
#
# Yield is set by callable bases, not chromosome length: chr21 is 18.6 % N, and the chr21 validation
# delivered 0.974 of the callable-base expectation. Do not read a shortfall against raw length as a bug.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
O=$IGI_WORK/release/wgs; W=$IGI_WORK/release/wgs_work
LB=(tumor normal)
i=$((SLURM_ARRAY_TASK_ID-1))
c=${CHROMS_ALL[$((i % 24))]}; j=$((i / 24)); l=${LB[$((j % 2))]}; d=${DATASETS[$((j / 2))]}
mkdir -p "$O/$d"
cd "$IGI_REPO/catalog"
# Tumour at 100x, normal at 30x. 30x tumour was a germline-grade depth carried into the tumour arm, and
# it cannot see the design's own subclones: at 30x the truncal clone yields 7-10 alt reads and every
# subclone under 5, so clone A at CCF 0.40 gives 4.2 and clone A1 at CCF 0.10 gives 0.7. The dataset exists
# to test subclonal neoantigen calling, so the depth has to support it. At 100x the subclones give 4-14 alt
# reads and the truncal clone 22-35. Normal stays at 30x, which is standard for a germline reference.
DEPTH=100; [ "$l" = normal ] && DEPTH=30
python3 run_wgs.py --design design.yaml --paths "$IGI_PATHS" --dataset "$d" --chrom "$c" --library "$l" \
  --depth $DEPTH \
  --catalog output/$d.snv_indel.tsv \
  --background output/$d.background.tsv \
  --extra-events output/$d.splice_causal.tsv \
  --out "$O/$d" --work "$W" --seed $((7000 + SLURM_ARRAY_TASK_ID))
echo "DONE $d $c $l $(date)"
