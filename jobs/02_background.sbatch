#!/bin/bash
#SBATCH -J igi_background
#SBATCH -c 24 --mem=32G -t 6:00:00
#SBATCH --array=1-2
# Everything that is derived from the designed catalog and consumed by the read builders. Runs AFTER
# 01_catalog; 10_wes, 20_wgs_illumina and 21_wgs_pacbio all read both of its outputs. ~2 min per dataset.
#
# run_splice_causal.py used to have no job at all. It was run by hand, and its output
# `<dataset>.splice_causal.tsv` is a tracked `--extra-events` input to three read arms -- so when the
# catalog moved on 2026-10-13 and the expressed table with it, nothing would have regenerated the splice
# variants and the reads would have carried the previous catalog's splice sites into the new one. That is
# the same shape as the MHC recovery of section 33: a step done by hand whose output nothing re-derives.
# It is first here because it only needs the catalog, while the background must reserve the catalog's
# spans; neither depends on the other.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
d=${DATASETS[$((SLURM_ARRAY_TASK_ID-1))]}
cd "$IGI_REPO/catalog"
python3 run_splice_causal.py --design design.yaml --paths "$IGI_PATHS" --dataset "$d" \
  --expressed output/"$d".expressed.tsv --out output/
python3 run_background.py --design design.yaml --paths "$IGI_PATHS" --dataset "$d" \
  --catalog output/"$d".snv_indel.tsv output/"$d".expressed.tsv --out output/
echo "DONE $d $(date)"
