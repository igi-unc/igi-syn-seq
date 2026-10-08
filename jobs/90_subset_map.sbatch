#!/bin/bash
#SBATCH -J igi_subset_map
#SBATCH -c 8 --mem=16G -t 2:00:00
#SBATCH --array=1-2
# The record -> chromosome map every chr1to6 filter task reads. One pass over the GTF and the transcript
# manifest; a record that cannot be resolved to a chromosome aborts the map rather than being dropped,
# because a silently unresolved source class would shrink the subset without saying so.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
d=${DATASETS[$((SLURM_ARRAY_TASK_ID-1))]}
GTF=$(python3 -c "import sys,yaml;print(yaml.safe_load(open(sys.argv[1]))['gtf'])" "$IGI_PATHS")
mkdir -p "$IGI_WORK/subset"
cd "$IGI_REPO/catalog"
python3 -u subset_chr1to6.py map --dataset "$d" --release "$IGI_WORK/release" \
  --catalog-out output --gtf "$GTF" --out-dir "$IGI_WORK/subset"
echo "DONE $d $(date)"
