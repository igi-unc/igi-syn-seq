#!/bin/bash
#SBATCH -J igi_subset
#SBATCH -c 8 --mem=24G -t 12:00:00
#SBATCH --array=1-12
# chr1to6 for the six whole-library arms (design section 12). The other ten sample types already have a
# chr1to6 form: WES, Illumina WGS, PacBio HiFi WGS and ONT WGS are simulated one chromosome at a time, so
# their subset is a selection of shards, and 10x TCR is kept whole.
#
# Twelve tasks, arm-major so the two datasets of one arm run side by side:
#   1,2 bulk RNA   3,4 10x GEX   5,6 Kinnex bulk   7,8 Kinnex sc   9,10 ONT bulk RNA   11,12 ONT sc RNA
#
# Run 90_subset_map.sbatch first: every task reads the record -> chromosome map it writes, and a task that
# cannot find its map fails rather than guessing.
#
# Memory is sized for the Kinnex arms, which hold the whole array table in memory (24.7 M molecules over
# 1.5 M arrays for the single-cell library, as array('i') rather than lists of ints -- the difference is
# 900 MB). The streaming arms need a few hundred MB.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
ARMS=(rna tenx_gex kinnex_bulk kinnex_sc ont_bulk_rna ont_sc_rna)
i=$((SLURM_ARRAY_TASK_ID-1))
arm=${ARMS[$((i/2))]}
d=${DATASETS[$((i%2))]}
M=$IGI_WORK/subset
R=$IGI_WORK/release
map=$M/$d.chr1to6_records.tsv.gz
[ -s "$map" ] || { echo "no record map at $map; run 90_subset_map.sbatch first" >&2; exit 1; }

# The samtools and pbindex commands are site specifics, so they come from the untracked paths file rather
# than being written here. Kinnex needs both: samtools to stream the segmented BAM through the filter and
# pbindex to rebuild the .pbi the design calls for.
read -r SAM PBI <<<"$(python3 - "$IGI_PATHS" <<'PY'
import sys, yaml
p = yaml.safe_load(open(sys.argv[1]))
print(p.get("samtools_cmd", "samtools {args}").replace(" ", "\x01"),
      p.get("pbindex_cmd", "").replace(" ", "\x01"))
PY
)"
SAM=${SAM//$'\x01'/ }; PBI=${PBI//$'\x01'/ }

echo "arm=$arm dataset=$d map=$(basename "$map")"
cd "$IGI_REPO/catalog"
python3 -u subset_chr1to6.py filter --arm "$arm" --dataset "$d" --release "$R" --map "$map" \
  --samtools-cmd "$SAM" --pbindex-cmd "$PBI"
echo "DONE $arm $d $(date)"
