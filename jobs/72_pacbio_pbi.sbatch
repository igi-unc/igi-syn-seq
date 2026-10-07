#!/bin/bash
#SBATCH -J igi_pbi
#SBATCH -c 4 --mem=16G -t 4:00:00
#SBATCH --array=1-4
# Build the .pbi for each merged PacBio library. samtools merge does not produce one for its output, so
# the delivered libraries had none; the merge job now does this inline, and this exists to backfill the
# libraries that were merged before that was true. ~10 min per 67 GB BAM, which is why it is a job and
# not a foreground command.
set -euo pipefail
source "${IGI_ENV:?}"
PBI=$(python3 -c "import yaml;print(yaml.safe_load(open('$IGI_PATHS')).get('pbindex_cmd',''))")
LB=(tumor normal)
i=$((SLURM_ARRAY_TASK_ID-1))
d=${DATASETS[$((i / 2))]}; l=${LB[$((i % 2))]}
B=$IGI_WORK/release/pacbio_merged/$d/${d}_${l}_hifi.bam
if [ -s "$B.pbi" ] && [ "$B.pbi" -nt "$B" ]; then echo "SKIP $B.pbi is current"; exit 0; fi
eval "${PBI/\{args\}/$B}"
echo "DONE $(basename $B).pbi $(stat -c%s "$B.pbi") bytes"
