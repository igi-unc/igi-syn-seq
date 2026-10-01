#!/bin/bash
# Submit a build stage with the site environment propagated. sbatch exports the caller's environment by
# default, which is how $IGI_ENV reaches the job.
#
#   export IGI_ENV=$PWD/jobs/env.sh
#   ./jobs/submit.sh 01_catalog.sbatch
#   ./jobs/submit.sh 10_wes.sbatch --dependency=afterok:<catalog jobid>
set -euo pipefail
: "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh first}"
[ -f "$IGI_ENV" ] || { echo "IGI_ENV=$IGI_ENV does not exist"; exit 1; }
source "$IGI_ENV"
script=$1; shift
d=$(dirname "${BASH_SOURCE[0]}")
[ -f "$d/$script" ] || { echo "no such stage: $script"; exit 1; }
mkdir -p "$IGI_WORK/logs"
exec sbatch -p "$IGI_PARTITION" -o "$IGI_WORK/logs/%x.%A_%a.out" "$@" "$d/$script"
