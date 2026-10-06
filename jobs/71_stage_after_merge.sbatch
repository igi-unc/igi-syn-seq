#!/bin/bash
#SBATCH -J igi_stage
#SBATCH -c 2 --mem=4G -t 1:00:00
set -euo pipefail
source "${IGI_ENV:?}"
IGI_RELEASE=$IGI_WORK/release IGI_RAFT=$IGI_RAFT_WS bash "$IGI_REPO/jobs/70_stage_inputs.sh"
