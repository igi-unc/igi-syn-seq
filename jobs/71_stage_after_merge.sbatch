#!/bin/bash
#SBATCH -J igi_stage
#SBATCH -c 2 --mem=4G -t 1:00:00
set -euo pipefail
source "${IGI_ENV:?}"
# IGI_RELEASE comes from env.sh (or the submitting environment, which sbatch propagates), so staging a
# candidate release is `IGI_RELEASE=... ./jobs/submit.sh 71_stage_after_merge.sbatch` and not an edit here.
# It is echoed because this stage REPLACES the symlinks the RAFT workspace serves, and the one thing worth
# seeing in the log afterwards is which root it pointed them at.
echo "staging $IGI_RELEASE -> $IGI_RAFT_WS/inputs"
IGI_RAFT=$IGI_RAFT_WS bash "$IGI_REPO/jobs/70_stage_inputs.sh"
