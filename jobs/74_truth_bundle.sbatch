#!/bin/bash
#SBATCH -J igi_truth
#SBATCH -c 4 --mem=32G -t 8:00:00
#SBATCH --array=1-2
# Assemble the truth bundle design section 11 specifies, and the RAFT manifests beside it.
#
# This is a STAGE rather than a command somebody remembers, for the reason the bundle was missing in the
# first place: every piece of the truth already existed -- catalog tables in the repo, read maps and
# molecule tables and cell rosters on scratch -- and nothing gathered them, so a benchmark user handed
# the staged release had nothing to score against. A deliverable that is only ever produced by hand is
# a deliverable that will be absent from the next release too.
#
# It writes into the RELEASE, not into the RAFT workspace. Section 2 puts the bundle at
# inputs/metadata/<dataset>/truth/, and 70_stage_inputs.sh is what puts it there -- the same split as
# every other deliverable, so a bundle can be built and read before anything the workspace serves
# changes. --manifest-dir is read only: it is where the staged manifests already are, and the bundle
# records their checksums.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
d=${DATASETS[$((SLURM_ARRAY_TASK_ID-1))]}
OUT=$IGI_RELEASE/truth/$d
MAN=${IGI_STAGED_METADATA:-$IGI_RAFT_WS/inputs/metadata/$d}
mkdir -p "$OUT"
cd "$IGI_REPO/catalog"

# The RAFT manifests go into the release too, for the same reason. --release both, because the subset
# mirrors the full release and a manifest for one is not a manifest for the other.
python3 make_lens_manifest.py --design design.yaml --release both \
  --out-dir "$IGI_RELEASE/manifests"

mopt=""
[ -d "$MAN" ] && mopt="--manifest-dir $MAN"
# No --allow-missing. A partial bundle exits non-zero and that is the point: it exists for assembling one
# mid-build, and a release-time run that needs it is a release-time run with a gap to report.
python3 make_truth_bundle.py --dataset "$d" --catalog-dir output --design design.yaml \
  --arms-bed "$IGI_ARMS_BED" --release-dir "$IGI_RELEASE" $mopt --out "$OUT"
echo "DONE $d truth bundle -> $OUT $(date)"
