#!/bin/bash
#SBATCH -J igi_drive
#SBATCH -c 1 --mem=2G -t 7-00:00:00
# A DRIVER MUST OUTLIVE ITS OWN CHAIN. This was 24:00:00, and 11579857 hit it at 2026-10-11
# 09:45:36 -- exactly 24 h -- four minutes after submitting 22_accept_wgs, which has a 40 h
# limit of its own and was 14 h into a ~21 h run. The driver died with the job it had just
# launched still going, so 14_accept and 15_accept_rna were never submitted at all and the
# chain simply stopped with nothing reporting that it had. Both stages had to be submitted by
# hand. The driver costs 1 CPU and 2 G to sit and wait, the partition's MaxTime is infinite,
# and the chain's own stages now total well over 24 h, so the limit was the cheapest thing
# here and the only one that bounded the whole build.
# Drive the rest of the release to completion, retrying failures instead of stalling on them.
#
# WHY THIS EXISTS. The chain was wired with sbatch `--dependency=afterok`, which is correct but brittle
# in one specific way: when a dependency fails, every downstream stage is left unsatisfiable and the
# pipeline stops where it broke. Unattended, that turns one OOM-killed chromosome into a lost night and a
# morning state of "stage 2 of 10 done, nothing retried".
#
# Today's evidence says that is the likely case rather than the unlucky one: the exome merge failed 4/4
# on its first attempt (a sourcing bug), the alphabet scan reported 96 COMPLETED having skipped 36 files,
# and the Illumina WGS arm had to be cancelled and resubmitted mid-flight. A clean unattended run is not
# this build's pattern.
#
# So this does what a human would: wait for a stage, look at what failed, resubmit JUST the failed tasks
# with more memory if they were OOM-killed, and only move on when the stage is genuinely clean. It gives
# up after RETRIES attempts and says so loudly rather than carrying a broken stage forward.
set -uo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
J=$(dirname "$IGI_ENV")
LOG=$IGI_WORK/logs/drive.$SLURM_JOB_ID.log
RETRIES=${RETRIES:-2}

say() { echo "[$(date '+%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

# Block until a job id has no task running, pending or completing.
wait_job() {
  local id=$1 n
  while :; do
    n=$(squeue -j "$id" -h -o %t 2>/dev/null | grep -cE 'R|PD|CG')
    [ "${n:-0}" -eq 0 ] && break
    sleep 60
  done
}

# JobID, NOT JobIDRaw. JobIDRaw prints the raw numeric id of each array task --
#   11578052|FAILED
# -- while JobID prints the array notation this needs:
#   11577784_1|FAILED
# The first version used JobIDRaw and would have resubmitted `--array=11578052`. Tested against a job
# with known failures before this was trusted with a night.
RETRYABLE='FAILED|OUT_OF_MEMORY|TIMEOUT|NODE_FAIL|DEADLINE'

# How many tasks did not succeed, for any reason. 0 means the stage is genuinely clean.
n_failed() {
  sacct -j "$1" -X -n --format=State -P 2>/dev/null \
    | grep -cE "$RETRYABLE|CANCELLED"
}

# A deliberate human cancel must not be fought. Counted as not-clean so the driver stops, but never
# resubmitted.
n_cancelled() {
  sacct -j "$1" -X -n --format=State -P 2>/dev/null | grep -c 'CANCELLED'
}

# Task numbers worth retrying, comma separated.
failed_tasks() {
  sacct -j "$1" -X -n --format=JobID,State -P 2>/dev/null \
    | grep -E "$RETRYABLE" | sed 's/|.*//' \
    | awk -F_ 'NF>1 {print $2}' | grep -E '^[0-9]+$' | sort -n -u | paste -sd, -
}

oom_tasks() {
  sacct -j "$1" -X -n --format=State -P 2>/dev/null | grep -c 'OUT_OF_MEMORY'
}

# run_stage <script> <array-or-empty> [extra sbatch args...]
# Returns 0 only when every task of the stage succeeded.
run_stage() {
  local script=$1 array=$2; shift 2
  local attempt=0 id spec bad mem_args=("$@")
  spec=$array
  # The driver's --mem is a FLOOR, not a ceiling. An sbatch command-line option BEATS an in-script
  # #SBATCH directive, so the flat --mem=32G below was silently LOWERING every stage that asks for
  # more: 14_accept requests 64G, 15_accept_rna 96G, 22_accept_wgs 64G. Nothing had failed from it,
  # which is exactly why it went unnoticed -- the exome acceptance peaks at 2.4 GB measured
  # (11559137_2, MaxRSS 2,379,535K) so 32G was never close to tight. 15_accept_rna has never run, so
  # its 96G is the one figure here that is still a guess rather than a measurement, and guessing LOW
  # on it costs an OOM, a retry and up to twelve hours of wall clock to learn nothing.
  local script_mem given_mem i
  script_mem=$(sed -n 's/^#SBATCH.*--mem=\([0-9][0-9]*\)G.*/\1/p' "$J/$script" | head -1)
  given_mem=${mem_args[*]}; given_mem=${given_mem##*--mem=}; given_mem=${given_mem%%[!0-9]*}
  if [ -n "$script_mem" ] && [ -n "$given_mem" ] && [ "$script_mem" -gt "$given_mem" ]; then
    say "$script asks ${script_mem}G and the driver offered ${given_mem}G; honouring the script"
    for i in "${!mem_args[@]}"; do
      case ${mem_args[$i]} in --mem=*) mem_args[$i]="--mem=${script_mem}G" ;; esac
    done
  fi
  while :; do
    attempt=$((attempt + 1))
    local args=()
    [ -n "$spec" ] && args+=("--array=$spec")
    args+=("${mem_args[@]}")
    id=$(IGI_ENV="$IGI_ENV" "$J/submit.sh" "$script" "${args[@]}" 2>&1 | tail -1 | awk '{print $NF}')
    if ! [[ $id =~ ^[0-9]+$ ]]; then
      say "FATAL $script: submit failed -> $id"
      return 1
    fi
    say "$script attempt $attempt submitted as $id ${args[*]}"
    wait_job "$id"
    local nf nc
    nf=$(n_failed "$id"); nc=$(n_cancelled "$id")
    if [ "${nf:-0}" -eq 0 ]; then
      say "$script OK ($id)"
      return 0
    fi
    if [ "${nc:-0}" -gt 0 ]; then
      say "FATAL $script ($id): $nc task(s) CANCELLED -- stopping rather than fighting a human"
      return 1
    fi
    bad=$(failed_tasks "$id")
    if [ -z "$bad" ]; then
      say "FATAL $script ($id): $nf task(s) failed but no task id parsed; not guessing"
      return 1
    fi
    local n_oom; n_oom=$(oom_tasks "$id")
    say "$script: failed tasks [$bad], $n_oom OOM"
    if [ "$attempt" -gt "$RETRIES" ]; then
      say "FATAL $script: giving up after $attempt attempts; failed [$bad]"
      return 1
    fi
    # Escalate memory only when the kernel said memory. Retrying a logic error with more RAM is
    # cargo-culting; retrying an OOM without more RAM is pointless.
    if [ "$n_oom" -gt 0 ]; then
      local cur=${mem_args[*]}; cur=${cur##*--mem=}; cur=${cur%%[!0-9]*}
      local next=$(( (${cur:-48} * 3) / 2 ))
      mem_args=("--mem=${next}G")
      say "$script: OOM -> retrying failed tasks at ${next}G"
    fi
    spec=$bad
  done
}

say "=== driver starting; waiting on the stages already in flight ==="
for w in ${WAIT_FOR:-}; do
  say "waiting for $w"
  wait_job "$w"
  bad=$(failed_tasks "$w")
  say "  $w finished; failed tasks [${bad:-none}]"
  if [ -n "$bad" ]; then
    say "  NOTE $w had failures; the stage that consumes it may fail and will be retried"
  fi
done

# Dependency order, from jobs/README.md. Each stage only starts once the previous is genuinely clean.
# The exome merge is driven here too, rather than left on afterok behind the exome rebuild: if one of
# the 96 exome tasks fails, an afterok merge never runs at all and nothing retries it.
run_stage 11_merge.sbatch          ""      --mem=32G || exit 1
run_stage 50_merge_longread.sbatch ""      --mem=32G || exit 1
run_stage 72_pacbio_pbi.sbatch     ""      --mem=16G || exit 1
run_stage 90_subset_map.sbatch     ""      --mem=32G || exit 1
run_stage 90_subset_chr1to6.sbatch ""      --mem=32G || exit 1
run_stage 91_merge_chr1to6.sbatch  ""      --mem=32G || exit 1
say "=== DATA COMPLETE: reads, BAMs and the chr1to6 release are built ==="

# Verification. These do not gate the data, so a failure here is reported and the rest still runs.
run_stage 92_scan_alphabet.sbatch      "" --mem=4G  || say "WARN alphabet scan failed"
run_stage 93_align_qc.sbatch           "" --mem=24G || say "WARN align qc failed"
run_stage 94_consumer_preflight.sbatch "" --mem=48G || say "WARN consumer preflight failed"
run_stage 22_accept_wgs.sbatch         "" --mem=32G || say "WARN wgs acceptance failed"
run_stage 14_accept.sbatch             "" --mem=32G || say "WARN exome acceptance failed"
run_stage 15_accept_rna.sbatch         "" --mem=32G || say "WARN rna acceptance failed"
say "=== DRIVER DONE ==="
