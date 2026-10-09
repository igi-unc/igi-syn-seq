# Shared, site-independent shell helpers for the job scripts. Sourced next to $IGI_ENV, which carries the
# site paths; nothing here may reference a path or a cluster.

# GNU sort sizes its default buffer from the machine's PHYSICAL memory and never reads the cgroup limit a
# SLURM step runs under. A bare `sort` in a 16 GB allocation on a 503 GB node is therefore OOM-killed, and
# a large --mem request does not bound the buffer -- it only has to exceed whatever sort picks, which is a
# property of the node. igi_rna task 1 died this way on 2026-10-08. Derive the buffer from the allocation
# so a step's footprint belongs to the job.
#
# Keep in step with sort_buffer_mb() in catalog/igi_catalog/readnames.py, which does the same for the
# Python caller.
sort_buffer_mb() {
  local frac_num=1 frac_den=4 lo=256 hi=8192 total=""
  if [[ "${SLURM_MEM_PER_NODE:-}" =~ ^[0-9]+$ ]]; then
    total=$SLURM_MEM_PER_NODE
  elif [[ "${SLURM_MEM_PER_CPU:-}" =~ ^[0-9]+$ ]] && [[ "${SLURM_CPUS_ON_NODE:-}" =~ ^[0-9]+$ ]]; then
    total=$((SLURM_MEM_PER_CPU * SLURM_CPUS_ON_NODE))
  fi
  if [ -z "$total" ] || [ "$total" -le 0 ]; then echo "$lo"; return; fi
  local buf=$((total * frac_num / frac_den))
  [ "$buf" -lt "$lo" ] && buf=$lo
  [ "$buf" -gt "$hi" ] && buf=$hi
  echo "$buf"
}

# Printed as a single `-S` argument, e.g. `sort $(sort_buffer_arg) -u ...`.
sort_buffer_arg() { echo "-S $(sort_buffer_mb)M"; }
