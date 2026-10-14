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

# Refuse to build into a path that leaves $IGI_RELEASE through a symlink.
#
# On 2026-10-14 the ds-02 rebuild was pointed at a candidate release root, $IGI_WORK/release-20261012,
# precisely so it would not touch the live staged release. It overwrote 361 live files anyway, two of
# them served to a consumer through ~/raft/inputs. The candidate root is not an independent tree: it was
# populated by 73_populate_release_root.sh, which LINKS the arms that are not being rebuilt back into
# $IGI_WORK/release so the candidate reads as a whole release. Those links do not care which direction
# you came from -- a builder opening $IGI_RELEASE/wes/<ds>/<file> for writing follows the link and lands
# in the live release.
#
# The check that was run beforehand asked whether anything in ~/raft/inputs pointed INTO the candidate
# (nothing did, 26,614 links). That is the wrong direction and it is the easy one to ask, because it is a
# question about the thing you are protecting rather than about the thing you are running. The only
# question that would have caught this is "is there a symlink in MY WRITE PATH", which is this function,
# and it costs one find per task.
#
# 73_populate_release_root.sh --except exists to leave the rebuilt arms unlinked. It depends on the
# operator naming the same arms they are about to rebuild; nothing checked that the two sets agreed, and
# on 2026-10-14 they did not. An assertion in the builder needs no such agreement: it is evaluated
# against the path actually about to be written.
assert_write_path() {
  local root d l t bad=""
  root=$(readlink -f "${IGI_RELEASE:?IGI_RELEASE must be set before asserting a write path}")
  for d in "$@"; do
    [ -e "$d" ] || continue
    while IFS= read -r l; do
      [ -n "$l" ] || continue
      t=$(readlink -f "$l" 2>/dev/null) || t="(broken link)"
      case "$t" in "$root"/*) ;; *) bad="$bad  $l -> $t"$'\n' ;; esac
    done < <(find "$d" -type l 2>/dev/null)
  done
  if [ -n "$bad" ]; then
    {
      echo "REFUSING TO BUILD: the write path contains symlink(s) that leave IGI_RELEASE=$root."
      echo "Writing through them would modify files outside this release root, which may be staged"
      echo "and delivered. Re-run 73_populate_release_root.sh with --except covering this arm, or"
      echo "replace the links with real files, then resubmit."
      printf '%s' "$bad" | head -10
      n=$(printf '%s' "$bad" | grep -c . || true); [ "${n:-0}" -gt 10 ] && echo "  ... and $((n - 10)) more"
    } >&2
    exit 1
  fi
}
