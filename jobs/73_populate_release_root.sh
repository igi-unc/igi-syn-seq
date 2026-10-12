#!/bin/bash
# Populate a CANDIDATE release root with the arms a rebuild does not touch, as symlinks into an existing
# root, so the merge, subset and acceptance stages see a complete release without 6.9 TB being copied.
#
#   ./jobs/73_populate_release_root.sh --from $IGI_WORK/release --to $IGI_WORK/release-20261012 \
#     --arm wes --arm wgs --arm pacbio --arm 'rna:*_full_rna_*' \
#     --except 'IGI-SYN-SEQ-02_chr*_normal*'
#
# WHY A CANDIDATE ROOT AT ALL. 70_stage_inputs.sh stages the release into the RAFT workspace as symlinks,
# so the deliverables have no copy of their own. A builder writing into the live release root changes what
# the workspace serves the instant it writes -- there is no staging step left to hold anything back, and
# nothing to compare the new library against because the old one is gone. Building beside it restores that.
#
# --arm <name>[:<glob>]  link $FROM/<name>/<ds>/ into $TO/<name>/<ds>/, all of it or just <glob>.
# --except <glob>        a basename pattern the REBUILD will write. Nothing matching it is linked, and the
#                        script then proves nothing matching it exists in the destination.
#
# THE --except CHECK IS THE POINT OF THIS SCRIPT. A shell redirect follows a symlink: `> out.gz` where
# out.gz links to the live release overwrites the live file, from a job that believes it is writing
# somewhere else entirely. The patterns are given rather than inferred because only the caller knows which
# array ids they are about to submit, and a wrong guess here is indistinguishable from success until the
# live release is already damaged.
set -euo pipefail

FROM=""; TO=""; ARMS=(); EXCEPT=(); CLI_DS=(); DATASETS=(${DATASETS[@]+"${DATASETS[@]}"})
DRY=${DRY_RUN:-0}
while [ $# -gt 0 ]; do
  case "$1" in
    --from) FROM=$2; shift 2 ;;
    --to) TO=$2; shift 2 ;;
    --arm) ARMS+=("$2"); shift 2 ;;
    --except) EXCEPT+=("$2"); shift 2 ;;
    --datasets) IFS=, read -r -a CLI_DS <<<"$2"; shift 2 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done
[ -n "$FROM" ] && [ -n "$TO" ] || { echo "--from and --to are required" >&2; exit 2; }
[ -d "$FROM" ] || { echo "--from $FROM is not a directory" >&2; exit 2; }
[ "${#ARMS[@]}" -gt 0 ] || { echo "at least one --arm is required" >&2; exit 2; }
# Resolved, not compared as strings: $IGI_WORK/release and $IGI_WORK/./release are the same directory and
# linking a root into itself would replace real files with links to themselves.
[ "$(readlink -f "$FROM")" != "$(readlink -f "$TO")" ] || {
  echo "--from and --to resolve to the same directory" >&2; exit 2; }

# Written out rather than nested defaults: "${CLI_DS[@]:-${DATASETS[@]:-a b}}" collapses the fallback
# into ONE element, so the first run looked for a directory named "IGI-SYN-SEQ-01 IGI-SYN-SEQ-02" and
# reported every arm absent -- which reads exactly like a correct run against an empty source.
if [ "${#CLI_DS[@]}" -gt 0 ] 2>/dev/null; then
  DS=("${CLI_DS[@]}")
elif [ "${#DATASETS[@]}" -gt 0 ] 2>/dev/null; then
  DS=("${DATASETS[@]}")
else
  DS=(IGI-SYN-SEQ-01 IGI-SYN-SEQ-02)
fi

excluded() {  # excluded <basename>
  local n=$1 g
  for g in ${EXCEPT[@]+"${EXCEPT[@]}"}; do
    # shellcheck disable=SC2254 -- $g is a glob on purpose
    case "$n" in $g) return 0 ;; esac
  done
  return 1
}

n_link=0; n_skip=0
for spec in "${ARMS[@]}"; do
  arm=${spec%%:*}
  pat="*"; [ "$spec" != "$arm" ] && pat=${spec#*:}
  for d in "${DS[@]}"; do
    src=$FROM/$arm/$d
    [ -d "$src" ] || { echo "  $arm/$d: absent in --from, nothing to link"; continue; }
    dst=$TO/$arm/$d
    [ "$DRY" = 1 ] || mkdir -p "$dst"
    a=0; b=0
    for f in "$src"/$pat; do
      [ -e "$f" ] || continue
      n=$(basename "$f")
      if excluded "$n"; then b=$((b+1)); continue; fi
      # -s to an ABSOLUTE target, and -n so a pre-existing link is replaced rather than followed into
      # its target directory. No -f on a regular file: if the destination already holds a real file the
      # link is refused, because that file is output somebody built and this script did not write it.
      if [ -e "$dst/$n" ] && [ ! -L "$dst/$n" ]; then
        echo "  REFUSING $arm/$d/$n: the destination holds a real file, not a link" >&2
        exit 1
      fi
      [ "$DRY" = 1 ] || ln -sfn "$(readlink -f "$f")" "$dst/$n"
      a=$((a+1))
    done
    n_link=$((n_link+a)); n_skip=$((n_skip+b))
    printf '  %-12s %-16s linked %4d   held back for the rebuild %4d\n' "$arm" "$d" "$a" "$b"
  done
done

echo
echo "  linked $n_link   held back $n_skip"
[ "$DRY" = 1 ] && { echo "  (dry run; nothing written)"; exit 0; }

# PROVE the two properties this script exists for, rather than reporting that it tried.
#
#   1. nothing the rebuild will write is reachable through a link, so no job can write into $FROM
#   2. every link resolves, so a merge reading the candidate root cannot find a dangling input
bad=0
if [ "${#EXCEPT[@]}" -gt 0 ]; then
  while IFS= read -r l; do
    if excluded "$(basename "$l")"; then
      echo "LEAK: $l matches a rebuild pattern and is a link into $FROM" >&2; bad=$((bad+1))
    fi
  done < <(find "$TO" -type l)
fi
n_dangling=$(find "$TO" -xtype l | wc -l)
[ "$n_dangling" -eq 0 ] || { echo "$n_dangling dangling link(s) in $TO" >&2; bad=$((bad+1)); }
[ "$bad" -eq 0 ] || { echo "the candidate root is not safe to build into" >&2; exit 1; }
echo "  verified: no link matches a rebuild pattern, no link dangles"
