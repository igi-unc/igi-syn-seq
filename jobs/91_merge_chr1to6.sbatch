#!/bin/bash
#SBATCH -J igi_merge_sub
#SBATCH -c 8 --mem=32G -t 24:00:00
#SBATCH --array=1-16
# Merge the chr1-chr6 shards of the four per-chromosome arms into ONE library each, so the chr1to6
# release has the same shape as the full one: one file per sample, not six.
#
# 16 tasks = 4 assays (wes, wgs, ont_wgs, pacbio) x 2 libraries x 2 datasets.
#
# WHY MERGE AT ALL, when the shards already exist and a FASTQ row may name several files: because a BAM
# row may not. lens-v2.0.0-dev's preflight_resolve_inputs gathers every file whose name starts with
# File_Prefix and then, for BAM input, fails hard on more than one -- "prefix matches N BAMs; expected
# one". PacBio HiFi WGS is BAM-delivered, so its six chromosome shards cannot be addressed by any single
# prefix. And with the release's chromosome-first naming (<ds>_chr1_tumor_wgs_R1) there is no prefix that
# selects one assay's six chromosomes without also selecting the other assays', so even the FASTQ arms
# have no working prefix until they are merged. One merged library per sample is also what design section 2
# means by the subset mirroring the full release.
#
# Rebuild when the output is MISSING OR OLDER THAN ITS NEWEST INPUT. `[ -s "$out" ]` alone is right for
# resuming an interrupted merge and wrong for every rebuild; it shipped stale data twice in one release.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
stale() {
  local out=$1; shift
  [ -s "$out" ] || return 0
  local f
  for f in "$@"; do
    [ -f "$f" ] || continue
    [ "$f" -nt "$out" ] && return 0
  done
  return 1
}

ASSAYS=(wes wgs ont_wgs pacbio); LB=(tumor normal)
i=$((SLURM_ARRAY_TASK_ID-1))
asy=${ASSAYS[$((i / 4))]}; j=$((i % 4)); l=${LB[$((j % 2))]}; d=${DATASETS[$((j / 2))]}
O=$IGI_WORK/release/$asy
M=$IGI_WORK/release/${asy}_chr1to6/$d
mkdir -p "$M"
echo "=== $asy $d $l chr1to6 ==="

# Six shards, explicitly. A glob of chr[1-6] would also take chr1 through chr6 only, but writing the list
# out means a missing chromosome is a visible failure rather than a quietly shorter library.
CH=(chr1 chr2 chr3 chr4 chr5 chr6)
shards() {  # shards <suffix-template with %s for the chromosome>
  local tmpl=$1 c f out=()
  for c in "${CH[@]}"; do
    f=$(printf "$tmpl" "$c")
    [ -f "$f" ] || { echo "MISSING shard: $f" >&2; return 1; }
    out+=("$f")
  done
  printf '%s\n' "${out[@]}"
}

case "$asy" in
  wes|wgs)
    sfx=""; [ "$asy" = wgs ] && sfx="_wgs"
    for m in 1 2; do
      out=$M/${d}_chr1to6_${l}${sfx}_R${m}.fastq.gz
      mapfile -t IN < <(shards "$O/$d/${d}_%s_${l}${sfx}_R${m}.fastq.gz")
      if stale "$out" "${IN[@]}"; then
        echo "  merging R${m} from ${#IN[@]} chromosomes"
        cat "${IN[@]}" > "$out"
      else
        echo "  R${m} is newer than every shard; keeping it"
      fi
    done
    n=$(zcat $M/${d}_chr1to6_${l}${sfx}_R1.fastq.gz | awk 'NR%4==1' | wc -l)
    u=$(zcat $M/${d}_chr1to6_${l}${sfx}_R1.fastq.gz | awk 'NR%4==1' | sort -T "$M" --parallel=8 $(sort_buffer_arg) -u | wc -l)
    echo "  reads $n, distinct names $u"
    [ "$n" -eq "$u" ] || { echo "FAIL: duplicate read names after merge"; exit 1; }
    for what in readmap records; do
      out=$M/${d}_chr1to6_${l}${sfx}_${what}.tsv; : > "$out"; first=1
      for f in $(shards "$O/$d/${d}_%s_${l}${sfx}_${what}.tsv.gz"); do
        if [ $first -eq 1 ]; then zcat "$f" >> "$out"; first=0; else zcat "$f" | tail -n +2 >> "$out"; fi
      done
      gzip -f "$out"
    done
    jo=$M/${d}_chr1to6_${l}${sfx}_junctions.tsv; : > "$jo"; first=1
    for c in "${CH[@]}"; do
      f=$O/$d/${d}_${c}_${l}${sfx}_junctions.tsv
      [ -f "$f" ] || continue
      if [ $first -eq 1 ]; then cat "$f" >> "$jo"; first=0; else tail -n +2 "$f" >> "$jo"; fi
    done
    echo "  junctions $(tail -n +2 "$jo" 2>/dev/null | wc -l)"
    ;;
  ont_wgs)
    out=$M/${d}_chr1to6_${l}_ont_wgs.fastq.gz
    mapfile -t IN < <(shards "$O/$d/${d}_%s_${l}_ont_wgs.fastq.gz")
    if stale "$out" "${IN[@]}"; then
      echo "  merging ${#IN[@]} chromosomes"
      cat "${IN[@]}" > "$out"
    else
      echo "  merged library is newer than every shard; keeping it"
    fi
    echo "  reads $(zcat "$out" | awk 'NR%4==1' | wc -l)"
    ;;
  pacbio)
    out=$M/${d}_chr1to6_${l}_hifi.bam
    mapfile -t IN < <(shards "$O/$d/${d}_%s_${l}_hifi.bam")
    if stale "$out" "${IN[@]}"; then
      rm -f "$out" "$out.pbi"
      printf '%s\n' "${IN[@]}" > $M/${l}.chr1to6.bamlist
      echo "  merging ${#IN[@]} BAMs"
      singularity exec -B "$IGI_BIND" "$IGI_SAMTOOLS_IMG" \
        samtools merge -@ 8 -b $M/${l}.chr1to6.bamlist -o "$out"
    fi
    echo "  reads $(singularity exec -B "$IGI_BIND" "$IGI_SAMTOOLS_IMG" samtools view -c "$out")"
    # samtools merge writes no .pbi, and the design calls for one beside every PacBio BAM.
    PBI=$(python3 -c "import yaml;print(yaml.safe_load(open('$IGI_PATHS')).get('pbindex_cmd',''))")
    if [ -n "$PBI" ] && stale "$out.pbi" "$out"; then
      eval "${PBI/\{args\}/$out}" && echo "  pbi $(stat -c%s "$out.pbi") bytes"
    fi
    ;;
esac
echo "DONE $asy $d $l chr1to6 $(date)"
