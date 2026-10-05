#!/bin/bash
#SBATCH -J igi_merge_lr
#SBATCH -c 8 --mem=32G -t 24:00:00
#SBATCH --array=1-12
# Merge the per-chromosome arms into per-library deliverables. 12 tasks = 3 assays x 2 libraries x 2
# datasets, so each runs independently and a failure costs one library rather than all of them.
#
#   wgs      per-chromosome FASTQ pairs -> one R1/R2 per library, plus readmap, records and junctions
#   pacbio   per-chromosome hifi BAMs   -> one BAM per library
#   ont_wgs  per-chromosome FASTQ       -> one FASTQ per library
#
# Read names are globally unique per chromosome by construction (NAME_SPACE_STRIDE plus the record prefix),
# and the FASTQ merges PROVE it by sorting rather than assuming it -- an earlier build had 667,090 of
# 1,000,832 record ids duplicated while the integrity check of the day passed.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
# Rebuild when the output is MISSING OR OLDER THAN ITS NEWEST INPUT, not merely when it is missing.
# The guards here used to be `[ -s "$out" ] || cat ...`, which is right for resuming an interrupted merge
# and wrong for every rebuild: after the copy-number fix re-simulated all 96 WGS chromosomes, the merge
# would have found the previous merged library in place, skipped, and reported success -- leaving the
# deliverable identical to the one the rebuild existed to replace.
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

ASSAYS=(wgs pacbio ont_wgs); LB=(tumor normal)
i=$((SLURM_ARRAY_TASK_ID-1))
asy=${ASSAYS[$((i / 4))]}; j=$((i % 4)); l=${LB[$((j % 2))]}; d=${DATASETS[$((j / 2))]}
O=$IGI_WORK/release/$asy
M=$IGI_WORK/release/${asy}_merged/$d
mkdir -p "$M"
echo "=== $asy $d $l ==="

case "$asy" in
  wgs)
    for m in 1 2; do
      out=$M/${d}_${l}_wgs_R${m}.fastq.gz
      if stale "$out" $O/$d/${d}_chr*_${l}_wgs_R${m}.fastq.gz; then
        echo "  merging R${m} ($(ls $O/$d/${d}_chr*_${l}_wgs_R${m}.fastq.gz | wc -l) chromosomes)"
        cat $O/$d/${d}_chr*_${l}_wgs_R${m}.fastq.gz > "$out"
      else
        echo "  R${m} is newer than every input; keeping it"
      fi
    done
    n=$(zcat $M/${d}_${l}_wgs_R1.fastq.gz | awk 'NR%4==1' | wc -l)
    u=$(zcat $M/${d}_${l}_wgs_R1.fastq.gz | awk 'NR%4==1' | sort -T "$M" --parallel=8 -u | wc -l)
    echo "  reads $n, distinct names $u"
    [ "$n" -eq "$u" ] || { echo "FAIL: duplicate read names"; exit 1; }
    for what in readmap records; do
      out=$M/${d}_${l}_wgs_${what}.tsv; : > "$out"; first=1
      for f in $O/$d/${d}_chr*_${l}_wgs_${what}.tsv.gz; do
        [ -f "$f" ] || continue
        if [ $first -eq 1 ]; then zcat "$f" >> "$out"; first=0; else zcat "$f" | tail -n +2 >> "$out"; fi
      done
      gzip -f "$out"
    done
    jo=$M/${d}_${l}_wgs_junctions.tsv; : > "$jo"; first=1
    for f in $O/$d/${d}_chr*_${l}_wgs_junctions.tsv; do
      [ -f "$f" ] || continue
      if [ $first -eq 1 ]; then cat "$f" >> "$jo"; first=0; else tail -n +2 "$f" >> "$jo"; fi
    done
    echo "  junctions $(tail -n +2 "$jo" 2>/dev/null | wc -l)"
    ;;
  pacbio)
    out=$M/${d}_${l}_hifi.bam
    if stale "$out" $O/$d/${d}_chr*_${l}_hifi.bam; then
      rm -f "$out" "$out.pbi"
      ls $O/$d/${d}_chr*_${l}_hifi.bam > $M/${l}.bamlist
      echo "  merging $(wc -l < $M/${l}.bamlist) BAMs"
      singularity exec -B "$IGI_BIND" "$IGI_SAMTOOLS_IMG" \
        samtools merge -@ 8 -b $M/${l}.bamlist -o "$out"
    fi
    echo "  reads $(singularity exec -B "$IGI_BIND" "$IGI_SAMTOOLS_IMG" samtools view -c "$out")"
    ;;
  ont_wgs)
    out=$M/${d}_${l}_ont_wgs.fastq.gz
    n=$(ls $O/$d/${d}_chr*_${l}_ont_wgs.fastq.gz 2>/dev/null | wc -l)
    [ "$n" -eq 24 ] || { echo "  only $n/24 chromosomes present; not merging a partial library"; exit 0; }
    if stale "$out" $O/$d/${d}_chr*_${l}_ont_wgs.fastq.gz; then
      echo "  merging $n chromosomes"
      cat $O/$d/${d}_chr*_${l}_ont_wgs.fastq.gz > "$out"
    else
      echo "  merged library is newer than every input; keeping it"
    fi
    echo "  reads $(zcat "$out" | awk 'NR%4==1' | wc -l)"
    ;;
esac
echo "DONE $asy $d $l $(date)"
