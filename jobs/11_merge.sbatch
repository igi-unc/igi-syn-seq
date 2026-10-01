#!/bin/bash
#SBATCH -J igi_merge
#SBATCH -c 8 --mem=32G -t 12:00:00
#SBATCH --array=1-4
# Merge the per-chromosome exome pieces into one library per dataset, and the truth maps alongside them.
# Read names are globally unique per chromosome run by construction (NAME_SPACE_STRIDE plus the record
# prefix); the check at the end PROVES it rather than assuming it. An earlier build had 667,090 of
# 1,000,832 record ids duplicated because ids restarted at 1 on each chromosome, and the integrity check
# of the day passed anyway -- it verified that names resolve, not that they are unique.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
O=$IGI_WORK/release/wes; M=$IGI_WORK/release/merged
LB=(tumor normal)
i=$((SLURM_ARRAY_TASK_ID-1)); d=${DATASETS[$((i/2))]}; l=${LB[$((i%2))]}
mkdir -p "$M/$d"
for m in 1 2; do
  out=$M/$d/${d}_${l}_R${m}.fastq.gz
  [ -s "$out" ] || cat $O/$d/${d}_chr*_${l}_R${m}.fastq.gz > "$out"
done
# Truth maps: first file whole, the rest without their header. `tail -n +2` reads to EOF so nothing gets
# SIGPIPE; `zcat | head` under `set -o pipefail` does, which killed this job the first time.
for what in readmap records; do
  out=$M/$d/${d}_${l}_${what}.tsv
  : > "$out"; first=1
  for f in $O/$d/${d}_chr*_${l}_${what}.tsv.gz; do
    if [ $first -eq 1 ]; then zcat "$f" >> "$out"; first=0; else zcat "$f" | tail -n +2 >> "$out"; fi
  done
  gzip -f "$out"
done
jout=$M/$d/${d}_${l}_junctions.tsv
: > "$jout"; first=1
for f in $O/$d/${d}_chr*_${l}_junctions.tsv; do
  [ -f "$f" ] || continue
  if [ $first -eq 1 ]; then cat "$f" >> "$jout"; first=0; else tail -n +2 "$f" >> "$jout"; fi
done
echo "--- $d $l ---"
n=$(zcat $M/$d/${d}_${l}_R1.fastq.gz | awk 'NR%4==1' | wc -l)
u=$(zcat $M/$d/${d}_${l}_R1.fastq.gz | awk 'NR%4==1' | sort -T "$M" --parallel=8 -u | wc -l)
echo "reads $n, distinct names $u"
[ "$n" -eq "$u" ] || { echo "FAIL: duplicate read names after merge"; exit 1; }
echo "records $(zcat $M/$d/${d}_${l}_records.tsv.gz | tail -n +2 | wc -l), junctions $(tail -n +2 $M/$d/${d}_${l}_junctions.tsv | wc -l)"
echo "DONE $d $l $(date)"
