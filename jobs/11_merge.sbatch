#!/bin/bash
#SBATCH -J igi_merge
#SBATCH -c 8 --mem=32G -t 12:00:00
#SBATCH --array=1-4
# Merge the per-chromosome exome pieces into one library per dataset, and the truth maps alongside them.
# Read names are globally unique per chromosome run by construction (NAME_SPACE_STRIDE plus the record
# prefix); the check at the end PROVES it rather than assuming it. An earlier build had 667,090 of
# 1,000,832 record ids duplicated because ids restarted at 1 on each chromosome, and the integrity check
# of the day passed anyway -- it verified that names resolve, not that they are unique.
# An output that merely EXISTS is not an output that is current. Guarding with `[ ! -s "$x" ]` skipped
# work whenever a previous run had left something behind, and that silently graded or shipped stale data
# twice in one release: the long-read merge would have kept pre-rebuild libraries, and the WGS acceptance
# scored a 2026-10-05 alignment of a 2026-10-06 library in 35 minutes instead of 21 hours. Compare mtimes.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
# NOT $BASH_SOURCE: sbatch copies the script to /var/spool/slurmd/job<N>/slurm_script, so
# dirname $BASH_SOURCE is the spool directory and lib.sh is not there. This failed all four tasks
# of the exome merge with "lib.sh: No such file or directory" -- introduced by the sort-buffer fix
# (0b89bd5) into three merge stages at once, none of which ran again until now. IGI_ENV is an
# absolute path by contract (submit.sh refuses without it), so its directory is the jobs/ dir.
source "$(dirname "$IGI_ENV")/lib.sh"
O=$IGI_WORK/release/wes; M=$IGI_WORK/release/merged
LB=(tumor normal)
i=$((SLURM_ARRAY_TASK_ID-1)); d=${DATASETS[$((i/2))]}; l=${LB[$((i%2))]}
mkdir -p "$M/$d"
# INTERLEAVED, not concatenated. `cat` of the per-chromosome pieces made every delivered library
# chromosome-block ordered -- ds-01 tumour WES began with chr10 and its first 3.4 M pairs were chr10 and
# chr11 only. Reads within a chromosome are already shuffled so no aligner and no acceptance result
# depends on this, but the consumer preflight takes "the first 5 %" of the library, and on the exome that
# slice contained no chr6 reads at all: the razers3 and OptiType preflight, whose whole purpose is the
# MHC, tested a sample with no MHC in it. Anyone running `zcat | head` saw the same distortion.
#
# Both mates are advanced from the same chosen chromosome in one pass, so R1 record n stays mated to R2
# record n. Cutting the mates independently is the §20.5 failure and would break every paired consumer.
out1=$M/$d/${d}_${l}_R1.fastq.gz
out2=$M/$d/${d}_${l}_R2.fastq.gz
newest=$(ls -t $O/$d/${d}_chr*_${l}_R1.fastq.gz 2>/dev/null | head -1)
if [ ! -s "$out1" ] || [ ! -s "$out2" ] || { [ -n "$newest" ] && [ "$newest" -nt "$out1" ]; }; then
  in1=$(ls $O/$d/${d}_chr*_${l}_R1.fastq.gz | paste -sd,)
  in2=$(ls $O/$d/${d}_chr*_${l}_R2.fastq.gz | paste -sd,)
  python3 "$IGI_REPO/catalog/interleave_fastq.py" \
    --in-r1 "$in1" --in-r2 "$in2" --out-r1 "$out1" --out-r2 "$out2" --threads 8
fi
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
u=$(zcat $M/$d/${d}_${l}_R1.fastq.gz | awk 'NR%4==1' | sort -T "$M" --parallel=8 $(sort_buffer_arg) -u | wc -l)
echo "reads $n, distinct names $u"
[ "$n" -eq "$u" ] || { echo "FAIL: duplicate read names after merge"; exit 1; }
echo "records $(zcat $M/$d/${d}_${l}_records.tsv.gz | tail -n +2 | wc -l), junctions $(tail -n +2 $M/$d/${d}_${l}_junctions.tsv | wc -l)"
echo "DONE $d $l $(date)"
