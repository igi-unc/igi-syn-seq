#!/bin/bash
#SBATCH -J igi_align
#SBATCH -c 24 --mem=64G -t 48:00:00
#SBATCH --array=1-4
# Align the merged exomes so acceptance runs on the release itself rather than a slice (D7 condition 1).
# An output that merely EXISTS is not an output that is current. Guarding with `[ ! -s "$x" ]` skipped
# work whenever a previous run had left something behind, and that silently graded or shipped stale data
# twice in one release: the long-read merge would have kept pre-rebuild libraries, and the WGS acceptance
# scored a 2026-10-05 alignment of a 2026-10-06 library in 35 minutes instead of 21 hours. Compare mtimes.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
M=$IGI_RELEASE/merged; A=$IGI_RELEASE/bams
LB=(tumor normal)
i=$((SLURM_ARRAY_TASK_ID-1)); d=${DATASETS[$((i/2))]}; l=${LB[$((i%2))]}
mkdir -p "$A/$d"
source "$(dirname "${IGI_ENV}")/lib.sh"
assert_write_path "$A/$d"
for e in amb ann bwt pac sa; do
  [ -f "$IGI_FASTA_NOALT.$e" ] || { echo "missing bwa index $IGI_FASTA_NOALT.$e"; exit 1; }
done
# $R1 and $R2 are NAMED here. The staleness test below used $R1 while nothing in this script ever set
# it, so under `set -u` the stage worked only when the BAM was absent -- `[ ! -s "$B" ]` short-circuits
# the `||` -- and died with "R1: unbound variable" on exactly the rebuild the test exists to serve. The
# comment above about mtimes was therefore describing a comparison that could not run.
R1=$M/$d/${d}_${l}_R1.fastq.gz
R2=$M/$d/${d}_${l}_R2.fastq.gz
for f in "$R1" "$R2"; do
  [ -s "$f" ] || { echo "missing merged library $f; run 11_merge.sbatch first"; exit 1; }
done
B=$A/$d/${d}_${l}.bam
if [ ! -s "$B" ] || [ "$R1" -nt "$B" ]; then
  rm -f "$B" "$B.bai"
  # bwa's stderr is kept: suppressing it once turned a missing index into an unreadable
  # "failed to read header".
  singularity exec -B "$IGI_BIND" "$IGI_BWA_IMG" bash -c "
    bwa mem -t 24 -R '@RG\tID:${d}_${l}\tSM:${d}_${l}\tPL:ILLUMINA\tLB:${d}_${l}' \
      $IGI_FASTA_NOALT $R1 $R2 2> $A/$d/${l}.bwa.err \
      | samtools sort -@6 -m 3G -o $B -" || { tail -20 "$A/$d/${l}.bwa.err"; exit 1; }
  singularity exec -B "$IGI_BIND" "$IGI_BWA_IMG" samtools index "$B"
fi
singularity exec -B "$IGI_BIND" "$IGI_BWA_IMG" samtools flagstat "$B" | head -6
echo "DONE $d $l $(date)"
