#!/bin/bash
#SBATCH -J igi_align_wgs_t
#SBATCH -c 16 --mem=64G -t 24:00:00
#SBATCH --array=1-2
# Align the merged TUMOUR short-read WGS library.
#
# This was inside 22_accept_wgs until 2026-10-19: the acceptance stage built the very BAM it then graded.
# Two costs, both paid. First, wall clock -- the tumour alignment is ~15 h of the stage's 21 h, so a
# verification stage held a 40 h limit and could not overlap the ~9 h of merges and subsets that ran
# ahead of it; splitting it out lets it run beside 23_align_wgs_normal instead and leaves acceptance at
# its check time. Second, and worse, a stage that produces its own input has no provenance to check: it
# can never report "the BAM I graded is older than the library it came from", because it is the thing
# that would have rebuilt it.
#
# The thread counts are DELIBERATELY unchanged from the version inside 22_accept_wgs (-t 16, sort -@4
# -m 3G against -c 16). bwa mem batches reads per thread and mem_pestat derives its insert-size model
# per batch, so a different -t can shift a handful of pairing decisions. There is no reason to accept
# that churn in a BAM that an acceptance stage then grades to 4 sigma.
#
# Mirrors 23_align_wgs_normal: same read-group convention, same mtime staleness guard, same refusal to
# merge a library that the merge stage owns.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
source "$(dirname "$IGI_ENV")/lib.sh"
O=$IGI_RELEASE/wgs; M=$IGI_RELEASE/wgs_merged; A=$IGI_RELEASE/wgs_bams
d=${DATASETS[$((SLURM_ARRAY_TASK_ID-1))]}
mkdir -p "$M/$d" "$A/$d"
assert_write_path "$A/$d" "$M/$d"
for e in amb ann bwt pac sa; do
  [ -f "$IGI_FASTA_NOALT.$e" ] || { echo "missing bwa index $IGI_FASTA_NOALT.$e"; exit 1; }
done
# Use 50_merge_longread.sbatch's output rather than merging again. Both scripts wrote the same libraries
# into this directory under different names and duplicated 234 GB: this one as ${d}_tumor_R1.fastq.gz and
# the merge job as ${d}_tumor_wgs_R1.fastq.gz. Fall back to merging only if the merge job has not run.
R1=$M/$d/${d}_tumor_wgs_R1.fastq.gz
R2=$M/$d/${d}_tumor_wgs_R2.fastq.gz
newest=$(ls -t $O/$d/${d}_chr*_tumor_wgs_R1.fastq.gz 2>/dev/null | head -1)
if [ ! -s "$R1" ] || { [ -n "$newest" ] && [ "$newest" -nt "$R1" ]; }; then
  echo "  merge job output absent or older than its chromosomes; merging here instead"
  # interleaved, not cat: see 11_merge.sbatch for why block order broke the consumer preflight
  python3 "$IGI_REPO/catalog/interleave_fastq.py" \
    --in-r1 "$(ls $O/$d/${d}_chr*_tumor_wgs_R1.fastq.gz | paste -sd,)" \
    --in-r2 "$(ls $O/$d/${d}_chr*_tumor_wgs_R2.fastq.gz | paste -sd,)" \
    --out-r1 "$R1" --out-r2 "$R2" --threads 8
fi
# Rebuild when the BAM is MISSING OR OLDER THAN ITS LIBRARY, not merely when it is missing: a rebuild
# that finds a stale BAM in place and skips is how a release gets graded on data it no longer contains.
B=$A/$d/${d}_tumor.bam
if [ ! -s "$B" ] || [ "$R1" -nt "$B" ]; then
  rm -f "$B" "$B.bai"
  singularity exec -B "$IGI_BIND" "$IGI_BWA_IMG" bash -c "
    bwa mem -t 16 -R '@RG\tID:${d}_tumor_wgs\tSM:${d}_tumor\tPL:ILLUMINA\tLB:${d}_tumor_wgs' \
      $IGI_FASTA_NOALT $R1 $R2 2> $A/$d/bwa.err \
      | samtools sort -@4 -m 3G -o $B -" || { tail -20 "$A/$d/bwa.err"; exit 1; }
  singularity exec -B "$IGI_BIND" "$IGI_BWA_IMG" samtools index "$B"
else
  echo "  $B is newer than its library; keeping it"
fi
singularity exec -B "$IGI_BIND" "$IGI_BWA_IMG" samtools flagstat "$B" | head -6
echo "DONE $d tumor wgs $(date)"
