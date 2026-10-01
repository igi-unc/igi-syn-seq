#!/bin/bash
#SBATCH -J igi_accept_wgs
#SBATCH -c 16 --mem=64G -t 24:00:00
#SBATCH --array=1-2
# Acceptance for the short-read WGS arm. The depth and allele-fraction checks need an aligned BAM, so this
# aligns the merged library first.
#
# --assay wgs changes three things, each because the assay genuinely differs rather than to make the check
# pass: depth is measured over non-N reference instead of capture intervals (chr21 is 18.6 % N, and an
# N-straddling window otherwise reads as a deletion); the GC curve is not applied, since it is capture
# efficiency and the builder uses a flat one; and every designed breakpoint on the chromosome counts as
# placeable, which is a stronger expectation than the exome's, not a weaker one.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
O=$IGI_WORK/release/wgs; M=$IGI_WORK/release/wgs_merged
A=$IGI_WORK/release/wgs_bams; R=$IGI_WORK/release/wgs_acceptance
d=${DATASETS[$((SLURM_ARRAY_TASK_ID-1))]}
mkdir -p "$M/$d" "$A/$d" "$R/$d"
for m in 1 2; do
  out=$M/$d/${d}_tumor_R${m}.fastq.gz
  [ -s "$out" ] || cat $O/$d/${d}_chr*_tumor_R${m}.fastq.gz > "$out"
done
B=$A/$d/${d}_tumor.bam
if [ ! -s "$B" ]; then
  singularity exec -B "$IGI_BIND" "$IGI_BWA_IMG" bash -c "
    bwa mem -t 16 -R '@RG\tID:${d}_tumor_wgs\tSM:${d}_tumor\tPL:ILLUMINA\tLB:${d}_tumor_wgs' \
      $IGI_FASTA_NOALT $M/$d/${d}_tumor_R1.fastq.gz $M/$d/${d}_tumor_R2.fastq.gz 2> $A/$d/bwa.err \
      | samtools sort -@4 -m 3G -o $B -" || { tail -20 "$A/$d/bwa.err"; exit 1; }
  singularity exec -B "$IGI_BIND" "$IGI_BWA_IMG" samtools index "$B"
fi
cd "$IGI_REPO/catalog"
for c in chr1 chr2 chr5 chr6 chr8 chr10 chr13 chr17 chr21 chrX; do
  python3 acceptance.py --dataset "$d" --chrom $c --catalog-dir output --design design.yaml \
    --assay wgs --bam "$B" --reference "$IGI_FASTA_NOALT" \
    --arms-bed "$IGI_ARMS_BED" --paths "$IGI_PATHS" \
    --junctions "$O/$d/${d}_${c}_tumor_junctions.tsv" \
    --work "$R/$d/$c" --report "$R/$d/${d}_${c}.json" 2>&1 | sed "s/^/[$d $c] /" || true
done
echo "DONE $d $(date)"
