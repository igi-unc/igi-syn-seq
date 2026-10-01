#!/bin/bash
#SBATCH -J igi_gcbias
#SBATCH -c 16 --mem=48G -t 12:00:00
# Fit capture efficiency against interval GC from a real exome, producing
# catalog/resources/ipisrc044_wes_gc_bias.json. Giving every capture interval the same depth is not a
# small simplification: measured coverage runs from 0.40 of baseline at 30-35 % GC to 1.25 at 60-65 %,
# and 28 % of captured bases sit in the 30-40 % band. A flat-depth exome hands a caller far more evidence
# than it would really have at a quarter of its targets.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
W=$IGI_WORK/gcbias; mkdir -p "$W"
SRC=$IGI_REF/inputs/fastqs/IPISRC044
N=12000000   # read pairs; enough for a stable per-interval mean over the exome
# `zcat | head` kills zcat with SIGPIPE once head has what it needs, which `set -o pipefail` turns into a
# failed pipeline. Subsampling runs with pipefail off; the record-count check afterwards is the real test.
set +o pipefail
for m in 1 2; do
  out=$W/wes_${m}.fq.gz
  [ -s "$out" ] || zcat $SRC/BG003082_WES-normal_${m}.fastq.gz | head -n $((N*4)) | gzip -1 > "$out"
done
set -o pipefail
for m in 1 2; do
  n=$(zcat "$W/wes_${m}.fq.gz" | wc -l)
  [ $((n % 4)) -eq 0 ] || { echo "wes_${m}.fq.gz: $n lines, not a whole number of records"; exit 1; }
done
for e in amb ann bwt pac sa; do
  [ -f "$IGI_FASTA_NOALT.$e" ] || { echo "missing bwa index $IGI_FASTA_NOALT.$e"; exit 1; }
done
if [ ! -s "$W/wes.bam" ]; then
  singularity exec -B "$IGI_BIND" "$IGI_BWA_IMG" bash -c "
    bwa mem -t 16 $IGI_FASTA_NOALT $W/wes_1.fq.gz $W/wes_2.fq.gz 2> $W/bwa.err       | samtools sort -@4 -m 2G -o $W/wes.bam -" || { tail -20 "$W/bwa.err"; exit 1; }
  singularity exec -B "$IGI_BIND" "$IGI_BWA_IMG" samtools index "$W/wes.bam"
fi
singularity exec -B "$IGI_BIND" "$IGI_BWA_IMG" samtools flagstat "$W/wes.bam" | head -5
# Summed base depth per capture interval, then the fit.
[ -s "$W/bedcov.txt" ] || singularity exec -B "$IGI_BIND" "$IGI_SAMTOOLS_IMG"   samtools bedcov "$IGI_EXOME_BED" "$W/wes.bam" > "$W/bedcov.txt"
python3 "$IGI_REPO/jobs/measure/fit_gc_bias.py" --bedcov "$W/bedcov.txt"   --reference "$IGI_FASTA_NOALT"   --source "IPISRC044 BostonGene WES normal, 12M read pairs, bwa mem to GRCh38 noalt"   --out "$IGI_REPO/catalog/resources/ipisrc044_wes_gc_bias.json"
echo "DONE $(date)"
