#!/bin/bash
#SBATCH -J igi_artprof
#SBATCH -c 16 --mem=32G -t 8:00:00
# Fit the Illumina per-cycle quality profile from real reads, producing
# catalog/resources/ipisrc044_novaseqx_R{1,2}.profile.txt. ART's stock HS25 profile is HiSeq 2500
# chemistry from 2013 with continuous quality scores; current instruments bin qualities into a handful of
# levels, so stock HS25 produces quality strings unlike anything a caller sees today. The fitted profile
# carries the four levels these instruments use, Q3/Q13/Q25/Q41.
#
# The design asks for a profile from IPISRC044 WES, but this patient has no 150 bp exome: the BostonGene
# exomes are 101 bp (normal) and 140 bp (tumour) while the design specifies PE150. A per-cycle quality
# profile is a property of the sequencer run rather than the capture, and capture effects are modelled
# separately by the GC curve, so the profile is fitted from the 150 bp Personalis blood-normal WGS
# (NovaSeq X, instrument LH00499) and PE150 is kept as designed.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
W=$IGI_WORK/artprof; mkdir -p "$W/reads" "$W/out"
SRC=$IGI_REF/inputs/fastqs/IPISRC044
N=4000000   # read pairs; per-cycle quality distributions converge well below this
set +o pipefail
for m in 1 2; do
  out=$W/reads/ipisrc044_${m}.fq.gz
  [ -s "$out" ] || zcat $SRC/DNA_PSNLDx240603252_AACAAGGC-ATCTCGCT_L001_R${m}_001.fastq.gz     | head -n $((N*4)) | gzip -1 > "$out"
done
set -o pipefail
for m in 1 2; do
  n=$(zcat "$W/reads/ipisrc044_${m}.fq.gz" | wc -l)
  [ $((n % 4)) -eq 0 ] || { echo "ipisrc044_${m}.fq.gz: $n lines, not whole records"; exit 1; }
  echo "ipisrc044_${m}.fq.gz: $((n/4)) reads"
done
cd "$W/out"
singularity exec -B "$IGI_BIND" "$IGI_ART_IMG" art_profiler_illumina ipisrc044_novaseqx "$W/reads" fq.gz 16
cp "$W/out/ipisrc044_novaseqx_R1.txt" "$IGI_REPO/catalog/resources/ipisrc044_novaseqx_R1.profile.txt"
cp "$W/out/ipisrc044_novaseqx_R2.txt" "$IGI_REPO/catalog/resources/ipisrc044_novaseqx_R2.profile.txt"
echo "DONE $(date)"
