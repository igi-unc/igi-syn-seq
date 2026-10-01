#!/bin/bash
#SBATCH -J igi_accept
#SBATCH -c 8 --mem=64G -t 24:00:00
#SBATCH --array=1-2
# Acceptance on the full release (D7 condition 1), per chromosome over the merged BAM so the copy-number
# and allele-fraction checks see every designed event. The chromosome list is the set carrying the designed
# CNAs and LOH events plus neutral controls.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
M=$IGI_WORK/release/merged; A=$IGI_WORK/release/bams; O=$IGI_WORK/release/acceptance
W_WES=$IGI_WORK/release/wes
d=${DATASETS[$((SLURM_ARRAY_TASK_ID-1))]}
mkdir -p "$O/$d"
cd "$IGI_REPO/catalog"

# The seed must be the one THIS chromosome was built with, because the off-target windows are derived from
# it. 10_wes.sbatch uses 5000 + array task id, and for the tumour library that is
# 5001 + chromosome index + 48 * dataset index. Passing a single --seed 5000 for every chromosome, which is
# what the first release check did, made the junction check test membership in windows the library was never
# built from: on dataset 01 chr13 the expected set {SV-0086, SV-0091, SV-0096} and the placed set {SV-0036}
# were completely disjoint, and the check passed anyway because it compared with >=.
#
# --build-meta is preferred where the builder recorded its own seed; it removes this arithmetic entirely.
k=0; [ "$d" = "${DATASETS[1]}" ] && k=1
for c in chr1 chr2 chr5 chr6 chr8 chr10 chr13 chr17 chr21 chrX; do
  ci=-1
  for n in "${!CHROMS_ALL[@]}"; do [ "${CHROMS_ALL[$n]}" = "$c" ] && ci=$n; done
  [ $ci -ge 0 ] || { echo "$c is not in CHROMS_ALL"; exit 1; }
  seed=$((5001 + ci + 48 * k))
  meta=$W_WES/$d/${d}_${c}_tumor.json
  mopt=""
  if [ -s "$meta" ] && grep -q '"seed"' "$meta"; then mopt="--build-meta $meta"; fi
  python3 acceptance.py --dataset "$d" --chrom "$c" --catalog-dir output --design design.yaml \
    --bam "$A/$d/${d}_tumor.bam" --reference "$IGI_FASTA_NOALT" \
    --arms-bed "$IGI_ARMS_BED" --capture-bed "$IGI_EXOME_BED" \
    --junctions "$M/$d/${d}_tumor_junctions.tsv" --paths "$IGI_PATHS" \
    --seed "$seed" $mopt \
    --work "$O/$d/$c" --report "$O/$d/${d}_${c}.json" 2>&1 | sed "s/^/[$d $c seed=$seed] /" || true
done

# Read-level checks once over the whole merged library rather than per chromosome.
python3 acceptance.py --dataset "$d" --catalog-dir output --design design.yaml \
  --r1 "$M/$d/${d}_tumor_R1.fastq.gz" \
  --read-map "$M/$d/${d}_tumor_readmap.tsv.gz" --record-map "$M/$d/${d}_tumor_records.tsv.gz" \
  --work "$O/$d" --report "$O/$d/${d}_reads.json" 2>&1 | sed "s/^/[$d reads] /" || true
echo "DONE $d $(date)"
