#!/bin/bash
#SBATCH -J igi_accept
#SBATCH -c 8 --mem=64G -t 24:00:00
#SBATCH --array=1-2
# Acceptance on the full release (D7 condition 1), per chromosome over the merged BAM so the copy-number
# and allele-fraction checks see every designed event. The chromosome list is the set carrying the
# designed CNAs and LOH events plus two neutral controls.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
M=$IGI_WORK/release/merged; A=$IGI_WORK/release/bams; O=$IGI_WORK/release/acceptance
d=${DATASETS[$((SLURM_ARRAY_TASK_ID-1))]}
mkdir -p "$O/$d"
cd "$IGI_REPO/catalog"
for c in chr1 chr2 chr5 chr6 chr8 chr10 chr13 chr17 chrX; do
  python3 acceptance.py --dataset "$d" --chrom $c --catalog-dir output --design design.yaml     --bam "$A/$d/${d}_tumor.bam" --reference "$IGI_FASTA_NOALT"     --arms-bed "$IGI_ARMS_BED" --capture-bed "$IGI_EXOME_BED"     --junctions "$M/$d/${d}_tumor_junctions.tsv" --paths "$IGI_PATHS" --seed 5000     --work "$O/$d/$c" --report "$O/$d/${d}_${c}.json" 2>&1 | sed "s/^/[$d $c] /" || true
done
# Read-level checks once over the whole merged library rather than per chromosome.
python3 acceptance.py --dataset "$d" --catalog-dir output --design design.yaml   --r1 "$M/$d/${d}_tumor_R1.fastq.gz"   --read-map "$M/$d/${d}_tumor_readmap.tsv.gz" --record-map "$M/$d/${d}_tumor_records.tsv.gz"   --work "$O/$d" --report "$O/$d/${d}_reads.json" 2>&1 | sed "s/^/[$d reads] /" || true
echo "DONE $d $(date)"
