#!/bin/bash
#SBATCH -J igi_accept_rna
#SBATCH -c 8 --mem=96G -t 12:00:00
#SBATCH --array=1-2
# D7 conditions 1 and 2 for the RNA arm, and the load-bearing measurement for the expression classes: no
# chr6 slice could exercise ERVs, splice isoforms or CTAs, so this is the first point at which all six
# designed classes reach an RNA library, and the first test of the 0.01 expression floor.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
O=$IGI_WORK/release/rna; M=$IGI_WORK/release/merged; A=$IGI_WORK/release/acceptance
d=${DATASETS[$((SLURM_ARRAY_TASK_ID-1))]}
CH=$(IFS=,; echo "${CHROMS_ALL[*]}")
mkdir -p "$A/$d"
cd "$IGI_REPO/catalog"
python3 acceptance.py --dataset "$d" --chrom "$CH" --catalog-dir output --design design.yaml \
  --rna-manifest "$O/$d/${d}_full_rna_transcripts.tsv" \
  --r1 "$O/$d/${d}_full_rna_R1.fastq.gz" \
  --work "$A/$d/rna" --report "$A/$d/${d}_rna_release.json"
python3 class_presence.py --dataset "$d" --catalog-dir output \
  --rna-manifest "$O/$d/${d}_full_rna_transcripts.tsv" \
  --records "$M/$d/${d}_tumor_records.tsv.gz" \
  --junctions "$M/$d/${d}_tumor_junctions.tsv" \
  --out "$A/$d/${d}_class_presence.json" > "$A/$d/${d}_class_presence.txt"
echo "DONE $d $(date)"
