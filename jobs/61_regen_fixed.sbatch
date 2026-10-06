#!/bin/bash
#SBATCH -J igi_regen
#SBATCH -c 48 --mem=48G -t 48:00:00
#SBATCH --array=1-6
# Regenerate the three arms whose defects are in the reads themselves, so no in-place repair can fix them.
#
#   tasks 1-2  ONT bulk RNA. Read length was 910/461, fitted from IPISRC044's 10x 5' SINGLE-CELL cDNA --
#              a 5'-biased, deliberately truncated library. The design calls this sample "cDNA, full
#              length", and at 910 bp most transcripts are not covered end to end, so the 30 designed
#              splice events had nothing to be detected by. Now 2091/843, measured on 400,000 reads of the
#              real HG002 Kinnex FLNC BAM, which is a genuine full-length cDNA library.
#   tasks 3-4  10x GEX. Read names were "@<dataset>:<n>", which no instrument produces; now the same
#              Illumina scheme WES, WGS and bulk RNA use. R1 qualities were a flat two-level string
#              (663 distinct in 20,000 reads); now drawn per cycle from the real four-level binned
#              IPISRC044 R1, which gives 2,947 against the real library's 2,819.
#   tasks 5-6  10x TCR. The same two naming and quality fixes, plus the junctions: CDR3s were random
#              peptides with internal cysteines and no relationship to the V and J genes the truth table
#              named. They are now recombined from the real GENCODE germline anchors, and the clonotype
#              table carries V/J/C, the nucleotide junction and the flagpost antigen link the design asks
#              for.
#
# Memory is 48G, not the 320G the ONT builder used to ask for: its measured peak was 29G, and the
# over-request is what left it queued behind the QOS cap instead of running.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
t=$SLURM_ARRAY_TASK_ID
d=${DATASETS[$(( (t-1) % 2 ))]}
cd "$IGI_REPO/catalog"
CH=$(IFS=,; echo "${CHROMS_ALL[*]}")

if [ "$t" -le 2 ]; then
  O=$IGI_WORK/release/ont_bulk_rna; W=$IGI_WORK/release/ont_bulk_rna_work
  mkdir -p "$O/$d"
  python3 -u run_ont.py --design design.yaml --paths "$IGI_PATHS" --dataset "$d" --assay bulk_rna \
    --chroms "$CH" --reads 20000000 --min-tpm 0.01 --jobs 48 --records-per-chunk 15000 --label full \
    --out "$O/$d" --work "$W" --seed $(( 9600 + ((t-1) % 2) + 1 ))

elif [ "$t" -le 4 ]; then
  O=$IGI_WORK/release/tenx_gex; W=$IGI_WORK/release/tenx_gex_work
  mkdir -p "$O/$d" "$W"
  python3 -u run_tenx_gex.py --design design.yaml --paths "$IGI_PATHS" --dataset "$d" \
    --chroms "$CH" --cells 4000 --reads-per-cell 40000 --label full \
    --out "$O/$d" --work "$W" --seed 9100

else
  O=$IGI_WORK/release/tenx_tcr; W=$IGI_WORK/release/tenx_tcr_work
  mkdir -p "$O/$d" "$W"
  python3 -u run_tenx_tcr.py --design design.yaml --paths "$IGI_PATHS" --dataset "$d" \
    --catalog-dir output --cells 4000 --reads-per-cell 5000 --molecules-per-cell 120 --label full \
    --out "$O/$d" --work "$W" --seed 9100
fi
echo "DONE task $t $d $(date)"
