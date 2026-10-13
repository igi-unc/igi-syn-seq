#!/bin/bash
#SBATCH -J igi_kinnex_sc
#SBATCH -c 48 --mem=80G -t 72:00:00
#SBATCH --array=1-2
# Kinnex single-cell RNA (PacBio MAS-seq, 10x 5' v2 cDNA), one task per dataset. Produces the pre-skera
# hifi_reads BAM, the segmented BAM from the real `skera split`, the cell roster and the clonotype table.
#
# Memory is 80G against a MEASURED peak of 56.5 GiB (2026-10-13, both tasks, 15 h 20). It ran at 64G and
# completed, at 88 % of the limit -- too thin to leave, since the next design change that raises the
# molecule count turns a 15 h run into an OOM at hour 14. Unlike the bulk arm, this one cannot defer the
# work: every segment carries a distinct barcode and UMI, so each molecule's sequence genuinely differs
# and there is nothing to apply lazily. See run_kinnex_sc.py's mol_objs list.
#
# The roster and the molecule pool are built from --seed scoped by dataset, and every single-cell assay must
# use the SAME --seed or a barcode stops meaning one cell across assays within a dataset. 10x GEX, 10x TCR
# and ONT single cell all take 9100.
#
# --mean-molecules is 10,600 rather than the 8,000 default because --segments is a POST-SKERA target and two
# measured losses sit in front of it: size selection keeps 69 % and skera returns 48.6 % of cDNAs as clean
# single-cDNA segments. 12 M segments therefore needs about 35.8 M molecules, which 4,000 cells only reach at
# ~10,600 captured molecules each. The first run took --segments as a molecule count and delivered 4.0 M
# segments, a third of spec.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
O=$IGI_RELEASE/kinnex_sc; W=$IGI_RELEASE/kinnex_sc_work
d=${DATASETS[$((SLURM_ARRAY_TASK_ID-1))]}
CH=$(IFS=,; echo "${CHROMS_ALL[*]}")
mkdir -p "$O/$d"
cd "$IGI_REPO/catalog"
python3 -u run_kinnex_sc.py --design design.yaml --paths "$IGI_PATHS" --dataset "$d" --chroms "$CH" \
  --cells 4000 --segments 12000000 --mean-molecules 10600 --min-tpm 0.01 --jobs 48 --label full \
  --out "$O/$d" --work "$W" --seed 9100
echo "DONE $d $(date)"
