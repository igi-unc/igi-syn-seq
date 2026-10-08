#!/bin/bash
#SBATCH -J igi_calib
#SBATCH -c 12 --mem=32G -t 6:00:00
#SBATCH --array=1-15
# Calibrate every long-read identity setting against the ALIGNMENT, not against the reads' own quality
# strings. See jobs/measure/alignment_identity.py for why the distinction matters; in short, Badread
# generates sequence from an error model and qualities from a separate qscore model, and for the PacBio
# model they disagree by a factor of five.
#
# Three settings are swept, four requests each:
#
#   pacbio_identity   drives PacBio HiFi WGS **and both Kinnex arms** (run_pacbio_wgs.py,
#                     run_kinnex_bulk.py, run_kinnex_sc.py all read this one key), target 0.99823
#   ont_wgs_identity  ONT genomic, target 0.98676
#   ont_identity      ONT cDNA -- single-cell AND bulk RNA, target 0.98220
#
# Reads are simulated from an unmodified reference window, so the alignment is ground truth: there are no
# germline or designed variants to mistake for error. One 10 Mb window at 2x is a few thousand reads and
# tens of millions of aligned bases, which pins the error rate far tighter than the 0.3-point offsets
# being measured.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
W=$IGI_WORK/calibrate; mkdir -p "$W"
MM="singularity exec -B $IGI_BIND $(dirname $IGI_SAMTOOLS_IMG)/quay.io_biocontainers_minimap2_2.28--h577a1d6_4.sif minimap2"
SAM="singularity exec -B $IGI_BIND $IGI_SAMTOOLS_IMG samtools"
REF=$IGI_FASTA_NOALT
BR=$(python3 -c "import yaml;print(yaml.safe_load(open('$IGI_PATHS'))['badread_cmd'])")

# setting | request | error/qscore model | length mean,sd | minimap2 preset
CASES=(
  "pacbio|99.4,99.8,0.5|pacbio2021|16689,4593|map-hifi"
  "pacbio|99.7,99.95,0.3|pacbio2021|16689,4593|map-hifi"
  "pacbio|99.85,99.99,0.2|pacbio2021|16689,4593|map-hifi"
  "pacbio|99.95,100,0.1|pacbio2021|16689,4593|map-hifi"
  "ont_wgs|99.17,99.6,1.5|nanopore2023|19117,15530|map-ont"
  "ont_wgs|99.4,99.8,1.2|nanopore2023|19117,15530|map-ont"
  "ont_wgs|99.6,99.9,1.0|nanopore2023|19117,15530|map-ont"
  "ont_wgs|99.8,99.95,0.8|nanopore2023|19117,15530|map-ont"
  "ont_cdna|98.41,99.6,1.5|nanopore2023|2091,843|map-ont"
  "ont_cdna|98.8,99.7,1.3|nanopore2023|2091,843|map-ont"
  "ont_cdna|99.1,99.8,1.1|nanopore2023|2091,843|map-ont"
  "ont_cdna|99.4,99.9,1.0|nanopore2023|2091,843|map-ont"
  # Tasks 13-14 test whether the PacBio model has an error FLOOR --identity cannot reach past. The
  # first four requests gave a constant ~3.0e-3 EXCESS over the nominal error rather than a constant
  # ratio, which is the signature of a floor; asking for perfect reads settles it.
  "pacbio|99.99,100,0.05|pacbio2021|16689,4593|map-hifi"
  "pacbio|100,100,0|pacbio2021|16689,4593|map-hifi"
  # Task 15 checks the floor at KINNEX read lengths. The floor was measured on 16.7 kb genomic reads and
  # both Kinnex arms share pacbio_identity, but their segments are 2 kb cDNA, so the per-base error of a
  # short read under the same model is worth confirming before committing 20-44 h per library to it.
  "kinnex|100,100,0|pacbio2021|2091,843|map-hifi"
)
IFS='|' read -r SET REQ MODEL LEN PRESET <<<"${CASES[$((SLURM_ARRAY_TASK_ID-1))]}"
TAG="${SET}_${REQ//,/-}"
echo "=== $SET request=$REQ model=$MODEL length=$LEN preset=$PRESET ==="

WIN=$W/win.fa
if [ ! -s "$WIN" ]; then $SAM faidx "$REF" chr1:20000000-30000000 > "$WIN"; fi
FQ=$W/$TAG.fq
if [ ! -s "$FQ" ]; then
  ARGS="simulate --reference $WIN --quantity 2x --length $LEN --identity $REQ \
    --error_model $MODEL --qscore_model $MODEL --seed 81 \
    --start_adapter_seq '' --end_adapter_seq '' --junk_reads 0 --random_reads 0 --chimeras 0"
  eval "${BR/\{args\}/$ARGS}" > "$FQ" 2>$W/$TAG.badread.err
fi
[ -s "$FQ" ] && [ "$(( $(wc -l < "$FQ") % 4 ))" -eq 0 ] \
  || { echo "badread produced a malformed FASTQ for $TAG"; exit 1; }

B=$W/$TAG.bam
if [ ! -s "$B" ] || [ "$FQ" -nt "$B" ]; then
  # Aligned against the WINDOW, not the whole genome: the reads come from it, so the window is the ground
  # truth, and a whole-genome index would let a read place off-window and show as error that is really
  # mismapping. --secondary=no keeps one alignment per read.
  $MM -ax "$PRESET" -t 12 --secondary=no "$WIN" "$FQ" 2>$W/$TAG.mm.err | $SAM sort -@4 -o "$B" -
  $SAM index "$B"
fi
cd "$IGI_REPO"
python3 -u jobs/measure/alignment_identity.py --bam "$B" --out "$W/$TAG.json" \
  --label "$SET request=$REQ" --max-reads 20000
echo "DONE $TAG"
