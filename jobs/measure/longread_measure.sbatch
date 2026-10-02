#!/bin/bash
#SBATCH -J igi_lr_measure
#SBATCH -c 4 --mem=16G -t 6:00:00
# Measure the real long-read sets that the pacbio_* and ont_* constants in the site paths file come from.
# Run this to recalibrate; then apply the model's identity offset (see measure_longread.py) before
# writing the value into the paths file, and re-measure the simulated output to confirm where it landed.
#
# GIAB's own HG002 ONT is all R9.4-era -- 2D reads from 2016 basecalled with Guppy v2/v3 -- and is not a
# usable R10.4.1 target, which is why the ONT WGS numbers come from ONT's open-data bucket instead
# (ont_wgs_bam in the site paths file). The ONT arms below are not written yet.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
M="python3 $IGI_REPO/jobs/measure/measure_longread.py --n 200000"
ST="singularity exec -B $IGI_BIND $IGI_SAMTOOLS_IMG samtools"
K=$IGI_HG002/Kinnex_fullLengthRNA_HG002_Revio

# HiFi genomic: pacbio_length_mean / pacbio_length_sd, and the accuracy target behind pacbio_identity
$M --samtools "$ST" --label "WGS HiFi (Revio genomic)" \
  --bam "$IGI_HG002/PacBio_HiFi-Revio_20231031/HG002_PacBio-Revio_m84039_230928_213653_s3.hifi_reads.bam"

# Post-skera segmented: kinnex_segment_length_mean / kinnex_segment_length_sd
$M --samtools "$ST" --label "Kinnex bulk, segmented (post-skera)" --bam "$K/1-Sreads/segmented.bam"
$M --samtools "$ST" --label "Kinnex bulk, FLNC" --bam "$K/2-FLNC/flnc.bam"

# Pre-segmentation, for the MAS array structure the Kinnex sc builder has to reproduce. The sc array holds
# 16 cDNAs (a median of 16 10x TSOs per read); only nine of its adapters are identified. See
# build-reference.md section 6.3.
$M --samtools "$ST" --label "Kinnex single cell, 0-CCS (pre-segmentation)" \
  --bam "$IGI_REF/inputs/bams/PacBio_Kinnex_scRNA_HG002_10x5p_2024/pacbio/0-CCS/m84039_240124_190648_s3.hifi_reads.bcM0002.bam"
echo "DONE $(date)"
