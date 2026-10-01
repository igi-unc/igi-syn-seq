# Site configuration for the IGI-SYN-SEQ job scripts. Copy to jobs/env.sh (untracked) and edit.
#
# Every job script sources this file and refers to these variables, so no tracked file contains an
# absolute path. jobs/env.sh is the only place a path appears outside config/sites/*.paths.yaml.

# This repository.
IGI_REPO=/path/to/igi-syn-seq

# Site paths file for the catalog and read builders.
IGI_PATHS=$IGI_REPO/config/sites/unc-lccc.paths.yaml

# Scratch root for build outputs, working directories and logs. Needs several TB.
IGI_WORK=/path/to/scratch/igi-syn

# Reference and container root (the RAFT workspace at this site).
IGI_REF=/path/to/raft

# Containers and references used by the alignment and measurement jobs.
IGI_BWA_IMG=$IGI_REF/imgs/michaelfranklin-bwasamtools-0.7.17-1.10.img
IGI_SAMTOOLS_IMG=$IGI_REF/imgs/quay.io-biocontainers-samtools-1.21--h96c455f_1.img
IGI_ART_IMG=$IGI_REF/imgs/quay.io_biocontainers_art_2016.06.05--h0704011_13.sif

# The noalt build, not the primary FASTA: the primary carries no bwa index, and alt contigs divert reads
# away from their primary locus, which is what buried the dataset 02 MHC in an earlier build.
IGI_FASTA_NOALT=$IGI_REF/references/homo_sapiens/fasta/Homo_sapiens.assembly38.noalt.fa
IGI_ARMS_BED=$IGI_REF/references/homo_sapiens/ucsc/hg38.arms.bed
IGI_EXOME_BED=$IGI_REF/references/homo_sapiens/beds/hg38_exome.bed

# Bind path for singularity (the filesystem holding $IGI_REF and $IGI_WORK).
IGI_BIND=/path/to/shared/filesystem

# SLURM partitions.
IGI_PARTITION=crunchnodes,allnodes

DATASETS=(IGI-SYN-SEQ-01 IGI-SYN-SEQ-02)
CHROMS_ALL=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)

# Public reference datasets used only by the measurement jobs.
IGI_HG002=/path/to/datasets/HG002
