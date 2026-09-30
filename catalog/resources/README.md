# Catalog resources

Reference material the catalog designer and read builders depend on. Provenance and licensing for each
is recorded here and in `docs/catalog-design-notes.md`.

## Expression baseline

`tcga_brca_basal.transcript_tpm.tsv.gz`, `tcga_brca_basal.samples.txt`

Per-transcript median TPM across 191 TCGA-BRCA basal-like tumours, from the UCSC Xena Toil recompute,
plus the contributing sample barcodes. Medians only, never per-sample values.

## COSMIC signatures

`cosmic_v3.4_sbs_subset.tsv`

Five COSMIC v3.4 SBS profiles (GRCh38) used by the passenger background model. Provenance in the file
header.

## ART quality profiles

`ipisrc044_novaseqx_R1.profile.txt`, `ipisrc044_novaseqx_R2.profile.txt`

Per-cycle base-quality profiles for `art_illumina -1/-2`, fitted with `art_profiler_illumina` from
4,000,000 read pairs of the IPISRC044 Personalis blood-normal WGS (NovaSeq X, instrument LH00499),
150 cycles.

The design asks for a profile from IPISRC044 WES, but this patient has no 150 bp exome: the BostonGene
exomes are 101 bp (normal) and 140 bp (tumour) while the design specifies PE150. A per-cycle quality
profile is a property of the sequencer run rather than the capture, and capture effects are modelled
separately, so the profile comes from the 150 bp WGS and PE150 is kept.

The fitted profile carries the four-level quality binning these instruments use, Q3/Q13/Q25/Q41, in
place of the continuous scores of ART's stock HS25 (HiSeq 2500) model. Some tools behave differently on
binned qualities; that is a property of current data, and the benchmark should show it rather than hide
it behind a 2013 error model.
