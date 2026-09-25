# Catalog designer

Seeded, config-driven selection and specification of every designed event in IGI-SYN-SEQ-01 and -02.
No site paths live in this directory: the design is `design.yaml`, and file locations come from a
local `paths.yaml` (template: `paths.example.yaml`).

```bash
python3 -m igi_catalog.designer --design design.yaml --paths /local/paths.yaml \
    --dataset IGI-SYN-SEQ-01 --out output/            # add --scale 0.05 for a smoke test
```

Requirements: Python 3.9+, `pysam`, `pyyaml`; netMHCpan 4.1 reachable through the command template in
`paths.yaml` (a container command works); reference FASTA with `.fai`; GENCODE v37 GTF; the phased germline
VCF of the baseline individual; capture BED; UCSC `rmsk` (bgzip + tabix) and merged `genomicSuperDups` BEDs;
chromosome-arm BED; the expression baseline shipped in `resources/`.

## Outputs per dataset (`output/`)
| File | Content |
|---|---|
| `<ds>.snv_indel.tsv` | one row per designed SNV/indel: locus, alleles, gene, consequence, haplotype, clone, timing, clonality tier, expected tumor VAF, multiplicity, local copy number, gene TPM and tier, allelic-expression setting, binding tier with best peptide/allele/%rank, wild-type counterpart peptide and its rank plus agretopicity, context stratum and raw context features, flagpost and chr1to6 flags, wild-type/mutant haplotype and protein windows |
| `<ds>.snv_indel.diffcards.txt` | one card per event: reference, both germline haplotypes, mutant haplotype with the change marked, wild-type vs mutant protein window with the first changed residue marked |
| `<ds>.snv_indel.summary.json` | counts by class, tier and context; grid fill; chr1to6 fraction |
| `<ds>.fusions.tsv` | one row per designed fusion: partners and transcripts, exon junction, genomic breakpoints, DNA mechanism, frame, fusion protein length, junction window, clone and expected DNA VAF, expression and binding tiers, WES visibility |
| `<ds>.fusions.diffcards.txt` | one card per fusion: the junction shown against both parent transcripts |

## Validation

```bash
python3 validate_catalog.py --table output/<ds>.snv_indel.tsv --reference /path/GRCh38.fa
```

Eight self-consistency checks, non-zero exit on failure so it can gate a regression run: REF alleles agree
with the reference; the mutant window carries ALT at the marked offset; window lengths track the indel
size; expected VAF is a probability; median VAF decreases down the clonality ladder; a binding tier implies
a peptide and `na` implies none; the recorded tier matches the recorded %rank; and a wild-type counterpart
is a same-length, different peptide.

`resources/tcga_brca_basal.transcript_tpm.tsv.gz` is the tumor expression baseline: per-transcript median TPM
across 191 TCGA-BRCA basal-like tumors (PanCanAtlas `Subtype_mRNA = Basal`) from the UCSC Xena Toil kallisto
table; `resources/tcga_brca_basal.samples.txt` lists the samples.

Module map: `annotation` (GTF, CDS/protein, coordinate mapping, in-silico mutation), `genome`, `intervals`,
`expression`, `germline` (phased haplotype windows), `context` (strata), `clones` (clone tree, copy number,
expected VAF), `binding` (netMHCpan), `diffcards`, `snv_indel` (the designer), `designer` (CLI).
