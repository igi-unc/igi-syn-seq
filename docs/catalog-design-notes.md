# Catalog design notes

How the designed somatic events of IGI-SYN-SEQ-01 and -02 are chosen, what each truth-table column means,
and what the final catalogs contain. Code: `catalog/` (module map in `catalog/README.md`). Design
parameters: `catalog/design.yaml`. Final tables: `catalog/output/`.

## 1. Principles

- **Seeded and reproducible.** One seed in `design.yaml` drives every random choice, keyed by dataset,
  so the two catalogs are independent draws from the same design.
- **Grid first, realism second.** Events are chosen to fill an evidence grid (clonality x expression x
  MHC binding x sequence context) rather than to mimic a mutational spectrum. Realism enters through the
  clone tree, the copy-number model, the expression baseline and the germline background, not through
  the choice of designed loci. Signature-driven passenger mutations are added separately at genome
  construction time and are not part of this catalog.
- **Every event carries its own evidence.** The truth table records the expected tumor VAF given purity,
  clone and local copy number; gene expression; the best neopeptide and its %rank; the context stratum;
  and a haplotype-aware diff card, so a missed call can be traced to a specific axis.

## 2. Data provenance and attribution

Both germline baselines are public datasets, and the published tables therefore contain derived sequence
from them: the diff cards and the `wt_hap_window` / `mut_hap_window` columns print each individual's
reconstructed haplotypes over a 41 bp window per event (about 71 kb per dataset, disclosing roughly 30
genotypes), which is what makes the `near_germline_het` stratum debuggable.

| Baseline | Source | Terms |
|---|---|---|
| HG002 / NA24385 | GIAB Ashkenazi trio son; Q100 v1.1 GRCh38 benchmark. HLA types from Chin et al. 2020, Nat Commun 11:4794, Supplementary Table 4 | consented public reference material |
| IPISRC044 | public tumor/normal dataset; germline called and phased here (see `baseline-references.md`) | public |

The expression baseline is derived from TCGA-BRCA RNA-seq, an open-access tier, via the UCSC Xena Toil
recompute. Cite: the TCGA Research Network (https://www.cancer.gov/tcga); Vivian et al. 2017, Nat
Biotechnol 35:314 (Toil); Goldman et al. 2020, Nat Biotechnol 38:675 (Xena). `resources/` holds only
per-transcript medians across 191 tumors plus the contributing sample barcodes, not per-sample values.

`resources/cosmic_v3.4_sbs_subset.tsv` holds five columns of the COSMIC Mutational Signatures v3.4 single
base substitution profiles for GRCh38 (SBS1, SBS2, SBS3, SBS5, SBS13), used by the passenger background
model. Cite: COSMIC Mutational Signatures v3.4 (https://cancer.sanger.ac.uk/signatures); Alexandrov et
al. 2020, Nature 578:94. The file was taken from the copy redistributed in AlexandrovLab's
SigProfilerAssignment (`COSMIC_v3.4_SBS_GRCh38.txt`, md5 `633b2e6708d168fb9f6e3942e8d69a1f`) and reduced
to the five signatures in use; the header records that provenance.

netMHCpan 4.1 is not redistributed here; the code calls whatever binary the site configuration points at,
under that user's own licence.

## 3. Inputs

| Input | Role |
|---|---|
| GENCODE v37 GTF | gene models; one representative CDS per protein-coding gene (MANE Select, else Ensembl canonical, else longest complete CDS) |
| GRCh38 reference | sequence, homopolymer runs, GC |
| phased germline VCF of the baseline individual | haplotype windows for diff cards; "near germline het" context; haplotype assignment of somatic events |
| TCGA-BRCA basal-like transcript TPM (191 tumors, Xena Toil kallisto) | gene expression tiers (`resources/`) |
| capture BED, UCSC RepeatMasker, UCSC segmental duplications, cytoband arms | context strata; arm-level copy-number regions |
| netMHCpan 4.1 | MHC-I binding of neopeptides against the patient's six class I alleles |

## 4. Method

### 4.1 Candidate generation
Random CDS positions are drawn from representative transcripts on chr1-22/X, excluding driver genes,
the chromothripsis arm and homozygous deletions. Each position is mutated in silico (all three
substitutions; for indels a deletion of 1-30 bp or a random insertion of the same lengths) and the
consequence is called from the translated CDS: missense, nonsense, synonymous, start/stop loss,
frameshift, in-frame insertion/deletion. Frameshifts are translated into the 3' UTR so the neo-ORF is
complete. Extra candidates are drawn inside LOH and amplified regions so those clonality tiers can be
filled.

### 4.2 Binding tier
Every 8-11mer that spans a changed residue (the whole neo-ORF for a frameshift) is scored with
netMHCpan 4.1 (`-BA`, eluted-ligand %rank) against the patient's six class I alleles. Peptides identical
to a wild-type peptide are discarded. The best remaining %rank sets the tier: strong <= 0.5, weak <= 2,
non > 2. Nonsense and synonymous events have no neopeptide and carry tier `na`.

Scoring the full 8-11mer range rather than 9-mers alone was a deliberate choice: 9-mers dominate class I
presentation but ignoring the other lengths would understate the epitope space a caller has to search.
The concern with using all lengths is that taking the best %rank over roughly 38 peptides instead of 9
could collapse the `non` tier, since more draws make a low rank more likely. Measured on 60 candidate
windows the spread stays usable: 52 % strong, 37 % weak, 12 % non, because overlapping peptides from one
window are correlated rather than independent draws.

### 4.3 Expression tier
Gene TPM is the sum of the baseline medians of the gene's transcripts (GENCODE v23 IDs mapped by
version-less ENST). Tiers: T0 < 0.5, T1 < 3, T10 < 30, T100 < 300, T1000 >= 300 TPM. A per-event
allelic-expression setting (balanced 0.5, silenced 0.1, dominant 0.9 mutant fraction) is drawn for
expressed genes and applied at transcriptome construction.

### 4.4 Context strata
Priority order: `deep_intronic` (outside the capture BED), `segdup`, `homopolymer` (run >= 6 touching
the site), `near_germline_het` (phased het within 30 bp), `exon_edge` (<= 10 bp from a CDS boundary),
`repeat_other`, `clean`. Raw features are kept as `ctx_*` columns. Phased pairs are two missense SNVs
<= 150 bp apart on the same haplotype in the same clone, cross-referenced by `paired_event`.

### 4.5 Clone assignment and expected VAF
Clones and CCFs come from `design.yaml`. The copy-number model keeps per-haplotype copy numbers per
clone (arm-level and focal events, inherited by descendants; WGD datasets start at 2+2). For an event on
haplotype h acquired in clone c:

    VAF = purity * sum_over_clones(excl_fraction * multiplicity) / (purity * mean_tumor_CN + 2 * (1 - purity))

Truncal events with `timing = pre_cna` carry multiplicity equal to the current copies of their
haplotype (this also encodes pre-WGD timing); post-CNA and subclonal events carry 1 where the haplotype
is retained and 0 where it is lost. Clonality tiers: `T_LOH` (truncal, in a one-haplotype region),
`T_amp` (truncal, haplotype copy number >= 4), `T_het`, and the subclone names.

### 4.6 Sites excluded by construction

A somatic allele whose reference span overlaps a germline variant is rejected at candidate generation.
Such a site cannot be written unambiguously in reference coordinates, and the haplotype edit would not
apply cleanly, so the truth row would claim a change the sequence does not contain. The first full run
produced exactly one such row before the check existed: a 21 bp deletion whose reference allele ran
through a nearby germline substitution, which the haplotype builder skipped silently. A skipped somatic
edit is now an error rather than a silent no-op, and `validate_catalog.py` checks window lengths against
the indel size so the class of failure cannot reappear unnoticed.

### 4.7 Diff cards
Every designed event gets a card showing the patient's own sequence before and after the change, so a
debugging session never has to reconstruct what was supposed to happen. For an SNV or indel the card
carries the 41-bp reference window, both germline haplotypes with the individual's phased variants
applied, the mutant haplotype with the change marked, and the wild-type against mutant protein window
with the first changed residue marked (for a frameshift, the whole neo-ORF start):

```
## IGI-SYN-SEQ-01-SNV-0001  TP53 R248Q missense  clone=T tier=T_LOH VAF=0.5385 expr=T10 bind=weak ctx=clean
  window    chr17:7674200-7674240  somatic C>T at 7674220 on hap0
  ref       TGATGGTGAGGATGGGCCTCCGGTTCATGCCGCCCATGCAG
  hap0      TGATGGTGAGGATGGGCCTCCGGTTCATGCCGCCCATGCAG
  hap1      TGATGGTGAGGATGGGCCTCCGGTTCATGCCGCCCATGCAG
  hap0+mut  TGATGGTGAGGATGGGCCTCTGGTTCATGCCGCCCATGCAG
                                ^
  protein   aa 236-260  first change at aa 248 (missense)
  WT        YMCNSSCMGGMNRRPILTIITLEDS
  MUT       YMCNSSCMGGMNQRPILTIITLEDS
                        ^
  peptide   NQRPILTII HLA-B*38:01 rank 1.472
```

A haplotype line that differs from `ref` marks a germline variant inside the window, which is what makes
the card useful for the `near_germline_het` stratum: the caller must get both changes right. A fusion card
instead shows the junction against both parent transcripts, so a wrong exon or frame is visible directly:

```
## IGI-SYN-SEQ-01-FUS-0001  SEC16A-NOTCH1 in_frame
  fusion    SEC16A(ex1-16) :: NOTCH1(ex18-end)  DUP  in-frame
  DNA       chr9:136456492(-) -> chr9:136510402(-)
  5' parent  ...TTGATCAGCCAGCTTGTGCAG | ATGGCTTCCCAGTTACGACTC...   (SEC16A continues)
  3' parent  ...ACATCGACGACTGCCGGCCCA | ACCCGTGTCACAACGGGGGCT...   (NOTCH1 continues)
  fusion     ...TTGATCAGCCAGCTTGTGCAG | ACCCGTGTCACAACGGGGGCT...
                                     ^ junction
```

The same windows are also columns of the truth table (`wt_hap_window`, `mut_hap_window`,
`wt_protein_window`, `mut_protein_window`, `junction_window`) so they can be joined against caller output
programmatically rather than read by eye.

### 4.8 Flagposts
Hotspot substitutions from `design.yaml` (TP53 R248Q/R273H, PIK3CA H1047R/E545K, KRAS G12D, BRAF V600E,
IDH1 R132H, NRAS Q61R, EGFR L858R, CTNNB1 S45F, AKT1 E17K, ESR1 Y537S) are resolved to genomic
positions in the representative transcript and placed truncally with pre-CNA timing.

## 5. Results

_Filled in from `catalog/output/*.summary.json` when the full runs complete._

## 6. Structural variants, viruses, CTAs, ERVs and splice variants

**Structural variants.** DEL, DUP, INV, TRA and mobile-element insertions across five log-spaced size
classes from 50 bp to 5 Mb. Half are placed deliberately inside coding sequence: an exon span is chosen,
the breakpoints are widened into the flanking introns, and the consequence is computed from how much CDS
the event removes or duplicates, giving whole-gene loss, in-frame or out-of-frame exon loss, exon
duplication and exon inversion. Each record carries the affected exons, the coding bases involved, and
whether either breakpoint falls in a capture region.

**Viruses.** Integrations carry a host junction, a haplotype and an expected DNA VAF, so they appear in
WGS and, when the site is captured, in WES. Episomes carry copy number and no junction. The trace EBV in
the normal sample is kept as a realistic negative that a caller should not report as tumour-specific.

**CTAs.** Drawn from the LENS cancer-testis gene list, spread across the five expression tiers, with some
restricted to a single clone so single-cell data can separate them. Each carries the normal-tissue 95th
percentile TPM from the pan-normal reference, so tumour restriction is checkable rather than assumed.

**ERVs.** Split into tumour-specific (absent from the normal panel), tumour-associated (low in normals and
raised in the tumour) and unexpressed negatives. Loci come from the curated ERV set and the full Hsap38
set. Where a locus has no annotated protein, the longest methionine-started ORF over the three forward
frames of the locus sequence is derived and scored, since an ERV with no peptide cannot be benchmarked as
an antigen.

**Splice variants.** Half are tumour-specific, each caused by a designed somatic base change at a
canonical donor or acceptor dinucleotide, recorded so the genomic cause and the RNA consequence can be
checked together; half are tumour-associated isoform switches with no genomic cause. Five mechanisms are
used: exon skip, intron retention, cryptic donor, cryptic acceptor and novel exon. An event is only placed
on an exon where its mechanism is realizable, because short exons and short introns make some mechanisms
impossible: three of the five fail outright on TP53, whose exon 3 is 22 bp. The novel isoform is
translated and compared against the reference protein, and peptides spanning the changed residues are
scored.

## 7. Gene fusions

Partners come from a list of cancer-relevant pairs, a list of adjacent same-strand pairs used for
read-through transcripts, and novel pairs drawn from expressed genes. For a chosen pair the designer
searches exon junctions for one with the requested frame, no stop codon before the junction, and a
usable protein, preferring junctions near the middle of both partners. The DNA mechanism follows from
the partners' positions and orientations: different chromosomes give a translocation, opposite strands
an inversion, and same strand a deletion or tandem duplication depending on order. Genomic breakpoints
are placed inside the introns flanking the chosen exons, which is what makes a subset visible in WES
(`wes_visible`) when a breakpoint falls inside a captured region. Read-through fusions get no DNA
breakpoint at all and no DNA VAF: they are the tumor-associated control that must be called from RNA
only. Junction neopeptides are the 9-mers spanning the junction codon after removing any peptide that
occurs in either parent protein.

## 8. Performance notes

The netMHCpan build available here scores about 40 peptide-allele pairs per second in one process and
forks a separate process per allele and per length, so scoring whole peptide windows serially would have
taken days per dataset. Three changes make the full 8-11mer range affordable:

- submit only the peptides that span a changed residue, as a peptide list (`-p`), instead of whole windows;
- treat each (allele, length) pair as its own concurrent job, 24 of them for six alleles and four lengths,
  so requesting more lengths adds parallel width rather than wall time;
- cache every (allele, peptide) result on disk, which makes re-runs and overlapping windows free.

Throughput scales with the cores available to the job, so the design run is given a 24-core allocation.

## 9. Known limitations

- Binding is class I only and netMHCpan only. The mhcflurry models bundled with LENS are an older layout
  that the 2.1.1 image does not load, so the planned second opinion is not yet wired in.
- Expression tiers are gene-level medians of a cohort, not the sample's own profile; the transcriptome
  step redistributes gene TPM over isoforms.
- Randomly oriented germline sites (4.7 % of IPISRC044 hets) are not yet excluded from flagpost placement.
- Some grid cells are structurally hard to fill: amplified regions contain only a few dozen genes, so
  `T_amp` combined with a rare expression tier has few candidates. The summary reports the fill rate per
  run rather than forcing these cells.
- The negative-control set (germline look-alikes, RNA-editing sites, pseudogene mismatches) is the one
  designed class still outstanding.
- Splice events disrupt donors only, never acceptors, because of how the causal site is chosen.
- The off-target coverage bands use one flat depth per band rather than a fitted decay curve, and the
  band depths are assumed rather than measured from a real capture library.
- Binding is scored for the whole candidate pool (about 12,000 missense candidates) although only some
  600 are placed, because a candidate's tier has to be known before it can be assigned to a grid cell.
  Scoring in waves until each cell fills would cut this several-fold and is the obvious next optimization.

## 10. Wild-type counterparts

For every event whose best mutant peptide keeps the wild-type reading frame (substitutions, start and
stop loss), the same window in the wild-type protein is scored against the same allele and recorded as
`wt_peptide`, `wt_peptide_rank_el` and `agretopicity` (wild-type over mutant affinity; above 1 means the
substitution genuinely improved binding). This distinguishes a novel epitope from one whose wild-type
version binds equally well, which a benchmark needs in order to judge false positives fairly. Indels and
frameshifts shift everything downstream out of register, so no counterpart is defined and the columns are
left blank rather than guessed.

## 11. Passenger background

`igi_catalog/background.py`, driven by `run_background.py`, writes `<dataset>.background.tsv`.

Mutations are drawn from a weighted mixture of COSMIC v3.4 SBS profiles: for IGI-SYN-SEQ-01
SBS3 0.60 / SBS1 0.15 / SBS5 0.15 / SBS13 0.10 at 1.0 per Mb; for IGI-SYN-SEQ-02
SBS2 0.25 / SBS13 0.25 / SBS3 0.20 / SBS1 0.15 / SBS5 0.15 at 4.0 per Mb.

A COSMIC profile is a distribution over the 96 pyrimidine-centred trinucleotide channels of *observed*
mutations, not a per-context rate. The sampler therefore draws a channel first and a genomic position
with the matching trinucleotide second, which reproduces the profile directly and needs no
context-frequency correction. Positions are rejected if they contain an N, overlap a germline variant,
or fall in a span already reserved by a designed event. Each mutation takes a clone by branch share
(`background_truncal_fraction`, 0.60 by default, with the remainder split between branches in proportion
to their cell fraction), a retained haplotype at that locus, and a timing that follows the clone. Short
indels are mixed in at `background_indel_fraction` (0.08).

Passengers that land in coding sequence are annotated with the same `CodingModel` the designed events
use, and protein-changing ones are scored through netMHCpan. A passenger missense in a well expressed
gene is a genuine neoepitope; a truth set that omitted it would score a caller's correct answer as a
false positive.

Measured on 3,817 SNVs across chr1-chr5 of IGI-SYN-SEQ-02: cosine similarity to the target mixture
0.9983, chi-square 91.9 on 87 degrees of freedom. Composition matched the configuration at 4.0 per Mb,
7.2 % indels, 57.6 % truncal and 1.3 % coding.

The background scores against its own netMHCpan cache file. The cache is rewritten whole on each flush,
so a background run sharing the catalog designer's cache would discard whatever the designer had added
since it started.

## 12. Designed germline alleles

`igi_catalog/germline_spikein.py`, driven by `run_germline_spikein.py`, writes the dataset's own
`germline.vcf.gz` plus a truth table and a diff card.

The allele is specified by coding position, not genomic coordinate. The script resolves it against the
representative transcript, requires the span to lie within one exon, builds the VCF record with an
anchor base and left-aligns it, then checks the codon number, the wild-type residue and the consequence
against what the design declared. A hardcoded coordinate would move silently under a different
annotation release; this fails loudly instead.

Both datasets carry BRCA1 c.68_69delAG, which resolves to chr17:43124027 ACT>A and yields a 38 aa
protein from 1,863 - p.Glu23ValfsTer17, the published consequence of the 185delAG founder allele. It is
phased onto the haplotype that the truncal 17q loss retains, so the tumour keeps only the pathogenic
copy. The baseline VCF is never read back through the `germline_vcf` key once the spiked file exists,
so re-running cannot double-spike or report its own allele as a pre-existing conflict.

## 13. Causal splice-site variants

`run_splice_causal.py` turns the `causal_variant` column of the expressed-class table into rows in the
column layout `genome_build.read_events` accepts, so that a tumour-specific splice junction has a
genomic cause in the sequence rather than only in the table. Each variant is checked against the
reference before it is written, and the reference dinucleotide is recorded.

One detail worth stating: the designer takes the donor dinucleotide downstream of the chosen exon, which
on a minus-strand gene reads as its reverse complement on the plus strand. All 15 tumour-specific events
in IGI-SYN-SEQ-02 sit at canonical donors (8 GT, 7 AC). Acceptor disruption is therefore not currently
exercised, and the substituted base is now chosen so the resulting dinucleotide is not recognised by
either spliceosome - the previous rule could turn a GT donor into AT, which is the U12 donor.
