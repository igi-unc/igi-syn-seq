# How IGI-SYN-SEQ is generated

Flowcharts of the generation process, from public inputs to the released FASTQ/BAM files and the truth
bundle. Solid boxes are implemented; dashed boxes are designed but not yet built. Details for each stage
live in `IGI-SYN-SEQ-design.md` (what is generated), `baseline-references.md` (stage 0) and
`catalog-design-notes.md` (stage 1).

## 1. End to end

```mermaid
flowchart TD
    subgraph IN["Public inputs"]
        GIAB["GIAB HG002<br/>Q100 v1.1 phased calls<br/>Revio HiFi WGS"]
        IPI["IPISRC044<br/>blood-normal WGS, tumour WES/RNA<br/>10x 5' GEX + TCR, ONT cDNA"]
        REF["GRCh38 + GENCODE v37<br/>capture BED, RepeatMasker, segdups"]
        EXP["TCGA-BRCA basal-like<br/>transcript TPM"]
        KIN["PacBio Kinnex bulk<br/>HG002 + GIAB MAS-seq"]
    end

    GIAB --> S0
    IPI --> S0
    S0["<b>Stage 0 - baselines</b><br/>phased diploid germline per dataset<br/>small variants + structural variants"]
    REF --> S1
    EXP --> S1
    S0 --> S1
    S1["<b>Stage 1 - catalog</b><br/>choose every designed event,<br/>score binding, assign clone and tier,<br/>emit truth tables + diff cards"]

    S1 -.-> S2
    S0 -.-> S2
    S2["Stage 2 - genomes<br/>per clone, per haplotype FASTA<br/>CNAs, SVs, viral contigs"]
    S2 -.-> S3
    EXP -.-> S3
    S3["Stage 3 - transcriptomes<br/>haplotype-aware transcripts,<br/>isoforms, fusions, ERVs, CTAs,<br/>per-cell expression"]
    KIN -.-> S4
    IPI -.-> S4
    S3 -.-> S4
    S2 -.-> S4
    S4["Stage 4 - simulate reads<br/>error models trained on real data"]
    S4 -.-> S5
    S5["Stage 5 - package<br/>BAM/FASTQ naming, pbindex, skera,<br/>LENS manifests"]
    S5 -.-> S6
    S6["Stage 6 - chr1to6 subset<br/>derived by alignment, not re-simulated"]
    S5 -.-> S7
    S6 -.-> S7
    S7["Stage 7 - validate<br/>align back, check VAF/expression vs truth,<br/>Cell Ranger, skera/lima/isoseq, LENS run"]

    S1 --> TRUTH[("Truth bundle")]
    S5 -.-> REL[("Released FASTQ / BAM")]

    classDef planned stroke-dasharray: 5 5
    class S2,S3,S4,S5,S6,S7,REL planned
```

## 2. Stage 0: how each germline baseline is built

The two datasets reach the same end state, a phased diploid germline, by different routes, because only
one of them comes with a benchmark call set.

```mermaid
flowchart TD
    subgraph A["IGI-SYN-SEQ-01 - HG002"]
        A1["GIAB Q100 v1.1 GRCh38<br/>assembly-derived, already phased"]
        A2["small variants<br/>chromosome-wide phase"]
        A3["structural variants<br/>records >= 50 bp"]
        A1 --> A2
        A1 --> A3 --> A4["SnpEff annotation"]
    end

    subgraph B["IGI-SYN-SEQ-02 - IPISRC044"]
        B1["blood-normal Illumina WGS"]
        B1 --> B2["DeepVariant WGS model<br/>genome-wide, 4.58 M PASS"]
        B2 --> B3["WhatsHap<br/>read-backed blocks<br/>60% of hets, median 128 bp"]
        B2 --> B4["SHAPEIT5 + 1000G panel<br/>90.7% of sites phased"]
        B3 --> B5["combine<br/>SHAPEIT5, else oriented WhatsHap block,<br/>else seeded random; PHASE_SOURCE recorded"]
        B4 --> B5
        B5 --> B6["one phase set per chromosome<br/>zero unphased hets"]
        B1 --> B7["Manta germline mode<br/>8,756 PASS SVs"] --> B8["SnpEff annotation"]
    end

    A2 --> OUT[("phased diploid germline<br/>+ annotated germline SVs")]
    A4 --> OUT
    B6 --> OUT
    B8 --> OUT
```

Note the asymmetry that drove the design: IPISRC044 has no long-read DNA. Its ONT files are long reads of
the 10x 5' single-cell cDNA libraries, so they carry phase only inside expressed exons and cannot support
long-read SV calling. That is why statistical phasing and a short-read SV caller are used for that dataset.

## 3. Stage 1: inside the catalog designer

```mermaid
flowchart TD
    G["GENCODE v37<br/>one representative CDS per gene"] --> C1
    R["GRCh38"] --> C1
    C1["sample random CDS positions<br/>skip drivers, chromothripsis arm,<br/>homozygous deletions"]
    C1 --> C2["mutate in silico<br/>translate, call consequence"]
    C2 --> C3{"consequence"}
    C3 -->|"missense, nonsense, synonymous,<br/>start/stop loss"| P1["candidate pool"]
    C3 -->|"frameshift, in-frame indel"| P1
    C3 -->|"no usable change"| X["discard"]

    P1 --> D1["expression tier<br/>TCGA-BRCA basal median TPM"]
    P1 --> D2["context stratum<br/>homopolymer, segdup, capture edge,<br/>near germline het, deep intronic"]
    P1 --> D3["binding<br/>8-11mers spanning the change,<br/>minus peptides present in wild type"]
    D3 --> D4["netMHCpan 4.1<br/>24 concurrent allele x length jobs<br/>disk cache on every pair"]
    D4 --> D5["best %rank -> strong / weak / non<br/>+ wild-type counterpart and agretopicity"]

    D1 --> S["select into the evidence grid<br/>clonality x expression x binding<br/>plus context strata and phased pairs"]
    D2 --> S
    D5 --> S
    CLN["clone tree + per-haplotype CNAs<br/>purity, WGD"] --> S
    S --> PL["place: haplotype, clone, timing<br/>expected VAF from purity, CCF and local CN"]
    GL["phased germline VCF"] --> PL
    PL --> O1[("truth table<br/>one row per event")]
    PL --> O2[("diff cards<br/>haplotypes + protein windows")]
    PL --> O3[("summary<br/>tier counts, grid fill")]
```

## 4. What a single somatic variant has to look like downstream

Every designed event must show up consistently in each assay that can physically see it. This is the rule
stage 3 and stage 4 implement, and the property stage 7 checks.

```mermaid
flowchart LR
    V["one designed exonic SNV<br/>haplotype h, clone c"] --> DNA["tumour genome<br/>copies of h in clone c"]
    V --> RNA0["tumour transcriptome<br/>weighted by expression tier<br/>and allelic-expression setting"]

    DNA --> WES["tumour WES<br/>VAF from purity, CCF, local CN"]
    DNA --> WGS["PacBio HiFi WGS<br/>same VAF, phased with nearby germline"]
    DNA --> NRM["normal WES / WGS<br/>absent"]
    RNA0 --> BRNA["bulk RNA<br/>present if expressed"]
    RNA0 --> KIN["Kinnex bulk<br/>full-length isoform"]
    RNA0 --> SC["Kinnex single-cell + 10x 5' GEX<br/>only in cells of clone c"]

    WES --> T[("truth row:<br/>expected VAF, tier labels,<br/>peptide, diff card")]
    WGS --> T
    NRM --> T
    BRNA --> T
    KIN --> T
    SC --> T
```

A deep intronic variant reaches only the WGS box; a synonymous variant reaches the RNA boxes but produces
no peptide; a read-through fusion has no DNA box at all. Those asymmetries are the point: they are what
distinguishes a caller that is reading evidence from one that is guessing.
