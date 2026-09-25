# igi-syn-seq

Generator and documentation for **IGI-SYN-SEQ** (datasets IGI-SYN-SEQ-01, HG002 baseline, and IGI-SYN-SEQ-02, IPISRC044 baseline), two fully synthetic, known-truth multi-assay
tumor/normal datasets (Illumina WES/RNA/10x 5' GEX + TCR, PacBio HiFi WGS, Kinnex bulk and
single-cell) built to benchmark LENS across every antigen class.

- `docs/pipeline-diagrams.md`: flowcharts of the whole generation process, the two baseline routes,
  the catalog designer, and how one variant must appear across the assays.
- `docs/IGI-SYN-SEQ-design.md`: design specification (baselines, TNBC tumor architecture,
  variant catalog and evidence tiers, assays, truth bundle, chr1to6 subset).
- `docs/data-sources.md`: external inputs and where they live.
- `docs/baseline-references.md`: how the phased germline baselines are built (DeepVariant, WhatsHap,
  SHAPEIT5, Manta, GIAB Q100), written with placeholders so it can be followed on any cluster.
- `docs/catalog-design-notes.md`: how the designed somatic events are chosen and scored, what every
  truth-table column means, and the data provenance and attribution.
- `catalog/`: the catalog designer (stage 1) as a path-free Python package, plus the design parameters,
  the expression baseline, and the generated truth tables and diff cards in `catalog/output/`.

Stages 2 onward (clone genomes, transcriptomes, read simulation, packaging, subsetting, validation) are
specified in the design document and the diagrams but not yet implemented.
