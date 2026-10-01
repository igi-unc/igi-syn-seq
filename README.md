# igi-syn-seq

Generator and documentation for **IGI-SYN-SEQ** (datasets IGI-SYN-SEQ-01, HG002 baseline, and IGI-SYN-SEQ-02, IPISRC044 baseline), two fully synthetic, known-truth multi-assay
tumor/normal datasets built to benchmark LENS across every antigen class. Fifteen sample types per
dataset: Illumina WES, bulk RNA, WGS, 10x 5' GEX and TCR; PacBio HiFi WGS, Kinnex bulk and single-cell;
and ONT WGS, bulk RNA and single-cell RNA.

- `docs/pipeline-diagrams.md`: flowcharts of the whole generation process, the two baseline routes,
  the catalog designer, and how one variant must appear across the assays.
- `docs/IGI-SYN-SEQ-design.md`: design specification (baselines, TNBC tumor architecture,
  variant catalog and evidence tiers, assays, truth bundle, chr1to6 subset).
- `docs/data-sources.md`: external inputs and where they live.
- `docs/baseline-references.md`: how the phased germline baselines are built (DeepVariant, WhatsHap,
  SHAPEIT5, Manta, GIAB Q100), written with placeholders so it can be followed on any cluster.
- `docs/catalog-design-notes.md`: how the designed somatic events are chosen and scored, what every
  truth-table column means, and the data provenance and attribution.
- `docs/build-reference.md`: **how each deliverable is actually produced** — which code runs, what it
  reads, what it costs, which constants are fitted from what, and what a given edit forces you to
  regenerate. Start here to correct one aspect without rebuilding everything.
- `catalog/`: the catalog designer and the read builders as a path-free Python package, plus the design
  parameters, the expression baseline, the fitted simulator resources, and the generated truth tables and
  diff cards in `catalog/output/`.
- `jobs/`: the SLURM array jobs that run the build, and the measurement jobs that fit the simulator
  constants. `jobs/README.md` gives the run order.

Status: the exomes and bulk RNA are built at release scale and pass acceptance (40/40 + 5/5 on both
datasets). Illumina and PacBio WGS are validated on chr21 against real-data targets. The Kinnex, 10x and
ONT builders are specified but not written; see `docs/build-reference.md` §1 for per-sample-type status.
