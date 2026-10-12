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

Status (2026-10-12, verified against the release tree): **all fifteen sample types are built at release
scale for both datasets.** Acceptance coverage is much thinner than that sounds -- exomes pass 44/44,
bulk RNA 5/5, Illumina WGS 40/40, and the remaining eleven sample types have no acceptance stage at all.
`docs/build-reference.md` §1 gives the per-sample-type picture in two separate columns, built and
accepted, because they are not the same claim.

Known defects, both found 2026-10-11 and fixed in the generator but **not yet rebuilt**: all four
IGI-SYN-SEQ-02 normal libraries were written at 69 % of their designed depth (§24), and the 10x TCR
libraries were sequenced over GEX's 5' window so Cell Ranger pairs only 38-54 % of cells (§25). Both are
read-level defects with correct underlying truth.
