#!/usr/bin/env python3
"""Assemble the truth bundle a benchmark user scores against.

    make_truth_bundle.py --dataset IGI-SYN-SEQ-01 --catalog-dir output --design design.yaml \
                         --arms-bed <arms.bed> --release-dir <release root> \
                         --manifest-dir <staged metadata dir> --out <.../truth>

WHY THIS EXISTS. Design §2 and §11 put a truth bundle at `inputs/metadata/<dataset>/truth/`, listing a
somatic VCF, CNA segments, SV and fusion tables, expression tables, cells and clonotypes. That directory
held only the four manifests and a README: **the bundle was never assembled.** The truth itself all
exists, and that is what made the gap easy to miss -- catalog tables in the repo under
`catalog/output/`, per-arm read maps and molecule tables and cell rosters on scratch. Scattered across
two trees and three naming conventions is not a deliverable. A benchmark user handed the staged release
has nothing to score against.

WHAT IT REFUSES TO DO. It does not write a partial bundle quietly. Every component §11 names is either
written or reported as missing, and a missing component is a non-zero exit. The release has already been
taught this lesson four times (§20.1, §22, §23, §24): a stage that produces something and reports success
without saying what it left out is worse than one that fails, because the gap then has to be rediscovered
from the outside. `--allow-missing` exists for assembling a bundle mid-build and names what it skipped,
in the bundle itself as well as on stdout.

IDENTIFIERS. The sample column of every VCF written here is the SYNTHETIC dataset name. IGI-SYN-SEQ-02's
germline baseline is called from a real donor and the baseline VCF carries that donor's cohort accession
as its sample name; copying it through would put an upstream human identifier into the one file a
benchmark user is certain to open. `check_no_upstream_identifiers` asserts it afterwards rather than
trusting that the rename happened, and the forbidden strings come from the design, so adding a dataset
does not quietly opt out.
"""
import argparse
import csv
import datetime
import glob
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Files whose job IS to record provenance, and which therefore name the baselines on purpose. design.yaml
# is the design document; excluding it is not a loophole, it is the distinction between "this release was
# built from IPISRC044" (a fact that belongs in the record) and an upstream accession appearing as a
# sample name in data a tool will read.
PROVENANCE_FILES = ("design.yaml", "MANIFEST.json", "README.md")


def rows(path):
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            yield r


def vcf_info(r):
    """INFO for one somatic record, carrying what §11 asks for."""
    def num(k, default=""):
        v = (r.get(k) or "").strip()
        return v if v else default
    bits = [
        f"EVENT={r['event_id']}",
        f"CLASS={num('class')}",
        f"SUBCLASS={num('subclass')}",
        f"GENE={num('gene')}",
        f"TRANSCRIPT={num('transcript')}",
        f"CSQ={num('consequence')}",
        f"AACHANGE={num('aa_change')}",
        f"CLONE={num('clone')}",
        f"CCF={num('ccf')}",
        f"TIMING={num('timing')}",
        f"HAP={num('haplotype')}",
        f"CLONALITY={num('clonality_tier')}",
        f"EXPVAF={num('expected_vaf_tumor')}",
        f"MULT={num('multiplicity')}",
        f"TUMORCN={num('tumor_cn_at_locus')}",
        f"TPM={num('gene_tpm')}",
        f"EXPTIER={num('expression_tier')}",
        f"AE={num('allelic_expression')}",
        f"BINDTIER={num('binding_tier')}",
        f"RANKEL={num('best_rank_el')}",
        f"ALLELE={num('best_allele')}",
        f"PEPTIDE={num('best_peptide')}",
        f"CONTEXT={num('context')}",
    ]
    if (r.get("flagpost") or "").lower() in ("1", "true", "yes"):
        bits.append("FLAGPOST")
    if (r.get("chr1to6") or "").lower() in ("1", "true", "yes"):
        bits.append("CHR1TO6")
    return ";".join(b for b in bits if not b.endswith("="))


VCF_HEADER_LINES = [
    ('EVENT', '1', 'String', 'designed event id in the catalog tables'),
    ('CLASS', '1', 'String', 'designed event class'),
    ('SUBCLASS', '1', 'String', 'designed event subclass'),
    ('GENE', '1', 'String', 'gene symbol'),
    ('TRANSCRIPT', '1', 'String', 'affected transcript'),
    ('CSQ', '1', 'String', 'consequence'),
    ('AACHANGE', '1', 'String', 'amino-acid change'),
    ('CLONE', '1', 'String', 'clone carrying the variant'),
    ('CCF', '1', 'Float', 'cancer cell fraction of that clone'),
    ('TIMING', '1', 'String', 'pre_cna or post_cna relative to the local copy-number change'),
    ('HAP', '1', 'Integer', 'haplotype the variant sits on'),
    ('CLONALITY', '1', 'String', 'clonality tier'),
    ('EXPVAF', '1', 'Float', 'expected tumour VAF from purity, copy number and multiplicity'),
    ('MULT', '1', 'Integer', 'variant copies at the locus'),
    ('TUMORCN', '1', 'String', 'total tumour copy number at the locus'),
    ('TPM', '1', 'Float', 'designed gene TPM'),
    ('EXPTIER', '1', 'String', 'expression tier'),
    ('AE', '1', 'String', 'allelic expression class'),
    ('BINDTIER', '1', 'String', 'binding tier'),
    ('RANKEL', '1', 'Float', 'best netMHCpan %rank_EL across the patient alleles'),
    ('ALLELE', '1', 'String', 'best-binding allele'),
    ('PEPTIDE', '1', 'String', 'best-binding neopeptide'),
    ('CONTEXT', '1', 'String', 'sequence-context class'),
]


def write_somatic_vcf(dataset, catalog_dir, out, notes):
    """somatic.vcf.gz from the designed SNV/indel table, sorted and bgzip-ready."""
    src = os.path.join(catalog_dir, f"{dataset}.snv_indel.tsv")
    if not os.path.exists(src):
        notes.append(f"somatic.vcf.gz: {os.path.basename(src)} does not exist")
        return None
    recs = []
    for r in rows(src):
        if not (r.get("chrom") and r.get("pos") and r.get("ref") and r.get("alt")):
            continue
        recs.append((r["chrom"], int(r["pos"]), r))
    # Karyotypic order, which is what every consumer expects and what the reference is in.
    order = {f"chr{c}": i for i, c in enumerate(list(range(1, 23)) + ["X", "Y", "M"])}
    recs.sort(key=lambda t: (order.get(t[0], 99), t[1]))

    path = os.path.join(out, "somatic.vcf")
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(f"##fileDate={datetime.date.today().strftime('%Y%m%d')}\n")
        fh.write(f"##source=igi-syn-seq make_truth_bundle.py\n")
        fh.write("##FILTER=<ID=PASS,Description=\"designed, present in the tumour by construction\">\n")
        for k, n, t, d in VCF_HEADER_LINES:
            fh.write(f'##INFO=<ID={k},Number={n},Type={t},Description="{d}">\n')
        fh.write('##INFO=<ID=FLAGPOST,Number=0,Type=Flag,Description="designed positive control">\n')
        fh.write('##INFO=<ID=CHR1TO6,Number=0,Type=Flag,Description="present in the chr1to6 subset">\n')
        # The sample name is the SYNTHETIC dataset, never the upstream donor the baseline came from.
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + dataset + "\n")
        for chrom, pos, r in recs:
            vaf = (r.get("expected_vaf_tumor") or ".").strip() or "."
            fh.write(f"{chrom}\t{pos}\t{r['event_id']}\t{r['ref']}\t{r['alt']}\t.\tPASS\t"
                     f"{vcf_info(r)}\tGT:AF\t0/1:{vaf}\n")
    gz = path + ".gz"
    with open(path, "rb") as fi, gzip.open(gz, "wb") as fo:
        shutil.copyfileobj(fi, fo)
    os.remove(path)
    return len(recs)


def write_cna_seg(dataset, design, arms_bed, out, notes):
    """cna.seg: per clone, per haplotype copy number, from the design rather than from a caller."""
    import yaml
    from igi_catalog.clones import CloneModel
    cfg = yaml.safe_load(open(design))["datasets"][dataset]
    cm = CloneModel(cfg, arms_bed)
    path = os.path.join(out, "cna.seg")
    n = 0
    with open(path, "w") as fh:
        fh.write("clone\tchrom\tstart\tend\tcn_hap0\tcn_hap1\tlabel\tbaseline_ploidy\n")
        base = ",".join(str(b) for b in cm.base)
        for clone, evs in (cfg.get("cna") or {}).items():
            for ev in evs:
                c, s, e = cm.region(ev["region"])
                cn = ev["cn"]
                fh.write(f"{clone}\t{c}\t{s}\t{e}\t{cn[0]}\t{cn[1]}\t"
                         f"{ev.get('label', '')}\t{base}\n")
                n += 1
    if n == 0:
        notes.append("cna.seg: the design declares no copy-number events for this dataset")
    return n


def copy_tables(dataset, catalog_dir, out, notes):
    """The designed-event tables §11 names, each under its specified bundle name."""
    want = {
        "svs.tsv": "sv.tsv",
        "fusions.tsv": "fusions.tsv",
        "viruses.tsv": "viral_integrations.tsv",
        "expressed.tsv": "expression.tumor_bulk.tsv",
        "hla_loh.tsv": "hla_loh.tsv",
        "splice_causal.tsv": "splice_events.tsv",
        "background.tsv": "passenger_background.tsv",
        "germline_spikein.tsv": "germline_spikein.tsv",
    }
    written = {}
    for src_suffix, dst in want.items():
        src = os.path.join(catalog_dir, f"{dataset}.{src_suffix}")
        if not os.path.exists(src):
            notes.append(f"{dst}: {os.path.basename(src)} does not exist")
            continue
        shutil.copy2(src, os.path.join(out, dst))
        with open(src) as fh:
            written[dst] = max(0, sum(1 for _ in fh) - 1)
    return written


def copy_single_cell(dataset, release_dir, out, notes):
    """cells.tsv and clonotypes.tsv, which live with the libraries that used them, not in the catalog."""
    found = {}
    for name in ("cells", "clonotypes"):
        hits = sorted(glob.glob(os.path.join(release_dir, "tenx_gex", dataset,
                                             f"{dataset}_full_{name}.tsv")) +
                      glob.glob(os.path.join(release_dir, "tenx_gex", dataset,
                                             f"{dataset}_*_{name}.tsv")))
        if not hits:
            notes.append(f"{name}.tsv: no roster found under {release_dir}/tenx_gex/{dataset}")
            continue
        # prefer the full release over a chr1to6 roster
        pick = next((h for h in hits if "chr1to6" not in os.path.basename(h)), hits[0])
        shutil.copy2(pick, os.path.join(out, f"{name}.tsv"))
        with open(pick) as fh:
            found[f"{name}.tsv"] = max(0, sum(1 for _ in fh) - 1)
    return found


def copy_manifests(dataset, manifest_dir, out, notes):
    """The LENS manifests, so the bundle is self-contained."""
    if not manifest_dir or not os.path.isdir(manifest_dir):
        notes.append("manifests/: no --manifest-dir given or it does not exist")
        return 0
    dst = os.path.join(out, "manifests")
    os.makedirs(dst, exist_ok=True)
    n = 0
    for f in sorted(glob.glob(os.path.join(manifest_dir, "*manifest*"))):
        shutil.copy2(f, os.path.join(dst, os.path.basename(f)))
        n += 1
    if n == 0:
        notes.append(f"manifests/: nothing matching *manifest* in {manifest_dir}")
    return n


def check_no_upstream_identifiers(dataset, design, out):
    """No file in the bundle may carry another dataset's baseline sample name.

    Asserted rather than assumed. The ds-02 baseline VCF's sample name is its real donor's cohort
    accession, and the rename to the synthetic name is exactly the kind of step that gets done for the
    file someone remembered and missed for the file they did not.
    """
    import yaml
    design_doc = yaml.safe_load(open(design))
    # The policy is the design's own: `publish_haplotype_windows` is owner decision D12, recording which
    # baselines are consented reference material and which are not. HG002 is GIAB public material and is
    # free to name; IPISRC044 is marked false, "public dataset, but not consented for redistribution".
    # Deriving the forbidden set from that flag rather than from a guess here means adding a baseline
    # cannot quietly opt out of the check, and it keeps one decision in one place.
    forbidden = {name for name, ok in (design_doc.get("publish_haplotype_windows") or {}).items()
                 if not ok and isinstance(name, str) and len(name) > 3}
    if not forbidden:
        return [], sorted(forbidden)

    hits = []
    for root, _dirs, files in os.walk(out):
        for f in files:
            if f in PROVENANCE_FILES:
                continue
            p = os.path.join(root, f)
            opener = gzip.open if f.endswith(".gz") else open
            try:
                with opener(p, "rt", errors="replace") as fh:
                    for i, line in enumerate(fh):
                        for bad in forbidden:
                            if bad in line:
                                hits.append((os.path.relpath(p, out), i + 1, bad))
                                break
                        if len(hits) > 50:
                            break
            except OSError:
                continue
    return hits, sorted(forbidden)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--catalog-dir", default="output")
    ap.add_argument("--design", default="design.yaml")
    ap.add_argument("--arms-bed", required=True)
    ap.add_argument("--release-dir", required=True,
                    help="the release root, for the cell and clonotype rosters")
    ap.add_argument("--manifest-dir", default="",
                    help="the staged metadata directory holding the LENS manifests")
    ap.add_argument("--out", required=True, help="the truth/ directory to write")
    ap.add_argument("--allow-missing", action="store_true",
                    help="write the bundle and exit 0 even when components are missing; what was "
                         "skipped is still named, in MANIFEST.json as well as on stdout")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    notes = []

    n_vcf = write_somatic_vcf(a.dataset, a.catalog_dir, a.out, notes)
    n_cna = write_cna_seg(a.dataset, a.design, a.arms_bed, a.out, notes)
    tables = copy_tables(a.dataset, a.catalog_dir, a.out, notes)
    sc = copy_single_cell(a.dataset, a.release_dir, a.out, notes)
    n_man = copy_manifests(a.dataset, a.manifest_dir, a.out, notes)

    shutil.copy2(a.design, os.path.join(a.out, "design.yaml"))
    try:
        ver = subprocess.run(["git", "-C", os.path.dirname(os.path.abspath(__file__)),
                              "describe", "--always", "--dirty"],
                             capture_output=True, text=True).stdout.strip()
    except OSError:
        ver = ""
    if not ver:
        notes.append("generator version: git describe produced nothing")

    import yaml
    seed = yaml.safe_load(open(a.design)).get("seed")
    with open(os.path.join(a.out, "seed.txt"), "w") as fh:
        fh.write(f"{seed}\n")

    hits, forbidden = check_no_upstream_identifiers(a.dataset, a.design, a.out)

    checksums = {}
    for root, _dirs, files in os.walk(a.out):
        for f in sorted(files):
            if f == "MANIFEST.json":
                continue
            p = os.path.join(root, f)
            h = hashlib.sha256()
            with open(p, "rb") as fh:
                for blk in iter(lambda: fh.read(1 << 20), b""):
                    h.update(blk)
            checksums[os.path.relpath(p, a.out)] = h.hexdigest()[:16]

    manifest = {
        "dataset": a.dataset,
        "assembled": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
        "generator_version": ver,
        "design_seed": seed,
        "somatic_vcf_records": n_vcf,
        "cna_segments": n_cna,
        "tables": tables,
        "single_cell": sc,
        "manifests": n_man,
        "missing": notes,
        "upstream_identifier_scan": {
            "forbidden_strings": forbidden,
            "hits": [{"file": f, "line": l, "string": b} for f, l, b in hits],
        },
        "files": checksums,
    }
    json.dump(manifest, open(os.path.join(a.out, "MANIFEST.json"), "w"), indent=2, sort_keys=True)

    print(f"[{a.dataset}] truth bundle -> {a.out}")
    print(f"  somatic.vcf.gz   {n_vcf if n_vcf is not None else 'NOT WRITTEN'} records")
    print(f"  cna.seg          {n_cna} segments")
    for k, v in sorted(tables.items()):
        print(f"  {k:<32} {v} rows")
    for k, v in sorted(sc.items()):
        print(f"  {k:<32} {v} rows")
    print(f"  manifests/       {n_man} file(s)")
    print(f"  {len(checksums)} file(s), checksummed in MANIFEST.json")

    if hits:
        print(f"\nUPSTREAM IDENTIFIER IN THE BUNDLE: {len(hits)} occurrence(s)", flush=True)
        for f, l, b in hits[:10]:
            print(f"  {f}:{l} contains {b!r}")
        print("a synthetic release must not carry the identifier of the donor its baseline came from")
        return 2

    if notes:
        print(f"\n{len(notes)} component(s) of design §11 could not be assembled:", flush=True)
        for n in notes:
            print(f"  {n}")
        if not a.allow_missing:
            print("this bundle is incomplete; that is a failure, not a pass "
                  "(pass --allow-missing to accept it mid-build)")
            return 1
        print("accepted under --allow-missing; the gaps are recorded in MANIFEST.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
