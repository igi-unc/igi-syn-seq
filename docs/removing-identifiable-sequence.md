# Removing identifiable sequence from the repository and its history

Owner decision D12. Reconstructed haplotype windows disclose a real individual's genotype over the
window. HG002 is consented GIAB reference material and its windows are published; IPISRC044 is public
data but not consented for redistribution of reconstructed haplotypes, so its windows are withheld.

This is written down because it is a procedure, not a one-off: a recipe that is not in the repository is
not a recipe, and the same steps apply if another non-consented baseline is ever added.

## What the rule is

`publish_haplotype_windows` in `catalog/design.yaml` is keyed by **baseline**, not by dataset, so the rule
follows the individual rather than the label a dataset happens to carry.

Four columns of `<dataset>.snv_indel.tsv` are affected, along with the corresponding lines of
`<dataset>.snv_indel.diffcards.txt`:

| Column | Why it is identifiable |
|---|---|
| `wt_hap_window` | the individual's own haplotype over a 41 bp window |
| `mut_hap_window` | the same window with the designed somatic change applied |
| `wt_protein_window` | the individual's own protein over the window |
| `mut_protein_window` | built on the individual's germline background |

Where the flag is false these read `withheld`. The window **coordinates** and the designed change stay,
so the truth table remains complete as a specification; only the individual's sequence goes. Whoever runs
the recipe regenerates the columns locally from their own copy of the baseline VCF, so nothing is lost.

## Part 1: stop publishing it (done)

Set the baseline to `false` in `publish_haplotype_windows` and regenerate the catalog. To strip tables
that already exist without a full regeneration, replace the four column values with `withheld` and drop
diff-card lines beginning `ref`, `hap0`, `hap1`, `WT` or `MUT`.

Verify:

    # every value withheld, no sequence lines left
    python3 - <<'PY'
    import csv, collections, re, os
    rows = list(csv.DictReader(open("catalog/output/<dataset>.snv_indel.tsv"), delimiter="\t"))
    for c in ("wt_hap_window", "mut_hap_window", "wt_protein_window", "mut_protein_window"):
        assert not [r[c] for r in rows if r[c] not in ("withheld", "")], c
    d = "catalog/output/<dataset>.snv_indel.diffcards.txt"
    assert not [l for l in open(d) if re.match(r"\s+(ref|hap0|hap1|WT|MUT)\s", l)]
    PY

## Part 2: remove it from history (requires a maintainer)

The data was committed before the rule existed and has been pushed, so it is reachable from earlier
commits and from GitHub. Removing it is a history rewrite and a force-push, which is why it is not
automated here: it changes commits other people may hold.

**Before anything, back up all refs.** A bundle plus a filesystem copy of the repository:

    BK=/path/to/backup-$(date +%Y%m%d-%H%M%S); mkdir -p "$BK"
    git bundle create "$BK/all-refs.bundle" --all
    cp -a <repo> "$BK/repo-copy"
    git rev-parse HEAD origin/main > "$BK/HEAD-before-rewrite.txt"

Then rewrite the affected blobs. `git filter-repo` is the supported tool; it is a single Python file and
is not bundled with git. Rewrite **contents** rather than deleting the files, so the non-identifiable
columns survive in history:

    git filter-repo --blob-callback '
      # for each historical version of the affected files, replace the four column values with
      # "withheld" and drop diff-card sequence lines, exactly as Part 1 does to the working tree
    '

Afterwards, confirm nothing remains reachable:

    for c in $(git log --format=%H --all -- catalog/output/<dataset>.snv_indel.tsv); do
      git show "$c:catalog/output/<dataset>.snv_indel.tsv" |
        awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="wt_hap_window")k=i}
                    NR>1 && $k!="withheld" && $k!=""{n++} END{if(n)print "'"$c"'", n}'
    done

Empty output means clean.

Then `git push --force-with-lease`.

## Part 3: what a force-push does not fix

State these rather than assume the rewrite is the end of it:

- **Anyone who already cloned or forked keeps the data**, and their history diverges from the rewritten
  one. There is no way to reach those copies.
- **GitHub keeps unreferenced objects addressable by commit SHA** until it garbage-collects. A rewrite
  alone does not purge them; ask GitHub Support to run a GC on the repository, citing the commits.
- Open pull requests and forks can hold references that keep the old objects alive; close or delete them
  first.
- Anyone who needs the withheld columns regenerates them locally. That is the intended workflow, not a
  degraded one.
