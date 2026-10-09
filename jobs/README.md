# Build jobs

The scripts that actually produce IGI-SYN-SEQ. `docs/build-reference.md` explains what each stage does
and what to edit to change a behaviour; this file is how to run them.

These are SLURM array jobs over the Python drivers in `catalog/`. The design specification describes a
Nextflow pipeline — that is the intended end state and does not exist yet.

## Setup

```bash
cp jobs/env.example.sh jobs/env.sh     # untracked; edit the paths
export IGI_ENV=$PWD/jobs/env.sh        # sbatch propagates the environment, which is how jobs find it
```

No tracked file contains an absolute path. `jobs/env.sh` and `config/sites/*.paths.yaml` are the only
two places a path belongs.

## Order

```bash
./jobs/submit.sh 01_catalog.sbatch                 # designed catalog      ~1 h x 2
./jobs/submit.sh 02_background.sbatch              # passenger background  ~2 min x 2
./jobs/submit.sh 10_wes.sbatch                     # exomes, 96 tasks      <16 h
./jobs/submit.sh 11_merge.sbatch                   # merge + prove name uniqueness
./jobs/submit.sh 12_rna.sbatch                     # bulk RNA, 16 GB      ~2 h x 2
./jobs/submit.sh 13_align.sbatch                   # bwa mem for acceptance
./jobs/submit.sh 14_accept.sbatch                  # acceptance, DNA
./jobs/submit.sh 15_accept_rna.sbatch              # acceptance, RNA + class presence
./jobs/submit.sh 20_wgs_illumina.sbatch            # short-read WGS, 96 tasks
./jobs/submit.sh 21_wgs_pacbio.sbatch              # HiFi WGS, 96 tasks
./jobs/submit.sh 22_accept_wgs.sbatch              # acceptance, short-read WGS
./jobs/submit.sh 30_kinnex_bulk.sbatch             # Kinnex bulk RNA, pre- and post-skera BAMs
./jobs/submit.sh 31_kinnex_sc.sbatch               # Kinnex single cell + cell roster + clonotypes
./jobs/submit.sh 32_tenx_gex.sbatch                # 10x 5' v2 gene expression (Cell Ranger naming)
./jobs/submit.sh 33_tenx_tcr.sbatch                # 10x 5' V(D)J (TCR)
./jobs/submit.sh 40_ont_wgs.sbatch                 # ONT WGS, 96 tasks
./jobs/submit.sh 41_ont_bulk_rna.sbatch            # ONT bulk cDNA
./jobs/submit.sh 42_ont_sc_rna.sbatch              # ONT single-cell cDNA
```

Then the stages that turn per-chromosome output into deliverables, derive the chr1to6 release, and stage
the result. These were missing from this list for a while, which made the post-rebuild sequence something
you had to reconstruct:

```bash
./jobs/submit.sh 50_merge_longread.sbatch          # merge per-chromosome long reads, 12 tasks
./jobs/submit.sh 72_pacbio_pbi.sbatch              # .pbi per merged PacBio library (samtools merge omits it)
./jobs/submit.sh 90_subset_map.sbatch              # record -> chromosome map, 2 tasks; read by every filter
./jobs/submit.sh 90_subset_chr1to6.sbatch          # chr1to6 for the six whole-library arms, 12 tasks
./jobs/submit.sh 91_merge_chr1to6.sbatch           # merge the chr1-6 shards of the four per-chromosome arms
./jobs/70_stage_inputs.sh                          # symlink the release into the RAFT workspace (not sbatch)
python3 catalog/make_lens_manifest.py --release both
python3 catalog/qa_release.py --release <release> --out qa.json --samtools "<cmd>"
python3 catalog/tests/test_guards.py               # the delivery guards; no cluster, no data, 0.2 s
```

`90_subset_map` must be re-run before `90_subset_chr1to6` whenever an arm has been rebuilt, because the
map is derived from that dataset's regenerated transcript manifest. `70_stage_inputs.sh` holds back any
chr1to6 deliverable older than what it was derived from, so a stale subset cannot be staged as current --
see build-reference §15.3 for the incident that made that necessary.

## Re-running one arm

Stages `10`-`21` and the Kinnex/10x/ONT builders rebuild unconditionally; none of them skips existing
output. So a rebuild is a resubmission of the same array, and the seeds derive from
`SLURM_ARRAY_TASK_ID`, which means the same task ids reproduce the same library. The 96-task arrays index
`i = TASK-1`, `c = CHROMS_ALL[i%24]`, `j = i/24`, `l = LB[j%2]`, `d = DATASETS[j/2]`, so **ds-01 is tasks
1-48 and ds-02 is 49-96**; the 2-task arrays use task 1 for ds-01 and task 2 for ds-02.

After rebuilding any arm, the chr1to6 subsets derived from it are stale and the merges above have to be
re-run.

## Shared shell helpers

`jobs/lib.sh` is sourced next to `$IGI_ENV` and holds site-independent helpers -- currently the GNU sort
buffer, which has to be derived from the step's allocation because sort sizes its default buffer from the
machine's physical memory and never reads the cgroup limit (build-reference §16).

`02` must follow `01`: the background model reserves the positions the designed catalog occupies so the
two can never collide. `10`, `12`, `20` and `21` are independent of each other and can run together; all
four depend on `01` and `02`. Chain them with `--dependency=afterok:<jobid>`.

Stages `10`-`21` are per-chromosome, so a partial failure is re-run by resubmitting the same array with
`--array=<failed ids>`. They do **not** skip work that already exists -- a re-run rebuilds and overwrites,
which is what makes a resubmission clean rather than leaving a mixture of old and new output. An earlier
version of this file claimed they skip; they do not.

## Measurement jobs

`measure/` fits the constants the builders use. These are run once and re-run only to recalibrate; their
outputs are committed under `catalog/resources/`.

```bash
./jobs/submit.sh measure/art_profile.sbatch        # Illumina per-cycle quality profile
./jobs/submit.sh measure/gc_bias.sbatch            # capture efficiency vs interval GC
./jobs/submit.sh measure/longread_measure.sbatch   # PacBio and Kinnex length/accuracy
```

Two more are run directly rather than through sbatch, because each needs a real BAM streamed through
`samtools view` and takes minutes rather than hours:

```bash
# Kinnex MAS adapters and array structure, decoded out of skera's own tags
samtools view segmented.bam > seg.sam
python3 jobs/measure/fit_mas_adapters.py --sam seg.sam \
  --out-fasta catalog/resources/kinnex_mas8_adapters.fasta \
  --out-profile catalog/resources/kinnex_mas8_profile.json

# Kinnex size selection: real segment lengths against the simulated molecule lengths
awk -F'\t' '$1!~/^@/{print length($10)}' seg.sam > seglen.txt
python3 jobs/measure/fit_kinnex_size_selection.py --real-lengths seglen.txt \
  --manifest <dataset>_full_rna_transcripts.tsv \
  --out catalog/resources/kinnex_size_selection.json
```

`measure/fit_gc_bias.py` and `measure/measure_longread.py` can also be run directly, and the second is
how a simulated long-read output is checked against its target.

Every single-cell assay (`31`, `32`, and the TCR and ONT single-cell arms when written) must run with the
**same seed, 9100**. The cell roster and the molecule pool are derived from it, so a different seed gives a
different roster and a barcode stops meaning one cell across assays -- the one property the single-cell
design exists to provide.

## Two things that will bite

**Read-name uniqueness is proved, not assumed.** `11_merge.sbatch` sorts the merged read names and fails
the job if the count of distinct names differs from the count of reads. An earlier build had 667,090 of
1,000,832 record ids duplicated because ids restarted at 1 on each chromosome, and the integrity check of
the day passed anyway — it verified that names *resolve*, not that they are *unique*. Do not remove that
check to save wall time.

**`zcat | head` fails under `set -o pipefail`.** `head` exits once satisfied, `zcat` takes SIGPIPE, and
pipefail turns that into a failed pipeline. This killed two jobs. Where subsampling is unavoidable the
scripts turn pipefail off around it and verify the record count afterwards; elsewhere they use
`tail -n +2`, which reads to EOF.
