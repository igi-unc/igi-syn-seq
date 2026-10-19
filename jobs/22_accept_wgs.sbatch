#!/bin/bash
#SBATCH -J igi_accept_wgs
#SBATCH -c 16 --mem=64G -t 8:00:00
#SBATCH --array=1-2
# WALL CLOCK: this stage no longer aligns anything. Its two full runs took 21:09 and 21:42 when the
# tumour alignment was inside it; with that moved to 22_align_wgs_tumor.sbatch, the measured cost is the
# checks alone -- 40:34 for ds-01 and 39:39/40:00 for ds-02 on 2026-10-16. The limit below is 8 h, a
# ~12x margin on those, rather than the 40 h it needed as a builder.
#
# RETRY: when a check FAILS this exits IGI_EXIT_VERDICT (see lib.sh), not 1, so 99_drive does not retry
# it. A verdict is deterministic; on 2026-10-16 the driver re-ran this stage three times across two
# hours to be told the same thing about ds-02 chr5 three times.
#
# STALENESS: the guard below compares the BAM against the library it came from and FAILS rather than
# rebuilding. A rebuild that finds stale work in place and skips is how a release gets graded on data it
# no longer contains -- the long-read merge would have kept pre-rebuild libraries, and this stage once
# scored a 2026-10-05 alignment of a 2026-10-06 library in 35 minutes instead of 21 hours.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
source "$(dirname "$IGI_ENV")/lib.sh"
O=$IGI_RELEASE/wgs; M=$IGI_RELEASE/wgs_merged
A=$IGI_RELEASE/wgs_bams; R=$IGI_RELEASE/wgs_acceptance
d=${DATASETS[$((SLURM_ARRAY_TASK_ID-1))]}
mkdir -p "$R/$d"
# The tumour BAM is 22_align_wgs_tumor.sbatch's output, not this stage's. It used to be built right here,
# which cost 15 h of this stage's 21 h and, worse, left the grader unable to say anything about the
# provenance of what it graded: a stage that would rebuild a stale BAM can never report one. Now it can,
# and does -- a BAM older than its library is a FAILURE here rather than something quietly refreshed.
R1=$M/$d/${d}_tumor_wgs_R1.fastq.gz
B=$A/$d/${d}_tumor.bam
[ -s "$R1" ] || { echo "missing merged tumour library $R1; run 50_merge_longread.sbatch first"; exit 1; }
[ -s "$B" ]  || { echo "missing tumour BAM $B; run 22_align_wgs_tumor.sbatch first"; exit 1; }
[ -s "$B.bai" ] || { echo "missing index $B.bai; run 22_align_wgs_tumor.sbatch first"; exit 1; }
if [ "$R1" -nt "$B" ]; then
  echo "STALE BAM: $B is older than the library it came from ($R1)."
  echo "  Grading this would describe a library the release no longer contains."
  echo "  Re-run 22_align_wgs_tumor.sbatch for $d, then this stage."
  exit 1
fi
cd "$IGI_REPO/catalog"
# One list, used BOTH to run the loop and to name the reports the guard expects, so the two can
# never drift apart.
ACCEPT_CHROMS=(chr1 chr2 chr5 chr6 chr8 chr10 chr13 chr17 chr21 chrX)
for c in "${ACCEPT_CHROMS[@]}"; do
  # The junction file is per chromosome and named with the _wgs_ infix. It is only passed when present:
  # builds made before run_wgs.py emitted one have no junction truth, and the check must then be absent
  # rather than silently reading nothing.
  j=$O/$d/${d}_${c}_tumor_wgs_junctions.tsv
  jopt=""; [ -s "$j" ] && jopt="--junctions $j"
  python3 acceptance.py --dataset "$d" --chrom $c --catalog-dir output --design design.yaml \
    --assay wgs --bam "$B" --reference "$IGI_FASTA_NOALT" \
    --normal-bam "$A/$d/${d}_normal.bam" \
    --designed-normal-depth 30 --designed-tumor-depth 100 \
    `# 23_align_wgs_normal.sbatch writes this BAM, and must have run before this stage. Until it did,` \
    `# wgs_merged held the normal's FASTQs and wgs_bams held only the tumour, so acceptance never` \
    `# opened a normal library at all -- which is how the 20.7x-against-30x defect survived a passing` \
    `# 40/40. A missing normal BAM is reported as a FAILURE here, not a skip, for the same reason.` \
    --arms-bed "$IGI_ARMS_BED" --paths "$IGI_PATHS" $jopt \
    --work "$R/$d/$c" --report "$R/$d/${d}_${c}.json" 2>&1 | sed "s/^/[$d $c] /" || true
done
# A zero-check run exits 2 and prints "no checks ran". The per-chromosome loop uses `|| true` so one bad
# chromosome does not abandon the rest, which means that exit code is swallowed -- so look for it in the
# reports instead, where it is a report with no checks at all.
# This scans $R, the acceptance reports -- not $O, which holds the builder's 48 per-chromosome meta files.
# Pointed at $O it found no "checks" key in any of them and would have failed every single run.
# Check the reports THIS STAGE OWNS, named from the same list the loop ran, rather than globbing the
# directory. A glob swept every .json in it, and anything else writing into the same directory is swept with it:
# its ${d}_class_presence.json is a summary document with no `checks` key at all, so the glob called a
# correct file an empty report and failed a stage whose own 44 checks had every one passed. The comment
# above records the SAME mistake once already -- the glob pointed at the builder's meta files -- which
# is the argument for naming the expected set rather than trusting what is in the directory.
#
# Naming it is also strictly stronger: a glob over what EXISTS cannot notice a report that was never
# written, and a chromosome whose acceptance died before writing anything is exactly the case the
# `|| true` above makes invisible.
expect=()
for c in "${ACCEPT_CHROMS[@]}"; do expect+=("$R/$d/${d}_${c}.json"); done
# The design against what the generator WROTE, with no aligner in the path (section 23). This is the
# check that can call a release defect a release defect: `allele fractions by tier` above compares the
# design against a pileup of an alignment, so it also answers whether a consumer can recover what was
# written, and on ds-02 it went red for that reason on a correct release. One invocation per dataset
# rather than per chromosome, because the record counts are one table for the whole library. No
# `|| true` here: set -e is wanted, and truth_af.py exits 1 when the generator and the design disagree
# and 2 when it measured nothing. Its report is deliberately NOT in the list below -- it has no `checks`
# key, and sweeping a document whose schema belongs to another stage is the mistake section 22 records.
python3 truth_af.py --dataset "$d" --catalog-dir output --design design.yaml \
  --records "$M/$d/${d}_tumor_wgs_records.tsv.gz" --readmap "$M/$d/${d}_tumor_wgs_readmap.tsv.gz" \
  --counts "$R/$d/${d}_record_counts.tsv" --report "$R/$d/${d}_truth_af.json"
[ -s "$R/$d/${d}_truth_af.json" ] || { echo "MISSING REPORT: truth_af wrote nothing"; exit 1; }

# Three different things can be wrong with a report, and until 2026-10-11 this guard only looked at two
# of them: the report is absent, it ran no checks, or it ran them and they FAILED. The third is the only
# one the stage exists to find, and it was the one not asked about -- the WGS acceptance exited 0 and
# SLURM said COMPLETED with two failed allele-fraction checks on ds-02 sitting in the reports. The
# `|| true` above swallows acceptance.py's exit 1 for the same reason it swallows exit 2, so a failed
# check has to be read back out of the report, exactly like a zero-check run.
bad=0; nfail=0; flagged=()
for r in "${expect[@]}"; do
  if [ ! -s "$r" ]; then echo "MISSING REPORT: $r was never written"; bad=$((bad+1)); continue; fi
  read -r n f < <(python3 -c \
    "import json;d=json.load(open('$r'));print(len(d['checks']),d.get('n_fail',0))" 2>/dev/null \
    || echo "0 0")
  if [ "${n:-0}" -eq 0 ]; then echo "EMPTY REPORT: $r ran no checks"; bad=$((bad+1)); continue; fi
  if [ "${f:-0}" -gt 0 ]; then nfail=$((nfail+f)); flagged+=("$(basename "$r") x$f"); fi
done
[ "$bad" -eq 0 ] || { echo "$bad of ${#expect[@]} report(s) missing or ran no checks; that is a failure, not a pass"; exit 1; }
if [ "$nfail" -gt 0 ]; then
  verdict_fail "FAILED CHECKS: $nfail check(s) failed across ${#flagged[@]} report(s): ${flagged[*]}"
fi
echo "DONE $d: ${#expect[@]} reports, every one with checks, 0 failed $(date)"
