#!/bin/bash
#SBATCH -J igi_accept_wgs
#SBATCH -c 16 --mem=64G -t 40:00:00
#SBATCH --array=1-2
# WALL CLOCK: the two full runs of this stage took 21:09 and 21:42 against a 24 h limit -- an 11 %
# margin. The junction fix restores the coverage of every junction-bearing interval (each was
# contributing ~1/4,000th of its sequence), so the merged WGS libraries this now aligns are ~2 %
# larger, eating a third of what was left. A TIMEOUT here is expensive in a way an over-request is
# not: the driver retries a timeout at the SAME limit (memory escalates only on OUT_OF_MEMORY), so
# it would burn ~21 h twice more and then give up. 40 h costs nothing if the job finishes in 22.
# Acceptance for the short-read WGS arm. The depth and allele-fraction checks need an aligned BAM, so this
# aligns the merged library first.
#
# --assay wgs changes three things, each because the assay genuinely differs rather than to make the check
# pass: depth is measured over non-N reference instead of capture intervals (chr21 is 18.6 % N, and an
# N-straddling window otherwise reads as a deletion); the GC curve is not applied, since it is capture
# efficiency and the builder uses a flat one; and every designed breakpoint on the chromosome counts as
# placeable, which is a stronger expectation than the exome's, not a weaker one.
# An output that merely EXISTS is not an output that is current. Guarding with `[ ! -s "$x" ]` skipped
# work whenever a previous run had left something behind, and that silently graded or shipped stale data
# twice in one release: the long-read merge would have kept pre-rebuild libraries, and the WGS acceptance
# scored a 2026-10-05 alignment of a 2026-10-06 library in 35 minutes instead of 21 hours. Compare mtimes.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
O=$IGI_WORK/release/wgs; M=$IGI_WORK/release/wgs_merged
A=$IGI_WORK/release/wgs_bams; R=$IGI_WORK/release/wgs_acceptance
d=${DATASETS[$((SLURM_ARRAY_TASK_ID-1))]}
mkdir -p "$M/$d" "$A/$d" "$R/$d"
# Use 50_merge_longread.sbatch's output rather than merging again. Both scripts wrote the same libraries
# into this directory under different names and duplicated 234 GB: this one as ${d}_tumor_R1.fastq.gz and
# the merge job as ${d}_tumor_wgs_R1.fastq.gz. Fall back to merging only if the merge job has not run.
R1=$M/$d/${d}_tumor_wgs_R1.fastq.gz
R2=$M/$d/${d}_tumor_wgs_R2.fastq.gz
newest=$(ls -t $O/$d/${d}_chr*_tumor_wgs_R1.fastq.gz 2>/dev/null | head -1)
if [ ! -s "$R1" ] || { [ -n "$newest" ] && [ "$newest" -nt "$R1" ]; }; then
  echo "  merge job output absent or older than its chromosomes; merging here instead"
  for m in 1 2; do cat $O/$d/${d}_chr*_tumor_wgs_R${m}.fastq.gz > $M/$d/${d}_tumor_wgs_R${m}.fastq.gz; done
fi
B=$A/$d/${d}_tumor.bam
if [ ! -s "$B" ] || [ "$R1" -nt "$B" ]; then
  rm -f "$B" "$B.bai"
  singularity exec -B "$IGI_BIND" "$IGI_BWA_IMG" bash -c "
    bwa mem -t 16 -R '@RG\tID:${d}_tumor_wgs\tSM:${d}_tumor\tPL:ILLUMINA\tLB:${d}_tumor_wgs' \
      $IGI_FASTA_NOALT $R1 $R2 2> $A/$d/bwa.err \
      | samtools sort -@4 -m 3G -o $B -" || { tail -20 "$A/$d/bwa.err"; exit 1; }
  singularity exec -B "$IGI_BIND" "$IGI_BWA_IMG" samtools index "$B"
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
bad=0
for r in "${expect[@]}"; do
  if [ ! -s "$r" ]; then echo "MISSING REPORT: $r was never written"; bad=$((bad+1)); continue; fi
  n=$(python3 -c "import json;print(len(json.load(open('$r'))['checks']))" 2>/dev/null || echo 0)
  [ "$n" -eq 0 ] && { echo "EMPTY REPORT: $r ran no checks"; bad=$((bad+1)); }
done
[ "$bad" -eq 0 ] || { echo "$bad of ${#expect[@]} report(s) missing or ran no checks; that is a failure, not a pass"; exit 1; }
echo "DONE $d: ${#expect[@]} reports, every one with checks $(date)"
