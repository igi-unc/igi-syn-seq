#!/bin/bash
#SBATCH -J igi_accept_wgs
#SBATCH -c 16 --mem=64G -t 24:00:00
#SBATCH --array=1-2
# Acceptance for the short-read WGS arm. The depth and allele-fraction checks need an aligned BAM, so this
# aligns the merged library first.
#
# --assay wgs changes three things, each because the assay genuinely differs rather than to make the check
# pass: depth is measured over non-N reference instead of capture intervals (chr21 is 18.6 % N, and an
# N-straddling window otherwise reads as a deletion); the GC curve is not applied, since it is capture
# efficiency and the builder uses a flat one; and every designed breakpoint on the chromosome counts as
# placeable, which is a stronger expectation than the exome's, not a weaker one.
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
if [ ! -s "$R1" ]; then
  echo "  merge job output absent; merging here instead"
  for m in 1 2; do cat $O/$d/${d}_chr*_tumor_wgs_R${m}.fastq.gz > $M/$d/${d}_tumor_wgs_R${m}.fastq.gz; done
fi
B=$A/$d/${d}_tumor.bam
if [ ! -s "$B" ]; then
  singularity exec -B "$IGI_BIND" "$IGI_BWA_IMG" bash -c "
    bwa mem -t 16 -R '@RG\tID:${d}_tumor_wgs\tSM:${d}_tumor\tPL:ILLUMINA\tLB:${d}_tumor_wgs' \
      $IGI_FASTA_NOALT $R1 $R2 2> $A/$d/bwa.err \
      | samtools sort -@4 -m 3G -o $B -" || { tail -20 "$A/$d/bwa.err"; exit 1; }
  singularity exec -B "$IGI_BIND" "$IGI_BWA_IMG" samtools index "$B"
fi
cd "$IGI_REPO/catalog"
for c in chr1 chr2 chr5 chr6 chr8 chr10 chr13 chr17 chr21 chrX; do
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
empty=0
for r in "$R/$d"/*.json; do
  [ -s "$r" ] || continue
  n=$(python3 -c "import json;print(len(json.load(open('$r'))['checks']))" 2>/dev/null || echo 0)
  [ "$n" -eq 0 ] && { echo "EMPTY REPORT: $r ran no checks"; empty=$((empty+1)); }
done
[ "$empty" -eq 0 ] || { echo "$empty report(s) ran no checks; that is a failure, not a pass"; exit 1; }
echo "DONE $d $(date)"
