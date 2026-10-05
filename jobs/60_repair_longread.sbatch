#!/bin/bash
#SBATCH -J igi_repair
#SBATCH -c 12 --mem=24G -t 12:00:00
#SBATCH --array=1-15
# Repair delivered long-read libraries in place. Four review findings, all in the packaging rather than in
# the reads, so none of them needs the simulator run again:
#
#   tasks 1-7    ONT: strip Badread's description from every FASTQ header into a sidecar map. The header
#                named the source reference, strand and position for a genomic read and the source
#                molecule for an RNA one, so every read carried its own answer.
#   tasks 8-11   HiFi WGS: rewrite the `samtools import` BAMs as Revio BAMs -- Revio read names,
#                zm/np/ec/rq/qs/qe, an @RG with PU/PM and an eight-hex-digit ID, and a .pbi.
#   tasks 12-15  Kinnex: redraw per-read np/ec and add the @RG PU that skera reads to name segments, then
#                re-run skera so the segment names carry the movie instead of starting with a slash.
#
# Memory is 24G because every task is a streaming pass; the 320G the builders asked for was never needed
# and only kept them queued behind the QOS cap.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
R=$IGI_WORK/release
cd "$IGI_REPO/catalog"
PBI=$(python3 -c "import yaml,sys; print(yaml.safe_load(open('$IGI_PATHS')).get('pbindex_cmd',''))")
SKERA=$(python3 -c "import yaml,sys; print(yaml.safe_load(open('$IGI_PATHS')).get('skera_cmd',''))")
ADAPT=$IGI_REPO/catalog/resources/kinnex_mas8_adapters.fasta
t=$SLURM_ARRAY_TASK_ID

ont=(
  "$R/ont_wgs_merged/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01_tumor_ont_wgs.fastq.gz"
  "$R/ont_wgs_merged/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01_normal_ont_wgs.fastq.gz"
  "$R/ont_wgs_merged/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02_tumor_ont_wgs.fastq.gz"
  "$R/ont_wgs_merged/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02_normal_ont_wgs.fastq.gz"
  "$R/ont_sc_rna/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01_full_ont_sc_rna.fastq.gz"
  "$R/ont_sc_rna/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02_full_ont_sc_rna.fastq.gz"
  "$R/ont_bulk_rna/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02_full_ont_bulk_rna.fastq.gz"
)
hifi=(
  "IGI-SYN-SEQ-01 tumor  $R/pacbio_merged/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01_tumor_hifi.bam"
  "IGI-SYN-SEQ-01 normal $R/pacbio_merged/IGI-SYN-SEQ-01/IGI-SYN-SEQ-01_normal_hifi.bam"
  "IGI-SYN-SEQ-02 tumor  $R/pacbio_merged/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02_tumor_hifi.bam"
  "IGI-SYN-SEQ-02 normal $R/pacbio_merged/IGI-SYN-SEQ-02/IGI-SYN-SEQ-02_normal_hifi.bam"
)
kinnex=(
  "IGI-SYN-SEQ-01 kinnex_bulk $R/kinnex_bulk/IGI-SYN-SEQ-01"
  "IGI-SYN-SEQ-02 kinnex_bulk $R/kinnex_bulk/IGI-SYN-SEQ-02"
  "IGI-SYN-SEQ-01 kinnex_sc   $R/kinnex_sc/IGI-SYN-SEQ-01"
  "IGI-SYN-SEQ-02 kinnex_sc   $R/kinnex_sc/IGI-SYN-SEQ-02"
)

if [ "$t" -le 7 ]; then
  python3 -u repair_longread.py --what ont --path "${ont[$((t-1))]}" --threads 10

elif [ "$t" -le 11 ]; then
  read -r ds lib bam <<<"${hifi[$((t-8))]}"
  python3 -u repair_longread.py --what hifi --path "$bam" \
    --sample "$ds" --library "${ds}_${lib}_hifi" --seed 5000 --pbindex-cmd "$PBI"

else
  read -r ds lib dir <<<"${kinnex[$((t-12))]}"
  pre="$dir/${ds}_full_${lib}_hifi_reads.bam"
  seg="$dir/${ds}_full_${lib}_segmented.bam"
  python3 -u repair_longread.py --what kinnex --path "$pre" \
    --sample "$ds" --library "${ds}_${lib}" --seed 9100 --pbindex-cmd "$PBI"
  # skera rebuilds the segmented BAM from the repaired array reads. Written beside the old one and only
  # swapped in once it exists, is non-trivial and carries a movie in its segment names: a half-written
  # 29 GB segmented BAM that has already replaced its source cannot be recovered without re-running the
  # whole Kinnex arm. The temporary name keeps the .bam suffix because skera derives its own
  # `.non_passing.bam` from it and rejects an output that is not a BAM.
  new="${seg%.bam}.rebuild.bam"
  rm -f "$new" "$new.pbi" "${new%.bam}.non_passing.bam"
  eval "${SKERA/\{args\}/split $pre $ADAPT $new}"
  read -r n first <<<"$(python3 -c 'import pysam,sys; f=pysam.AlignmentFile(sys.argv[1],"rb",check_sq=False); n=0; fi=""
for r in f.fetch(until_eof=True):
    fi = fi or r.query_name
    n += 1
print(n, fi)' "$new")"
  [ "${n:-0}" -gt 1000 ] || { echo "skera produced only ${n:-0} segments; refusing to replace $seg"; exit 1; }
  case "$first" in
    m*/*/ccs/*) : ;;
    *) echo "segment name '$first' still has no movie; refusing to replace"; exit 1 ;;
  esac
  rm -f "$seg" "$seg.pbi" "${seg%.bam}.non_passing.bam"
  mv "$new" "$seg"
  if [ -e "$new.pbi" ]; then mv "$new.pbi" "$seg.pbi"; fi
  if [ -e "${new%.bam}.non_passing.bam" ]; then mv "${new%.bam}.non_passing.bam" "${seg%.bam}.non_passing.bam"; fi
  echo "segmented rebuilt: $n segments, first name $first"
fi
echo "DONE task $t $(date)"
