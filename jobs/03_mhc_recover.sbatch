#!/bin/bash
#SBATCH -J igi_mhc_recover
#SBATCH -c 24 --mem=96G -t 2-00:00:00
# Owner decision D3, recover step. The IPISRC044 blood normal was aligned alt-aware, so its HLA reads went to
# the 525 HLA-* contigs and the primary chr6 MHC was left at 7-11x -- too shallow to genotype the class I
# loci, which is where the designed HLA LOH has to be visible. Realign the reads that could belong in the MHC
# to a reference with no alt and no HLA contigs, so they have nowhere to go but chr6, then recall.
#
# Measured outcome on IPISRC044 (job 11552481, 2026-09-28, 1h26m):
#   2,580 of 3,366 contigs kept; 547,704 pairs and 8,680 singletons extracted
#   HLA-A depth 7-11x -> 40.8x, HLA-B -> 48.9x
#   7,532 PASS records in chr6:28-34 Mb, 3,991 heterozygous
#   per locus: HLA-A 32 PASS / 0 het, HLA-B 93/64, HLA-C 105/100
#
# HLA-A has no heterozygous calls because IPISRC044 is homozygous there, A*01:01/A*01:01, which design.yaml
# records. That is the correct answer, not a failure of the recovery, and it means HLA-A LOH cannot be
# detected by allelic imbalance in IGI-SYN-SEQ-02. IGI-SYN-SEQ-01 is heterozygous at all three class I loci
# and covers that case.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"
W=$IGI_WORK/hla_recover; mkdir -p "$W"
DV=$IGI_REF/imgs/google-deepvariant-1.8.0.img
FA=$IGI_FASTA_PRIMARY
PRIM=$IGI_FASTA_NOALT
ln -sf "$IGI_RECOVER_BAM" "$W/n.bam"; ln -sf "$IGI_RECOVER_BAM.bai" "$W/n.bam.bai"

# 1. Primary-only reference: keep the primary chromosomes and the unplaced/random/decoy contigs, drop every
#    *_alt and HLA-* contig so their reads have nowhere else to go.
if [ ! -s "$PRIM.bwt" ]; then
  awk '$1 !~ /_alt$/ && $1 !~ /^HLA-/ {print $1}' "$FA.fai" > "$W/keep.txt"
  echo "keeping $(wc -l < "$W/keep.txt") contigs of $(wc -l < "$FA.fai")"
  singularity exec -B "$IGI_BIND" "$IGI_BWA_IMG" bash -c "xargs samtools faidx $FA < $W/keep.txt > $PRIM"
  singularity exec -B "$IGI_BIND" "$IGI_BWA_IMG" bash -c "samtools faidx $PRIM; bwa index $PRIM"
fi
echo "reference ready $(date)"

# 2. Reads that could belong in the MHC: those parked on HLA and alt contigs, those already on chr6 MHC, and
#    the unplaced unmapped pool.
grep -E "^HLA-|_alt$" "$FA.fai" | cut -f1 > "$W/divert.txt"
singularity exec -B "$IGI_BIND" "$IGI_BWA_IMG" bash -c "
  samtools view -@ 8 -b $W/n.bam \$(tr '\n' ' ' < $W/divert.txt) chr6:28000000-34000000 '*' > $W/cand.bam
  samtools collate -@ 8 -u -O $W/cand.bam $W/collate_tmp \
    | samtools fastq -@ 4 -1 $W/c_1.fq -2 $W/c_2.fq -s $W/c_s.fq -0 /dev/null -n -
"
echo "extracted: $(( $(wc -l < "$W/c_1.fq") / 4 )) pairs, $(( $(wc -l < "$W/c_s.fq") / 4 )) singletons $(date)"

# 3. Realign to the primary-only reference.
singularity exec -B "$IGI_BIND" "$IGI_BWA_IMG" bash -c "
  bwa mem -t 20 -R '@RG\tID:recover\tSM:IPISRC044\tPL:ILLUMINA' $PRIM $W/c_1.fq $W/c_2.fq 2>/dev/null \
    | samtools sort -@ 4 -m 2G -o $W/mhc_realigned.bam -
  samtools index $W/mhc_realigned.bam
  for r in 'HLA-A chr6:29940000-29945000' 'HLA-B chr6:31353000-31358000' 'HLA-C chr6:31268000-31272000'; do
    g=\${r%% *}; reg=\${r##* }
    samtools depth -r \$reg $W/mhc_realigned.bam \
      | awk -v g=\$g '{s+=\$3;n++} END{printf \"%s mean depth after realignment: %.1f\n\", g, (n?s/n:0)}'
  done
"

# 4. Recall the MHC.
printf 'chr6\t28000000\t34000000\n' > "$W/mhc.bed"
V=$W/IPISRC044.mhc.recovered.vcf.gz
singularity exec -B "$IGI_BIND" "$DV" /opt/deepvariant/bin/run_deepvariant \
  --model_type WGS --ref "$PRIM" --reads "$W/mhc_realigned.bam" --regions "$W/mhc.bed" \
  --output_vcf "$V" --num_shards 20 --intermediate_results_dir "$W/dv_tmp"

# 5. Verify, in the bcftools image. The first version of this ran the block in the bwa/samtools image, which
#    has no bcftools, so every count printed "command not found" followed by 0 -- and a reader of the log
#    would conclude the recovery had produced nothing when it had in fact produced 7,532 PASS records. The
#    counts are asserted now rather than merely printed.
BCF="singularity exec -B $IGI_BIND $IGI_BCFTOOLS_IMG bcftools"
$BCF --version | head -1
tot=$($BCF view -H -f PASS -r chr6:28000000-34000000 "$V" | wc -l)
het=$($BCF view -H -f PASS -g het -r chr6:28000000-34000000 "$V" | wc -l)
echo "recovered in chr6:28-34Mb: $tot PASS, $het heterozygous"
for r in "HLA-A chr6:29942000-29945000" "HLA-B chr6:31353000-31358000" "HLA-C chr6:31268000-31272000"; do
  g=${r%% *}; reg=${r##* }
  printf "  %-6s %s PASS, %s het\n" "$g" \
    "$($BCF view -H -f PASS -r "$reg" "$V" | wc -l)" "$($BCF view -H -f PASS -g het -r "$reg" "$V" | wc -l)"
done
# A recovery that genotyped nothing is a failure, and must not be reported as a completed job.
[ "$tot" -gt 1000 ] || { echo "FAIL: only $tot PASS records in the MHC; the recovery did not work"; exit 1; }
[ "$het" -gt 500 ] || { echo "FAIL: only $het heterozygous calls; phasing will not be possible"; exit 1; }
rm -rf "$W/dv_tmp" "$W/cand.bam" "$W/c_1.fq" "$W/c_2.fq" "$W/c_s.fq"
echo "DONE $(date)"
