#!/bin/bash
# Stage the release into the RAFT workspace as SYMLINKS, per the layout in design §2:
#
#   inputs/fastqs/<ds>/<full|chr1to6>/   Illumina + ONT FASTQs
#   inputs/bams/<ds>/<full|chr1to6>/     PacBio unaligned BAMs (+ .pbi) and Kinnex segmented BAMs
#
# Symlinks, not copies: the full release is ~1.3 TB and lives on scratch, which the input manifest marks
# `durable: no`. A symlink makes that visible rather than hiding it behind a copy that looks durable and
# is not. If the release is ever moved off scratch, re-run this and the links follow.
#
# Only the FINAL deliverables are linked. A link to a library that is mid-rebuild would be worse than a
# missing one, because a RAFT run would consume it and succeed.
#
# chr1to6 is handled in two halves, because the design handles it in two halves:
#
#   per-chromosome arms (WES, WGS, PacBio, ONT WGS)   the builder already emits one file per chromosome,
#                                                     so the subset is a plain selection of chr1-chr6
#   whole-library arms (bulk RNA, 10x, Kinnex, ONT RNA)
#                                                     design §12: "derived from the full release by
#                                                     alignment, not re-simulated" -- keep reads whose
#                                                     alignment overlaps chr1to6, plus viral and unmapped
#                                                     reads and mates. That filter is not built yet, so
#                                                     these are reported as pending rather than linked.
#   10x TCR                                           design §12 keeps it COMPLETE in chr1to6, because
#                                                     TRA/TRB/TRG lie outside the subset, so chr1to6
#                                                     links the full library.
set -euo pipefail
: "${IGI_RELEASE:?export IGI_RELEASE=/path/to/release}"
: "${IGI_RAFT:?export IGI_RAFT=/path/to/raft/workspace}"
DATASETS=(${IGI_DATASETS:-IGI-SYN-SEQ-01 IGI-SYN-SEQ-02})
DRY=${DRY_RUN:-0}

n_link=0; n_skip=0; n_pend=0
# A merged library that is OLDER than its own per-chromosome inputs is mid-rebuild, and linking it would
# be worse than leaving it out: a RAFT run would consume the superseded library and succeed. Right now the
# three WGS arms are exactly in that state -- every chromosome was re-simulated for the copy-number fix
# and the merge has not run yet -- so they are held back by this check rather than by a date I typed in.
fresher_than_inputs() {
  local out="$1"; shift
  [ -e "$out" ] || return 1
  local f
  for f in "$@"; do
    [ -e "$f" ] || continue
    [ "$f" -nt "$out" ] && return 1
  done
  return 0
}

link() {  # link <src> <dest_dir>
  local src="$1" dir="$2"
  if [ ! -e "$src" ]; then n_skip=$((n_skip+1)); echo "    skip (absent): $(basename "$src")"; return; fi
  if [ "$DRY" = 1 ]; then echo "    would link $(basename "$src")"; n_link=$((n_link+1)); return; fi
  mkdir -p "$dir"
  ln -sfn "$src" "$dir/$(basename "$src")"
  n_link=$((n_link+1))
}
pending() { n_pend=$((n_pend+1)); echo "    PENDING (needs the chr1to6 alignment filter): $1"; }

for d in "${DATASETS[@]}"; do
  FQ_FULL=$IGI_RAFT/inputs/fastqs/$d/full
  FQ_SUB=$IGI_RAFT/inputs/fastqs/$d/chr1to6
  BAM_FULL=$IGI_RAFT/inputs/bams/$d/full
  BAM_SUB=$IGI_RAFT/inputs/bams/$d/chr1to6
  echo "=== $d: full ==="
  for l in tumor normal; do
    for m in 1 2; do link "$IGI_RELEASE/merged/$d/${d}_${l}_R${m}.fastq.gz" "$FQ_FULL"; done
    if fresher_than_inputs "$IGI_RELEASE/wgs_merged/$d/${d}_${l}_wgs_R1.fastq.gz" \
         $IGI_RELEASE/wgs/$d/${d}_chr*_${l}_wgs_R1.fastq.gz; then
      for m in 1 2; do link "$IGI_RELEASE/wgs_merged/$d/${d}_${l}_wgs_R${m}.fastq.gz" "$FQ_FULL"; done
    else
      echo "    HELD (older than its chromosomes; merge pending): ${d}_${l}_wgs_R1/R2"; n_pend=$((n_pend+1))
    fi
    if fresher_than_inputs "$IGI_RELEASE/ont_wgs_merged/$d/${d}_${l}_ont_wgs.fastq.gz" \
         $IGI_RELEASE/ont_wgs/$d/${d}_chr*_${l}_ont_wgs.fastq.gz; then
      link "$IGI_RELEASE/ont_wgs_merged/$d/${d}_${l}_ont_wgs.fastq.gz" "$FQ_FULL"
    else
      echo "    HELD (older than its chromosomes; merge pending): ${d}_${l}_ont_wgs"; n_pend=$((n_pend+1))
    fi
    if fresher_than_inputs "$IGI_RELEASE/pacbio_merged/$d/${d}_${l}_hifi.bam" \
         $IGI_RELEASE/pacbio/$d/${d}_chr*_${l}_hifi.bam; then
      link "$IGI_RELEASE/pacbio_merged/$d/${d}_${l}_hifi.bam" "$BAM_FULL"
      link "$IGI_RELEASE/pacbio_merged/$d/${d}_${l}_hifi.bam.pbi" "$BAM_FULL"
    else
      echo "    HELD (older than its chromosomes; merge pending): ${d}_${l}_hifi.bam"; n_pend=$((n_pend+1))
    fi
  done
  for m in 1 2; do link "$IGI_RELEASE/rna/$d/${d}_full_rna_R${m}.fastq.gz" "$FQ_FULL"; done
  for r in R1 R2; do
    link "$IGI_RELEASE/tenx_gex/$d/${d}-GEX_S1_L001_${r}_001.fastq.gz" "$FQ_FULL"
    link "$IGI_RELEASE/tenx_tcr/$d/${d}-TCR_S1_L001_${r}_001.fastq.gz" "$FQ_FULL"
  done
  link "$IGI_RELEASE/ont_bulk_rna/$d/${d}_full_ont_bulk_rna.fastq.gz" "$FQ_FULL"
  link "$IGI_RELEASE/ont_sc_rna/$d/${d}_full_ont_sc_rna.fastq.gz" "$FQ_FULL"
  for a in kinnex_bulk kinnex_sc; do
    link "$IGI_RELEASE/$a/$d/${d}_full_${a}_segmented.bam" "$BAM_FULL"
    link "$IGI_RELEASE/$a/$d/${d}_full_${a}_segmented.bam.pbi" "$BAM_FULL"
  done

  echo "=== $d: chr1to6 (per-chromosome arms are a plain selection) ==="
  for c in 1 2 3 4 5 6; do
    for l in tumor normal; do
      for m in 1 2; do
        link "$IGI_RELEASE/wes/$d/${d}_chr${c}_${l}_R${m}.fastq.gz" "$FQ_SUB"
        link "$IGI_RELEASE/wgs/$d/${d}_chr${c}_${l}_wgs_R${m}.fastq.gz" "$FQ_SUB"
      done
      link "$IGI_RELEASE/ont_wgs/$d/${d}_chr${c}_${l}_ont_wgs.fastq.gz" "$FQ_SUB"
      link "$IGI_RELEASE/pacbio/$d/${d}_chr${c}_${l}_hifi.bam" "$BAM_SUB"
      link "$IGI_RELEASE/pacbio/$d/${d}_chr${c}_${l}_hifi.bam.pbi" "$BAM_SUB"
    done
  done
  # kept complete by design §12, so chr1to6 points at the full library
  for r in R1 R2; do link "$IGI_RELEASE/tenx_tcr/$d/${d}-TCR_S1_L001_${r}_001.fastq.gz" "$FQ_SUB"; done
  echo "=== $d: chr1to6 (whole-library arms) ==="
  for a in "bulk RNA" "10x GEX" "Kinnex bulk" "Kinnex sc" "ONT bulk RNA" "ONT sc RNA"; do pending "$a"; done
done

echo
echo "  linked $n_link   absent $n_skip   pending the chr1to6 filter $n_pend"
