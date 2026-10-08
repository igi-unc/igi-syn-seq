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
# chr1to6 is built in two halves, because the design builds it in two halves:
#
#   per-chromosome arms (WES, WGS, PacBio, ONT WGS)   the builder emits one file per chromosome, and
#                                                     91_merge_chr1to6.sbatch merges chr1-chr6 into one
#                                                     library per sample. The shards are NOT linked: a
#                                                     BAM row may match only one file, and with
#                                                     chromosome-first names no prefix selects one
#                                                     assay's six chromosomes anyway.
#   whole-library arms (bulk RNA, 10x GEX, Kinnex x2, ONT RNA x2)
#                                                     catalog/subset_chr1to6.py, driven by
#                                                     90_subset_chr1to6.sbatch: each read is kept or
#                                                     dropped by the transcript record its truth sidecar
#                                                     names. See that file on why selection is by truth
#                                                     map rather than by re-aligning.
#   10x TCR                                           design §12 keeps it COMPLETE in chr1to6, because
#                                                     TRA/TRB/TRG lie outside the subset. It stays linked
#                                                     under full/ ONLY: a second link under the same name
#                                                     made "<ds>-TCR_S1_L001" match four FASTQs, and the
#                                                     chr1to6 manifest points at the full/ copy, which is
#                                                     the same file.
#
# Every chr1to6 name carries a `chr1to6` label that the full name does not, and no File_Prefix is a
# prefix of another. That is a requirement, not tidiness: preflight_resolve_inputs walks every search
# directory recursively and matches basenames, so full/ and chr1to6/ are one namespace.
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

link_as() {  # link <src> <dest_dir> <new_basename>
  local src="$1" dir="$2" name="$3"
  if [ ! -e "$src" ]; then n_skip=$((n_skip+1)); echo "    skip (absent): $(basename "$src")"; return; fi
  if [ "$DRY" = 1 ]; then echo "    would link $name"; n_link=$((n_link+1)); return; fi
  mkdir -p "$dir"; ln -sfn "$src" "$dir/$name"; n_link=$((n_link+1))
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
    # Linked with a _wes_ infix the release does not have. RAFT resolves a sample by gathering every
    # FASTQ matching File_Prefix, and the release names WES as <ds>_<lib>_R1.fastq.gz -- so the prefix
    # <ds>_tumor would also match <ds>_tumor_wgs_R1 and <ds>_tumor_ont_wgs, and one manifest row would
    # claim three assays. A symlink's name is free to differ from its target.
    for m in 1 2; do
      link_as "$IGI_RELEASE/merged/$d/${d}_${l}_R${m}.fastq.gz" "$FQ_FULL" "${d}_${l}_wes_R${m}.fastq.gz"
    done
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

  # The chr1to6 directories are rebuilt from scratch each run, because this script is the only thing that
  # writes them and an earlier layout's links would otherwise survive forever. Only symlinks are removed,
  # and the ${VAR:?} guards mean an unset path deletes nothing.
  echo "=== $d: chr1to6 ==="
  if [ "$DRY" != 1 ]; then
    for dir in "${FQ_SUB:?}" "${BAM_SUB:?}"; do
      [ -d "$dir" ] && find "$dir" -maxdepth 1 -type l -delete
    done
  fi
  # EVERY chr1to6 deliverable is derived, so every one of them can be stale, and a stale derived library
  # is the most dangerous thing in this tree: a RAFT run consumes it and succeeds. Each is therefore
  # linked only if it is newer than what it was derived FROM -- the chr1-6 shards for the merged arms,
  # the full library for the filtered ones. This is not hypothetical: the `*`-allele fix invalidated
  # every ds-01 library while that dataset's subsets were still being written from the old reads.
  sub_link() {  # sub_link <src> <dest_dir> <rename|-> -- <inputs...>
    local src="$1" dir="$2" name="$3"; shift 4
    if fresher_than_inputs "$src" "$@"; then
      if [ "$name" = "-" ]; then link "$src" "$dir"; else link_as "$src" "$dir" "$name"; fi
    else
      echo "    HELD (older than what it is derived from): $(basename "$src")"; n_pend=$((n_pend+1))
    fi
  }
  for l in tumor normal; do
    for m in 1 2; do
      # WES keeps the _wes_ infix the release does not have, for the same reason the full release does.
      sub_link "$IGI_RELEASE/wes_chr1to6/$d/${d}_chr1to6_${l}_R${m}.fastq.gz" "$FQ_SUB" \
               "${d}_chr1to6_${l}_wes_R${m}.fastq.gz" -- \
               $IGI_RELEASE/wes/$d/${d}_chr[1-6]_${l}_R${m}.fastq.gz
      sub_link "$IGI_RELEASE/wgs_chr1to6/$d/${d}_chr1to6_${l}_wgs_R${m}.fastq.gz" "$FQ_SUB" - -- \
               $IGI_RELEASE/wgs/$d/${d}_chr[1-6]_${l}_wgs_R${m}.fastq.gz
    done
    sub_link "$IGI_RELEASE/ont_wgs_chr1to6/$d/${d}_chr1to6_${l}_ont_wgs.fastq.gz" "$FQ_SUB" - -- \
             $IGI_RELEASE/ont_wgs/$d/${d}_chr[1-6]_${l}_ont_wgs.fastq.gz
    for x in "" ".pbi"; do
      sub_link "$IGI_RELEASE/pacbio_chr1to6/$d/${d}_chr1to6_${l}_hifi.bam$x" "$BAM_SUB" - -- \
               $IGI_RELEASE/pacbio/$d/${d}_chr[1-6]_${l}_hifi.bam
    done
  done
  for m in 1 2; do
    sub_link "$IGI_RELEASE/rna/$d/${d}_chr1to6_rna_R${m}.fastq.gz" "$FQ_SUB" - -- \
             "$IGI_RELEASE/rna/$d/${d}_full_rna_R${m}.fastq.gz"
    sub_link "$IGI_RELEASE/tenx_gex/$d/${d}-chr1to6-GEX_S1_L001_R${m}_001.fastq.gz" "$FQ_SUB" - -- \
             "$IGI_RELEASE/tenx_gex/$d/${d}-GEX_S1_L001_R${m}_001.fastq.gz"
  done
  for a in ont_bulk_rna ont_sc_rna; do
    sub_link "$IGI_RELEASE/$a/$d/${d}_chr1to6_${a}.fastq.gz" "$FQ_SUB" - -- \
             "$IGI_RELEASE/$a/$d/${d}_full_${a}.fastq.gz"
  done
  for a in kinnex_bulk kinnex_sc; do
    for x in "" ".pbi"; do
      sub_link "$IGI_RELEASE/$a/$d/${d}_chr1to6_${a}_segmented.bam$x" "$BAM_SUB" - -- \
               "$IGI_RELEASE/$a/$d/${d}_full_${a}_segmented.bam"
    done
  done
  echo "    10x TCR: kept complete per design §12; the chr1to6 manifest points at the full/ copy"
done

echo
echo "  linked $n_link   absent $n_skip   pending the chr1to6 filter $n_pend"
