#!/bin/bash
#SBATCH -J igi_alphabet
#SBATCH -c 2 --mem=4G -t 8:00:00
#SBATCH --array=1-96
# Count every non-ACGTN character in every delivered base, across the whole release.
#
# WHY THIS IS A JOB AND NOT A FLAG ON qa_release.py: the release is ~1.3 TB over 612 files, and the
# offending characters sit at roughly ONE READ PER MILLION in the exome -- 94 positions in GRCh38 are
# ambiguous, and a targeted assay covers few of them. The first version of the alphabet check sampled
# 100,000 reads per file and reported `none` for the exome library that aborted OptiType's razers3 with
# `seqan::ParseError: Unexpected character 'Y' found`. A sample can witness contamination; it cannot
# establish absence below its own resolution. So the scan is complete, and completeness costs a pass
# over every base -- about 40 min at this width, against ~33 h if it ran serially.
#
# IT SCANS THE PER-CHROMOSOME UNITS, NOT THE MERGED LIBRARIES. Those units are exactly the array tasks
# that build them, so a dirty file names the task to re-run: the ambiguous positions are confined to
# particular chromosomes (none at all on chr4, 5, 8, 11, 14, 15, 18, 19, 20 or chrY), and rebuilding
# only the chromosomes that can carry a code is the difference between a targeted repair and
# regenerating two whole datasets.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh before submitting}"

R=$IGI_WORK/release
OUT=$IGI_WORK/release/alphabet_scan
mkdir -p "$OUT"

# ONE deterministic list, built by ONE task. The first version let whichever task arrived first build it
# via tmp+`mv`, which looked atomic and was -- but a task that already had the file OPEN when another's
# `mv` replaced the inode lost its NFS file handle:
#
#   [1] ..._kinnex_bulk_segmented.bam -> 5610724   *=11
#   slurm_script: line 48: read: read error: 0: Stale file handle
#   SCAN TASK 1 DONE
#
# A `while read` loop treats a failed read as end of input, so the loop exited normally, the script
# printed DONE and SLURM recorded COMPLETED -- with 6 of 7 files unscanned. Six tasks did that and 36
# files went unexamined behind a clean exit status. That is the same defect this whole scan exists to
# catch: a check that reports success without having done its work.
#
# So: task 1 builds the list and nobody else writes it, every task reads it ONCE into memory so a
# replaced inode cannot truncate the loop mid-way, and the task FAILS LOUDLY if it did not process
# everything it was assigned.
LIST=$OUT/scanlist.txt
if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
  find "$R"/wes "$R"/wgs "$R"/ont_wgs "$R"/pacbio "$R"/ont_bulk_rna "$R"/ont_sc_rna "$R"/rna \
       "$R"/tenx_gex "$R"/tenx_tcr "$R"/kinnex_bulk "$R"/kinnex_sc -maxdepth 2 \
       \( -name '*.fastq.gz' -o -name '*_hifi.bam' -o -name '*_segmented.bam' \) 2>/dev/null \
    | sort > "$LIST.building"
  mv "$LIST.building" "$LIST"
else
  for _ in $(seq 1 600); do [ -s "$LIST" ] && break; sleep 2; done
  [ -s "$LIST" ] || { echo "FATAL: $LIST never appeared; task 1 must run"; exit 1; }
fi
# Read the whole list up front. This is the fix for the stale handle: no descriptor stays open across the
# scan, so nothing that happens to the file afterwards can silently end the loop.
mapfile -t FILES < "$LIST"
[ "${#FILES[@]}" -gt 0 ] || { echo "FATAL: empty scan list"; exit 1; }

RES=$OUT/task_${SLURM_ARRAY_TASK_ID}.tsv
: > "$RES"

# The regex pre-test matters: gsub rebuilds the string on every read it touches, and skipping it on the
# clean reads took one exome library from 105 s to 11 s for the same counts.
# `ph` is set for BAM only: in SAM a SEQ of exactly `*` means "no sequence stored", which is not a base.
# Counting it as a symbolic allele was a false positive -- ds-02's Kinnex BAM reported `*=7` with zero `*`
# embedded in any sequence. Reported separately, as EMPTY.
TALLY='{ n++; if (ph && $0 == "*") { e++; next }
         if ($0 ~ /[^ACGTNacgtn]/) { s = $0; gsub(/[ACGTNacgtn]/, "", s);
         for (i = 1; i <= length(s); i++) c[substr(s, i, 1)]++ } }
       END { printf "%d", n; if (e) printf "\tEMPTY=%d", e;
             for (k in c) printf "\t%s=%d", k, c[k]; printf "\n" }'

i=0
assigned=0
done_count=0
for f in "${FILES[@]}"; do
  i=$((i + 1))
  [ $(( (i - 1) % 96 + 1 )) -eq "$SLURM_ARRAY_TASK_ID" ] || continue
  assigned=$((assigned + 1))
  if [[ $f == *.bam ]]; then
    # No `head` anywhere in the pipeline, so pipefail can stay on: "nothing was read" must never be
    # indistinguishable from "nothing offended". `< /dev/null` because singularity drains stdin.
    out=$(set -o pipefail
          singularity exec -B "$IGI_BIND" "$IGI_SAMTOOLS_IMG" samtools view "$f" < /dev/null \
            | cut -f10 | awk -v ph=1 "$TALLY") || out="ERROR"
  else
    out=$(set -o pipefail; pigz -dc "$f" < /dev/null | awk 'NR % 4 == 2' | awk -v ph=0 "$TALLY") \
      || out="ERROR"
  fi
  printf '%s\t%s\n' "$f" "$out" >> "$RES"
  done_count=$((done_count + 1))
  echo "[$i] $(basename "$f") -> $out"
done

# Completeness is part of the result. Exiting 0 having scanned a fraction of the assignment is what the
# stale-handle bug did, and it cost 36 files behind 96 COMPLETED tasks.
if [ "$done_count" -ne "$assigned" ]; then
  echo "FATAL: assigned $assigned files, completed $done_count"
  exit 1
fi
echo "SCAN TASK $SLURM_ARRAY_TASK_ID DONE ($done_count/$assigned files)"
