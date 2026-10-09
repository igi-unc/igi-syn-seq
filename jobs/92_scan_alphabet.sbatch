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

# One deterministic list, so task N always takes the same files and a re-run of one task is meaningful.
LIST=$OUT/scanlist.txt
if [ ! -s "$LIST" ]; then
  find "$R"/wes "$R"/wgs "$R"/ont_wgs "$R"/pacbio "$R"/ont_bulk_rna "$R"/ont_sc_rna "$R"/rna \
       "$R"/tenx_gex "$R"/tenx_tcr "$R"/kinnex_bulk "$R"/kinnex_sc -maxdepth 2 \
       \( -name '*.fastq.gz' -o -name '*_hifi.bam' -o -name '*_segmented.bam' \) 2>/dev/null \
    | sort > "$LIST.$SLURM_ARRAY_TASK_ID.tmp"
  mv -n "$LIST.$SLURM_ARRAY_TASK_ID.tmp" "$LIST"
  rm -f "${LIST:?}.$SLURM_ARRAY_TASK_ID.tmp"
fi

RES=$OUT/task_${SLURM_ARRAY_TASK_ID}.tsv
: > "$RES"

# The regex pre-test matters: gsub rebuilds the string on every read it touches, and skipping it on the
# clean reads took one exome library from 105 s to 11 s for the same counts.
TALLY='{ n++; if ($0 ~ /[^ACGTNacgtn]/) { s = $0; gsub(/[ACGTNacgtn]/, "", s);
         for (i = 1; i <= length(s); i++) c[substr(s, i, 1)]++ } }
       END { printf "%d", n; for (k in c) printf "\t%s=%d", k, c[k]; printf "\n" }'

i=0
while IFS= read -r f; do
  i=$((i + 1))
  [ $(( (i - 1) % 96 + 1 )) -eq "$SLURM_ARRAY_TASK_ID" ] || continue
  if [[ $f == *.bam ]]; then
    # No `head` anywhere in the pipeline, so pipefail can stay on: "nothing was read" must never be
    # indistinguishable from "nothing offended".
    out=$(set -o pipefail
          singularity exec -B "$IGI_BIND" "$IGI_SAMTOOLS_IMG" samtools view "$f" \
            | cut -f10 | awk "$TALLY") || out="ERROR"
  else
    out=$(set -o pipefail; pigz -dc "$f" | awk 'NR % 4 == 2' | awk "$TALLY") || out="ERROR"
  fi
  printf '%s\t%s\n' "$f" "$out" >> "$RES"
  echo "[$i] $(basename "$f") -> $out"
done < "$LIST"
echo "SCAN TASK $SLURM_ARRAY_TASK_ID DONE"
