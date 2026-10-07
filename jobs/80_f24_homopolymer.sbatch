#!/bin/bash
#SBATCH -J igi_f24
#SBATCH -c 16 --mem=48G -t 12:00:00
#SBATCH --array=1-2
# F24: does Badread put HiFi indel errors where real HiFi puts them?
#
#   task 1  SIMULATED, from the plain reference. Badread is run with the release's exact PacBio
#           parameters over unmodified GRCh38, so every indel in the alignment is a simulator error by
#           construction and nothing has to be excluded. Simulating from a derived chromosome instead
#           would mix designed and germline indels into the measurement.
#   task 2  REAL HG002 Revio HiFi, the same movie the length and identity models were fitted from, with
#           indels within 10 bp of a GIAB Q100 truth indel excluded. Without that exclusion the
#           measurement is mostly germline biology: HiFi indel error is ~1e-5/base and human germline
#           indels ~5e-4/base, fifty times larger, and germline indels are themselves
#           homopolymer-enriched, which would flatter the simulator.
set -euo pipefail
DS=/datastore/nextgenout5/share/labs/Vincent_Lab/datasets
source "${IGI_ENV:?}"
W=$IGI_WORK/f24; mkdir -p "$W"
MM="singularity exec -B /datastore $(dirname $IGI_SAMTOOLS_IMG)/quay.io_biocontainers_minimap2_2.28--h577a1d6_4.sif minimap2"
SAM="singularity exec -B /datastore $IGI_SAMTOOLS_IMG samtools"
REF=$IGI_FASTA_NOALT
t=$SLURM_ARRAY_TASK_ID

if [ "$t" = 1 ]; then
  LABEL=simulated-from-reference
  FQ=$W/sim.fq
  if [ ! -s "$FQ" ]; then
    # one 40 Mb window of chr1 is plenty for 20,000 reads at 16.7 kb, and keeps the reference small
    $SAM faidx "$REF" chr1:20000000-60000000 > $W/win.fa
    BR=$(python3 -c "import yaml;print(yaml.safe_load(open('$IGI_PATHS'))['badread_cmd'])")
    IDENT=$(python3 -c "import yaml;print(yaml.safe_load(open('$IGI_PATHS')).get('pacbio_identity','99.4,99.8,0.5'))")
    LM=$(python3 -c "import yaml;print(yaml.safe_load(open('$IGI_PATHS')).get('pacbio_length_mean',16689))")
    LS=$(python3 -c "import yaml;print(yaml.safe_load(open('$IGI_PATHS')).get('pacbio_length_sd',4593))")
    ARGS="simulate --reference $W/win.fa --quantity 10x --length $LM,$LS --identity $IDENT \
      --error_model pacbio2021 --qscore_model pacbio2021 --seed 24 \
      --start_adapter_seq '' --end_adapter_seq '' --junk_reads 0 --random_reads 0 --chimeras 0"
    eval "${BR/\{args\}/$ARGS}" > "$FQ" 2>$W/badread.err
  fi
  EXTRA=""
else
  LABEL=real-HG002-Revio
  FQ=$W/real.fq
  if [ ! -s "$FQ" ]; then
    R=/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/HG002/PacBio_HiFi-Revio_20231031/HG002_PacBio-Revio_m84039_230928_213653_s3.hifi_reads.bam
    $SAM view "$R" | head -40000 | awk '{print "@"$1"\n"$10"\n+\n"$11}' > "$FQ"
  fi
  GIAB=$(python3 -c "import yaml;print(yaml.safe_load(open('$IGI_PATHS'))['germline_vcf_baseline']['HG002'])")
  EXTRA="--exclude-vcf $GIAB"
fi

B=$W/${LABEL}.bam
if [ ! -s "$B" ] || [ "$FQ" -nt "$B" ]; then
  $MM -ax map-hifi -t 16 --secondary=no "$REF" "$FQ" 2>$W/${LABEL}.mm.err \
    | $SAM sort -@4 -o "$B" -
  $SAM index "$B"
fi
cd "$IGI_REPO"
python3 -u jobs/measure/homopolymer_indel_rate.py --bam "$B" --ref "$REF" \
  --out "$W/${LABEL}.f24.json" --label "$LABEL" --min-run 6 --max-reads 20000 $EXTRA
echo "DONE $LABEL"
