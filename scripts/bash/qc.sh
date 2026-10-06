#!/bin/bash

# ============================================================
# 16S rRNA DATASET - RAW READ QC
# Project: PRJNA759642
# ============================================================


# ============================================================
# 1. ENVIRONMENT
# ============================================================

conda activate rnaseq

# ============================================================
# 2. DIRECTORIES
# ============================================================

RAW_DIR="../data/raw_fastq_16S_main_cohort"
FASTQC_DIR="../results/fastqc"
MULTIQC_DIR="../results/multiqc"

mkdir -p "$FASTQC_DIR"
mkdir -p "$MULTIQC_DIR"

RUNS=(
SRR15702528 SRR15702529 SRR15702530 SRR15702531 SRR15702532 SRR15702534
SRR15702535 SRR15702536 SRR15702537 SRR15702538 SRR15702540 SRR15702542
SRR15702544 SRR15702546 SRR15702548 SRR15702552 SRR15702554 SRR15702556
SRR15702558 SRR15702559 SRR15702560 SRR15702561 SRR15702562 SRR15702563
SRR15702564 SRR15702566 SRR15702567 SRR15702568 SRR15702628 SRR15702629
SRR15702630 SRR15702631 SRR15702632 SRR15702633 SRR15702634 SRR15702635
SRR15702636 SRR15702638 SRR15702639 SRR15702640 SRR15702641 SRR15702642
)

# ============================================================
# 3. WORKING DIRECTORY
# ============================================================

echo "Working directory: $(pwd)"
echo "Raw data directory: $RAW_DIR"

# ============================================================
# 5. CHECK INPUT FILES
# ============================================================

echo "Checking FASTQ files..."

for RUN in "${RUNS[@]}"; do

  if ! ls "$RAW_DIR/${RUN}"_*.fastq.gz >/dev/null 2>&1; then
    echo "ERROR: FASTQ files not found for $RUN"
    exit 1
  fi

done

echo "All FASTQ files found."


# ============================================================
# 6. FASTQC
# ============================================================

echo "Running FastQC..."

for RUN in "${RUNS[@]}"; do

  echo "Processing $RUN..."

  fastqc \
    -t 4 \
    -o "$FASTQC_DIR" \
    "$RAW_DIR/${RUN}"_*.fastq.gz

done


# ============================================================
# 7. MULTIQC
# ============================================================

echo "Running MultiQC..."

multiqc \
  "$FASTQC_DIR" \
  -o "$MULTIQC_DIR"


# ============================================================
# 8. COMPLETION
# ============================================================

echo "QC completed successfully."
echo "FastQC reports: $FASTQC_DIR"
echo "MultiQC report: $MULTIQC_DIR"