#!/bin/bash

# ================================
# Genotype Quality and Depth Filtering Script
# ================================

# --------- Files and Directories ---------
OUTPUT_VCF_DIR="results"
INPUT_VCF="${OUTPUT_VCF_DIR}/impact_NFE_and_AF_filtered.vcf.gz"
OUTPUT_VCF="${OUTPUT_VCF_DIR}/GQ_nonref_filtered.vcf"

# --------- Logging ---------
LOG_FILE="${OUTPUT_VCF_DIR}/gq_filter_$(date '+%Y%m%d_%H%M%S').log"
mkdir -p "$OUTPUT_VCF_DIR"

echo "Starting genotype quality and non-reference filtering..." | tee -a "$LOG_FILE"
echo "Input VCF:  $INPUT_VCF" | tee -a "$LOG_FILE"
echo "Output VCF: $OUTPUT_VCF" | tee -a "$LOG_FILE"

# --------- bcftools Filtering ---------
bcftools view \
  -i 'COUNT(GT!="0/0" && GT!="./." && FMT/GQ >= 20 && FMT/DP >= 15) > 0' \
  -Oz -o "$OUTPUT_VCF" "$INPUT_VCF" \
  2>&1 | tee -a "$LOG_FILE"

# --------- Final Status ---------
if [[ $? -eq 0 ]]; then
  echo "Genotype quality filtering completed successfully." | tee -a "$LOG_FILE"
else
  echo "An error occurred. Check $LOG_FILE for details." | tee -a "$LOG_FILE"
fi
