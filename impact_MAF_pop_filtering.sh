#!/bin/bash

# ================================
# Impact + Allele Frequency + NFE Population Filtering Script
# ================================

# --------- Files and Directories ---------
OUTPUT_VCF_DIR="results"
INPUT_VCF="${OUTPUT_VCF_DIR}/pathogenicity_scored.vcf.gz"
OUTPUT_VCF="${OUTPUT_VCF_DIR}/impact_NFE_and_AF_filtered.vcf"
VEP_DATA_DIR="vep_data"
VEP_IMAGE="vep.sif"
FILTER_VEP="/opt/vep/src/ensembl-vep/filter_vep"

# --------- Logging ---------
LOG_FILE="${OUTPUT_VCF_DIR}/impact_af_nfe_filter_$(date '+%Y%m%d_%H%M%S').log"
mkdir -p "$OUTPUT_VCF_DIR"

echo "Starting impact, allele frequency, and population-based filtering..." | tee -a "$LOG_FILE"
echo "Input VCF:  $INPUT_VCF" | tee -a "$LOG_FILE"
echo "Output VCF: $OUTPUT_VCF" | tee -a "$LOG_FILE"

# --------- Run Filter_VEP ---------
singularity exec --bind "${VEP_DATA_DIR}:/opt/vep/.vep" "$VEP_IMAGE" \
  "$FILTER_VEP" \
  -i "$INPUT_VCF" \
  -o "$OUTPUT_VCF" \
  --format vcf \
  --force_overwrite \
  --filter "(IMPACT is HIGH or IMPACT is MODERATE) and (gnomADe_NFE_AF <= 0.01 or not gnomADe_NFE_AF)" \
  2>&1 | tee -a "$LOG_FILE"

# --------- Exit Status Check ---------
if [[ $? -eq 0 ]]; then
  echo "Impact, AF, and population filtering completed successfully." | tee -a "$LOG_FILE"
else
  echo "An error occurred during impact/population filtering. Check log: $LOG_FILE" | tee -a "$LOG_FILE"
fi
