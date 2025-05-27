#!/bin/bash

# =======================================
# Final Step: Filter Variants in Known FTD Genes
# =======================================

# --------- Files and Directories ---------
OUTPUT_VCF_DIR="results"
INPUT_VCF="${OUTPUT_VCF_DIR}/GQ_nonref_filtered.vcf"
OUTPUT_VCF="${OUTPUT_VCF_DIR}/candidate_variants.vcf"
GENE_LIST="FTD_genes_HP_0002145.txt"
SIF_IMAGE="vep.sif"
VEP_DATA_DIR="vep_data"
FILTER_VEP="/opt/vep/src/ensembl-vep/filter_vep"

# --------- Logging ---------
LOG_FILE="${OUTPUT_VCF_DIR}/ftd_filter_$(date '+%Y%m%d_%H%M%S').log"
mkdir -p "$OUTPUT_VCF_DIR"

echo "Starting final filtering step for FTD-associated genes..." | tee -a "$LOG_FILE"
echo "Input VCF:  $INPUT_VCF" | tee -a "$LOG_FILE"
echo "Gene List:  $GENE_LIST" | tee -a "$LOG_FILE"
echo "Output VCF: $OUTPUT_VCF" | tee -a "$LOG_FILE"

# --------- Run filter_vep ---------
singularity exec \
  --bind "${VEP_DATA_DIR}:/opt/vep/.vep" \
  "$SIF_IMAGE" \
  "$FILTER_VEP" \
  -i "$INPUT_VCF" \
  -o "$OUTPUT_VCF" \
  --format vcf \
  --only_matched \
  --force_overwrite \
  --filter "SYMBOL in $GENE_LIST" \
  2>&1 | tee -a "$LOG_FILE"

# --------- Final Status ---------
if [[ $? -eq 0 ]]; then
  echo "FTD gene filtering completed successfully." | tee -a "$LOG_FILE"
else
  echo "An error occurred. Check $LOG_FILE for details." | tee -a "$LOG_FILE"
fi