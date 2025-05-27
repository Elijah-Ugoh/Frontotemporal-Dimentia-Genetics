#!/bin/bash

# ========================
# Pathogenicity Filtering Script using VEP
# ========================

# --------- Files and Directories ---------
VEP_DATA_DIR="vep_data"
VEP_IMAGE="vep.sif"
FILTER_VEP="/opt/vep/src/ensembl-vep/filter_vep"

INPUT_VCF="results/consequence_filtered.vcf.gz"
OUTPUT_VCF="results/pathogenicity_scored.vcf"

# --------- Logging ---------
LOG_FILE="results/pathogenicity_filter_$(date '+%Y%m%d_%H%M%S').log"
mkdir -p results

echo "Starting pathogenicity filtering..." | tee -a "$LOG_FILE"
echo "Input VCF:  $INPUT_VCF" | tee -a "$LOG_FILE"
echo "Output VCF: $OUTPUT_VCF" | tee -a "$LOG_FILE"

# --------- Filtering Expression ---------
FILTER_EXPR="((CADD_PHRED ne '.' and CADD_PHRED > 20) or \
           ((SIFT_score ne '.' and SIFT_score < 0.05) or SIFT_pred is deleterious) and \
           (Polyphen2_HDIV_score ne '.' and (Polyphen2_pred is possibly_damaging or Polyphen2_pred is probably_damaging or Polyphen2_HDIV_score > 0.8)) or \
           (MetaLR_score ne '.' and MetaLR_score > 0.5) or \
           (GERP++_RS ne '.' and GERP++_RS > 4) or \
           (M-CAP_score ne '.' and M-CAP_score > 0.025) or \
           (DANN_score ne '.' and DANN_score > 0.96) or \
           (SpliceAI_DS_AG ne '.' and SpliceAI_DS_AG > 0.5) or \
           (SpliceAI_DS_AL ne '.' and SpliceAI_DS_AL > 0.5) or \
           (SpliceAI_DS_DG ne '.' and SpliceAI_DS_DG > 0.5) or \
           (SpliceAI_DS_DL ne '.' and SpliceAI_DS_DL > 0.5) or \
           (ada_score ne '.' and ada_score > 0.6) or \
           (rf_score ne '.' and rf_score > 0.6))"

# --------- Run Filter_VEP ---------
singularity exec --bind "${VEP_DATA_DIR}:/opt/vep/.vep" "$VEP_IMAGE" \
  "$FILTER_VEP" \
  -i "$INPUT_VCF" \
  -o "$OUTPUT_VCF" \
  --filter "$FILTER_EXPR" \
  --format vcf \
  --force_overwrite \
  2>&1 | tee -a "$LOG_FILE"

# --------- Exit Status Check ---------
if [[ $? -eq 0 ]]; then
  echo "Pathogenicity scoring completed successfully." | tee -a "$LOG_FILE"
else
  echo "An error occurred during pathogenicity filtering. Check the log: $LOG_FILE" | tee -a "$LOG_FILE"
fi