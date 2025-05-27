#!/bin/bash

# === Files and Directories ===
OUTPUT_VCF_DIR="results"
INPUT_VCF="${OUTPUT_VCF_DIR}/output.annotated.vcf.gz"

# Log file for this step
LOG_FILE="${OUTPUT_VCF_DIR}/consequence_filter.log"
mkdir -p "$OUTPUT_VCF_DIR"

echo "Starting consequence-based variant filtering..." | tee -a "$LOG_FILE"

# Run consequence filtering using VEP's filter_vep tool
singularity exec \
  --bind vep_data/:/opt/vep/.vep \
  vep.sif \
  /opt/vep/src/ensembl-vep/filter_vep \
  -i ${OUTPUT_VCF_DIR}/output.annotated.vcf.gz \
  -o ${OUTPUT_VCF_DIR}/consequence_filtered.vcf \
  --format vcf \
  --force_overwrite \
  --only_matched \
  -filter "Consequence is stop_gained or Consequence is stop_lost or Consequence is start_lost or Consequence is frameshift_variant or Consequence is splice_donor_variant or Consequence is splice_acceptor_variant or Consequence is missense_variant or Consequence is splice_region_variant or Consequence is coding_sequence_variant or Consequence is regulatory_region_variant or Consequence is inframe_deletion or Consequence is inframe_insertion or Consequence is protein_altering_variant" \
  2>&1 | tee -a "$LOG_FILE"

echo "Consequence filtering completed." | tee -a "$LOG_FILE"
