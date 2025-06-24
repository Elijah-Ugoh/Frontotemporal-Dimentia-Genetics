#!/bin/bash

# === Files and Directories ===
OUTPUT_VCF_DIR="results"
mkdir -p "$OUTPUT_VCF_DIR"
SIF_IMAGE="vep.sif"
VEP_DATA_DIR="vep_data"
INPUT_VCF="vcf/landqvist.postVQSR.unfiltered.vcf"
PLUGIN_DIR="Plugins"

# Plugin Files
SPLICEAI_SNV="Plugins/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz"
SPLICEAI_INDEL="Plugins/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz"
DBNSFP="Plugins/dbNSFP/dbNSFP5.1a_grch38-003.gz"
DBSCSNV="Plugins/dbscSNV/dbscSNV1.1_GRCh38.txt.gz"
CADD="Plugins/CADD/whole_genome_SNVs.tsv.gz"

# === Logging Setup ===
LOG_FILE="${OUTPUT_VCF_DIR}/annotation_$(date '+%Y%m%d_%H%M%S').log"

echo "Starting VEP annotation..." | tee -a "$LOG_FILE"
echo "Using input VCF: $INPUT_VCF" | tee -a "$LOG_FILE"

# === Run VEP Annotation ===
singularity exec \
  --bind $(pwd)/${VEP_DATA_DIR}:/opt/vep/.vep \
  ${SIF_IMAGE} \
  vep --dir /opt/vep/.vep \
      --cache --offline \
      --species homo_sapiens \
      --assembly GRCh38 \
      --fasta /opt/vep/.vep/homo_sapiens/113_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
      --format vcf \
      -i ${INPUT_VCF} \
      -o ${OUTPUT_VCF_DIR}/output.annotated.vcf.gz \
      --vcf \
      --compress_output bgzip \
      --force_overwrite \
      --hgvs \
      --symbol \
      --numbers \
      --canonical \
      --ccds \
      --biotype \
      --uniprot \
      --tsl \
      --appris \
      --gene_phenotype \
      --af \
      --af_gnomad \
      --max_af \
      --variant_class \
      --sift b \
      --polyphen b \
      --plugin SpliceAI,snv=${SPLICEAI_SNV},indel=${SPLICEAI_INDEL},cutoff=0.5 \
      --plugin dbNSFP,consequence=ALL,${DBNSFP},GERP++_NR,GERP++_RS,CADD_phred,CADD_raw,AlphaMissense_score,ClinPred_pred,ClinPred_score,Polyphen2_HDIV_pred,Polyphen2_HDIV_score,Polyphen2_HVAR_pred,Polyphen2_HVAR_score,MetaLR_pred,SIFT_score,SIFT_pred,MutationTaster_score,MutationTaster_pred,M-CAP_score,M-CAP_pred,clinvar_trait,Ensembl_geneid,LIST-S2_score,LIST-S2_score,DANN_score \
      --plugin dbscSNV,${DBSCSNV} \
      --plugin CADD,snv=${CADD} \
      --plugin NMD \
      --dir_plugins ${PLUGIN_DIR} \
      --fork 16 \
      2>&1 | tee -a "$LOG_FILE"

echo "Variant annotation completed!" | tee -a "$LOG_FILE"
