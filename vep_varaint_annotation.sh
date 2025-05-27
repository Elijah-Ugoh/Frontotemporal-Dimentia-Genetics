#!/bin/bash
#SBATCH --job-name=vep_annot
#SBATCH --output=vep_annot_%j.out
#SBATCH --error=vep_annot_%j.err
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=16
#SBATCH --mem=240G
#SBATCH --account=csens2024-3-26
#SBATCH --time=24:00:00
#SBATCH --partition=sens-xl

# Description:
# SLURM script for full variant annotation using Ensembl VEP (v113) in a Singularity container,
# including plugin annotations, consequence filtering, pathogenicity scoring, and final gene-level filtering.
# Final output: bgzipped, indexed VCF with candidate variants.

# ml GCC/13.3.0 Singularity  # Optional. This script uses a locally provided compilation of singularity

set -euo pipefail

# Ensure results are always copied back to the submission dir at any point even if the pipeline fails
cleanup() {
  echo "Job ending. Now copying results..."
  cp -rp results $SLURM_SUBMIT_DIR || echo "Failed to copy output dir"
}
trap cleanup EXIT

# write to standard output
cat $0

start_time=$(date +%s) 

# copy input files and directories to temp dir
date +"%T" && cp -rp vep_data vcf Plugins singularity $NAISS_TMP && date +"%T" 
cp -p vep.sif annotation.sh consequence_filtering.sh impact_MAF_pop_filtering.sh pathogenicity_scoring.sh ftd_gene_filter.sh genotype_quality_filter.sh FTD_genes_HP_0002145.txt $NAISS_TMP

echo "All input files and directories successfully copied to $NAISS_TMP"
cd $NAISS_TMP

# Annotation using Ensembl VEP
# Load necessary modules
ml GCC/12.3.0 BCFtools/1.18
ml GCC/10.2.0 tabixpp/1.1.0 

# Define output directory for all results
mkdir results && export OUTPUT_VCF_DIR="results"

# The vcf file is compressed, indexed with tabix, and merges information for all samples. 
# It also outputs a bgzip compressed file.
echo "==== Starting VEP Annotation ===="
sh annotation.sh

date +"%T" 

# index annotated output 
tabix -p vcf ${OUTPUT_VCF_DIR}/output.annotated.vcf.gz

# consequence filteration
# Filter the annotation output based on selected VEP consequence filters 
sh consequence_filtering.sh

# Compress and tabix index the output
bgzip ${OUTPUT_VCF_DIR}/consequence_filtered.vcf
tabix -p vcf ${OUTPUT_VCF_DIR}/consequence_filtered.vcf.gz

# pathogenicity scoring/prediction using CADD, SIFT, and Polyphen
sh pathogenicity_scoring.sh


# Compress and tabix index
bgzip ${OUTPUT_VCF_DIR}/pathogenicity_scored.vcf
tabix -p vcf ${OUTPUT_VCF_DIR}/pathogenicity_scored.vcf.gz

# Impact (HIGH or MODERATE) + allele frequency (< 0.01) + NFE population filtering
sh impact_MAF_pop_filtering.sh

# Compress and tabix index
bgzip ${OUTPUT_VCF_DIR}/impact_NFE_and_AF_filtered.vcf
tabix -p vcf ${OUTPUT_VCF_DIR}/impact_NFE_and_AF_filtered.vcf.gz

# Remove varaints with poor genotyping calls
sh genotype_quality_filter.sh

# Filter variants that are close to known FTD genes
sh ftd_gene_filter.sh

# Compress and tabix index the final out
bgzip ${OUTPUT_VCF_DIR}/candidate_variants.vcf
tabix -p vcf ${OUTPUT_VCF_DIR}/candidate_variants.vcf.gz

end_time=$(date +%s)

runtime=$((end_time - start_time))

echo "Total script run time is: $((runtime / 60)) minutes, $((runtime % 60)) seconds" 
