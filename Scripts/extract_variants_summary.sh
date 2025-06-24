#!/bin/bash
# Written by Elijah C. Ugoh
# Date: 21-06-2025
# This bash script uses bcftools and awk to parse the final filtered vcf file after VEP annotation and filteration
# It checks for non-refernce genotypes (GT != 0/0, 0|0, or ./.) and extracts useful variant annotation data (CHRM, POS, etc.), incluing the sample IDs carrying each candidate variant
# This information is then saved in a tab-separated .tsv file for easy reading
# NB: Input vcf can either be compressed or not, and must be located in the same dir as this script (otherwise, modify the input file path accordingly)
    # BCFtools module must be loaded for script to work
    # The fields of interest can be modified

ml GCC/12.3.0 BCFtools/1.18 # load necessary modules

# Extract sample names from the VCF
sample_names=($(bcftools query -l final_filtered.vcf.gz))

# Use bcftools to query CHROM, POS, REF, ALT, per-sample GTs, and CSQ field
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\t%INFO/CSQ\n' final_filtered.vcf.gz | awk -v OFS="\t" -v samples="${sample_names[*]}" '
BEGIN {
    split(samples, sample_arr, " ");
    print "CHROM", "POS", "REF", "ALT", "rsid", "Consequence", "Impact", "Gene", "HGVSc", "HGVSp", "gnomADe_AF", "SIFT", "PolyPhen", "clinvar_clnsig", "CADD_PHRED", "Sample_IDs"
}
{
    # Iterate through each transcript entry in the CSQ
    num_samples = length(sample_arr);
    split($NF, transcripts, ",");  # CSQ field is the last column
    for (t in transcripts) {
        split(transcripts[t], info, "|");
        
        # Assign basic variant info from bcftools query and extract annotated functional info from the CSQ field for each variant
        CHROM = $1;
        POS = $2;
        REF = $3;
        ALT = $4;
        rsid = info[18];
        Consequence     = info[2];
        Impact          = info[3];
        Gene            = info[4];
        HGVSc           = info[11];
        HGVSp           = info[12];
        gnomADe_AF      = info[38];
        SIFT            = info[34];
        PolyPhen        = info[35];
        clinvar_clnsig  = info[85];
        CADD_PHRED      = info[89];

        Sample_IDs = ""; # Initialize an empty string to hold the names of samples carrying the variants

        for (i = 1; i <= num_samples; i++) {
            GT = $(4 + i);   #  check each sample genotype satrting from column 5 onward. 
            if (GT != "0/0" && GT != "0|0" && GT != "./.") { # Only non-reference genotype is appended to the Sample_IDs list above
                Sample_IDs = Sample_IDs sample_arr[i] ",";
            }
        }
        # Remove trailing comma (just added above) after the last sample ID
        sub(/,$/, "", Sample_IDs);
        
       # Return all the values all each extracted fields
        print CHROM, POS, REF, ALT, rsid, Consequence, Impact, Gene, HGVSc, HGVSp, gnomADe_AF, SIFT, PolyPhen, clinvar_clnsig, CADD_PHRED, Sample_IDs;
    }
}' > final_variant_sample_report.tsv # write output to a file
