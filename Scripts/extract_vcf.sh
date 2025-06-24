#!/bin/bash

# Set thresholds
THRESHOLD_EXPANSION=30
THRESHOLD_PATHOGENIC=100

# Header
echo "SampleID CHR POS LocusId Genotype CI REPCN InRepeatReads Coverage Flag" > C9ORF72_vcf_summary.txt

# Loop through VCFs
for vcf in STR_analysis_results/*.vcf.gz; do
  sample=$(basename "$vcf" .vcf.gz)

  # Extract line containing C9ORF72 info
  line=$(zgrep "C9ORF72" "$vcf")
  
  if [ -z "$line" ]; then
    echo "$sample NA NA NA NA NA NA NA NA NoCall" >> C9ORF72_vcf_summary.txt
    continue
  fi

  chr=$(echo "$line" | cut -f1)
  pos=$(echo "$line" | cut -f2)
  info=$(echo "$line" | cut -f8)
  format=$(echo "$line" | cut -f9)
  values=$(echo "$line" | cut -f10)

  # Parse INFO
  locusid=$(echo "$info" | sed -n 's/.*VARID=\([^;]*\).*/\1/p')

  # Extract fields using FORMAT and corresponding VALUES
  gt=$(echo "$values" | awk -F':' '{print $1}')
  repci=$(echo "$values" | awk -F':' '{print $4}')
  repcn=$(echo "$values" | awk -F':' '{print $3}')
  adir=$(echo "$values" | awk -F':' '{print $7}')
  lc=$(echo "$values" | awk -F':' '{print $8}')

  # Determine flag based on largest allele
  allele1=$(echo "$gt" | cut -d'/' -f1)
  allele2=$(echo "$gt" | cut -d'/' -f2)
  allele1_val=$(echo "$repcn" | cut -d'/' -f1)
  allele2_val=$(echo "$repcn" | cut -d'/' -f2)
  max_allele=$(echo -e "$allele1_val\n$allele2_val" | sort -nr | head -n1)

  if [[ "$max_allele" -ge "$THRESHOLD_PATHOGENIC" ]]; then
    flag="Pathogenic Expansion"
  elif [[ "$max_allele" -ge "$THRESHOLD_EXPANSION" ]]; then
    flag="Intermediate Expansion"
  else
    flag="Normal"
  fi

  # Output
  echo "$sample $chr $pos $locusid $gt $repci $repcn \"$adir\" $lc $flag" >> C9ORF72_vcf_summary.txt
done
