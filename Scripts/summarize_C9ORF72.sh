#!/bin/bash

# Set thresholds
THRESHOLD_EXPANSION=30
THRESHOLD_PATHOGENIC=100

# Header
echo "SampleID AlleleCount LocusId Genotype CI InRepeatReads Coverage Flag" > C9ORF72_json_summary.txt

# Loop through JSONs
for json in STR_analysis_results/*.json; do
  sample=$(basename "$json" .json)

  # Extract values with jq
  allele_count=$(jq -r '.LocusResults.C9ORF72.AlleleCount' "$json")
  locusid=$(jq -r '.LocusResults.C9ORF72.LocusId' "$json")
  genotype=$(jq -r '.LocusResults.C9ORF72.Variants.C9ORF72.Genotype' "$json")
  ci=$(jq -r '.LocusResults.C9ORF72.Variants.C9ORF72.GenotypeConfidenceInterval' "$json")
  inrepeat=$(jq -r '.LocusResults.C9ORF72.Variants.C9ORF72.CountsOfInrepeatReads' "$json")
  coverage=$(jq -r '.LocusResults.C9ORF72.Coverage' "$json")

  # Calculate flag
  allele1=$(echo "$genotype" | cut -d'/' -f1)
  allele2=$(echo "$genotype" | cut -d'/' -f2)
  max_allele=$(echo -e "$allele1\n$allele2" | sort -nr | head -n1)

  if [ "$max_allele" -ge "$THRESHOLD_PATHOGENIC" ]; then
    flag="Pathogenic Expansion"
  elif [ "$max_allele" -ge "$THRESHOLD_EXPANSION" ]; then
    flag="Intermediate Expansion"
  else
    flag="Normal"
  fi

  # Write to CSV
  echo "$sample $allele_count $locusid $genotype $ci \"$inrepeat\" $coverage $flag" >> C9ORF72_json_summary.txt
done
