# Investigating the Genetic Causes of Fronto-temporal Dimentia (FTD)

This is a short pipeline that analyzes short reads whole-genome sequencing data from Fronto-temporal Dimentia (FTD) patients, starting with quality control to remove poor quality reads, mapping to a reference geneome, variant calling, annotation, and candidate variant filtering. 

## Starting with Installing GATK
The containerized version is used and pulled from Docker hub (This is necessary since we want to run GATK in a Docker container within a cluster).
- First, Docker Engine is installed and started by following the guide [here](https://docs.docker.com/engine/install/ubuntu/).
  - Confirm Docker Engine was properly installed and ready for use: ```docker --version```.
- Next, GATK is installed by following the install guide [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035889991--How-to-Run-GATK-in-a-Docker-container) for the containerized version of the toolkit.
- Convert the pulled Docker image to Singularity for use on the cluster. ```singularity build gatk_4.6.2.0.sif docker://broadinstitute/gatk:4.6.2.0``` (requires Singualrity to be installed).
  - ```sudo apt update && sudo apt install singularity-container```. 
  - Confirm all is fine afterwards ```singularity exec gatk_4.6.2.0.sif gatk --help```

NB: The version of GATK used can be changed. At the time of this work, v/4.6.2.0 was used. 

### Pre-processing of BAM files
- The variant discovery process starts with a prelimary preprocessing of the the data. The raw FASTQ file per sample has been cleaned and the sequence reads mapped to the reference genome to produce files in BAM format and sorted by coordinates (.bai files). 
- Next, GATK best practices recommend marking duplicates to mitigate biases introduced by data generation steps such as PCR amplification. This has been performed on the data.

- Finally, we re-calibrate the base quality scores (BQSR), as the variant calling process relies heavily on the quality scores assigned to the individual base calls in each sequence read to produce the best quality results.
- - GATK recommends using common, validated polymorphic sites for BQSR to mask true variants (so GATK doesn’t mistakenly treat them as errors).

For hg38, we used Homo_sapiens_assembly38.dbsnp138.vcf.gz from dbSNP and Mills + 1000G indels, both downloadable from [Genomics Public Data](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/).

```bash
mkdir GATK_Analysis
mkdir GATK_Anlayis/BaseRecalibrator
mkdir GATK_Anlayis/BaseRecalibrator/recal_data
```
## Generate VCF Files



## SNV and Indel Variant Discovery

### Annotation in Ensemble Variant Effect Predictor (VEP)
[VEP](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0974-4) determines the effect of variants (SNPs, insertions, deletions, CNVs or structural variants) on genes, transcripts, and protein sequence, as well as regulatory regions.

The command line version  of VEP is used locally, as it is faster and extendable. [Installation instructions](https://www.ensembl.org/info/docs/tools/vep/index.html) are below:

- Since a HPC is used, we install a conternarized (Singularity) version of VEP. [Singularity](https://sylabs.io/singularity/) is first (if not already installed) installed following instructions [here](https://github.com/sylabs/singularity/blob/main/INSTALL.md): 
- VEP is then installed with Singularity using the VEP Docker image:

```bash
singularity pull --name vep.sif docker://ensemblorg/ensembl-vep

# Once done, create a directory to store to store the VEP cache and fasta file 
mkdir vep_data
# Install VEP in the dir. This directory is pointed to when hen running VEP via Singularity
singularity exec vep.sif INSTALL.pl -c vep_data -a cf -s homo_sapiens -y GRCh38

# Now, VEP is ready for use via Singularity. 
```

VEP requires a VCF file for variant annotation. 

After installations, run ```singularity exec vep.sif vep --dir vep_data/ --help``` to see usage option.

On the HPC, make a new dir where analysis input and output files wil be stored.

```bash
cd gr14
cd home/elugoh/Documents/
mkdir genetic_analysis
```

The input vcf file, the ```vep.sif```, ```Singularity``` compilation, and ```vep_data``` containing VEP cache and fasta files are all moved to the ```genetic_analysis``` dir on the HPC from where VEP can be run using a SLURM job script. 

```bash
# Sample command line options for running VEP

SIF_IMAGE="vep.sif"
VEP_DATA="vep_data"
INPUT_VCF="vcf/landqvist.postVQSR.unfiltered.vcf.gz"
OUTPUT_VCF="annotation_output_dir/"

# Run VEP with Singularity
singularity exec \
  --bind $(pwd)/${VEP_DATA_DIR}:/opt/vep/.vep \
  ${SIF_IMAGE} \
  vep --dir /opt/vep/.vep \
      --cache --offline \
      --species homo_sapiens \
      --assembly GRCh38 \
      --format vcf \
      --input_file ${INPUT_VCF} \
      --output_file ${OUTPUT_VCF}/output.vcf.gz \
      --vcf \
      --force_overwrite \
      --hgvs \
```

It is ideal to use a compressed vcf file for processing speed. If compressed, it should be tabix-indexed as well.

```bash
tabix -p vcf vcf/landqvist.postVQSR.unfiltered.vcf.gz # Allows genomic coordinates to be read quickly and translated into file offsets in the vcf file
```

In the actual annotation, we make use of multiple [VEP plugins](https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html) to optimize the annotations. 
Specifically, we use: CADD, dbNSFP, dbscSNV, NMD, SpliceAI. These plugins all have execution scripts, which are alaready installed in the containerized VEP download, but the required  plugin files must be downloaded also. 

Summary of all plugins and usage:
------------------------------------
| Plugin | Needs external data file             | Function | 
|--------|--------------------------------------|----------| 
| dbNSFP | ✅ Download the dbNSFP tabix-indexed file ([dbNSFP5.1a](https://sites.google.com/site/jpopgen/dbNSFP)) | Pathogenicity prediction by retrieving data for missense variants from the tabix-indexed dbNSFP file |	
| CADD | ✅ Requires CADD score files ([All possible SNVs and indels in GRCh38](https://cadd.gs.washington.edu/download)) | Pathogenicity prediction by retrieving CADD scores for variants from one or more tabix-indexed CADD data files | 
|SpliceAI | ✅ Requires precomputed SpliceAI scores for indels and SNVs from [Illumina BaseSpace](https://basespace.illumina.com/analyses/194103939/files) | Splicing variants predictions by retrieving pre-calculated annotations from SpliceAI to predict splice junctions from an arbitrary pre-mRNA transcript sequence. Annotates genetic variants with their predicted effect on splicing |
| dbscSNV | ✅ Needs [dbscSNV1.1.gz](http://www.liulab.science/dbscsnv.html) indexed by tabix | Splicing predictions by retrieving data for splicing variants from the tabix-indexed dbscSNV file |
| NMD	| ❌ Does not require external data, calculates from transcript models |  Transcript annotation by predicting if a variant allows the transcript escape nonsense-mediated mRNA decay based on certain rules |
------------------------------------


- Upload these plugin files to the cluster in a separate dir, ```Plugins```. 
- In case of less local disk space, download the plugins using an external drive, mount the drive, and then transfer the files directly to the HPC.

```bash
# Navigate to the /mnt dir on wsl
sudo mkdir /mnt/d
sudo mount -t drvfs d: /mnt/d -o uid=$(id -u $USER),gid=$(id -g $USER),metadata
sudo umount mnt/d # after use 
```

NB: For dbNSFP, since the tabixed file for ```dbNSFP51a_grch38-003.gz``` was not provide, we create it as follows: 

```bash
ml GCC/10.2.0 tabixpp/1.1.0 
tabix -s 1 -b 2 -e 2 dbNSFP51a_grch38-003.gz

# Tip: check the file field headers to see which pathogencity prediction parameters to use, as using the "ALL" option generates an unnecessarily large amount of data

zcat dbNSFP51a_grch38-003.gz | head -n 1 | tr "\t" "\n" | nl
```
Aso for dbscSNV, the process is similar after downloading the plugon file:

```bash
unzip dbscSNV1.1.zip
head -n1 dbscSNV1.1.chr1 > h

cat dbscSNV1.1.chr* | grep -v ^chr | sort -k5,5V -k6,6n | cat h - | awk '$5 != "."' | bgzip -c > dbscSNV1.1_GRCh38.txt.gz # the file is bgzipped and sorted by chromosome and position

tabix -s 5 -b 6 -e 6 -c c dbscSNV1.1_GRCh38.txt.gz # tabix requires the processed file
```
### 1. Annotation and variant filtering for pathogenicity
Once the files are downloaded and ready, the consolidated vep annotation command is run in the ```annotation.sh``` script. This produces a VCF file with rich functional data. 

Using the ```vep_filter``` function, the annotation output is subjected to rigorous filtering steps that systematically prioritizes candidate variants by retaining only rare, high-quality variants with strong predicted deleteriousness or clinical significance, based on standard thresholds for pathogenicity scores (including SIFT, PolyPhen, CADD, ClinVar, SpliceAI, and others), consequence, genotype quality, and population allele frequency (<1% and <5%).

This is done using a SLURM job script. This pipeiline automatically uses the ```annotation.sh```,  ```consequence_filtering.sh```, ```pathogenicity_scoring.sh```, ```impact_MAF_pop_filtering.sh```, ```genotype_quality_filter.sh```, ```ftd_gene_filter.sh``` scripts to process the input VCF file.

The final output is a bgzipped and tabix-indexed VCF file containing candidate variants, which have also been filtered based on [known FTD genes](https://hpo.jax.org/browse/term/HP:0002145). 

```bash
# The file must be processed as so after download before used
cut - f 2 genes_for_HP_0002145 > FTD_genes.txt | sed -i 's/ //g' 
```

```bash
sbatch vep_varaint_annotation.sh
```

### 2. Applied pathogenicity prediction thresholds
Summary pathogenicity thresholds applied.

| **Predictor**           | **Field (from dbNSFP plugin)**        | **Interpretation**                          | **Recommended Filter**                            |
|-------------------------|---------------------------------------|----------------------------------------------|---------------------------------------------------|
| **SIFT**                | `SIFT_pred` / `SIFT_score`            | `deleterious` or score < 0.05                | `SIFT_pred is deleterious` or `SIFT_score < 0.05` |
| **PolyPhen2**           | `Polyphen2_HDIV_pred` / `_score`      | `probably_damaging` (best), `possibly_damaging` | `Polyphen2_HDIV_score > 0.8`                      |
| **CADD**                | `CADD_phred`, `CADD_raw`              | Higher is worse; > 20 is deleterious         | `CADD_phred > 20`                                 |
| **GERP++**              | `GERP++_RS`                            | Conservation score; > 4 is conserved         | `GERP++_RS > 4`                                   |
| **M-CAP**               | `M-CAP_score`, `M-CAP_pred`           | Pathogenic if score > 0.025                  | `M-CAP_score > 0.025`                             |
| **MetaLR**              | `MetaLR_pred`, `MetaLR_score`         | `D` (damaging), score > 0.5 recommended      | `MetaLR_score > 0.5`                              |
| **MutationTaster**      | `MutationTaster_pred`, `_score`       | `A` or `D` = disease-causing                 | `MutationTaster_pred is A or D`                  |
| **DANN**                | `DANN_score`                          | > 0.9 suggested for deleterious              | `DANN_score > 0.9`                                |
| **dbscSNV (splice)**    | `ada_score`, `rf_score`               | > 0.6 or > 0.8 often used                    | `ada_score > 0.6 and rf_score > 0.6`              |
| **SpliceAI**            | `SpliceAI_DS_*`                       | Score > 0.5 recommended                      | `SpliceAI_DS_AG > 0.5 or SpliceAI_DS_AL > 0.5` etc. |

Nominated candidate variants are saved in a txt file and compared in databases like VarSome, ClinVar, UniProt, and gnomAD browser.

```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ\n' final_filtered.vcf.gz | awk -F '\t' '
BEGIN {
    OFS = "\t";
    print "CHROM", "POS", "REF", "ALT", "rsid", "Consequence", "Impact", "Gene", "HGVSc", "HGVSp", "gnomADe_AF", "SIFT", "PolyPhen", "clinvar_clnsig", "CADD_PHRED"
}
{
    split($5, transcripts, ",");
    for (i in transcripts) {
        split(transcripts[i], info, "|");
        CHROM = $1;
        POS = $2;
        rsid = info[18];
        REF = $3;
        ALT = $4;
        Consequence     = info[2];
        Impact          = info[3];
        Gene            = info[4]; HGVSc = info[11]; HGVSp = info[12]; gnomADe_AF = info[38];
        SIFT            = info[34];
        PolyPhen        = info[35];
        clinvar_clnsig  = info[85]; CADD_PHRED = info[89];
        print CHROM, POS, REF, ALT, rsid, Consequence, Impact, Gene, HGVSc, HGVSp, gnomADe_AF, SIFT, PolyPhen, clinvar_clnsig, CADD_PHRED;
    }
}' | uniq > gene_list.txt
```

This process can also be repeated using less stringent filters, such as a MAF of 5% instead of 1% and 0.3 for SpliceAI scores, etc. 
The annotation can also be streamlined by supplying the --pick or --flag_pick flag in ```annotation.sh``` to tell VEP to report only the most deleterious or clinically significant transcript.  


## Short Tandem Repeats Analysis
ExpansionHunter is used to genotype the WGS data for STRs. This tool is designed for targeted genotyping of short tandem repeats and flanking variants. It operates by performing a targeted search through a BAM/CRAM file for reads that span, flank, and are fully contained in each repeat.

ExpansionHunter requires the binary alignment/map (BAM) file. 

Below, we install the ExpansionHunter tool by extracting the compield binary. 

```bash
# Download the compiled executable from GitHub into the project dir
wget https://github.com/Illumina/ExpansionHunter/releases/download/v5.0.0/ExpansionHunter-v5.0.0-linux_x86_64.tar.gz
tar -xzvf ExpansionHunter-v5.0.0-linux_x86_64.tar.gz # Extract the files

# ExpansionHunter is now ready to use
```
The ```str_analysis.sh``` script is used to analyse the samples for repeat expansion units. 
```bash
# Finding Repeat Expansions.
# All the sorted bam files are supplied as inputs. Exapnsion hunter returns a .vcf, .json, and a re-alinged bam file for each sample. 
for bam in Data/bam/*.bam; do
    sample_name=$(basename "$bam" .dedup.sorted.bam)

    ./ExpansionHunter-v5.0.0-linux_x86_64/bin/ExpansionHunter \
        --reads "$bam" \
        --reference Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
        --variant-catalog ExpansionHunter-v5.0.0-linux_x86_64/variant_catalog/hg38/variant_catalog.json \
        --output-prefix "STR_analysis_results/${sample_name}" \
        --analysis-mode streaming \
        --log-level info
done
```
Summary information for nominated STR variants can be extracted from either the .json or .vcf outputs and saved to a txt file.

```bash
./extract_vcf.sh
./summarize_C9ORF72.sh
```
### Inspecting the Read Alignments in Genomic Regions Containing Repeat Expansions
The graph realigned reads output generated by ExpansionHunter can be further visualized in REViewer to  inspect the quality of the alignments of the reads used to genotype the repeat regions. REViewer was also developed by the authors of ExpansionHunter is a tool for visualizing alignments of reads in regions containing tandem repeats. 

REViewer requires:
- the BAMlet with graph-realigned reads generated by ExpansionHunter
- the corresponding variant catalog used for the genotyping
- VCF file generated by ExpansionHunter, and 
- a FASTA file with reference genome (hg38)

```bash
# The Linux binary for the latest release of REViewer is download thus
wget https://github.com/Illumina/REViewer/releases/download/v0.2.7/REViewer-v0.2.7-linux_x86_64.gz
gunzip REViewer-v0.2.7-linux_x86_64.gz # Decompress the file
```
Usage:

Sort and index the  the BAMlet generated by ExpansionHunter:
```bash
for bam in STR_analysis_results/*.bam; do
  sorted_bam="${bam%.bam}_sorted.bam"
  samtools sort "$bam" -o "$sorted_bam"
  samtools index "$sorted_bam"
done
```
Because of contig naming inconsitency with different reference genome databses, REViewer may return an error. To prevent this, a compatible hg38 database if downloaded from GATK's Google Cloud Bucket of [useful human genomic resources](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/). 

```bash
wget https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict?inv=1&invt=Ab0XXQ
wget https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta?inv=1&invt=Ab0XXQ
wget https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai?inv=1&invt=Ab0XXQ
```

REViewer can then be run to generate .svg files, which can be viewed in a browser.
```bash
# standing in the genetic_analysis dir

for realigned_bam in STR_analysis_results/*_realigned_sorted.bam; do
  sample_name=$(basename "$realigned_bam" _realigned_sorted.bam)
  ./REViewer-v0.2.7-linux_x86_64 \
    --reads "$realigned_bam" \
    --vcf "STR_analysis_results/${sample_name}.vcf" \
    --reference resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
    --catalog ExpansionHunter-v5.0.0-linux_x86_64/variant_catalog/hg38/variant_catalog.json \
    --locus C9ORF72 \
    --output-prefix "REV_outputs/${sample_name}"
done
```