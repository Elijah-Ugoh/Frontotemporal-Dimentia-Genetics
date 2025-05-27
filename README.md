# Investigating the Genetic Causes of Fronto-temporal Dimentia (FTD)

This is a short pipeline that analyzes short reads whole-genome sequencing data from Fronto-temporal Dimentia (FTD) patients, starting with quality control to remove poor quality reads, mapping to a reference geneome, variant calling, annotation, and candidate variant filtering. 


## Annotation in Ensemble Variant Effect Predictor (VEP)
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
Specifically, we use: CADD, dbNSFP, dbscSNV, NMD, SpliceAI. These plugins all have execution scipts, which are alaready installed in the containerized VEP download, but the required  plugin files must be downloaded also. 

Summary of all pluggins and usage:
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
sudo mount -t drvfs h: /mnt/h -o uid=$(id -u $USER),gid=$(id -g $USER),metadata
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

Candidate variants and saved in a txt file and compared in databases like VarSome, ClinVar, UniProt, and gnomAD browser.

This process can also be repeated using less stringent filters, such as a MAF of 5% instead of 1% and 0.3 for SpliceAI scores, etc. 



To read: [ACMG Criteria](https://www.nature.com/articles/gim201530)