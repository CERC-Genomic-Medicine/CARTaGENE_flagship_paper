# TOPMed imputation preparation
This folder contains the scripts and documentation used to prepare the dataset for imputation, specifically split into two overlapping batches of 25,000 samples, fill-info tags and split into individuals chromosomes. This is performed solely in preparation to imputation on TOPMed data. This step is reliant on the [Harmonization step](https://github.com/CERC-Genomic-Medicine/CARTaGENE_flagship_paper/tree/main/Genotype_processing/7_Pre_TOPMed/).

## Parameter
### Input files
Bed file - Merged CaG data v.1.1 in plink .bed format (including .bim .fam companion file)  
  
### Used Software
- [bcftools/1.19](https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2)
- [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020)

### Output files : 
- VCFs - CaG dataset in VCF format **per chromosome** and **in two batches**
- index - index of vcfs

## Steps
### Step 1 -  Annotate
Annote variants
1) Convert plink .bed format to vcf format using [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020)
2) Annotate the vcf file using [bcftools/1.19](https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2)'s Plugin fill-tags(https://samtools.github.io/bcftools/howtos/plugin.fill-tags.html)
3) index annotated vcf using [bcftools/1.19](https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2)'s index function

### Step 2 - Split into batches and chromosomes
1) create batches of 25,000 samples (randomly selected)
2) Separate vcf into batches and chromosome using [bcftools/1.19](https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2)'s -s and -r options respectively

