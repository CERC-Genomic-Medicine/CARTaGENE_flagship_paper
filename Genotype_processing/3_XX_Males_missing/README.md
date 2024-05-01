# XX males missing
This folder contains the code to address the male heterozygotes in position which should be haploid (i.e. non-PAR)

## Parameter
### Input files
Bedfiles - CaG data v.1.1 in plink .bed format (including .bim .fam companion file) **per Genotyping arrays** issued from the [ Reference Alignement](https://github.com/CERC-Genomic-Medicine/CARTaGENE_flagship_paper/tree/main/Genotype_processing/2_Reference_Alignement)  

### Software
- [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020)

### Output files :
Bedfile - CaG data v.1.1 in plink .bed format (including .bim .fam companion file) **per Genotyping arrays**

## Step
### Step 1 - Set HH-missing
1) Split X chromosome with [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020)'s --split-x 'b38' option
2) Set HH to missint with [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020)'s --set-hh-missing option
3) Merge XY back into X with [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020)'s --merge-x option

### Step 2 - Verification
1) Verfication X final region
2) Verification all male status of heterozygote haploid X chromosome and Y chromosome
