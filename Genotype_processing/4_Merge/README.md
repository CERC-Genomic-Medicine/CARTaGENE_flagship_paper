# Merge
This folder contains the scripts and documentation used to Merge the files from the several genotyping arrays into a single CaG v.1.1 dataset. This folder relies on the output of several previous step. See [previous step](https://github.com/CERC-Genomic-Medicine/CARTaGENE_flagship_paper/tree/main/Genotype_processing/3_XX_Males_missing). 


## Parameter
### Input files
Bedfiles - CaG data v.1.1 in plink .bed format (including .bim .fam companion file) **per Genotyping arrays** issued from the [ hh_missing](This folder contains the scripts and documentation used to liftover CaG v.1.1 dataset. This folder relies on the output of several previous step. See [previous step](https://github.com/CERC-Genomic-Medicine/CARTaGENE_flagship_paper/tree/main/Genotype_processing/3_XX_Males_missing)  

### Used Software
- [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020)

### Output files :  
Bed file - CaG data v.1.1 in plink .bed format (including .bim .fam companion file)  
  
## Steps  
While this folder contains multiple scripts.  
### Step 1 - Find & remove duplicated samples (by name)  
1) find IID with count > 1
2) Order files per SNP count (smaller to larger)
3) remove IID from the less covered array to the most covered. Remove IID only in one array, with [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020)'s --remove option
### Step 2 - Find common ID  
1) Find SNP id occuring the same amount of time as the number of array
2) Filter genotyping files for shared variants with [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020)'s extract option with [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020)'s extract option
### Step 3 - Merge Data Set  
1) Merge bed file using [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020)'s --bmerge function

