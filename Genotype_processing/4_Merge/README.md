# Liftover
This folder contains the scripts and documentation used to Merge the files from the several genotyping arrays into a single CaG v.1.1 dataset. This folder relies on the output of several previous step. See [previous step](https://github.com/CERC-Genomic-Medicine/CARTaGENE_flagship_paper/tree/main/Genotype_processing/3_XX_Males_missing). 


## Parameter
### Input files
Bedfiles - CaG data v.1.1 in plink .bed format (including .bim .fam companion file) **per Genotyping arrays** issued from the [ hh_missing](This folder contains the scripts and documentation used to liftover CaG v.1.1 dataset. This folder relies on the output of several previous step. See [previous step](https://github.com/CERC-Genomic-Medicine/CARTaGENE_flagship_paper/tree/main/Genotype_processing/3_XX_Males_missing)  

### Used Software
- [python 3.11.5](https://www.python.org/downloads/release/python-3115/)
  - argparse
  - pysam-0.22.0
- [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020)

### Output files :  
Bedfile - CaG data v.1.1 in plink .bed format (including .bim .fam companion file)  
  
## Steps  
While this folder contains multiple scripts.  
### Step 1 - Find common ID  
1) Genotype_processing/4_Merge/Find_shared_variants.py  
### Step 2 - Find duplicated samples (by name)  
1) Genotype_processing/4_Merge/Analysis_Duplicated_samples.py  
### Step 3 - Merge Data Set  
Genotype_processing/4_Merge/Merge_script.sh  
1) remove Y chromosome from shared variants
2) Filter genotyping files for shared variants with [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020)'s extract option
3) Remove duplicated samples from 760 array with [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020)'s --remove option
4) Merge bed file using [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020)'s --bmerge function

