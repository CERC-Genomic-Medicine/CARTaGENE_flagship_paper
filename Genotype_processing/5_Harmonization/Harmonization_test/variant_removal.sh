#!/bin/bash
## Parameters
### Input Files
CaG="Path/to/[file.bim]" #Merged CaG data v.1.1 in plink .bed format (including .bim .fam companion file)
LRT_fail="Path/to/[file]" # List of variant with LRT pvalue below threshold


plink -bfile ${CaG%.*} --keep-allele-order  --output-chr chrMT --exclude ${LRT_fail} --make-bed --out ${CaG%.*}_AF_filtered
plink -bfile ${CaG%.*} --a1-allele ${CaG} 6 2  --output-chr chrMT --exclude ${LRT_fail} --recode vcf-iid bgz --out ${CaG%.*}_AF_filtered
