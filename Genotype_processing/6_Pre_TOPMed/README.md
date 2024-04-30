# preparation for TOPMed imputation

## Goal :
 - create 2 batch of partially overlapping randomly selected samples
 - create 1 vcf file per Chromosome per said batch

## Input files :
- CARTaGENE_hg38_shared_TOPMed_filtered


## Output files :
- CARTaGENE_hg38_shared_filtered_chr{1..22}_batch{1..2}.vcf.gz
- CARTaGENE_hg38_shared_filtered_chr{1..22}_batch{1..2}.vcf.gz.csi
- CARTaGENE_hg38_shared_filtered_chrX_batch{1..2}.vcf.gz
- CARTaGENE_hg38_shared_filtered_chrX_batch{1..2}.vcf.gz.csi


## Used Software :
-plink/1.9
-bcftools (w/plugins)
