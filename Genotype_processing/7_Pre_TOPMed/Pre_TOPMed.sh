#!/bin/bash

plink --bfile CARTaGENE_hg38_shared_TOPMed_filtered  --a1-allele CARTaGENE_hg38_shared_TOPMed_filtered.bim 6 2 --output-chr chrMT --export vcf-iid bgz --out CARTaGENE_hg38_shared_filtered_temp ; done
bcftools index CARTaGENE_hg38_shared_filtered_temp.vcf.gz
bcftools +fill-tags CARTaGENE_hg38_shared_filtered_temp.vcf.gz -Oz -o CARTaGENE_hg38_shared_filtered.vcf.gz
rm *temp.vcf.gz
bcftools index CARTaGENE_hg38_shared_filtered.vcf.gz
bcftools query -l CARTaGENE_hg38_shared_filtered.vcf.gz > samples
shuf samples > samples_random.txt
head -n 25000 samples_random.txt > batch1.txt
tail -n 25000 samples_random.txt > batch2.txt

for  i in chr{1..22} chrX ; do 
  bcftools view -r $i CARTaGENE_hg38_shared_filtered.vcf.gz -S batch1.txt  -Oz -o CARTaGENE_hg38_shared_filtered_"$i"_batch1.vcf.gz; 
  bcftools view -r $i CARTaGENE_hg38_shared_filtered.vcf.gz -S batch2.txt  -Oz -o CARTaGENE_hg38_shared_filtered_"$i"_batch2.vcf.gz; 
  bcftools index CARTaGENE_hg38_shared_filtered_"$i"_batch1.vcf.gz ; 
  bcftools index CARTaGENE_hg38_shared_filtered_"$i"_batch2.vcf.gz ; 
done
