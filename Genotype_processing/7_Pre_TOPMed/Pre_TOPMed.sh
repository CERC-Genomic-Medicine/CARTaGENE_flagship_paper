#!/bin/bash
# Parameters
## File
### Input file
CaG="Path/to/file/*.bim"        # path to CAG PLINK binary format's bim file with AF concordant with TOPMed, set should contain a .bed .bim .fam file.
### Variables
Batch_length=25000              # Length of batches (limited by TOPMed imputation server)

plink --bfile ${CaG%.*}  --a1-allele ${CaG} 6 2 --output-chr chrMT --export vcf-iid bgz --out ${CaG%.*}_temp ; done
bcftools index ${CaG%.*}_temp.vcf.gz
bcftools +fill-tags ${CaG%.*}_temp.vcf.gz -Oz -o ${CaG%.*}_filtered.vcf.gz
rm *temp.vcf.gz
bcftools index ${CaG%.*}_filtered.vcf.gz
bcftools query -l ${CaG%.*}_filtered.vcf.gz > samples
shuf samples > samples_random.txt
head -n ${Batch_length} samples_random.txt > batch1.txt
tail -n ${Batch_length} samples_random.txt > batch2.txt

for chromosome in chr{1..22} chrX ; do 
  bcftools view -r ${chromosome} ${CaG%.*}_filtered.vcf.gz -S batch1.txt  -Oz -o ${CaG%.*}_filtered_${chromosome}_batch1.vcf.gz; 
  bcftools view -r ${chromosome} ${CaG%.*}_filtered.vcf.gz -S batch2.txt  -Oz -o ${CaG%.*}_filtered_${chromosome}_batch2.vcf.gz; 
  bcftools index ${CaG%.*}_filtered_${chromosome}_batch1.vcf.gz ; 
  bcftools index ${CaG%.*}_filtered_${chromosome}_batch2.vcf.gz ; 
done
