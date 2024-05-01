#!/bin/bash

# Script
## Declare Variables

path=“/path/to/plink_files/*bim”  # path to CAG PLINK binary format's bim file, each set should contain a .bed .bim .fam file.

## Remove chrY
grep -v chrY shared.txt > shared_variants.txt

## Merge 
for bim in ${path} ; do
        plink --bfile ${bim%.*} --extract shared.txt  --keep-allele-order --make-bed --output-chr chrMT --out ${bim%.*}_shared
done

## 760 genotype array was selected since it had less genotyping position initially.

mv 760_hg38_shared.bed 760_hg38_shared_dup.bed
mv 760_hg38_shared.bim 760_hg38_shared_dup.bim
mv 760_hg38_shared.fam 760_hg38_shared_dup.fam

plink --bfile 760_hg38_shared_dup --remove exclude_dup_sample.txt --keep-allele-order --make-bed --output-chr chrMT --out 760_hg38_shared

files=($(ls ${path}))
base_file="${files[0]%.*}_shared"
unset files[0] # pop out the first file

output="merge_temp"

# Loop through the remaining files and merge them one by one
for file in "${files[@]}"
do
    # Determine the next output name
    next_output="${output}_next"

    # Run plink to merge the current base with the next file
    plink --bfile $base_file --keep-allele-order --output-chr chrMT --bmerge ${file%.*}_shared --out $next_output

    # Update base_file to the new merged file for the next iteration
    base_file=$next_output
    output=$next_output
done

mv ${output}.bed CARTaGENE_hg38_shared.bed
mv ${output}.bim CARTaGENE_hg38_shared.bim
mv ${output}.fam CARTaGENE_hg38_shared.fam

rm *temp*
