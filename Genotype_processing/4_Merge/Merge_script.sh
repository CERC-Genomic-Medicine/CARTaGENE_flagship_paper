#!/bin/bash

# Script
## Declare Variables

declare -a arrays=("17k" "5300" "4224" "760" "archi")
path=''
## Find common ID

comm -1 -2 <(cut -f 2 ${path}/17k_hg38.bim | sort) <(cut -f 2 ${path}/archi_hg38.bim | sort) > temp1.txt
comm -1 -2 <(cut -f 2 ${path}/4224_hg38.bim| sort) <(cut -f 2 ${path}/5300_hg38.bim | sort) > temp2.txt
comm -1 -2 <(cut -f 2 ${path}/760_hg38.bim | sort) <( sort temp1.txt) > temp3.txt
comm -1 -2 <(sort temp3.txt) <(sort temp2.txt) > temp4.txt
grep -v 'ID ' temp4.txt > tmp

## Remove chrY
grep -v chrY tmp > shared.txt
## Test
wc -l shared.txt # 481,223
## Analysis of duplicated samples

Python3 Analysis_Duplicated_samples.py

## Merge 
for array in "${arrays[@]}" ; do
        plink --bfile ${path}/${array}_hg38 --extract shared.txt  --keep-allele-order --make-bed --output-chr chrMT --out ${array}_hg38_shared
done

for array in "${arrays[@]}" ; do 
        plink --bfile ${path}/${array}_hg38 --keep exclude_dup_sample.txt --missing --out ${array}.total
        plink --bfile ${path}/${array}_hg38_shared --keep exclude_dup_sample.txt --missing --out ${array}.dup
done


mv 760_hg38_shared.bed 760_hg38_shared_dup.bed
mv 760_hg38_shared.bim 760_hg38_shared_dup.bim
mv 760_hg38_shared.fam 760_hg38_shared_dup.fam

plink --bfile 760_hg38_shared_dup --remove exclude_dup_sample.txt --keep-allele-order --make-bed --output-chr chrMT --out 760_hg38_shared

plink --bfile 17k_hg38_shared --keep-allele-order --output-chr chrMT --bmerge 5300_hg38_shared --out merge_temp
plink --bfile merge_temp --keep-allele-order --output-chr chrMT  --bmerge 4224_hg38_shared --out merge_temp_2
plink --bfile merge_temp_2 --keep-allele-order --output-chr chrMT --bmerge archi_hg38_shared --out merge_temp_3
plink --bfile merge_temp_3 --keep-allele-order --output-chr chrMT --bmerge 760_hg38_shared --out CARTaGENE_hg38_shared
rm *temp*
