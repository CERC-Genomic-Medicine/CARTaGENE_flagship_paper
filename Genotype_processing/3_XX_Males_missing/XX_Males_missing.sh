#!/bin/bash

## Declare Variable
declare -a arrays=("17k" "5300" "4224" "760" "archi")
path='Mooser_433651_genotypes_hg38_ref'

## Set HH-missing
for array in "${arrays[@]}"
do
        plink --bfile ${path}/${array}_hg38_renamed --split-x 'b38'  --keep-allele-order --make-bed --output-chr chrMT --out ${array}_hg38_Xsplit
        plink --bfile ${array}_hg38_Xsplit --keep-allele-order --set-hh-missing --make-bed --output-chr chrMT --out ${array}_hg38_Xsplit_temp
        plink --bfile ${array}_hg38_Xsplit_temp --merge-x --keep-allele-order --make-bed --output-chr chrMT --out ${array}_hg38
done

## Verfication X region

cat *_hg38_Xsplit_temp.hh | grep 'chrX' | cut -f 3 | sort -V | uniq ## head and tail
cat *_hg38_Xsplit_temp.hh | grep 'chrY' | cut -f 3 | sort -V | uniq ## head and tail

## Verficiation X male Y female

cat ${path}/*_hg38_renamed.fam | cut -d " " -f 1,5 | grep 2$ | cut -d " " -f 1 | sort | uniq > female.txt
cat ${path}/*_hg38_renamed.fam | cut -d " " -f 1,5 | grep 1$ | cut -d " " -f 1 | sort | uniq > male.txt

cat *_hg38_Xsplit_temp.hh | grep 'chrX' | cut -f 1 | sort | uniq > X_ind.hh
cat *_hg38_Xsplit_temp.hh | grep 'chrY' | cut -f 1 | sort | uniq > Y_ind.hh

grep -v -f male.txt X_ind.hh ## expected : no output
grep -f female.txt X_ind.hh ## expected : no output
grep -f female.txt Y_ind.hh ## expected : no output

####### clean up ####
for i in archi 4224 5300 760 17k
do
rm "$i"_hg38_Xsplit_temp.*
rm "$i"_hg38_Xsplit*
done
