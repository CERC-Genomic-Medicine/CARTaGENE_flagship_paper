#!/bin/bash

## Declare Variable

path=“/path/to/plink_files/*bim”  # path to CAG PLINK binary format's bim file, each set should contain a .bed .bim .fam file.

## Set HH-missing
for bim in ${path}
do
        plink --bfile ${bim%.*} --split-x 'b38'  --keep-allele-order --make-bed --output-chr chrMT --out ${bim%.*}_Xsplit
        plink --bfile ${bim%.*}_Xsplit --keep-allele-order --set-hh-missing --make-bed --output-chr chrMT --out ${bim%.*}_Xsplit_temp
        plink --bfile ${bim%.*}_Xsplit_temp --merge-x --keep-allele-order --make-bed --output-chr chrMT --out ${bim%.*}_hh
done

## Verfication X region

cat *_Xsplit_temp.hh | grep 'chrX' | cut -f 3 | sort -V | uniq ## head and tail
cat *_Xsplit_temp.hh | grep 'chrY' | cut -f 3 | sort -V | uniq ## head and tail

## Verficiation X male Y female
for bim in ${path}
do
        cat ${bim%.*}.fam | cut -d " " -f 1,5 | grep 2$ | cut -d " " -f 1 | sort | uniq >> female.txt
        cat ${bim%.*}.fam | cut -d " " -f 1,5 | grep 1$ | cut -d " " -f 1 | sort | uniq >> male.txt
done

cat *)Xsplit_temp.hh | grep 'chrX' | cut -f 1 | sort | uniq > X_ind.hh
cat *_Xsplit_temp.hh | grep 'chrY' | cut -f 1 | sort | uniq > Y_ind.hh

grep -v -f male.txt X_ind.hh ## expected : no output
grep -f female.txt X_ind.hh ## expected : no output
grep -f female.txt Y_ind.hh ## expected : no output

####### clean up ####
for bim in ${path}
do
        rm ${bim%.*}_Xsplit_temp.*
        rm ${bim%.*}_hg38_Xsplit*
done
