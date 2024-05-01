#!/bin/bash

## Declaring variable

GRCh38='Homo_sapiens.GRCh38.fa'
path=“/path/to/plink_files/*bim”  # path to build CAG PLINK binary format's bim file (hg38 build version), each set should contain a .bed .bim .fam file. 


## Alignement

for bim in ${path} ; do
        python3 plink2reference.py -b ${bim%} -f ${GRCh38} -o ${bim%.*}                             #Produce documents needed for correct strand-flip and REF/ALT alteration
        plink --bfile ${bim%.*}  --exclude ${bim%.*}.remove.txt --make-bed --output-chr chrMT --out temp     #Removes impossible to adjust SNPs
        plink --bfile temp --flip ${bim%.*}.strand_flip.txt --make-bed --output-chr chrMT --out temp2                            # Strand filp
        plink --bfile temp2 --a1-allele ${bim%.*}.force_a1.txt --make-bed --output-chr chrMT --out ${bim%.*}_aligned         # Forces the REF/ALT Designation
        rm temp*
done

## Deduplication

for bim in ${path} ; do
        plink --bfile ${bim%.*}_aligned  --list-duplicate-vars                                                      #Identifies duplicate variant Position:REF:ALT
        plink --bfile ${bim%.*}_aligned  --freq counts
        python3 plink_dupvar_prioritize.py -d plink.dupvar -c plink.frq.counts -o ${bim%.*}_exclude_dup.txt              #Identifies the duplicate variant with less missingness
        plink --bfile ${bim%.*}_hg38_aligned --exclude ${bim%.*}_exclude_dup.txt --output-chr chrMT --keep-allele-order --make-bed --out ${bim%.*}_aligned_dedup # Removes duplicates
        rm plink.dupvar plink.frq.counts
done

## RENAMING SNPs

for bim in ${path} ; do
        awk '{print $1":"$4":"$5":"$6 , $2}' ${bim%.*}_aligned_dedup.bim > rename.txt                                                   # Create list of old var ID and new var id (chr:pos:REF:ALT)
        plink --bfile ${bim%.*}_aligned_dedup --update-name rename.txt 1 2  --make-bed --keep-allele-order --output-chr chrMT --out ${bim%.*}_renamed     # renames
        rm rename.txt
done

## Verification
for bim in ${path} ; do
        echo ${bim%.*} ; grep chrM ${bim%.*}_renamed.bim | wc -l ; grep chrX ${bim%.*}_renamed.bim | wc -l ; 
        grep chrY ${bim%.*}_renamed.bim | wc -l ; grep chrXY ${bim%.*}_renamed.bim | wc -l ; 
        wc -l ${bim%.*}_renamed.bim
done
