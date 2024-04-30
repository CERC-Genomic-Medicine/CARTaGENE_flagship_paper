#!/bin/bash

## Declaring variable

GRCh38='Homo_sapiens.GRCh38.fa'
declare -a arrays=("17k" "5300" "4224" "760" "archi")
path='' # Obscured for security

## Alignement

for array in "${arrays[@]}" ; do
        python3 plink2reference.py -b ${path}/${array}_Hg38_Xmerged.bim -f ${GRCh38} -o ${array}                                #Produce documents needed for correct strand-flip and REF/ALT alteration
        plink --bfile ${path}/${array}_Hg38_Xmerged  --exclude ${array}.remove.txt --make-bed --output-chr chrMT --out temp     #Removes impossible to adjust SNPs
        plink --bfile temp --flip ${array}.strand_flip.txt --make-bed --output-chr chrMT --out temp2                            # Strand filp
        plink --bfile temp2 --a1-allele ${array}.force_a1.txt --make-bed --output-chr chrMT --out ${array}_hg38_aligned         # Forces the REF/ALT Designation
        rm temp*
done

## Deduplication

for array in "${arrays[@]}" ; do
        plink --bfile ${array}_hg38_aligned  --list-duplicate-vars                                                      #Identifies duplicate variant Position:REF:ALT
        plink --bfile ${array}_hg38_aligned  --freq counts
        python3 plink_dupvar_prioritize.py -d plink.dupvar -c plink.frq.counts -o ${array}_exclude_dup.txt              #Identifies the duplicate variant with less missingness
        plink --bfile ${array}_hg38_aligned --exclude ${array}_exclude_dup.txt --output-chr chrMT --keep-allele-order --make-bed --out ${array}_hg38_aligned_dedup # Removes duplicates
        rm plink.dupvar plink.frq.counts
done

## RENAMING SNPs

for array in "${arrays[@]}" ; do
        awk '{print $1":"$4":"$5":"$6 , $2}' ${array}_hg38_aligned_dedup.bim > rename.txt                                                   # Create list of old var ID and new var id (chr:pos:REF:ALT)
        plink --bfile ${array}_hg38_aligned_dedup --update-name rename.txt 1 2  --make-bed --keep-allele-order --output-chr chrMT --out ${array}_hg38_renamed     # renames
        rm rename.txt
done

## Verification
for array in "${arrays[@]}" ; do
        echo ${array} ; grep chrM ${array}_hg38_renamed.bim | wc -l ; grep chrX ${array}_hg38_renamed.bim | wc -l ; 
        grep chrY ${array}_hg38_renamed.bim | wc -l ; grep chrXY ${array}_hg38_renamed.bim | wc -l ; 
        wc -l ${array}_hg38_renamed.bim
done
