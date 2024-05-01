#!/bin/bash

path=“/path/to/plink_files/*bim”  # path to CAG PLINK binary format's bim file, each set should contain a .bed .bim .fam file.

mkdir logs

# Step 1 : Merge pseudo-autosomal regions (XY to X)
for  bim in ${path} archi; do 
        plink -bfile ${bim%.*} --merge-x --make-bed --output-chr MT --out ${bim%.*}_XYmerged
done

# Step 2 Liftover 
for bim in *_XYmerged.bim; do 
        awk '{print "chr" $1, $4-1, $4, $2}' ${bim} > ${bim%.*}_bedfile           ## Produces a bed file of all variants
        liftOver ${bim%.*}_bedfile hg19ToHg38.over.chain.gz ${bim%.*}_mapfile ${bim%.*}_unmappedfile   ## Produces Map file and unmapped variant file
        grep ^chr[0-9A-Za-z]*_ ${bim%.*}_mapfile | cut -f 4 > ${bim%.*}_excludefile    ## Identifies Alternate contig mod A-Z a-z
        grep -v '^#' ${bim%.*}_unmappedfile | cut -f 4 >> ${bim%.*}_excludefile              ## Total list of variant to be excluded
        grep -v ${bim%.*}_excludefile ${bim%.*}_mapfile > ${bim%.*}_mapfile_final          ## Remove excluded variant

        plink --bfile ${bim%.*}_XYmerged \
        --exclude ${bim%.*}_excludefile \
        --make-bed \
        --not-chr MT \
        --out tmp_${bim%.*}

        plink --bfile tmp_${bim%.*} \
        --update-chr ${bim%.*}_mapfile_final 1 4 \
        --update-map ${bim%.*}_mapfile_final 3 4 \
        --make-bed \
        --output-chr chrMT \
        --out ${bim%.*}_Hg38_Xmerged

        mv *.log logs/
        rm tmp* # Clean Up.
done

## Verification Step ##

for bim in *_Hg38_Xmerged.bim; do
        echo $i ; grep ^chrM ${bim}| wc -l ; grep ^chrX ${bim}| wc -l ; grep ^chrY ${bim} | wc -l ; grep ^chrXY ${bim} | wc -l ;
done
 ## check No mitochondria (chrM or chrMT) and no Pseudoautosomal regions coded chrXY

for bim in *_Hg38_Xmerged.bim; do
        grep rs2341354 ${bim} ;
done
# produces postion corresponding to dbsnp Hg38
