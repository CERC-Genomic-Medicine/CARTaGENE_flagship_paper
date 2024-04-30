#!/bin/bash

plink -bfile CARTaGENE_hg38_shared --keep-allele-order  --output-chr chrMT --exclude Variants_failing_ancestry_LRTtest.txt --make-bed --out CARTaGENE_hg38_shared_AF_filtered
plink -bfile CARTaGENE_hg38_shared --a1-allele CARTaGENE_hg38_shared.bim 6 2  --output-chr chrMT --exclude Variants_failing_ancestry_LRTtest.txt --recode vcf-iid bgz --out CARTaGENE_hg38_shared_AF_filtered
