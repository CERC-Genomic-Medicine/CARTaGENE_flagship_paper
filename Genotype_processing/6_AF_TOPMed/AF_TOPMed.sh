#!/bin/bash

plink --bfile CARTaGENE_hg38_Harmonized  --keep-allele-order --freq --out CARTaGENE_hg38_Harmonized
python3 plink_freq_vs_topmed.py -s CARTaGENE_hg38_Harmonized.frq -tALL.BRAVO_TOPMed_Freeze_8.vcf.gz -m 0.2 -o Diff_02.txt
plink --bfile CARTaGENE_hg38_Harmonized --exclude Diff_02.txt.af_diff.txt --keep-allele-order --make-bed --out CARTaGENE_hg38_shared_TOPMed_filtered
plink -bfile CARTaGENE_hg38_shared_TOPMed_filtered --a1-allele CARTaGENE_hg38_shared_TOPMed_filtered.bim 6 2  --output-chr chrMT --recode vcf-iid bgz --out CARTaGENE_hg38_shared_TOPMed_filtered
