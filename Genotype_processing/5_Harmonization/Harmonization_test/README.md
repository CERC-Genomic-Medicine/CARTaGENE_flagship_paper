# Harmonization Test

## Goal 
- Compare ancestry aware allele frequency across arrays

## Input :
- CARTaGENE_hg38_shared (bed files)
- ${array}_hg38_shared.fam (include arrays "17k" "5300" "4224" "760" "archi")
- script compare_ancestry_adjusted_af.py from CERC github

## Workflow

1. Preparation of Covariate / explanatory variables
 - File_prep_1.sh
 - File_prep_2.py
2. Harmonisation test (testing)
 - Harmonization_test.sh
3. Determination at adjusted p-value 0.05
 - Hamonization_2.py
