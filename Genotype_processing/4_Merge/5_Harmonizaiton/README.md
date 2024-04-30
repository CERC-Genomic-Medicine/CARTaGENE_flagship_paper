# Harmonization

## Goal
- Identify variants which the array contribute significantly to the allele frequency (thus likely a batch effect) accounting for Genetic ancestry.
- Filter said variants.

## Order
1) PCA_Projection
2) Harmonization_test
3) Remove Variants (Remove_Variants.sh

## Software

Plink 1.9 (1.9b_6.21-x86_64)
plink2 (PLINK v2.00a3LM)

bcftools 1.19
trace (Laser v.1.03)
vcf2geno (Laser v.1.03)
