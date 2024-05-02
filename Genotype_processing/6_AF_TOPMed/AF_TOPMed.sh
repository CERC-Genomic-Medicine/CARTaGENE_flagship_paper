#!/bin/bash
CaG="Path/to/file.bim" # path to a CAG (harmonized and merged) PLINK binary format's bim file, each should contain a set of .bed .bim .fam file.
Bravo="Path/to/file.vcf" # path to TOPMed's freeze 8 
Threshold=0.2 # Frequency concordance threshold

plink --bfile %{CaG%.*}  --keep-allele-order --freq --out %{CaG%.*} # Make Allele frequency
python3 plink_freq_vs_topmed.py -s %{CaG%.*}.frq -t ${Bravo} -m 0.2 -o Diff_${Threshold} # Test frequency concordance
plink --bfile %{CaG%.*} --exclude Diff_${Threshold}.af_diff.txt --keep-allele-order --make-bed --out %{CaG%.*}_TOPMed_filtered # Make Bed
plink -bfile %{CaG%.*}_filtered --a1-allele %{CaG%.*}_TOPMed_filtered.bim 6 2  --output-chr chrMT --recode vcf-iid bgz --out %{CaG%.*}_TOPMed_filtered # make VCF
