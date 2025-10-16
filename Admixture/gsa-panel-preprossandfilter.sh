#!/bin/bash

module load plink/1.9b_6.21-x86_64
module load r

# working directory
# cd Morocco_Cartagene

# Assign plink extension
#plink=data/gsa_cartagene/gsa_merged_hg38
plink=$1
base=$(basename $plink)
#keepbase=NAfrica_FC
keepbase=$2

# 1. Get allele frequencies before any preprocessing step
# plink --bfile ${plink} --freq --missing --out allele_freqs/${base}
# plink --bfile ${plink} --family --freq --missing --out allele_freqs/${base}

# 2. Filter individuals from gsa - Morocco and sample of frenchCAD

# Get individuals to filter
# python scripts/extract_metadata.py MOROCCO "FRENCH CANADIAN" 30 ${plink}.${keepbase}.keep
# Run ipynb using filter 

# Filter individuals from list
plink --bfile ${plink} --keep <(cut -f2,3 ${plink}.${keepbase}.keep) --make-bed --out ${plink}.${keepbase}.tmp

# Change Pop ID of individuals
paste <(cut -f2,3 ${plink}.${keepbase}.keep) <(cut -f1,2 ${plink}.${keepbase}.keep) > ${plink}.${keepbase}.updateIDs.tmp
plink --bfile ${plink}.${keepbase}.tmp --update-ids ${plink}.${keepbase}.updateIDs.tmp --make-bed --out ${plink}.${keepbase}

# 3. Get allele frequencies after filtering step
plink --bfile ${plink}.${keepbase} --freq --missing --out allele_freqs/${base}.${keepbase}
plink --bfile ${plink}.${keepbase} --family --freq --missing --out allele_freqs/${base}.${keepbase}

rm *.tmp*



