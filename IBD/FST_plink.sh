#!/bin/bash

# Load necessary modules
module load plink/2.00a3.6  # PLINK for genetic data processing
module load r               # R for statistical analysis

# Convert VCF file to PLINK binary format (BED, BIM, FAM)
# Exclude specific individuals from group3_small_cluster_to_exclude.txt
plink2 --vcf /lustre07/scratch/justinp/NextFlow/CaG_genotyping/Data/CARTaGENE_hg38_shared.phased.maf.hwe.missingness.chr_all.vcf.gz \
       --exclude group3_small_cluster_to_exclude.txt \
       --make-bed \
       --out CARTaGENE_hg38_shared.phased.maf.hwe.missingness.chr_all

# Convert VCF to PLINK binary format without exclusions
plink2 --vcf /lustre07/scratch/justinp/NextFlow/CaG_genotyping/Data/CARTaGENE_hg38_shared.phased.maf.hwe.missingness.chr_all.vcf.gz \
       --make-bed \
       --out CARTaGENE_hg38_shared.phased.maf.hwe.missingness.chr_all

# Adjust FAM file format for compatibility
cut -f 2 CARTaGENE_hg38_shared.phased.maf.hwe.missingness.chr_all.fam > tmp1
cut -f 3-6 CARTaGENE_hg38_shared.phased.maf.hwe.missingness.chr_all.fam > tmp2
paste tmp1 tmp1 tmp2 > CARTaGENE_hg38_shared.phased.maf.hwe.missingness.chr_all.fam
rm tmp1 tmp2  # Clean up temporary files

# Compute pairwise FST between clusters for group 3
iplink2 --bfile CARTaGENE_hg38_shared.phased.maf.hwe.missingness.chr_all \
        --within group3_cluster_info.txt \
        --fst CATPHENO method=hudson \
        --exclude group3_small_cluster_to_exclude.txt \
        --out pairwise_FST_CAG_cluster3

# Prepare cluster information for group 4
cut -d" " -f1 v2_clusters_hap_IBD_louvain_iteration4.txt > tmp1
cut -d" " -f2 v2_clusters_hap_IBD_louvain_iteration4.txt > tmp2
paste tmp1 tmp1 tmp2 > group4_cluster_info.txt
rm tmp1 tmp2  # Clean up temporary files

# Compute pairwise FST between clusters for group 4
plink2 --bfile CARTaGENE_hg38_shared.phased.maf.hwe.missingness.chr_all \
       --within group4_cluster_info.txt \
       --fst CATPHENO method=hudson \
       --out pairwise_FST_CAG_cluster4

# Clean up output files by removing unwanted character 'C'
sed -i 's/C//g' pairwise_FST_CAG_cluster3.fst.summary

# Run R script to merge FST results
Rscript FST_merge_clusters.R
