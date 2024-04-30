# Subsetting WGS data to genotyping array positions

## Step 1. Download array positions
Run `download_array_positions.sh` bash script to download positions of variants on Illumna GSA v2 array.

## Step 2. Generate VCF for genotyping array positions
Run `subset_wgs.nf` [Nextflow](https://www.nextflow.io/) pipeline to combine single-sample VCF files from WGS into multi-sample VCF with genotyping array positions. The pipeline assumes that single-sample WGS-based VCFs were generated using GATK DRAGEN and variants were hard-filtered based on GATK DRAGEN's hard-filter thresholds. The pipeline also checks that males are coded as haploids in the chrX PAR region. The pipeline keeps only non-monomorphic SNVs with missingness <10%.

The pipeline uses [bcftools](https://samtools.github.io/bcftools/bcftools.html).
