#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --account=def-vmooser
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-141

module load bcftools

LIST=sample_list
FILEDIR='./VCF_gatk_hard_filter'
date
echo "Starting task $SLURM_ARRAY_TASK_ID"
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ./${LIST} )
VCF=$FILEDIR/$SAMPLE.gatk_hard_filter.vcf.gz

echo "Sample: $SAMPLE // VCF: $VCF"

bcftools view -m2 -M2 -v snps -R positions.txt $VCF -Oz -o SUBSET/$SAMPLE.gatk.SNV_Illumina_760K.vcf.gz
