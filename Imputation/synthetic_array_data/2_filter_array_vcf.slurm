#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --account=def-vmooser
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-141

module load bcftools

LIST=sample_list
FILEDIR='../array_pos_vcfs'
echo "Starting task $SLURM_ARRAY_TASK_ID"
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ./${LIST} )
VCF=$FILEDIR/$SAMPLE.*.vcf.gz

echo "Sample: $SAMPLE // VCF:" && ls $VCF

bcftools view $VCF -Ou | bcftools  +setGT -Ou - -- -t q -n . -e 'FILTER="PASS"' | bcftools annotate -x INFO,^FORMAT/GT -Oz -o $SAMPLE.array.setmissing.vcf.gz
bcftools index -t  $SAMPLE.array.setmissing.vcf.gz
