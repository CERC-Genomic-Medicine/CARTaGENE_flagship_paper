#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --account=def-vmooser
#SBATCH --mem-per-cpu=8G

module load bcftools

LIST=array_vcfs
OUT=merged.array-snps.qc.vcf.gz

bcftools merge -l $LIST -Ou | \
bcftools +fill-tags -Ou -- -t AN,AC,AF,NS,F_MISSING | \
# set sites where missing for all samples to reference 0/0
bcftools +setGT -Ou - -- -t q -n 0 -i 'F_MISSING=1' | \
bcftools view -e 'F_MISSING>0.5' -m 2 -M 2 -v snps -Ou | \
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o $OUT

bcftools index -t $OUT
