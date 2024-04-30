# Goal :
- Create a PCA projection of samples on Reference HGDPK1

## REFERENCE Prep

# Create Shared position file

sed -i 's/:/\t/g' shared.txt
awk '{print $1, $2-1, $2}' shared.txt | sed 's/ /\t/g' - > shared.bed

# Process and filter gnomad
# As a Job where SLURM_ARRAY_TASK_ID = {1..22}

bcftools view -fPASS -q0.05 -Q0.95 -i 'F_MISSING < 0.001' -S /project/rrg-vmooser/shared/HGDP_1KG/HGDP_1KG.PassedQC.id -R shared.bed [Original]/gnomad.genomes.v3.1.2.hgdp_tgp.chr${SLURM_ARRAY_TASK_ID}.vcf.bgz -Oz -o gnomad.genomes.v3.1.2.hgdp_tgp.chr${SLURM_ARRAY_TASK_ID}.vcf.gz --threads 5

# Process and filter gnomad

bcftools concat gnomad.genomes.v3.1.2.hgdp_tgp.chr{1..22}.vcf.gz -n -Oz -o gnomad.genomes.v3.1.2.hgdp_tgp.vcf.gz
plink2 --vcf gnomad.genomes.v3.1.2.hgdp_tgp.vcf.gz --indep-pairwise 1000 100 0.9 --out gnomad.genomes.v3.1.2.hgdp_tgp.LD_prune --set-all-var-ids '@:#:\$r:\$a' --new-id-max-allele-len 145 --threads 5
wget https://github.com/CERC-Genomic-Medicine/Regenie_nextflow/blob/3a7f3e3bb5e5a9cbc8a075b54917c423237cdc28/util/Low_complexity_regions/HG38.Anderson2010.bed.gz
plink2 --vcf gnomad.genomes.v3.1.2.hgdp_tgp.vcf.gz --set-all-var-ids 'chr@:#:$r:$a' --new-id-max-allele-len 145 --threads 5  --exclude gnomad.genomes.v3.1.2.hgdp_tgp.LD_prune.prune.out --make-pgen --out tmp
plink2 --pfile tmp --exclude bed0 HG38.Anderson2010.bed.gz --export vcf bgz id-paste=iid --out gnomad.genomes.v3.1.2.hgdp_tgp.LD_prune
bcftools index 
vcf2geno --inVcf gnomad.genomes.v3.1.2.hgdp_tgp.LD_prune.vcf.gz --out gnomad.genomes.v3.1.2.hgdp_tgp.LD_prune

## Samples Prep

plink -bfile CARTaGENE_hg38_shared --a1-allele CARTaGENE_hg38_shared.bim 6 2  --autosome --output-chr MT --recode vcf-iid bgz --out CARTaGENE_hg38_shared_autosome
bcftools query -l CARTaGENE_hg38_shared_autosome | split -l 3000 -d --additional-suffix=.lst
for i in {0..9} ; do
  vcf2geno --inVcf CARTaGENE_hg38_shared_autosome.vcf.gz --peopleIncludeFile x0"$i".lst --out CARTaGENE_hg38_shared_"$i"
done

## PCA projection

trace -s CARTaGENE_hg38_shared_0.geno -g 1000genomes_unrelate.geno -o CARTaGENE_"$i" -k 20 -K 20 ### produces CARTaGENE_0.RefPC.coord used later

for i in {1..9} ; do
trace -s CARTaGENE_hg38_shared_"$i".geno -g 1000genomes_unrelate.geno -c CARTaGENE_0.RefPC.coord -o CARTaGENE_"$i" -k 20 -K 20
done

cat CARTaGENE_0.ProPC.coord > CARTaGENE.ProPC.coord
for i in {1..9} ; do sed 1d CARTaGENE_"$i".ProPC.coord >>CARTaGENE.ProPC.coord ; done
