# Parameters
## File
### Input file
HGDP1K="Path/to/file/*vcf.bgz"  # Files from gnomads HGDP+1K dataset
CaG="Path/to/file/*.bim"        # path to  Merge CAG PLINK binary format's bim file, set should contain a .bed .bim .fam file.
### Filtering Ressources file
LongLD="Path/to/file"           # Path to long range linkage regions as defined in Anderson et al. 2010 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3861250/)
PassedQC="Path/to/file"         # Path to gnomad's HGDP_1KG.PassedQC.id
Shared="Path/to/file"           # Path to file of CAG genotyping array shared variants  
## Variables
### General
Threads=5                # Multi Processing   
ind_trace_batch=3000     # Number of individual per batch for trace
### Filter HGDP + 1KG
min_allelefreq=0.05      # Minimal allele frequncing
max_allelefreq=0.95      # Minimal allele frequncing
MISSING_CUTOFF=0.001     # Missingness cut-off   
### Prunning HGDP + 1KG 
window=1000
step=100
Rsq=0.9

# Step 1 Create Shared position file

sed -i 's/:/\t/g' ${Shared}
awk '{print $1, $2-1, $2}' shared.txt | sed 's/ /\t/g' - > shared.bed

# step 2 Process and filter gnomad

for vcf in ${HGDP1K} # Filters each gnomAD file (can be performed in // on slurm with minimal modification
do
  bcftools view -fPASS -q${min_allelefreq} -Q${max_allelefreq} -i "F_MISSING < ${MISSING_CUTOFF}" -S ${PassedQC} -R shared.bed ${vcf} -Oz -o ${vcf%.*}.gz --threads ${Threads}
done

bcftools concat ${vcf%.*}.gz -n -Oz -o gnomad.genomes.v3.1.2.hgdp_tgp.vcf.gz  # Concatenate files
plink2 --vcf gnomad.genomes.v3.1.2.hgdp_tgp.vcf.gz --indep-pairwise ${window} ${step} ${Rsq} --out gnomad.genomes.v3.1.2.hgdp_tgp.LD_prune --set-all-var-ids '@:#:\$r:\$a' --new-id-max-allele-len 145 --threads 5 # ID var to be pruned
plink2 --vcf gnomad.genomes.v3.1.2.hgdp_tgp.vcf.gz --set-all-var-ids 'chr@:#:$r:$a' --new-id-max-allele-len 145 --threads ${Threads}  --exclude gnomad.genomes.v3.1.2.hgdp_tgp.LD_prune.prune.out --make-pgen --out tmp # Remove pruned variants
plink2 --pfile tmp --exclude bed0 ${LongLD} --export vcf bgz id-paste=iid --out gnomad.genomes.v3.1.2.hgdp_tgp.LD_prune # remove long range linkage regions as defined in Anderson et al. 2010 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3861250/)
bcftools index gnomad.genomes.v3.1.2.hgdp_tgp.LD_prune.vcf.gz
vcf2geno --inVcf gnomad.genomes.v3.1.2.hgdp_tgp.LD_prune.vcf.gz --out gnomad.genomes.v3.1.2.hgdp_tgp.LD_prune # convert to geno/site format

## Samples Prep

plink -bfile ${CaG%.*} --a1-allele ${CaG} 6 2  --autosome --output-chr MT --recode vcf-iid bgz --out ${CaG%.*}_autosome 
bcftools query -l ${CaG%.*}_autosome  | split -l ind_trace_batch -d --additional-suffix=.lst

for list in *.lst ; do
  vcf2geno --inVcf ${CaG%.*}_autosome _autosome.vcf.gz --peopleIncludeFile $list --out CaG_${list%.*}
done

## PCA projection

files=($(ls CaG_*.geno))
base_file="${files[0]%.*}"
unset files[0] # pop out the first file

trace -s ${base_file}.geno -g gnomad.genomes.v3.1.2.hgdp_tgp.LD_prune.geno -o ${base_file} -k 20 -K 20 ### produces RefPC.coord used later

for file in ${files} ; do
trace -s ${file} -g gnomad.genomes.v3.1.2.hgdp_tgp.LD_prune.geno -c ${base_file}.RefPC.coord -o ${file%.*} -k 20 -K 20
done

cat ${base_file} > CARTaGENE.ProPC.coord
for file in ${files} ; do sed 1d ${file%.*}.ProPC.coord >>CARTaGENE.ProPC.coord ; done
