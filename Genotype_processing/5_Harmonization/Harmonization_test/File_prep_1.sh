# Script
## CARTaGENE unrelated file
plink --bfile ../Mooser_433651_genotypes_hg38_Merged/CARTaGENE_hg38_shared --keep-allele-order --output-chr chrMT --export vcf bgz id-paste=iid --out CARTaGENE_hg38_shared

plink2 --vcf CARTaGENE_hg38_shared.vcf.gz --king-cutoff 0.0442 --output-chr chrMT --make-pgen --out CARTaGENE_hg38_shared_unrelated

plink2 --pfile CARTaGENE_hg38_shared_unrelated --output-chr chrMT --export vcf bgz id-paste=iid --out CARTaGENE_hg38_shared_unrelated

## Label files
declare -a arrays=("17k" "5300" "4224" "760" "archi")

for array in "${arrays[@]}"; do 
  cut -f 1 -d ' ' ${array}_hg38_shared.fam | sed "s/$/\t${array}/g" >> label;
done

cut -f 1,5 ../Mooser_433651_genotypes_hg38_Merged/CARTaGENE_hg38_shared.fam | grep 2$ | cut -f 1 > female.filter
cut -f 1,5 ../Mooser_433651_genotypes_hg38_Merged/CARTaGENE_hg38_shared.fam | grep 1$ | cut -f 1 > male.filter
cut -f 1  CARTaGENE_hg38_shared_unrelated.king.cutoff.in.id >  CARTaGENE_hg38_shared_unrelated.king.cutoff.in.iid

grep -f CARTaGENE_hg38_shared_unrelated.king.cutoff.in.iid label > label_unrelated
grep -f male.filter label_unrelated > label_unrelated_male
grep -f female.filter label_unrelated > label_unrelated_female
