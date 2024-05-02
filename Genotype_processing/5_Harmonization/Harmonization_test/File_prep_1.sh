# Parameters
## File
### Input file
CaG="Path/to/file/*.bim"        # path to  Merge CAG PLINK binary format's bim file, set should contain a .bed .bim .fam file.
CaG_array="Path/to/file/*.bim"  # path to  not merged CAG PLINK binary format's bim file, set should contain a .bed .bim .fam file.
## Variable
king_cutoff=0.0442

## Create files of unrelated CaG individuals

plink --bfile ${CaG%.*} --keep-allele-order --output-chr chrMT --export vcf bgz id-paste=iid --out ${CaG%.*}

plink2 --vcf ${CaG%.*}.vcf.gz --king-cutoff ${king_cutoff} --output-chr chrMT --make-pgen --out ${CaG%.*}_unrelated

plink2 --pfile ${CaG%.*}_unrelated --output-chr chrMT --export vcf bgz id-paste=iid --out ${CaG%.*}_unrelated

## Created Label files

for bim in CaG_array; do 
  cut -f 1 -d ' ' ${CaG_array%.*}.fam | sed "s/$/\t${CaG_array%.*}.fam/g" >> label;
done

cut -f 1,5 ${CaG%.*}.fam | grep 2$ | cut -f 1 > female.filter
cut -f 1,5 ${CaG%.*}.fam | grep 1$ | cut -f 1 > male.filter
cut -f 1  ${CaG%.*}_unrelated.king.cutoff.in.id >  ${CaG%.*}_unrelated.king.cutoff.in.iid

grep -f ${CaG%.*}_unrelated.king.cutoff.in.iid label > label_unrelated
grep -f male.filter label_unrelated > label_unrelated_male
grep -f female.filter label_unrelated > label_unrelated_female
