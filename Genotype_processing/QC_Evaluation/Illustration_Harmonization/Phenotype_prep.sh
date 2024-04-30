#!/bin/bash

grep 0$ CARTaGENE_samples_phase.txt | cut -f 2 > phase_1.txt
grep 1$ CARTaGENE_samples_phase.txt | cut -f 2 > phase_2.txt
for i in {1..2}; do for j in archi 760 5300 4224 17k; do comm -1 -2 <(cut -f 2 -d " " "$j"_hg38_shared.fam | sort) <(sort phase_"$i".txt) > "$j"_phase"$i".samples; done; done

awk '{print $1 "\t" $1 "\t" "1"}' 4224_phase1.samples > 4Kphase1vs5Kphase1.pheno_temp
awk '{print $1 "\t" $1 "\t" "2"}' 5300_phase1.samples >> 4Kphase1vs5Kphase1.pheno_temp

awk '{print $1 "\t" $1 "\t" "1"}' 4224_phase1.samples > 4Kphase1vsarchiphase1.pheno_temp
awk '{print $1 "\t" $1 "\t" "2"}' archi_phase1.samples >> 4Kphase1vsarchiphase1.pheno_temp

awk '{print $1 "\t" $1 "\t" "1"}' 17k_phase1.samples > 17Kphase1vs5Kphase1.pheno_temp
awk '{print $1 "\t" $1 "\t" "2"}' 5300_phase1.samples >> 17Kphase1vs5Kphase1.pheno_temp

awk '{print $1 "\t" $1 "\t" "1"}' archi_phase1.samples > archiphase1vs760phase1.pheno_temp
awk '{print $1 "\t" $1 "\t" "2"}' 760_phase1.samples >> archiphase1vs760phase1.pheno_temp

awk '{print $1 "\t" $1 "\t" "1"}' 17k_phase1.samples > 17Kphase1vs17Kphase2.pheno_temp
awk '{print $1 "\t" $1 "\t" "2"}' 17k_phase2.samples >> 17Kphase1vs17Kphase2.pheno_temp

awk '{print $1 "\t" $1 "\t" "1"}' 5300_phase1.samples > 5Kphase1vs5Kphase2.pheno_temp
awk '{print $1 "\t" $1 "\t" "2"}' 5300_phase2.samples >> 5Kphase1vs5Kphase2.pheno_temp

for i in 4Kphase1vs5Kphase1 4Kphase1vsarchiphase1  17Kphase1vs5Kphase1 archiphase1vs760phase1 17Kphase1vs17Kphase2 5Kphase1vs5Kphase2; do grep -f CARTaGENE_hg38_shared_unrelated.king.cutoff.in.id "$i".pheno_temp | sed "1i #FID\tIID\t$i" > $i.pheno; done

plink2 --pfile CARTaGENE_hg38_shared_unrelated --pca 5 --out CARTaGENE
cut -f 1,2,5  CARTaGENE_hg38_shared.fam | grep -f CARTaGENE_hg38_shared_unrelated.king.cutoff.in.id | sed 's/1$/0/g' | sed 's/2$/1/g' | sed '1i FID\tIID\tSEX'  > SEX_ID

python3 covar_PCA_prep.py

