#!/bin/bash

for i in {1..22}
do
  python compare_ancestry_adjusted_af.py -v CARTaGENE_hg38_shared_unrelated.vcf.gz -l label.unrelated -r chr"$i" -p CARTaGENE_projection.PCA -k 4 -o AF_ftest_4PC_chr$i
done
python compare_ancestry_adjusted_af.py -v CARTaGENE_hg38_shared_unrelated.vcf.gz -l label.unrelated -r chrX:10001-2781479 -p CARTaGENE_projection.PCA -k 4 -o AF_ftest_4PC_PAR1
python compare_ancestry_adjusted_af.py -v CARTaGENE_hg38_shared_unrelated.vcf.gz -l label.unrelated -r chrX:155701383-156030895 -p CARTaGENE_projection.PCA -k 4 -o AF_ftest_4PC_PAR2
python compare_ancestry_adjusted_af.py -v CARTaGENE_hg38_shared_unrelated.vcf.gz -l label_male.unrelated -r chrX:2781479-155701383  -p CARTaGENE_projection.PCA -k 4 -o AF_ftest_4PC_Xmale
python compare_ancestry_adjusted_af.py -v CARTaGENE_hg38_shared_unrelated.vcf.gz -l label_female.unrelated -r chrX:2781479-155701383  -p CARTaGENE_projection.PCA -k 4 -o AF_ftest_4PC_Xfemale

## concat results
cat AF_ftest_4PC_chr1 > AF_ftest_4PC
for i in chr{2..22} PAR1 PAR2 Xfemale Xmale
do sed 1d AF_ftest_4PC_"$i" >> AF_ftest_4PC
done
