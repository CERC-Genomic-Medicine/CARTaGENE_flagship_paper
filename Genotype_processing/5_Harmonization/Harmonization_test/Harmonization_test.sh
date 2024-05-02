# Parameters
## File
### Input file
CaG_unrelated="Path/to/file/[file].vcf.gz"        # path to  unrelated CAG VCF file
Label="Path/to/file/[file]"                       # Array label of unrelated inviduals (arrays)
Label_male="Path/to/file/[file]"                  # Array label of unrelated male inviduals 
Label_female="Path/to/file/[file]"                # Array label of unrelated female inviduals
CaG_pca="Path/to/file/[file].PCA"
## Variables
PAR1='chrX:10001-2781479'       # PAR1 hg38
PAR2='chrX:155701383-156030895' # PAR2 hg38
nPAR='chrX:2781479-155701383'   # non-par hg38

for i in {1..22}
do
  python compare_ancestry_adjusted_af.py -v ${CaG_unrelated} -l ${Label} -r chr"$i" -p ${CaG_pca} -k 4 -o AF_ftest_4PC_chr$i
done
python compare_ancestry_adjusted_af.py -v ${CaG_unrelated} -l ${Label} -r ${PAR1} -p ${CaG_pca} -k 4 -o AF_ftest_4PC_PAR1
python compare_ancestry_adjusted_af.py -v ${CaG_unrelated} -l ${Label} -r ${PAR2} -p ${CaG_pca} -k 4 -o AF_ftest_4PC_PAR2
python compare_ancestry_adjusted_af.py -v ${CaG_unrelated} -l ${Label_male} -r ${nPAR}  -p ${CaG_pca} -k 4 -o AF_ftest_4PC_Xmale
python compare_ancestry_adjusted_af.py -v ${CaG_unrelated} -l ${Label_female} -r ${nPAR}  -p ${CaG_pca} -k 4 -o AF_ftest_4PC_Xfemale

## concat results
cat AF_ftest_4PC_chr1 > AF_ftest_4PC
for i in chr{2..22} PAR1 PAR2 Xfemale Xmale
do 
sed 1d AF_ftest_4PC_"$i" >> AF_ftest_4PC
done
