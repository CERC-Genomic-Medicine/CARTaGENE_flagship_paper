One step before performing the imputation, is to compare the allele frequency of variants in the genotyping array data with the reference panel.

The first step is to extract allele frequencies from genotyping array VCFs and reference panel VCFs using the following commands. We need to make sure that AF column is present in VCF files and if not we can use bcftools +fill-tag plugin to fill in the allele frequency column.

```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' /path/to/genotyping_array_vcf > genotyping_array_AF.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' /path/to/reference_vcf > reference_AF.txt
```

Then we can calculate the sample size of each dataset using the following commands:

```bash
bcftools query -l path/to/genotyping_array_vcf | wc
bcftools query -l /path/to/reference_vcf | wc
```

Then we run the fisher exact python program using the following command: 

```bash
python Fisher_exact.py -ig genotyping_array_AF.txt -ir reference_AF.txt -ng genotyping_array_sample_size -nr referene_sample_size -o path/to/output_txt_file
```

Then we filter out variants from the genotyping array data which had significant differences in allele frequencies
(significance p-value threshold of 1e-50).

This codes can be parallelized based on the chromosomes if needed to increase the speed.