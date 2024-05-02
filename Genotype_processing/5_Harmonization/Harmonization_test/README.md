# Harmonization Test
This folder contains the scripts and documentation used to identify the list of variant failling the likelihood ratio test comparing two linear regression models G ~ m + PC1 + ... + PC4 + A2 + ... + A5 and G ~ m + PC1 + ... + PC4, where PCs are the Principal components of this projection and A are the arrays as binairy variables. The failling variant are considered to be non-harmonious across arrays. This test is performed using the of CaG merge samples' projection on [gnomAD's HGDP+1KG dataset](https://gnomad.broadinstitute.org/downloads#v3-hgdp-1kg).This process relies on the output of the [array merger step](Genotype_processing/5_Harmonization/PCA_Projection/PCA_Projection.sh) and some of the files produced in the [array merger step](https://github.com/CERC-Genomic-Medicine/CARTaGENE_flagship_paper/tree/main/Genotype_processing/4_Merge).

## Parameter
### Input files
PCA Coordinate file - File of the coordinate of CaG samples within the PCA space of gnomAD's HGDP+1KG.
bed file - CaG data v.1.1 in plink .bed format (including .bim .fam companion file) **per array**
bed file - **Merged**  CaG data v.1.1 in plink .bed format (including .bim .fam companion file)
### Scripts
[compare_ancestry_adjusted_af.py]([https://github.com/CERC-Genomic-Medicine/scripts/compare_ancestry_adjusted_af.py](https://github.com/CERC-Genomic-Medicine/scripts/blob/master/compare_ancestry_adjusted_af.py))  
### Software
- [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020)
- [PLINK v2.00a3LM 64-bit (22 Mar 2022)](www.cog-genomics.org/plink/2.0/)
- [python 3.11.5](https://www.python.org/downloads/release/python-3115/)
  - pandas 1.4.4
  - argparse
  - pysam-0.22.0
  - numpy-1.26.4
  - scipy-1.11.2
  - statsmodels.api-0.15.0
### Analysis parameters
- Pval_threshold : 0.05 - P-value threshold of likelihood test

## Steps
### Step 1 - Variable preparation (Covariate / explanatory)
 - File_prep_1.sh
 - File_prep.py
1) Convert CaG files to unrelated CaG individuals with [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020), and [PLINK v2.00a3LM 64-bit (22 Mar 2022)](www.cog-genomics.org/plink/2.0/)'s --king-cutoff
- Step Parameter
  - 0.0442 kinship coefficients cut off
2) Create Genotyping array label files for unrelated individuals, along with sex-specific files
3) Convert header of projection Coordinate (File_prep.py)

### Step 2 - Testing
 - Harmonization_test.sh
1) test each autosome chromosome, PAR regions, and non-PAR regions (for each sex).
2) Concatenate the results

### Step 3 - Identify non-harmonious variants
 - Hamonization_2.py
1) Identify variant with P-value < 0.05

### Step 4 - Filter non-harmonious variants
variant_removal.sh
1) Filter out previously identified variant using PLINK v2.00a3LM 64-bit (22 Mar 2022)](www.cog-genomics.org/plink/2.0/)'s exclude
