# TOPMed filter
This folder contains the scripts and documentation used to identify and remove the variants that display different allele frequency between the CaG data and TOPMed's data. This is performed in preparation to imputation on TOPMed data. This step is reliant on the [Harmonization step](https://github.com/CERC-Genomic-Medicine/CARTaGENE_flagship_paper/tree/main/Genotype_processing/5_Harmonization/Harmonization_test).
  
## Parameters 
### Input file
bed files - **Merged and harmonized**  CaG data v.1.1 in plink .bed format (including .bim .fam companion file)
- TOPMed's freeze v8 dataset 
### Scripts
[plink_freq_vs_topmed.py](https://github.com/CERC-Genomic-Medicine/scripts/blob/master/plink_freq_vs_topmed.py)
### Software
- [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020)
- [python 3.11.5](https://www.python.org/downloads/release/python-3115/)
   - argparse
   - pandas 1.4.4
   - pysam-0.22.0
   - matplotlib-3.7.2

## Steps
### Step 1 - identify and remove variants with different allele frequency

1) Calculate CaG allele frequency with [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020) 's --freq function
2) identify discordant variants with [plink_freq_vs_topmed.py](https://github.com/CERC-Genomic-Medicine/scripts/blob/master/plink_freq_vs_topmed.py)
   - Step Parameter
     - Discordant threshold 0.2
3) Exlude identified variants with [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020) 's --exclude option
