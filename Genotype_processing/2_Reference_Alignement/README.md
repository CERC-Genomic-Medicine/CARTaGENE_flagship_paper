# Reference Alignement
This folder contains the scripts and documentation used to align the CaG v.1.1 dataset to the reference. This folder's scripts relies on a hg38 dataset, this was created in [liftOver](https://github.com/CERC-Genomic-Medicine/CARTaGENE_flagship_paper/tree/main/Genotype_processing/1_LiftOver). 

## Parameter
### Input files
Bedfiles - CaG data v.1.1 in plink .bed format (including .bim .fam companion file) **per Genotyping arrays** issued from the [LiftOver step](https://github.com/CERC-Genomic-Medicine/CARTaGENE_flagship_paper/tree/main/Genotype_processing/1_LiftOver)  

[Homo_sapiens.GRCh38.fa](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000001405.28/)

### Output files :
Bedfile - aligned and renamed CaG data v.1.1 in plink .bed format (including .bim .fam companion file) per Genotyping arrays

### Used Software :
- [python 3.11.5](https://www.python.org/downloads/release/python-3115/)
  - argparse
  - pysam-0.22.0
- [plink2reference.py](https://github.com/CERC-Genomic-Medicine/scripts/plink2reference.py)
- [plink_dupvar_prioritize.py](https://github.com/CERC-Genomic-Medicine/scripts/dupvar_prioritize.py)
- [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020)

## Steps
### Step 1 - Alignement
Alignement to the reference, transitionning form plink1.9's binary format's normal encoding, i.g. Major/Minor Allele, to REF/ALT. To achieve this :
1) identify elements to remove, strand-flips, and force allele, using [plink2reference.py](https://github.com/CERC-Genomic-Medicine/scripts/plink2reference.py)
2) remove, strand-flips, and force allele allele, using [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020)'s --exclude --flip --a1-allele respectively.


### Step 2 - Deduplication
Redundant genotyping were present in all arrays (i.e genotyping same variant twice or more). To address these redundancies :
1) duplicated elements were identified using [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020) --list-duplicate-vars.
2) calculate frequency with [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/)'s --freq counts
3) Identify duplicate with lower coverage using [plink_dupvar_prioritize.py](https://github.com/CERC-Genomic-Medicine/scripts/dupvar_prioritize.py)
4) Remove duplicate with lower coverage with using [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/)'s -- exclude

### Step 3 - Renaming variants
change variants naming system to chr:pos:REF:ALT
1) Rename [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/)'s --update-name

### Step 4 - Verification (Optional)
1) Verify the absence of M and XY chromosome.
