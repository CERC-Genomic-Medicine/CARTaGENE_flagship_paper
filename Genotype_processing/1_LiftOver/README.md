# Liftover
This folder contains the scripts and documentation used to liftover CaG v.1.1 dataset. This folder relies on input dataset with up-to-date samples composition (i.e. individual with withdrawn consent removed) this was created in [Preprocesssing](https://github.com/CERC-Genomic-Medicine/CARTaGENE_flagship_paper/tree/main/Genotype_processing/0_preprocessing). 

## Parameters
### Software used
- [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020)  
- [LiftOver](https://genome-store.ucsc.edu/)

### Input files
Bedfiles - CaG data v.1.1 in plink .bed format (including .bim .fam companion file) **per Genotyping arrays** issued from the [Preprocesssing step](https://github.com/CERC-Genomic-Medicine/CARTaGENE_flagship_paper/tree/main/Genotype_processing/0_preprocessing)   
[hg19ToHg38.over.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz) from UCSC  

### Output file
Hg38 Bedfile - Hg38 CaG data v.1.1 in plink .bed format (including .bim .fam companion file) **per Genotyping arrays** 

## Steps
### Step 1 - Merge pseudo-autosomal regions  (XY to X)
In order to perform liftover, pseudo-autosomal region were merged into their chromosome X coordinates. To achieve this [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/)'s --merge-x option was used. 
  
### Step 2 - LiftOver
Liftover from hg19 (b37) to hg38 (b38) requiered the identification of corresponding coordinate for each variant accross build to do this, the [LiftOver](https://genome-store.ucsc.edu/) tool was used with the corresponding ressource file [hg19ToHg38.over.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz) for the desired builds. Genotyping files were updated using [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/), Unmapped variants and mapped to alternate contig were withdrawn.

### Step 3 - Validation (optional)
To validate liftover, the absence of mitochondria (chrM or chrMT) and Pseudoautosomal regions coded chrXY can be evaluated, along with the change in the position of variant rs2341354 (1:918573 (GRCh37) -> 1:983193 (GRCh38))
