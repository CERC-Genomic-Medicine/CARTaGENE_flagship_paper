
# Genotype processing 

## About

This folder contains the nextflow pipeline used to LiftOver, Merge and harmonize the CARTaGENE data (N=29,333 merged final version). This folder was conceived to documents the preparation of CARTaGENE (CaG) data for imputation. 

## Workflow order

```diff
- Important!
This pipeline uses plink version 1 and version 2, the two are called using plink and plink2 respectively. Importantly, they are not interchangeable.
PLINK v2.00a3LM 64-bit Intel is used to establish KING kinship coefficient for kinship cutoff in the 'harmonization.nf' pipeline. Otherwise plink v1.90b6.21 64-bit is used. 
```

### Step 1. Array file pre-process
The 'preprocess.nf' pipeline removes individuals with consent withdrawn, performs LiftOver from b37 to b38 and aligns to the reference. The pipeline is self-explanatory and includes the following steps:

1) Remove samples from individuals with withdrawn consent.
2) LiftOver variant positions from b37 to b38.
3) Aligns variants records to the b38 genome. In this step, the allele A1/A2 encoding is updated from Major/Minor to Reference/Alternate, removing palindromic variants, records that did not contain the reference allele and variants that are not SNPs.
4) Deduplicate variant by removing the variant copy with the highest missingness rate.
4) Set haploid heterozygote (i.e. male genotype in chromosome X non-PAR) to missing.

This pipeline uses :

- Python version 3.11 (with packages in Requierement.txt and bin folder)
- [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020) **Assumed in Path**
- [LiftOver](https://genome-store.ucsc.edu/) **Assumed in Path**

### Step 2. Merge

The 'merge.nf' pipeline deduplicates individuals with the same IID and performs the merger of all arrays plink files. To deduplicate, individuals were removed for the array with the smallest number of variants genotyped. To merge all files, only variants present in all arrays were kept. This pipeline was run using SLURM job scheduler on the Digital Research Alliance of Canada high performance compute clusters. The pipeline is self-explanatory and includes the following steps:

1) Removes duplicate individuals from arrays (with the smallest number of variants genotyped) until there is no shared individual IIDs
2) Establishes a list of shared variants between all arrays and keep only those shared variants for each array.
3) Merges all array on keeping only shared variants

This pipeline uses :

- [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020) **Assumed in Path**

### Step 3. Harmonization

The 'harmonization.nf' pipeline removes position which display a batch effect based on array. To assess the batch effect two models are evaluated using likelihood ratio test (LRT), one with only the projected principal components (representing genetic ancestry) and one with both the projected principal components and the arrays as explanatory variables. This pipeline was run using SLURM job scheduler on the Digital Research Alliance of Canada high performance compute clusters. The pipeline is self-explanatory and includes the following steps:

1) PCA projection onto a reference parallelized. This step is divided in two, one instance is run first to obtain a reference PCA, as to no recompute it. This step was performed on a panel of [gnomAD's HGDP + 1KG callset](https://gnomad.broadinstitute.org/downloads#v3-hgdp-1kg)
2) Establish list of unrelated individuals using plink2's king-cutoff.
3) Test each variant (parallelized by chromosome). This step consists of fitting two linear regression models G ~ PC_1 + ... + PC_4 + dataset_1 + ... + dataset_n and G ~ PC_1 + ... + PC_4, where G = {0, 1, 2} - are individual genotypes, PC_i - principal components, dataset_j - dataset indicator variable. This step also computes F-test P-value and Likelihood Ratio Test P-value comparing the two models. With the assumption of random sampling, the dataset_i variables should not provide significant improvements compared to the nested model.'
4) Filtering the input files for variant with p-value below the chosen threshold (0.05 / number of tests)
5) Formatting output file, this step converts from plink files to vcf, annotates the vcf, set the representation of males genotypes in the non-PAR to haploid.
6) Creates an index for the vcf.

This pipeline uses :
- Python version 3.11 (with packages in Requierement.txt and bin folder)
- [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020) **Assumed in Path**
- [PLINK v2.00a3LM 64-bit Intel ](https://www.cog-genomics.org/plink/2.0/) (22 Mar 2022) **Assumed in Path** 
- [Laser v.1.03](https://csg.sph.umich.edu/chaolong/LASER/) (LASER's vcf2geno and trace are used)
- [bcftools/1.19](https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2) **Assumed in Path**

### Step 4. TOPMed imputation preparation 

The 'topmed_prep.nf' pipeline removes position with more than 0.2 differences in allele frequency between TOPMed and our data (due to poor imputation capacity for those allele), the pipeline also split the data in two random overlapping batches of 25k individuals to satisfy the current limit of TOPMed. The pipeline was run using SLURM job scheduler on the Digital Research Alliance of Canada high performance compute clusters.  The pipeline is self-explanatory and includes the following steps:

1) Test the difference in allele frequency between bravo Freeze 8 and your file.
2) Remove variants with greater differences in allele frequency than the threshold.
3) Create two overlapping list of 25k individuals of representing all individuals. Precisely this step, randomize individuals and separates them into two lists:  the first 25k individuals and the last 25k individuals, with overlap between them. This is done to satisfy TOPMed's current limits.
4) Create two sets of files based on these lists.

This pipeline uses :

- Python version 3.11 (with packages in Requierement.txt and bin folder)
- [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020) **Assumed in Path**
- [bcftools/1.19](https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2) **Assumed in Path**

## Additional

Additionally, we provided script used for the creation of a reference panel in the desired format.

### Panel creation
The 'reference_panel/panel.nf' pipeline creates a reference panel with a maximal number of shared position between the target data and the references, while using standard data curation. This pipeline was run using SLURM job scheduler on the Digital Research Alliance of Canada high performance compute clusters. The pipeline is self-explanatory and includes the following steps:

1) The variants are filtered on shared variant position with the target, quality filters are also imposed, and variants are also selected based on their allele frequency.
2) The variants files are concatenated, LD-pruned and converted to geno/site format.

This pipeline uses :

- [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020) **Assumed in Path**
- [Laser v.1.03](https://csg.sph.umich.edu/chaolong/LASER/) (LASER's vcf2geno was used for vcf conversion to .geno/.site file format)
