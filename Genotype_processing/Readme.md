# Genotype Processing 

## About

This folder contains the Nextflow pipeline used to LiftOver, merge, and harmonize the CARTaGENE data (N=29,333 merged final version). This folder documents the preparation of CARTaGENE (CaG) data for imputation.

## Workflow Order

```diff
- Important!
These pipelines use PLINK version 1 and version 2. The two are called using the 'plink' and 'plink2' commands respectively. Importantly, they are not interchangeable.
PLINK v2.00a3LM 64-bit Intel is used to establish the KING kinship coefficient for kinship cutoff in the 'harmonization.nf' pipeline. Otherwise, PLINK v1.90b6.21 64-bit is used.
```

### Step 1. Array Data Preprocessing

The 'preprocess.nf' pipeline removes individuals with withdrawn consent, performs LiftOver from b37 to b38, and aligns to the reference. This pipeline was run on the Digital Research Alliance of Canada high-performance compute clusters. The pipeline includes the following steps:

1. Remove samples from individuals with withdrawn consent.
2. LiftOver variant positions from b37 to b38.
3. Align variant records to the b38 genome. In this step, the allele A1/A2 encoding is updated from Major/Minor to Reference/Alternate. Palindromic variants, records that did not contain the reference allele, and variants that are not SNPs are removed.
4. Deduplicate variants by removing the variant copy with the highest missingness rate.
5. Set haploid heterozygote (i.e., heterozygous male genotypes in chromosome X non-PAR) to missing.

This pipeline uses:

- Python version 3.11 (with packages in requirements.txt and scripts in the bin folder)
- [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020) **Assumed in Path**
- [LiftOver](https://genome-store.ucsc.edu/) **Assumed in Path**

### Step 2. Merge

The 'merge.nf' pipeline deduplicates individuals with the same IID (i.e., the PLINK file's within-family IDs), as in this dataset, this value is meant to be unique to each sample and performs the merger of all arrays' PLINK files. To deduplicate, individuals with shared IID were removed from the array with the smallest number of variants genotyped at the beginning of this step. To merge all files, only variants present in all arrays were retained. This pipeline was run using the SLURM job scheduler on the Digital Research Alliance of Canada high-performance compute clusters. The pipeline is self-explanatory and includes the following steps:

1. Remove duplicate individuals from arrays (with the smallest number of variants genotyped) until there are no shared individual IIDs.
2. Establish a list of shared variants between all arrays and keep only those shared variants for each array.
3. Merge all arrays, keeping only shared variants.

This pipeline uses:

- [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020) **Assumed in Path**

### Step 3. Harmonization

The 'harmonization.nf' pipeline removes positions that display a batch effect based on the array. To assess the batch effect, two models are evaluated using the likelihood ratio test (LRT): one with only the first four projected principal components (representing genetic ancestry) and one with both the projected principal components and the arrays as explanatory variables. This pipeline was run using the SLURM job scheduler on the Digital Research Alliance of Canada high-performance compute clusters. The pipeline is self-explanatory and includes the following steps:

1. Perform parallelized PCA projection onto a reference. This step is divided into two parts: first, an instance is run to obtain a reference PCA, so the second part, the other parallel instances, do not recompute it. This step was performed using a reference panel made from [gnomAD's HGDP + 1KG callset](https://gnomad.broadinstitute.org/downloads#v3-hgdp-1kg) (N=4,119).
2. Establish a list of unrelated individuals using PLINK2's king-cutoff.
3. Perform a likelihood ratio test on each variant (parallelized by chromosome).
4. Filter out the genotype files for variants with p-values below or equal to the chosen threshold (i.e. ≤ 0.05 / number of tests).
5. Format the output file: convert PLINK files to VCF, annotate the VCF, set the representation of male genotypes in the non-PAR to haploid.
6. Create an index for the VCF.

This pipeline uses:

- Python version 3.11 (with packages in requirements.txt and scripts in the bin folder)
- [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020) **Assumed in Path**
- [PLINK v2.00a3LM 64-bit Intel](https://www.cog-genomics.org/plink/2.0/) (22 Mar 2022) **Assumed in Path**
- [Laser v.1.03](https://csg.sph.umich.edu/chaolong/LASER/) (LASER's vcf2geno and trace are used)
- [bcftools/1.19](https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2) **Assumed in Path**
  - with 'fill-tags' and 'fixploidy' plugins

### Step 4. TOPMed Imputation Preparation

The 'topmed_prep.nf' pipeline removes positions with ≥ 20% difference in allele frequency between TOPMed and this dataset. The pipeline also splits the data into two random overlapping batches of 25,000 individuals to satisfy the current limit of TOPMed. The pipeline was run using the SLURM job scheduler on the Digital Research Alliance of Canada high-performance compute clusters. The pipeline is self-explanatory and includes the following steps:

1. Test the difference in allele frequency between [Bravo Freeze 8](https://topmed.nhlbi.nih.gov/topmed-whole-genome-sequencing-methods-freeze-8) and CaG data.
2. Remove variants with greater or equal differences in allele frequency than the threshold.
3. Create two overlapping lists of 25,000 individuals each, representing the entire cohort. Randomize the order of individuals and then divide them into two lists: the first 25,000 individuals and the last 25,000 individuals, ensuring an overlap between the lists. This procedure is done to satisfy TOPMed's current limits.
4. Create two sets of VCF files based on these lists.

This pipeline uses:

- Python version 3.11 (with packages in requirements.txt and scripts in the bin folder)
- [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020) **Assumed in Path**
- [bcftools/1.19](https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2) **Assumed in Path**

## Additional

Additionally, we provided the pipeline used for the creation of a reference panel in the desired format.

### Panel Creation

The 'reference_panel/panel.nf' pipeline creates a reference panel with a maximized number of shared positions between the target data and the references, while using standard data curation. This pipeline was run using the SLURM job scheduler on the Digital Research Alliance of Canada high-performance compute clusters. The pipeline is self-explanatory and includes the following steps:

1. Filter variants. Variants are filtered based on shared variant positions with the target dataset. Quality filters are applied, and variants are selected based on their allele frequency.
2. File concatenation. The filtered VCFs are concatenated.
3. Preform LD-Pruning. The VCF is linkage disequilibrium (LD) pruned.
4. Perform VCF -> GENO file conversion. The VCF is converted to GENO/SITE files format using vcf2geno.

This pipeline uses:

- [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020) **Assumed in Path**
- [Laser v.1.03](https://csg.sph.umich.edu/chaolong/LASER/) (LASER's vcf2geno is used for VCF conversion to .geno/.site file format)
