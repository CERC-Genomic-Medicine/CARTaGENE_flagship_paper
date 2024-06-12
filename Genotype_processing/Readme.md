# Genotype processing 

## About

This folder contains the nextflow pipeline used to LiftOver, Merge and harmonize the CARTaGENE data (N=29,333 merged final version). This folder was conceived to documents the preparation of CARTaGENE (CaG) data for imputation. 

<<<<<<< HEAD


=======
>>>>>>> 1b06d2172488ea895354d9e2dd0827d251912442
## Workflow order

```diff
- Important!
This pipeline uses plink version 1 and version 2, the two are called using plink and plink2 respectively. Importantly they are not interchangeable.
```

### Step 1. Array file pre-process
The 'preprocess.nf' pipeline removes individual with conscent withdrawn, performs LiftOver from b37 to b38 and aligns to the reference. The pipeline is self-explanatory and includes the following steps:

1) Remove samples from individuals with withdrawn consent.
2) LiftOver positions from b37 to b38.
3) Aligns to reference b38
4) Remove duplicate variant (based on chr/position/ref/alt)
4) set haploid heterozygote (male genotype in non-PAR) to missing

This pipeline uses :

- Python version 3.11 (with packages in Requierement.txt and bin folder)
- [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020) **Assumed in Path**
- [LiftOver](https://genome-store.ucsc.edu/) **Assumed in Path**

### Step 2. Merge

The 'merge.nf' pipeline deduplicate individual with the same ID and performs the merger of all array plink files. To deduplicate, individuals were removed for the array with the least variants genotyped. To merge all files, only variants present in all arrays were kept. The pipeline is self-explanatory and includes the following steps:

1) Removes individuals from arrays (in growing order of number of variant genotyped) until there is no shared individual names
2) Establishes a list of shared variants and keep only shared variants
3) Merges all array on keeping only shared variants

This pipeline uses :

- [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020) **Assumed in Path**

### Step 3. Harmonization

The 'harmonization.nf' pipeline removes position which display a batch effect based on array. To assess the batch effect two models are evaluated using LRT, one with only the projected principal components representing genetic ancestry and one with both the projected principal components and the Arrays as explanatory viariables. The pipeline is self-explanatory and includes the following steps:

1) PCA projection onto a reference parallelized. This step is divided in two, one instance is run first to obtain a reference PCA, as to no recompute it. This step was performed on a panel of [gnomAD's HGDP + 1KG callset](https://gnomad.broadinstitute.org/downloads#v3-hgdp-1kg)
2) Establish list of unrelated individuals using plink2's king-cutoff.
3) Test each variant (parallelized by chromosome). This step consists of fitting two linear regression models G ~ PC_1 + ... + PC_4 + dataset_1 + ... + dataset_n and G ~ PC_1 + ... + PC_4, where G = {0, 1, 2} - are individual genotypes, PC_i - principal components, dataset_j - dataset indicator variable. This step also computes F-test P-value and Likelihood Ratio Test P-value comparing the two models. With the assumption of random sampling, the dataset_i variables should not provide significant improvements compared to the nested model.'
4) Filtering the input files for variant with p-value below the chosen threshold (0.05 / number of test)
5) Formatting output file, this step converts from plink files to vcf, annotates the vcf, set the representation of males genotypes in the non-PAR to haploid.
6) Creates an index the vcf.

This pipeline uses :
- Python version 3.11 (with packages in Requierement.txt and bin folder)
- [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020) **Assumed in Path**
- [PLINK v2.00a3LM 64-bit Intel ](https://www.cog-genomics.org/plink/2.0/) (22 Mar 2022) **Assumed in Path ; Used to establish relationship based on King metric ** 
- [Laser v.1.03](https://csg.sph.umich.edu/chaolong/LASER/) (LASER's vcf2geno and trace are used)
- [bcftools/1.19](https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2) **Assumed in Path**

### Step 4. TOPMed imputation preparation 

The 'topmed_prep.nf' pipeline removes position with more than 0.2 differences in allele frequency between TOPMed and our data (due to poor imputation capacity for those allele), the pipeline also split the data in two random overlapping batches of 25k individuals to satisfy the current limit of TOPMed.  The pipeline is self-explanatory and includes the following steps:

1) Test the difference in allele frequency between bravo Freeze 8 and your file
2) Remove variants with greater differences than the threshold in allele frequency
3) create two overlapping list of individuals (randomized)
4) Create two sets of file based on this list.

This pipeline uses :

- Python version 3.11 (with packages in Requierement.txt and bin folder)
- [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020) **Assumed in Path**
- [bcftools/1.19](https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2) **Assumed in Path**

## Additionnal

Additionnaly we provided script used for the the creation of a reference pannel in the desired format.

### Panel creation
The 'reference_panel/panel.nf' pipeline creates a reference pannel with a maximal number of shared position between the target data and the references, while using standard data curation. The pipeline is self-explanatory and includes the following steps:

1) The variants are filtered on shared variant position with the target, quality filters are also imposed, and variant are also selected based on their allele frequency.
2) The variants filees are concatenated, ld-pruned, converted to geno/site format.

This pipeline uses :

- [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020) **Assumed in Path**
- [Laser v.1.03](https://csg.sph.umich.edu/chaolong/LASER/) (LASER's vcf2geno was used for vcf conversion to .geno/.site file format)

