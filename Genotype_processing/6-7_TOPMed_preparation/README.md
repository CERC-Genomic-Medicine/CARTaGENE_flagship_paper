# Genotype processing 
 
## About
This folder contains the files used to produce 2 x23 chromosome specific vcf reprensenting 2 overlapping batches of 25000 individuals. These vcfs The pipelines were run using SLURM job scheduler on the [Digital Research Alliance of Canada](https://alliancecan.ca/en) high performance compute clusters. The pipeline is a single script composed of multiple steps.


### Step 6 (following genotype array merger and harmonization)
- Test the difference in allele frequency between bravo Freeze 8 and your file
- Remove variants with greater differences than the threshold in allele frequency
- create two overlaping list of individuals (randomized)

### Step 7
- Create individual files per chromosome for each batches

## Necessary software

- Python version 3.11 (with packages in Requierement.txt)
- [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020) **Assumed in Path**
- [PLINK v2.00a3LM 64-bit (22 Mar 2022)](www.cog-genomics.org/plink/2.0/) **Assumed in Path**

## Parameters (used in CaG flagship)

- Reference ([TOPMed freeze 8](https://www.nature.com/articles/s41586-021-03205-y)
- Threshold (0.2)
