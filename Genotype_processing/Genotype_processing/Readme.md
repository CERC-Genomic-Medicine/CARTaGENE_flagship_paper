# Genotype processing 
 
## About
this folder contains the nextflow pipeline used to LiftOver Merge dans harmonize the CARTaGENE data (N=29,333 merged final version). The pipelines were run using SLURM job scheduler on the [Digital Research Alliance of Canada](https://alliancecan.ca/en) high performance compute clusters. The pipeline is a single script composed of multiple steps.



0) Prior to any data handling, samples for individuals who have withdrawn their consent must be removed
1) LiftOver from hg19 to Hg38 form [UCSC](https://hgdownload.soe.ucsc.edu/downloads.html)
2) Align genotyping Arrays to [Hg38](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000001405.28/) as (REF/ALT) 
3) Set haploid heterozygotes to missing
4) Deduplicating samples and merge
5) Test variants' agreement accross array using ancestry ajusted LRT
	5.1) Using [Lazer](https://csg.sph.umich.edu/chaolong/LASER/)'s trace program perform PCA projection
	5.2) Use the resulting Projection to ajust ancestry within test

## Necessary software

- Python version 3.11 (with packages in Requierement.txt)
- [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020) **Assumed in Path**
- [PLINK v2.00a3LM 64-bit (22 Mar 2022)](www.cog-genomics.org/plink/2.0/) **Assumed in Path**
- [LiftOver](https://genome-store.ucsc.edu/) **Assumed in Path**
- trace [Laser v.1.03](https://csg.sph.umich.edu/chaolong/LASER/)
- vcf2geno [Laser v.1.03](https://csg.sph.umich.edu/chaolong/LASER/)
- [bcftools/1.19](https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2) **Assumed in Path**


## Parameters (parameters used for CARTaGENE Flagship)

- **REF_projection**: Reference file for PCA projection. [gnomAD](https://gnomad.broadinstitute.org/downloads#v3-hgdp-1kg).

- **REF_ID**: Reference file for unrelated QCed individuals for PCA projection. [gnomAD](https://gnomad.broadinstitute.org/downloads#v3-hgdp-1kg).

- **chain**: liftOver chain file from hg19 to hg38. [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz).

- **genome**: Fasta reference file. [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000001405.28/).

- **K**: Number of principal components to consider for ancestry adjustment in LRT.
  ```plaintext
  K = 4
  ```
- **threshold**: Threshold in likelihood ratio test (LRT).
  ```plaintext
  threshold = "0.05"
  ```
- **Rsq**: LD Pruning Projection Reference.
  ```plaintext
  Rsq = 0.9
  ```
- **step**: LD Pruning step size.
  ```plaintext
  step = 100
  ```
- **window**: LD Pruning window size.
  ```plaintext
  window = 1000
  ```
- **max_allelefreq**: Maximum allele frequency for reference file filtering.
  ```plaintext
  max_allelefreq = 0.95
  ```
- **min_allelefreq**: Minimum allele frequency for reference file filtering.
  ```plaintext
  min_allelefrac = 0.05
  ```
- **MISSING_CUTOFF**: Missing data cutoff for reference file filtering.
  ```plaintext
  MISSING_CUTOFF = 0.001
  ```
- **King_cutOff**: Relatedness threshold to determine unrelated individuals for LRT.
  ```plaintext
  King_cutOff = "0.0442"
  ```

