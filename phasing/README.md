
# High-depth whole-genome sequencing (WGS) data phasing

This folder contains a collections of automated [Nextflow](https://www.nextflow.io/) pipelines used to phase SNVs, indels, and SVs called from CARTaGENE WGS data (N = 2,173).
The pipelines were run using SLURM job scheduler on the [Digital Research Alliance of Canada](https://alliancecan.ca/en) high performance compute clusters.

## Step 1. Preparing SNVs and indels for phasing.

The `prepare_SNVs_indels.nf` pipeline prepares VCF/BCF files with SNVs and indels for phasing. The pipeline is self-explanatory and includes the following steps:
1) Set individual GT values that are flagged as not 'PASS' to missing values. Note: input VCF/BCF files were already pre-filtered to contain only PASS variants.
2) For autosomal chromosomes and chrX PAR regions, force ploidy to 2.
3) For chrX non-PAR region, check ploidy for males and females; and set males as haploids (i.e. GT = 0 or GT = 1 instead of GT = 0/0 or GT = 1/1).
4) Recalculate AN, AC, and F_MISSING fields and remove monomorphic variants and variants with missigness >0.1 (recommended threshold by statistical phasing methods).
5) Left-align indels and deduplicate.

This pipeline uses [Vt](https://genome.sph.umich.edu/wiki/Vt) for left-aligning and [bcftools](https://samtools.github.io/bcftools/bcftools.html) for everything else.

## Step 2. Preparing SVs for phasing and merging them with SNVs and indels.

The `prepare_and_merge_SVs.nf` pipeline prepares VCF/BCF files with SVs for phasing and merges them with SNVs and indels from Step #1. The pipeline is self-explanatory and includes the following steps:
1) Apply checks, filters, and transformations on SVs similar to those from Step #1.
2) Fill out the missing REF field with the reference allele from the GRCh38 file.
3) Merge SVs with SNVs and indels. For chrX non-PAR region, check if ploidy in SVs matches to ploidy in SNVs and indels before merging.

This pipeline uses [bcftools](https://samtools.github.io/bcftools/bcftools.html).

## Step 3. Phasing.

The `phase_SNVs_indels_SVs.nf` pipeline phases BCF file with SNVs, indels, and SVs using [Shapeit v5.1.1](https://odelaneau.github.io/shapeit5/) and the [1000G+HGDP reference panel](https://doi.org/10.1101%2F2023.01.23.525248) (N~4,000).
In our benchmarking experiments, we confirmed that reference-based phasing imroved phasing accuracy. The pipeline is self-explanatory and contains all parameters used. Common (AF>0.001) and rare variants were phased in two steps following the [Shapeit v5.1.1](https://odelaneau.github.io/shapeit5/) recommendations. Because the WGS data was reasonably small (N = 2,173) we didn't chunk individual chromosomes when performing common variant phasing.

This pipeline uses [Shapeit v5.1.1](https://odelaneau.github.io/shapeit5/) and the corresponding resources; and [bcftools](https://samtools.github.io/bcftools/bcftools.html).

## (Optional) Benchmarking.

The `benchmark_beagle.nf`, `benchmark_phasing.nf`, and `benchmark_shapeit5.nf` automated pipelines were used to benchmark phasing accuracy with [Shapeit v5.1.1](https://odelaneau.github.io/shapeit5/), [Beagle v5.4](https://faculty.washington.edu/browning/beagle/beagle.html), and [Eagle v2.4.1](https://alkesgroup.broadinstitute.org/Eagle/). For benchmarking, we used trios from [1000 Genomes Project](Trios_1KG) and [Genome In The Bottle (GIAB)](Trios_GIAB). In summary, the pipelines append a trio's child to the study data and run statistical phasing on combined genotypes, one at a time. The child's phased genotypes are then compared to the trio-based phasing. This serves as a proxy for the phasing accuracy of the genotyoes of the other study individuals.
