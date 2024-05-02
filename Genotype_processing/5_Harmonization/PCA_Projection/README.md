# PCA projection (sub-process of the Harmonization
This folder contains the scripts and documentation used to create the PCA projection requiered for the Harmonization test, a likelihood ratio test comparing two linear regression models G ~ m + PC1 + ... + PC4 + A2 + ... + A5 and G ~ m + PC1 + ... + PC4, where PCs are the Principal components of this projection and A are the arrays as binairy variables. This projection is performed on [gnomAD's HGDP+1KG dataset](https://gnomad.broadinstitute.org/downloads#v3-hgdp-1kg).This process relies on publicly availlable data and the output of the [array merger step](https://github.com/CERC-Genomic-Medicine/CARTaGENE_flagship_paper/tree/main/Genotype_processing/4_Merge).

## Parameter
### Input files
Bed file - Merged CaG data v.1.1 in plink .bed format (including .bim .fam companion file)
Shared Variant files - List of shared variant in the CaG data in chr:pos:REF:ALT format.
Gnomad dataset 
- VCF files of [gnomAD's HGDP+1KG dataset v3.1.2](https://gnomad.broadinstitute.org/downloads#v3-hgdp-1kg)  
- list of QC passing individuals    

Long-Range LD bed file - long range regions as defined in [Anderson et al. 2010](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3861250/) availlable in hg38 [here](https://github.com/CERC-Genomic-Medicine/Regenie_nextflow/blob/3a7f3e3bb5e5a9cbc8a075b54917c423237cdc28/util/Low_complexity_regions/HG38.Anderson2010.bed.gz)  

### Used Software
- trace [Laser v.1.03](https://csg.sph.umich.edu/chaolong/LASER/)
- vcf2geno [Laser v.1.03](https://csg.sph.umich.edu/chaolong/LASER/)
- [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020)
- [PLINK v2.00a3LM 64-bit (22 Mar 2022)](www.cog-genomics.org/plink/2.0/)
- [bcftools/1.19](https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2)

### Output files :  
Coordinate file - File of the coordinate of CaG samples within the PCA space of gnomAD's HGDP+1KG.

## Steps
### Step 1 - Create Shared position file
1) this step converts list of variant to a bed file.

### Step 2 - Process Reference dataset
This step remove low quality/redundant genomic elements and sample with low quality sequencing from the dataset and convert it to the geno/site format used by laser.
1) Using bcftools' options : Remove non-passing quality check variant (-fPASS option), low variance alleles (using -q and -Q options), high missingness alleles ( -i option), QC failling samples (as identified by gnomAD, using -S option), and alleles in regions not shared with CaG dataset (using -R option).
   - Step Parameter - Allele selected
     - Allele frequency selected [0.05 - 0.95] (low variance alleles)
     - Allele passing quality check
     - Allele with low missingness, < 0.001.
     - Allele with CaG overlapping position
   - Step Parameters -  Samples  selected
     - Samples passing QC defined by gnomAD
2) Contatenate resulting file
3) Prune Allele based on independance using plink2 (--indep-pairwise funciton)
   - Step Parameter
     - window size 1000
     - step size 100
     - r^2 threshold 0.9
4) Remove long range Linkage Disequilibrium unsing plink2's exclude
   - Step Parameter
     - long range regions as defined in [Anderson et al. 2010](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3861250/) availlable in hg38 [here](https://github.com/CERC-Genomic-Medicine/Regenie_nextflow/blob/3a7f3e3bb5e5a9cbc8a075b54917c423237cdc28/util/Low_complexity_regions/HG38.Anderson2010.bed.gz)  
5) Convert to geno/site format used by laser, using vcf2geno

### Step 3 - Process CaG dataset
this step convert CaG data to the geno/site format used by laser.
1) Convert from bed file to vcf using [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/) (19 Oct 2020)
3) Convert to geno/site format used by laser, using vcf2geno

### Step 4 - Project CaG
This step uses laser's trace to create a PCA space for the reference and projects CaG onto it using overlapping alleles.
1) In paralele project CaG dataset onto HGDP+1KG dataset
   - Step Parameter - Trace Parameter
     - Number of PCs to compute (-k) 20
     - Number of dimension of the sample-specific PCA map to project (-K) 20
