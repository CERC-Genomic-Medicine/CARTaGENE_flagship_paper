# Pre-Processing
This folder contains the documentation and script used to create CaG input data compliant with consent status at the moment of writing. The bash companion of this file reflect the methods used on the [Digital Research Alliance of Canada](https://alliancecan.ca/en) high performance compute clusters. 

## Parameters
### Software used
PLINK v1.90b6.21 64-bit (19 Oct 2020) (availlable at https://www.cog-genomics.org/plink/)  

### Input files
{Bedfiles} - CaG data v.1.1 in plink .bed format (including .bim .fam companion file) **per Genotyping arrays**.  
{Removed_samples} -  list of consent retracted individuals **per Genotyping array**  

## Step 
### Pre-processing
The script named preprocessing.sh remove individuals with withdrawn consent. To this end this script uses :
- [PLINK v1.90b6.21 64-bit](https://www.cog-genomics.org/plink/)
