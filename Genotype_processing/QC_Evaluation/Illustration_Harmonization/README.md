# Illustration Harmonization
This folder contains the scripts and documentation used to prepare some analysis of the harmonization stpe of the genotype processing. These script relies on the whole process up to and including the [Harmonization step](https://github.com/CERC-Genomic-Medicine/CARTaGENE_flagship_paper/tree/main/Genotype_processing/5_Harmonization).

Result can be found at [Pheweb](https://cerc-genomic-medicine.ca/pheweb/cartagene/about) (in the next update)

## Parameters
### Software used
- [PLINK v2.00a3LM 64-bit (22 Mar 2022)](www.cog-genomics.org/plink/2.0/)
- [python 3.11.5](https://www.python.org/downloads/release/python-3115/)
- pandas 2.2.1
- numpy 1.22.2
- glob
- os
- subprocess
- matplotlib 3.5.1
### Files 
- Phase_file - File detailing the recruitement phase of each sample
CaG_unmerged - path to **FOLDER** containing the genotyping arrays PLINK binary format
Variants - Variant to be compared (mahattan with and without those variants)
CaG_merged_unfiltered - path to CaG unfiltered merged PLINK binary format (with bed bim fam files)
Unrelated_individuals - path to file containing the list of unrelated individuals (obtained in step 5 Harmonization)
### Other parameters
Array_list - ["archi","760","5300","4224","17k"]
PC - '5' 																# number of PCs as covariates to be used in association testing
threads - '5'

## Steps
### Step 1 Create plink2's input files for variants association analysis
Within Illustration_Harmonization.py
- read table IID -> Phase 
- create table IID - > Array
- Create Matrix Phase_Array vs Phase_Array
- Create PCs covariates using [PLINK v2.00a3LM 64-bit (22 Mar 2022)](www.cog-genomics.org/plink/2.0/)'s --PCA function on only unrelated individuals
### Step 2 perform variant association for each 'phenotypes'  
Within Illustration_Harmonization.py calling plink using subprocess
- Perform variant association for each 'phenotypes' using [PLINK v2.00a3LM 64-bit (22 Mar 2022)](www.cog-genomics.org/plink/2.0/)'s glm function on only unrelated individuals
### Step 3 Illustrate
Within Illustration_Harmonization.py
- Illustrate Manhattan plot wih all variant 
- Illustrate Manhattan plot wih without non harmonious variant 
- **Perform these illustration for all comparisons**
### Step 4 Illustrate Allele frequency comparison
Using frq_comparison.py 
- Parameters for autosome & PAR
  - i : all Autosome and PAR output of [Harmonization test](https://github.com/CERC-Genomic-Medicine/CARTaGENE_flagship_paper/tree/main/Genotype_processing/5_Harmonization/Harmonization_test)
  - n : number of test total (including non-PAR)
  - t : Autosomal chromosome and chromosome X PAR
  - A : "archi" "760" "5300" "4224" "17k"
  - out : output name
- Parameters for Male chromosome X PAR
  - i : Male PAR output of [Harmonization test](https://github.com/CERC-Genomic-Medicine/CARTaGENE_flagship_paper/tree/main/Genotype_processing/5_Harmonization/Harmonization_test)
  - n : number of test total (including non-PAR)
  - t : Male chromosome X PAR
  - A : "archi" "760" "5300" "4224" "17k"
  - out : output name
- Parameters for Female chromosome X PAR
  - i : Femal PAR output of [Harmonization test](https://github.com/CERC-Genomic-Medicine/CARTaGENE_flagship_paper/tree/main/Genotype_processing/5_Harmonization/Harmonization_test)
  - n : number of test total (including non-PAR)
  - t : Female chromosome X PAR
  - A : "archi" "760" "5300" "4224" "17k"
  - out : output name
