# Harmonization
This folder contains the scripts and documentation used to identify the list of variant failling the likelihood ratio test comparing two linear regression models G ~ m + PC1 + ... + PC4 + A2 + ... + A5 and G ~ m + PC1 + ... + PC4, where PCs are the Principal components of this projection and A are the arrays as binairy variables. The failling variant are considered to be non-harmonious across arrays. This test is performed using the of CaG merge samples' projection on [gnomAD's HGDP+1KG dataset](https://gnomad.broadinstitute.org/downloads#v3-hgdp-1kg).This process relies on the [array merger step](https://github.com/CERC-Genomic-Medicine/CARTaGENE_flagship_paper/tree/main/Genotype_processing/4_Merge).  

# Overview
This process is divided into two, PCA_Projection and Harmonization_test
## Steps
### Step 1 - PCA projection
To account for ancestry, samples are mapped onto the PCA space of gnomAD's HGDP+1KG Dataset .
### Step 2 - Harmonization Test
Perform the likelihood ratio test for all variant and filter the CaG file according to the threshold.
