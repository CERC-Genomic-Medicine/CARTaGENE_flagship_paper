# HLA analysis in CARTaGENE

1. HLA region extraction using `hla_extract.nf` pipeline at [github.com/pmccl027/hla_analysis/tree/main/extract_hla_region](https://github.com/pmccl027/hla_analysis/tree/main/extract_hla_region)
2. HLA allele calling via [HLA*LA](https://pubmed.ncbi.nlm.nih.gov/30942877) with `hla_la.nf` at [github.com/pmccl027/hla_analysis/tree/main/hla_calling](https://github.com/pmccl027/hla_analysis/tree/main/hla_calling)
3. HLA imputation panel generated and implemented using [HLA-TAPAS](https://github.com/immunogenomics/HLA-TAPAS) pipeline (Luo et al. 2021, NatGen)
4.  `HLA_Regenie/` | nextflow pipeline to implement HLA-specific phenome-wide association testing through Regenie
5. `Cartagene_HLA_figures.ipynb` | code used to generate S11 and S17 in paper
