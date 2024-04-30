# Goal :
- Align A1 to Reference (set strand to +, force REF/ALT designation  
- remove duplicate (keep variant with highest call rate)  
- Renames variant to chr:POS:REF:ALT convention
 **P.S. X chromosome is merged

## Input files :
- bed file set {Arrays}_Hg38_Xmerged
- Arrays Included : 17k 5300 4224 760 archi
- Homo_sapiens.GRCh38.fa (Software Ressource file)


## Output files :
- {Arrays}_Hg38_renamed.bed
- {Arrays}_Hg38_renamed.bim
- {Arrays}_Hg38_renamed.fam

## Used Software :
- plink2reference.py (from https://github.com/CERC-Genomic-Medicine/scripts)
- plink_dupvar_prioritize.py (from https://github.com/CERC-Genomic-Medicine/scripts)
- plink/1.9
