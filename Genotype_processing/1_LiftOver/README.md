# Liftover

## GOAL

- Change array genotyping position from b37 to b38
- intermediary goal : change encoding to MT
- intermediary goal : Merge XY into X chromosome 

## Software

- plink (v1.90b6.21)
- LiftOver (UCSC)

## Input

- path=''
- arrays ${path}/gsa.${i}.final.v1.1.bed ${path}/gsa.${i}.final.v1.1.fam ${path}/gsa.${i}.final.v1.1.bim
- ** 17K ${path}/gsa.17k.final.v1.1.hg19.bed ${path}/gsa.17k.final.v1.1.hg19.fam ${path}/gsa.17k.final.v1.1.hg19.bim
- hg19ToHg38.over.chain.gz (from UCSC)

### Arrays

- 4224 5300 760 archi 17k
