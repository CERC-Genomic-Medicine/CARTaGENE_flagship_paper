#!/usr/bin/env nextflow
/*
* AUTHOR: Vincent Chapdelaine
* VERSION: 3.0
* YEAR: 2024
*/

nextflow.enable.dsl = 2


// Pre-requisites/assumptions:
// 1. Script was designed for genotyping arrays.
// 2. Inputs are in plink bed files
// 3. Make sure that plink is installed/loaded on your system

// How to execute:
// nextflow run liftover.nf --beds_initial "/home/user/data/*.bed" --Consent "/path/to/consent.txt" --genome "/home/user/reference/GRCh38.fa" --chain "/home/user/reference/hg19ToHg38.over.chain.gz --OutDir /home/user/out"

params.beds_initial = "/path/to/data/*.bed" // Absolute path to genotyping array in plink bed file (with bim and fam files in the same folder)
params.Consent = '/path/to/file.txt'        // Absolute path to individual having withdrawn their conscent ( format : FID IID ; with header #FID IID )
params.chain = '/path/to/chain.gz'          // Absolute path to liftover chain file availlable at http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/
params.genome = '/path/to/data/*.fa'        // Absolute path to the genome of interest in fasta format
params.OutDir = 'path/to/Dir/'              // Absolute path to the output directory


process consent_remove {
  cache 'lenient'

  input:
  tuple path(bed), path(bim), path(fam)
  each path(Consented)  // Samples without consent

  output:
  tuple path("${bim.getBaseName()}_consented.bed"), path("${bim.getBaseName()}_consented.bim"), path("${bim.getBaseName()}_consented.fam"), emit: Preprocessed

  script:
  """
  plink -bfile ${bim.getBaseName()} --remove ${Consented} --merge-x 'no-fail' --make-bed --output-chr MT --out ${bim.getBaseName()}_consented
  """
}

process liftover {
  cache 'lenient'

  input:
  tuple path(bed), path(bim), path(fam)
  each path(chain)

  output:
  tuple path("${bim.getBaseName()}_Hg38_Xmerged.bed"), path("${bim.getBaseName()}_Hg38_Xmerged.bim"), path("${bim.getBaseName()}_Hg38_Xmerged.fam"), emit: Lifted

  script:
  """
  awk '{print "chr" \$1, \$4-1, \$4, \$2}' ${bim} > ${bim.getBaseName()}_bedfile           ## Produces a bed file of all variants
  liftOver ${bim.getBaseName()}_bedfile ${chain} ${bim.getBaseName()}_mapfile ${bim.getBaseName()}_unmappedfile   ## Produces Map file and unmapped variant file
  grep ^chr[0-9A-Za-z]*_ ${bim.getBaseName()}_mapfile | cut -f 4 > ${bim.getBaseName()}_excludefile    ## Identifies Alternate contig mod A-Z a-z
  grep -v '^#' ${bim.getBaseName()}_unmappedfile | cut -f 4 >> ${bim.getBaseName()}_excludefile              ## Total list of variant to be excluded
  grep -v ${bim.getBaseName()}_excludefile ${bim.getBaseName()}_mapfile > ${bim.getBaseName()}_mapfile_final          ## Remove excluded variant

  plink --bfile ${bim.getBaseName()} \
  --exclude ${bim.getBaseName()}_excludefile \
  --make-bed \
  --not-chr MT \
  --out tmp_${bim.getBaseName()}

  plink --bfile tmp_${bim.getBaseName()} \
  --update-chr ${bim.getBaseName()}_mapfile_final 1 4 \
  --update-map ${bim.getBaseName()}_mapfile_final 3 4 \
  --make-bed \
  --output-chr chrMT \
  --out ${bim.getBaseName()}_Hg38_Xmerged
  """
}

// Align to the proper reference & set HH to missing

process alignement_hhmissing {
  cache 'lenient'

  input:
  tuple path(bed), path(bim), path(fam)
  each path(fasta)

  output:
  tuple path("${bim.getBaseName()}_hh.bed"), path("${bim.getBaseName()}_hh.bim"), path("${bim.getBaseName()}_hh.fam"), emit: Aligned
  
  publishDir "${params.OutDir}/S2", pattern: "*_hh.bim", mode: "copy"
  publishDir "${params.OutDir}/S2", pattern: "*_hh.bed", mode: "copy"
  publishDir "${params.OutDir}/S2", pattern: "*_hh.fam", mode: "copy"

  
  script:
  """
  plink2reference.py -b ${bim} -f ${fasta} -o ${bim.getBaseName()}                             #Produce documents needed for correct strand-flip and REF/ALT alteration
  plink --bfile ${bim.getBaseName()}  --exclude ${bim.getBaseName()}.remove.txt --make-bed --output-chr chrMT --out temp     #Removes impossible to adjust SNPs
  plink --bfile temp --flip ${bim.getBaseName()}.strand_flip.txt --make-bed --output-chr chrMT --out temp2                            # Strand filp
  plink --bfile temp2 --a1-allele ${bim.getBaseName()}.force_a1.txt --make-bed --output-chr chrMT --out ${bim.getBaseName()}_aligned         # Forces the REF/ALT Designation

  plink --bfile ${bim.getBaseName()}_aligned  --list-duplicate-vars                                                      #Identifies duplicate variant Position:REF:ALT
  plink --bfile ${bim.getBaseName()}_aligned  --freq counts
  plink_dupvar_prioritize.py -d plink.dupvar -c plink.frq.counts -o ${bim.getBaseName()}_exclude_dup.txt              #Identifies the duplicate variant with less missingness
  plink --bfile ${bim.getBaseName()}_aligned --exclude ${bim.getBaseName()}_exclude_dup.txt --output-chr chrMT --keep-allele-order --make-bed --out ${bim.getBaseName()}_aligned_dedup # Removes duplicates

  awk '{print \$1":"\$4":"\$5":"\$6 , \$2}' ${bim.getBaseName()}_aligned_dedup.bim > rename.txt                                                   # Create list of old var ID and new var id (chr:pos:REF:ALT)
  plink --bfile ${bim.getBaseName()}_aligned_dedup --update-name rename.txt 1 2  --make-bed --keep-allele-order --output-chr chrMT --out ${bim.getBaseName()}_renamed     # renames

  plink --bfile ${bim.getBaseName()}_renamed --split-x 'b38'  --keep-allele-order --make-bed --output-chr chrMT --out ${bim.getBaseName()}_renamed_Xsplit
  plink --bfile ${bim.getBaseName()}_renamed_Xsplit --keep-allele-order --set-hh-missing --make-bed --output-chr chrMT --out ${bim.getBaseName()}_renamed_Xsplit_temp
  plink --bfile ${bim.getBaseName()}_renamed_Xsplit_temp --merge-x --keep-allele-order --make-bed --output-chr chrMT --out ${bim.getBaseName()}_hh
  """
}


workflow {
  //Input files
     G_Arrays = Channel.fromPath(params.beds_initial).map(f -> [f, file("${f.getParent()}/${f.getBaseName()}.bim"), file("${f.getParent()}/${f.getBaseName()}.fam")])
     Consent = Channel.fromPath(params.Consent)
  // Ressource file
     Chain = Channel.fromPath(params.chain)
     Genome = Channel.fromPath(params.genome)
  // Steps
    //Step 0
        S0 = Consent_remove(G_Arrays, Consent)
    //Step 1
        S1 = liftover(S0.Preprocessed, Chain)
    // Step 2-3
        S2 = alignement_hhmissing(S1.Lifted, Genome)

}
