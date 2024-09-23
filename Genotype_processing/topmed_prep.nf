#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


// Pre-requisites/assumptions:
// 1. Script assumes ploidy does not requiere correction
// 2. Make sure that bcftools is installed/loaded on your system

// minimal execution:
// nextflow run TOPMed_prep.nf --vcf "path/to/file.vcf.gz"  --ref "path/to/file.vcf.gz"  --name "CaG" --OutDir "path/to/dir"


//Files
params.vcf="path/to/file.vcf.gz"            // vcf of the cohort to be imputed
params.ref="path/to/file.vcf.gz"           // vcf of the TOPMed bravo freeze (assumes a .tbi index is present
params.OutDir="path/to/dir"          // output directory
params.name='CaG'                       // name of project


//Parameters
params.Batch_length=25000           // number of individual per batch (currently TOPMed server has a 25000 individual limit)
params.Threshold=0.2                 // Maximum difference in frequency


process testing {
  label 'testing'
  cache 'lenient'
  scratch false

  maxRetries 3
  scratch false
  errorStrategy 'retry'
  errorStrategy { task.exitStatus == 137 || task.exitStatus == 143 ? 'retry' : 'terminate' }
  
  executor "slurm"
  clusterOptions "--account=rrg-vmooser"

  cpus 1
  memory { 5.GB * task.attempt }
  time "4h"

  input:
  path(vcf)
  tuple path(Bravo), path(index) 

  output:
  tuple path("${params.name}.vcf.gz"), path("${params.name}.vcf.gz.csi"), emit: vcf
  path("batch*.txt"), emit : batchs
  
  publishDir "${params.OutDir}/Batches", pattern: "batch*.txt", mode: "copy"
  publishDir "${params.OutDir}/vcf", pattern: "${params.name}.vcf.*", mode: "copy"
  script:
  """
  plink --vcf ${vcf}  --keep-allele-order --freq --out ${vcf.getBaseName()} # Make Allele frequency
  plink_freq_vs_topmed.py -s ${vcf.getBaseName()}.frq -t ${Bravo} -m 0.2 -o Diff_${params.Threshold} # Test frequency concordance
  bcftools view ${vcf} --exclude ID=@Diff_${params.Threshold}.af_diff.txt -o ${params.name}_TEMP.vcf.gz -Oz 
  bcftools index ${params.name}_TEMP.vcf.gz
  bcftools +fill-tags ${params.name}_TEMP.vcf.gz -Oz -o ${params.name}.vcf.gz
  bcftools index ${params.name}.vcf.gz
  bcftools query -l ${params.name}.vcf.gz > samples
  shuf samples > samples_random.txt
  head -n ${params.Batch_length} samples_random.txt > batch1.txt
  tail -n ${params.Batch_length} samples_random.txt > batch2.txt
  """
}

process chromosome_split {
  label 'chromosome_split'
  cache 'lenient'
  scratch true

  label 'testing'
  cache 'lenient'
  scratch false

  maxRetries 3
  scratch false
  errorStrategy 'retry'
  errorStrategy { task.exitStatus == 137 || task.exitStatus == 143 ? 'retry' : 'terminate' }
  
  executor "slurm"
  clusterOptions "--account=rrg-vmooser"

  cpus 1
  memory { 5.GB * task.attempt }
  time "4h"


  input:
  tuple path(batch), path(vcf), path(index)

  output:
  tuple path("*.vcf.gz"), path("*.vcf.gz.csi"), emit : batch1
  
  publishDir "${params.OutDir}", pattern: "${params.name}_*", mode: "copy"
  script:
  """
  for chromosome in chr{1..22} chrX ; do 
  bcftools view -r \${chromosome} ${vcf} -S ${batch}  -Oz -o ${params.name}_\${chromosome}_${batch.getBaseName()}}.vcf.gz; 
  bcftools index ${params.name}_\${chromosome}_${batch.getBaseName()}.vcf.gz ; 
done
  """
}

workflow {
  //Input files
     VCF = Channel.fromPath(params.vcf)
     REF = Channel.fromPath(params.ref).map(f -> [f, file("${f}.tbi")])

     S1 = testing(VCF,REF)
     S1.batchs.flatten().combine(S1.vcf).view()
     chromosome_split(S1.batchs.flatten().combine(S1.vcf))
}
