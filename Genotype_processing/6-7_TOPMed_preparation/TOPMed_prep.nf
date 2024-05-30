#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process S1_Testing {
  label 'Testing'
  cache 'lenient'
  scratch false

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

process S2_Chromosome_split {
  label 'Chromosome_split'
  cache 'lenient'
  scratch true

  input:
  tuple path(batch), path(vcf), path(index)

  output:
  tuple path("*.vcf.gz"), path("*.vcf.gz.csi"), emit : batch1
  
  publishDir "${params.OutDir}/VCF_batches", pattern: "${params.name}_*", mode: "copy"
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

     S1 = S1_Testing(VCF,REF)
     S1.batchs.flatten().combine(S1.vcf).view()
     S2_Chromosome_split(S1.batchs.flatten().combine(S1.vcf))
}
