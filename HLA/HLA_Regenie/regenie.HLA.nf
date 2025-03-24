/*
* Author: Peyton McClelland <peyton.mcclelland@mail.mcgill.ca>
* Year: 2024
*/

process step1 {
  label 'step1'
  //debug true

  input:
  tuple path(step1_geno_pgen), path(step1_geno_psam), path(step1_geno_pvar), path(pheno), path(covar), val(phenoCol)
  
  output:
  tuple val(phenoCol), path("regenie_step1_out*.loco"), path("regenie_step1_out_pred.list"), emit: output
  path("*.log"), emit: step1_log

  publishDir "${params.out_dir}/step1/logs", pattern: "*.log", mode: "copy", saveAs: { fn -> "step1_out.${phenoCol}.log"}
  

  script:
  """
  echo step1_geno_pgen
  regenie \
  --step 1 \
  --pgen ${step1_geno_pgen.getBaseName()} \
  --covarFile $covar \
  --phenoFile $pheno \
  --phenoCol $phenoCol \
  --catCovarList ${params.cat_covar} \
  --qt \
  --bsize ${params.bsize} \
  --threads ${params.threads_step1} \
  --out regenie_step1_out
  """
}

process prepare_step2_geno {
  debug true

  input:
  path(step2_geno_input)

  output:
  tuple path("step2_geno_fileset.pgen"), path("step2_geno_fileset.pvar"), path("step2_geno_fileset.psam")

  script:
  """
  plink2 --vcf $step2_geno_input --make-pgen --double-id --out step2_geno_fileset
  """

}

process step2 {
  label 'step2'
  debug true

  input:
  tuple val(phenoCol), path(loco), path(pred), path(step2_geno_pgen), path(step2_geno_pvar), path(step2_geno_psam), path(pheno), path(covar)

  output:
  tuple val(phenoCol), path("*.regenie"), emit: pheno_result
  path("*.log"), emit: step2_log

  publishDir "${params.out_dir}/step2/logs", pattern: "*.log", mode: "copy", saveAs: { fn -> "step2_out.${phenoCol}.log"}
  publishDir "${params.out_dir}/step2/assoc_results", pattern: "*.regenie", mode: "copy"

  script:
  """
  echo $phenoCol
  regenie \
  --step 2 \
  --pgen ${step2_geno_pgen.getBaseName()} \
  --phenoFile $pheno \
  --phenoCol $phenoCol \
  --covarFile $covar \
  --catCovarList ${params.cat_covar} \
  --bsize ${params.bsize} \
  --qt \
  --firth --approx \
  --pThresh ${params.p_threshold} \
  --pred regenie_step1_out_pred.list \
  --out ${params.out_prefix}
  """
}

process pre_merge {
  debug true

  input:
  tuple val(phenoCol), path(regenie)

  output:
  path("*.regenie.to_merge.txt")

  script:
  """
  head -n1 ${regenie} | sed 's/^/PHENOTYPE /' > header.txt
  tail -n+2 ${regenie} | sed 's/^/${phenoCol} /'  > to_merge.no_header.txt
  cat header.txt to_merge.no_header.txt > ${regenie}.to_merge.txt

  """
} 

process merge_result {
  label 'step2'
  // debug true

  input:
  path("*.regenie.to_merge.txt")

  output:
  path("*.MERGED.regenie")

  publishDir "${params.out_dir}/step2/", pattern: "*MERGED.regenie", mode: "copy"

  script:
  """
  cat *.regenie.to_merge.txt | head -n1 > header.txt
  cat *.regenie.to_merge.txt | grep -vw "^PHENOTYPE" >> no_header.txt
  cat header.txt no_header.txt >> ${params.out_prefix}.MERGED.regenie
  """
}
workflow {

step1_geno = Channel.fromPath(params.step1_genotype) | map { file -> [ file, 
								       "${file.getParent()}/${file.baseName}.psam", 
								       "${file.getParent()}/${file.baseName}.pvar" ]}

pheno = Channel.fromPath(params.phenotype_file, checkIfExists:true) 
covar = Channel.fromPath(params.covariate_file, checkIfExists:true)

pheno_cols = Channel.from(file(params.phenotype_summary)).splitCsv(skip: 1 ).map { row -> [ row[0] ] } 

input_step1 = step1_geno.concat( pheno, covar ).collect().combine( pheno_cols )
s1_out = step1( input_step1 )

step2_geno_input = Channel.fromPath(params.step2_genotype)
step2_geno_fileset = prepare_step2_geno( step2_geno_input )

input_step2 = s1_out.output.combine(step2_geno_fileset.concat(pheno, covar).collect() )

s2_out = step2( input_step2 )

merge_s2_result = merge_result( pre_merge( s2_out.pheno_result ).collect() )
}
