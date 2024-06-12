#!/usr/bin/env nextflow
/*
* AUTHOR: Vincent Chapdelaine
* VERSION: 3.0
* YEAR: 2024
*/

nextflow.enable.dsl = 2

// Pre-requisites/assumptions:
// 1. Script was designed for genotyping arrays.
// 2. bed file originate frome multiple arrays
// 3. Make sure that plink is installed/loaded on your system


// Minimal usage : nextlfow run harmonization.nf --target_data "/path/to/data/target_data.bed"  --reference_pannel_vcfs "/path/to/pannel_*.vcf.gz" --reference_pannel_qc_id /path/to/data/pannel_qc.txt --OutDir 'path/to/Dir/'  --laser_dir "path/to/laser_directory/"


// Minimal Parameters

params.target_data_variants = "/path/to/data/target_data.bim" // Absolute path to genotyping array in plink bed file (with bim and fam files in the same folder)
params.reference_pannel_vcfs = "/path/to/pannel_*.vcf.gz"      // Absolute path Reference File for PCA projection (assumes the presence of .csi index file in the same directory)
params.reference_pannel_qc_id = '/path/to/data/pannel_qc.txt'                            // Absolute path to Reference File with individual passing QC (and unrelated)
params.OutDir = 'path/to/Dir/'                  // Absolute path to the output directory
params.laser_dir="path/to/laser_directory/"     // Absolute path to Umichigan's LASER program for vcf2geno


// Filtering Parameters

Rsq=0.9                   // LD Pruning Projection Reference
step=100                  // LD Pruning Projection Reference
window=1000               // LD Pruning Projection Reference
max_allelefreq=0.95       // Reference File filtering
min_allelefreq=0.05       // Reference File filtering
MISSING_CUTOFF=0.001      // Reference File filtering



process reference_filter_base {
  label 'reference_filter_base'
  cache 'lenient'
  scratch false
  memory { 32.GB * task.attempt }
  errorStrategy 'retry'
  errorStrategy { task.exitStatus == 137 || task.exitStatus == 143 ? 'retry' : 'terminate' }
  maxRetries 5
  cpus 5

  input:
  tuple path(vcf), path(index), path(bim)
  each path(pannel_samples)

  output:
  path("temp_*.vcf.gz"), emit : refrence_vcfs
  path("temp_*.vcf.gz.csi"), emit: refrence_index

  script:
  """
  cut -f 2 ${bim} | sed 's/:/\t/g' - | awk '{print \$1, \$2-1, \$2}' - | sed 's/ /\t/g' - > shared.bed
  bcftools view -f 'PASS' -q${params.min_allelefreq} -Q${params.max_allelefreq} -S ${pannel_samples} -i "F_MISSING < ${params.MISSING_CUTOFF}" -R shared.bed ${vcf} -Ou  --threads 5|  bcftools annotate -x INFO,^GT - -Oz -o temp_${vcf.getBaseName()}.gz --threads 5
  bcftools index temp_${vcf.getBaseName()}.gz
  """
}

// 1) Create Reference for PCA projection 2) First batch of PCA projection (to create and reuse REF PCA) 
process reference_localization_conversion {
  label 'reference_localization_conversion'
  scratch false
  cache 'lenient'
  memory { 32.GB * task.attempt }
  errorStrategy 'retry'
  errorStrategy { task.exitStatus == 137 || task.exitStatus == 143 ? 'retry' : 'terminate' }
  maxRetries 5
  cpus 5

  input:
  path(bim)
  path(refrence)
  path(index)


  output:
   path("pannel.geno"), emit: refrence
   path("pannel.site"), emit: refrence_site


  publishDir "${params.OutDir}", pattern: "pannel.geno", mode: "copy"
  publishDir "${params.OutDir}", pattern: "pannel.site", mode: "copy"
  script:
  """
  bcftools concat ${refrence} -n -Oz -o refrence.vcf.gz  # Concatenate files
  plink2 --vcf refrence.vcf.gz --indep-pairwise ${params.window} ${params.step} ${params.Rsq} --out refrence.LD_prune --set-all-var-ids '@:#:\$r:\$a' --new-id-max-allele-len 145 --threads $5# ID var to be pruned
  plink2 --vcf refrence.vcf.gz --set-all-var-ids 'chr@:#:\$r:\$a' --new-id-max-allele-len 145 --sort-vars --threads 5  --exclude refrence.LD_prune.prune.out --make-pgen --out tmp # Remove pruned variants
  plink2 --pfile tmp --export vcf bgz id-paste=iid --out refrence.LD_prune --threads 5
  bcftools index refrence.LD_prune.vcf.gz
  ${params.laser_dir}/vcf2geno/vcf2geno --inVcf refrence.LD_prune.vcf.gz --out pannel # convert to geno/site format

    """
}

workflow {
  //Input files
     reference_pannel_vcfs = Channel.fromPath(params.reference_pannel_vcfs).map(f -> [f, file("${f}.csi")])
     target_data_variants = Channel.fromPath(params.target_data_variants) 
     reference_pannel_qc_id = Channel.fromPath(params.reference_pannel_qc_id)
  // Ressource file
        base = reference_filter_base(reference_pannel_vcfs.combine(target_data_variants), reference_pannel_qc_id)
        reference = reference_localization_conversion(target_data_variants, reference_filter_base.refrence_vcfs.collect(sort:true), reference_filter_base.refrence_index.collect(sort:true))

}