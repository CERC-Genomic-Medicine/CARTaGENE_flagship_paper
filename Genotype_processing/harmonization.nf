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


// Minimal usage : nextlfow run harmonization.nf --unharmonized_data "/path/to/data/unharmonized_data.bed"  --reference_pca_projection_pannel "/path/to/pannel.geno" --list_array_origin /path/to/data/samples_array_list.txt --OutDir 'path/to/Dir/'  --laser_dir "path/to/laser_directory/"


// Minimal Parameters

params.unharmonized_data = "/path/to/data/unharmonized_data.bed" // Absolute path to genotyping array in plink bed file (with bim and fam files in the same folder)
params.reference_pca_projection_pannel = "/path/to/pannel.geno"      // Absolute path Reference File for PCA projection (assumes the presence of .site file in the same directory)
params.list_array_origin = "/path/to/data/samples_array_list.txt" // Absolute path to the sample - array correspondance (format sample\tArray)
params.OutDir = 'path/to/Dir/'                  // Absolute path to the output directory
params.laser_dir="path/to/laser_directory/"     // Absolute path to Umichigan's LASER program


// Parallelisation option
params.ind_trace_batch='100'      // N batches for projection [2,100]

// LRT Parameters
params.K=4                       // Number of PC to consider for Genetic Ancestry in LRT
params.threshold="0.05"          // Pvalue threshold in LRT
params.King_cutOff="0.0442"      // Relatedness Threshold to determine unrelated individuals for LRT


//Constant
  params.PAR1='chrX:10001-2781479'       # PAR1 hg38
  params.PAR2='chrX:155701383-156030895' # PAR2 hg38
  params.nPAR='chrX:2781479-155701383'   # non-par hg38
  params.nPAR_ploidy1='chrX 2781479 155701383 M 1' ## Ploidy reformating to haploid
  params.nPAR_ploidy2='chrX 1 10000 M 1'           ## Ploidy reformating to haploid

// 1) Create Reference for PCA projection 2) First batch of PCA projection (to create and reuse REF PCA) 
process first_pca_projection {
  label 'first_pca_projection'
  cache 'lenient'
  maxRetries 3
  scratch false
  errorStrategy 'retry'
  errorStrategy { task.exitStatus == 137 || task.exitStatus == 143 ? 'retry' : 'terminate' }
  
  executor "slurm"
  clusterOptions "--account=rrg-vmooser"

  cpus 12
  memory { 16.GB + 5.GB * task.attempt }
  time "4h"

  input:
  tuple path(bed), path(bim), path(fam), path(reference_geno), path(reference_site)


  output:
   path("*.lst"), emit : lists
   path("${bed.getBaseName()}_autosome.vcf.gz"), emit :vcf_autosome
   path("*.RefPC.coord"), emit: reference_pca
   path("*.ProPC.coord"), emit: coord
  script:
  """
  ## Convert bed to vcf.gz
  plink -bfile ${bed.getBaseName()} --a1-allele ${bim} 6 2  --autosome --output-chr MT --recode vcf-iid bgz --out ${bed.getBaseName()}_autosome 

  ## Create even lists of samples
  bcftools query -l ${bed.getBaseName()}_autosome.vcf.gz  | split -n r/${params.ind_trace_batch} -d --additional-suffix=.lst

  ## Create geno/site file for first set of samles
  ${params.laser_dir}/vcf2geno/vcf2geno --inVcf ${bed.getBaseName()}_autosome.vcf.gz --peopleIncludeFile x00.lst --out ${bed.getBaseName()}_x00.lst
  rm x00.lst

  ## PCA projection (1st unparralelised for the reference PCA to be reused)
  ${params.laser_dir}/trace -s ${bed.getBaseName()}_x00.lst.geno -g ${reference_geno} -o 00 -k ${params.K} -K ${params.K} -nt 12 ### produces RefPC.coord used later


    """
}

// Parallelize PCA projection by individuals
process parallel_pca_projection {
  label 'parallel_pca_projection'
  memory { 16.GB + 5.GB * task.attempt }
  errorStrategy 'retry'
  errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate' }
  maxRetries 5
  scratch false

  executor "slurm"
  clusterOptions "--account=rrg-vmooser"

  cpus 12
  time "5h"
  
  input: 
  tuple path(list), path(reference_geno), path(reference_site)
  each path(vcf_autosome)
  each path(ref_PCA)
  
  
  output:
  path("*.ProPC.coord"), emit : coord

  script:
  """
  ${params.laser_dir}/vcf2geno/vcf2geno  --inVcf ${vcf_autosome} --peopleIncludeFile ${list} --out ${vcf_autosome.getBaseName()}_${list.getBaseName()} 

  ${params.laser_dir}/trace -s ${vcf_autosome.getBaseName()}_${list.getBaseName()}.geno -c ${ref_PCA} -g ${reference_geno} -o ${vcf_autosome.getBaseName()}_${list.getBaseName()} -k ${params.K} -K ${params.K} -nt 12 ### produces RefPC.coord used later

  """
}

// Create unrelated individuals list
process unrelated {
  label 'unrelated'
  cache 'lenient'
    maxRetries 3
  scratch false
  errorStrategy 'retry'
  errorStrategy { task.exitStatus == 137 || task.exitStatus == 143 ? 'retry' : 'terminate' }

  executor "slurm"
  clusterOptions "--account=rrg-vmooser"

  cpus 12
  memory { 16.GB * task.attempt }
  time "2h"

  input :
  tuple path(bed), path(bim), path(fam)

  output :
  tuple path("*unrelated.vcf.gz"), path("*unrelated.vcf.gz.csi"), emit : unrelated
  path("*_tmp.king.cutoff.in.id"), emit: list_unrelated
  """
  plink2 --bfile ${bed.getBaseName()} --alt1-allele ${bim} 6 2 --output-chr chrMT --export vcf bgz id-paste=iid --out ${bed.getBaseName()} --threads 12

  plink2 --vcf ${bed.getBaseName()}.vcf.gz --king-cutoff ${params.King_cutOff} --output-chr chrMT --make-pgen --out ${bed.getBaseName()}_tmp --threads 12

  plink2 --pfile ${bed.getBaseName()}_tmp --alt1-allele ${bim} 6 2 --output-chr chrMT --keep ${bed.getBaseName()}_tmp.king.cutoff.in.id --export vcf bgz id-paste=iid --out ${bed.getBaseName()}_unrelated --threads 12

  bcftools index ${bed.getBaseName()}_unrelated.vcf.gz

  """
}

// LRT of Array 
process pca_test {
  label 'pca_test'
  cache 'lenient'
  maxRetries 3
  scratch false
  errorStrategy 'retry'
  errorStrategy { task.exitStatus == 137 || task.exitStatus == 143 ? 'retry' : 'terminate' }
  
  executor "slurm"
  clusterOptions "--account=rrg-vmooser"

  cpus 1
  memory { 4.GB * task.attempt }
  time "2h"

  input:
  tuple val(chrom), path(bed), path(bim), path(fam), path(vcf), path(index)
  each path(list_unrelated)
  each path(label)
  path(coord)
  each path(coord1)


  output:
  path("AF_ftest_${params.K}PC*"), emit : Tests

  publishDir "${params.OutDir}/S5", pattern: "AF_ftest_${params.K}PC", mode: "copy"
  publishDir "${params.OutDir}/S5", pattern: "cohort.pca", mode: "copy"
  script:
  
  """
  cat ${coord1} > ProPC.coord
  for file in ${coord} ; do sed 1d \$file >> ProPC.coord ; done
  ## Create files of unrelated CaG individuals

  ## Created Label files

  cut -f 1,5 ${bed.getBaseName()}.fam | grep 2\$ | cut -f 1 > female.filter
  cut -f 1,5 ${bed.getBaseName()}.fam | grep 1\$ | cut -f 1 > male.filter
  cut -f 1  ${list_unrelated} >  unrelated_king.txt

  grep -f unrelated_king.txt ${label} > label_unrelated
  grep -f male.filter label_unrelated > label_unrelated_male
  grep -f female.filter label_unrelated > label_unrelated_female

  sed "s/indivID/IID/g" ProPC.coord > cohort.PCA

  ## Regions to be tested on X chromosome (male female separated)


  if [ "${chrom}" -eq 23 ]; then
  compare_ancestry_adjusted_af.py -v ${vcf} -l label_unrelated -r \${params.PAR1} -p  cohort.PCA -k ${params.K} -o AF_ftest_${params.K}PC_PAR1
  compare_ancestry_adjusted_af.py -v ${vcf} -l label_unrelated -r \${params.PAR2} -p  cohort.PCA -k ${params.K} -o AF_ftest_${params.K}PC_PAR2
  compare_ancestry_adjusted_af.py -v ${vcf} -l label_unrelated_male -r \${params.nPAR}  -p  cohort.PCA -k ${params.K} -o AF_ftest_${params.K}PC_Xmale
  compare_ancestry_adjusted_af.py -v ${vcf} -l label_unrelated_female -r \${params.nPAR}  -p  cohort.PCA -k ${params.K} -o AF_ftest_${params.K}PC_Xfemale
  else
  compare_ancestry_adjusted_af.py -v ${vcf} -l label_unrelated -r chr${chrom} -p  cohort.PCA -k ${params.K} -o AF_ftest_${params.K}PC_chr${chrom}
  fi

  """
}

// Produce Harmonious vcf
process filtering {
  label 'filtering'
  cache 'lenient'

  maxRetries 3
  scratch false
  errorStrategy 'retry'
  errorStrategy { task.exitStatus == 137 || task.exitStatus == 143 ? 'retry' : 'terminate' }

  executor "slurm"
  clusterOptions "--account=rrg-vmooser"

  cpus 1
  memory { 32.GB * task.attempt }
  time "10h"
  
  input:
  tuple path(bed), path(bim), path(fam)
  path(AF_Test)


  output:
   tuple path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit : vcf
   path("Variants_failing_ancestry_LRTtest.txt")

  publishDir "${params.OutDir}/S5/VCF", pattern: "*_Harmonized.vcf.gz", mode: "copy"
  publishDir "${params.OutDir}/S5/Index", pattern: "*.tbi", mode: "copy"
  publishDir "${params.OutDir}/S5/Unharmonious_variant", pattern: "Variants_failing_ancestry_LRTtest.txt", mode: "copy"

  """
    ## concat results
  cat AF_ftest_${params.K}PC_chr1 > AF_ftest_${params.K}PC
  for i in chr{2..22} PAR1 PAR2 Xfemale Xmale
  do
  sed 1d AF_ftest_${params.K}PC_\${i} >> AF_ftest_${params.K}PC
  done

  echo ${params.nPAR_ploidy1} > ploidy.txt ## Male haploid for non PAR
  echo ${params.nPAR_ploidy2} >> ploidy.txt 
  cut -f 1,5 ${fam} | sed 's/2\$/F/g' | sed 's/1\$/M/g'  > samples_sex.txt 

  filtering.py -i AF_ftest_${params.K}PC -p ${params.threshold}
  plink -bfile ${bed.getBaseName()} --a1-allele ${bed.getBaseName()}.bim 6 2  --output-chr chrMT --exclude Variants_failing_ancestry_LRTtest.txt --recode vcf-iid bgz --out temp_Harmonized
  bcftools +fixploidy temp_Harmonized.vcf.gz  -- -p ploidy.txt -s samples_sex.txt > temp.vcf  ## Creates a non PAR ploidy haploid in males
  bcftools +fill-tags temp.vcf -Oz -o Harmonized.vcf.gz
  bcftools index -t Harmonized.vcf.gz
  rm temp*
  """
  }


workflow {
  //Input files
     unharmonized_data = Channel.fromPath(params.unharmonized_data).map(f -> [f, file("${f.getParent()}/${f.getBaseName()}.bim"), file("${f.getParent()}/${f.getBaseName()}.fam")])
     label = Channel.fromPath(params.list_array_origin)
     reference = Channel.fromPath(params.reference_pca_projection_pannel).map(f -> [f, file("${f.getParent()}/${f.getBaseName()}.site")])
     chromosome = channel.from( 1..23 )

  // PCA projection 1st
    S5_1 = first_pca_projection(unharmonized_data.combine(reference))
  // PCA projection (remainder)
    S5_2 = parallel_pca_projection(S5_1.lists.flatten().combine(reference),S5_1.vcf_autosome, S5_1.reference_pca)
  // Get Unrelated individuals
    unrelated = unrelated(unharmonized_data)
    ch1=unharmonized_data.concat(unrelated.unrelated).flatten().collect()
    chromosome.combine(ch1).flatten().collate( 6 ).set { result }
    result.view()
    S5_3 = pca_test(result, unrelated.list_unrelated,label,S5_2.coord.collect(), S5_1.coord)
  // Step 5.4
    S5_4 = filtering(unharmonized_data , S5_3.Tests.collect())
}
