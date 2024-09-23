#!/usr/bin/env nextflow
/*
* AUTHOR: Vincent Chapdelaine
* VERSION: 3.0
* YEAR: 2024
*/

nextflow.enable.dsl = 2


// Pre-requisites/assumptions:
// 1. Script was designed for genotyping arrays.
// 2. Overapping samples (by IDs) will be eliminated from the array with the least variants genotyped.
// 3. Script assumes the IIDs are meant to be unique
// 4. Script assumes the same reference build
// 5. Make sure that plink is installed/loaded on your system

// How to execute:
// nextflow run merge.nf --beds "/home/user/data/*.bed" --OutDir "/home/user/out"

params.beds = "/path/to/data/*.bed" // Absolute path to genotyping arrays (multiple) in plink bed file (with bim and fam files in the same folder)
params.OutDir = 'path/to/Dir/'      // Absolute path to the output directory

process merging {
  label 'merging'
  cache "lenient"
  
  executor "slurm"
  clusterOptions "--account=rrg-vmooser"

  cpus 5
  memory "16GB"
  time "2h"


  input:
  path(files)

  output:
  
  tuple path("unharmonized.bed"), path("unharmonized.bim"), path("unharmonized.fam"), emit: unharmonized
  path("samples_array_list.txt"), emit: Family_Arrays
  
  publishDir "${params.OutDir}/S4", pattern: "unharmonized*", mode: "copy"  
  publishDir "${params.OutDir}/S4", pattern: "samples_array_list.txt", mode: "copy"  
  
  script:
  """
  ## establish priority order in array deduplication
   list_files_order=(\$(wc -l *.bim | sort -n | awk '{print \$2}' | head -n -1))

  ## Focus on samples list
  for iter in \${!list_files_order[@]}; do
    list_files_order[iter]=\${list_files_order[iter]%.*}.fam
  done

  ## iterate over order list and drop shared samples until none remain
    for iter in \${!list_files_order[@]}
  do
    filename=\${list_files_order[iter]##*/}
    awk '{print \$1}' \$(ls \${list_files_order[@]}) | sort | uniq -c | awk '\$1 > 1 {print \$2,\$2}' > Samples_to_exclude.txt  ## Determine the shared samples in remainder
    if [ \$(wc -l < Samples_to_exclude.txt) -eq 0 ]; then # loop break if empty
        break
    else
      sed -i '1i#FID IID' Samples_to_exclude.txt                    ## Formating
      head Samples_to_exclude.txt
      plink --bfile \${list_files_order[iter]%.*} --remove Samples_to_exclude.txt --keep-allele-order --make-bed --output-chr chrMT --threads 5 --out \${filename%.*}_tmp 
      list_files_order[iter]=\${filename%.*}_tmp.fam  ## !! Replaces old file with new deduplicated file !! ##
    fi
  done

  ## Focus on variant files
  for iter in \${!list_files_order[@]}
  do
    list_files_order[iter]=\${list_files_order[iter]%.*}.bim
  done

  ## Establish list of shared variants
  nfile=\$(ls -lh \${list_files_order[@]} | wc -l) ## number of arrays
  ## Count variants occurance and if equals to number of arrays -> shared
  awk '{print \$2}' \$(ls \${list_files_order[@]}) | sort | uniq -c | awk -v n="\$nfile" '\$1 == n {print \$2}' | grep -v 'chrY' - > shared.txt
  
  ## Drop all non shared array
  for bim in \${list_files_order[@]} ; do
    out=\${bim##*/}
    plink --bfile \${bim%.*} --extract shared.txt  --keep-allele-order --make-bed --output-chr chrMT --out \${out%.*}_shared --threads 5
  done

  ## Update list to shared plink files
  for i in "\${!list_files_order[@]}";
  do
    list_files_order[\$i]="\${list_files_order[\$i]##*/}_shared"
  done

  ## keep data correspondance samples -> array (used for harmonization)
  for FAM in \${list_files_order[@]}; do
    cut -f 1 -d ' ' "\${FAM%.*}_shared".fam | sed "s/\$/\t\${FAM%.*}/g" >> samples_array_list.txt;
  done

  ## Iterate over files and merging each one
  files=(\${list_files_order[@]})
  base_file="\${files[0]%.*}_shared"
  unset files[0] # pop out the first file (reason : not merging with itself)
  output="merge_temp"
  for file in "\${files[@]}"
  do
  printf "\$base_file \n"
  printf "\${file%.*} \n"
    next_output="\${output}_next"
    plink --bfile \$base_file --keep-allele-order --output-chr chrMT --bmerge \${file%.*}_shared --threads 5 --out \$next_output
    base_file=\$next_output
    output=\$next_output
  done

  ## Create final product of this workflow
  mv \${output}.bed unharmonized.bed
  mv \${output}.bim unharmonized.bim
  mv \${output}.fam unharmonized.fam

  """
}

workflow {
  plink_files = Channel.fromPath(params.beds).map(f -> [f, file("${f.getParent()}/${f.getBaseName()}.bim"), file("${f.getParent()}/${f.getBaseName()}.fam")]).collect()
  plink_files.view()
  merging(plink_files)
}
