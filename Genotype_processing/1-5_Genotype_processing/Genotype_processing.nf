#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process S0_preprocess {
  label 'S0_preprocess'
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

process S1_LiftOver {
  label 'S1_LiftOver'
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
process S2_Alignement_S3_HHmissing {
  label 'S2_Alignement'
  cache 'lenient'

  input:
  tuple path(bed), path(bim), path(fam)
  each path(fasta)

  output:
  tuple path("${bim.getBaseName()}_hh.bed"), path("${bim.getBaseName()}_hh.bim"), path("${bim.getBaseName()}_hh.fam"), emit: Aligned
  
  publishDir "${params.OutDir}/S2", pattern: "*.bim", mode: "copy"
  
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

// Deduplicate samples & Merge
process S4_Merge {
  label 'S4_Merge'
  cache 'lenient'
  scratch true

  input:
  path(files)

  output:
  tuple path("${params.Merged_name}_unharmonized.bed"), path("${params.Merged_name}_unharmonized.bim"), path("${params.Merged_name}_unharmonized.fam"), emit: unharmonized
  path("*_shared.fam"), emit: Family_Arrays
  
  publishDir "${params.OutDir}/S4", pattern: "${params.Merged_name}_unharmonized.*}", mode: "copy"  
  
  script:
  """
  list_files_order=(\$(wc -l *.fam | sort -n | awk '{print \$2}' | head -n -1))
  for iter in \${!list_files_order[@]}
  do
    filename=\${list_files_order[iter]##*/}
    awk '{print \$1}' \$(ls \${list_files_order[@]}) | sort | uniq -c | awk '\$1 > 1 {print \$2,\$2}' > Samples_to_exclude.txt
    sed -i '1i#FID IID' Samples_to_exclude.txt
    head  Samples_to_exclude.txt
    if [ \$(wc -l < Samples_to_exclude.txt) -eq 0 ]; then # End if list is empty
        break
    else
      plink --bfile \${list_files_order[iter]%.*} --remove Samples_to_exclude.txt --keep-allele-order --make-bed --output-chr chrMT --threads ${params.threads_S4} --out \${filename%.*}
      list_files_order[iter]="\${list_files_order[iter]##*/}"
    fi
  done
  for iter in \${!list_files_order[@]}
  do
    list_files_order[iter]=\${list_files_order[iter]%.*}.bim
  done
  nfile=\$(ls -lh \${list_files_order[@]} | wc -l)
  awk '{print \$2}' \$(ls \${list_files_order[@]}) | sort | uniq -c | awk -v n="\$nfile" '\$1 == n {print \$2}' | grep -v 'chrY' - > shared.txt
  for bim in \${list_files_order[@]} ; do
    out=\${bim##*/}
    plink --bfile \${bim%.*} --extract shared.txt  --keep-allele-order --make-bed --output-chr chrMT --out \${out%.*}_shared --threads ${params.threads_S4}
  done
  for i in "\${!list_files_order[@]}";
  do
    list_files_order[\$i]="\${list_files_order[\$i]##*/}_shared"
  done
  files=(\${list_files_order[@]})
  base_file="\${files[0]%.*}_shared"
  unset files[0] # pop out the first file
  output="merge_temp"
  for file in "\${files[@]}"
  do
    next_output="\${output}_next"
    plink --bfile \$base_file --keep-allele-order --output-chr chrMT --bmerge \${file%.*}_shared --threads ${params.threads_S4} --out \$next_output
    base_file=\$next_output
    output=\$next_output
  done
  mv \${output}.bed ${params.Merged_name}_unharmonized.bed
  mv \${output}.bim ${params.Merged_name}_unharmonized.bim
  mv \${output}.fam ${params.Merged_name}_unharmonized.fam
  """
}
// Unfiltered Reference VCF to filtered Reference VCF
process REF_projection_prep {
  label 'REF_projection_prep'
  cache 'lenient'

  input:
  tuple path(vcf), path(index), path(bed), path(bim), path(fam)
  each path(HGDP1K_samples)

  output:
  path("temp_*.vcf.gz"), emit : REF_VCFs
  path("temp_*.vcf.gz.csi"), emit: REF_index

  script:
  """
  cut -f 2 ${bim} | sed 's/:/\t/g' - | awk '{print \$1, \$2-1, \$2}' - | sed 's/ /\t/g' - > shared.bed
  bcftools view -f 'PASS' -q${params.min_allelefreq} -Q${params.max_allelefreq} -S ${HGDP1K_samples} -i "F_MISSING < ${params.MISSING_CUTOFF}" -R shared.bed ${vcf}  --threads ${params.threads_S5} |  bcftools annotate -x INFO,^GT - -Oz -o temp_${vcf.getBaseName()}.gz --threads ${params.threads_S5}
  bcftools index temp_${vcf.getBaseName()}.gz
  """
}

// 1) Create Reference for PCA projection 2) First batch of PCA projection (to create and reuse REF PCA) 
process S5_PCAprojection_prep {
  label 'S5_PCAprojection_prep'
  cache 'lenient'
  memory { 32.GB * task.attempt }
  errorStrategy 'retry'
  errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate' }
  maxRetries 5

  input:
  tuple path(bed), path(bim), path(fam)
  path(REF)
  path(index)


  output:
   path("*.lst"), emit : lists
   path("${bed.getBaseName()}_autosome.vcf.gz"), emit :vcf_autosome
   path("REF.geno"), emit: REF
   path("REF.site"), emit: REF_site
   path("*.RefPC.coord"), emit: REF_PCA
   path("*.ProPC.coord"), emit: coord
  script:
  """
  echo ${REF}
  bcftools concat ${REF} -n -Oz -o REF.vcf.gz  # Concatenate files
  plink2 --vcf REF.vcf.gz --indep-pairwise ${params.window} ${params.step} ${params.Rsq} --out REF.LD_prune --set-all-var-ids '@:#:\$r:\$a' --new-id-max-allele-len 145 --threads ${params.threads_S5} # ID var to be pruned
  plink2 --vcf REF.vcf.gz --set-all-var-ids 'chr@:#:\$r:\$a' --new-id-max-allele-len 145 --sort-vars --threads ${params.threads_S5}  --exclude REF.LD_prune.prune.out --make-pgen --out tmp # Remove pruned variants
  plink2 --pfile tmp --export vcf bgz id-paste=iid --out REF.LD_prune --threads ${params.threads_S5}
  bcftools index REF.LD_prune.vcf.gz
  ${params.vcf2geno_exec} --inVcf REF.LD_prune.vcf.gz --out REF # convert to geno/site format

  ## Samples Prep

  plink -bfile ${bed.getBaseName()} --a1-allele ${bim} 6 2  --autosome --output-chr MT --recode vcf-iid bgz --out ${bed.getBaseName()}_autosome 
  bcftools query -l ${bed.getBaseName()}_autosome.vcf.gz  | split -n r/${params.ind_trace_batch} -d --additional-suffix=.lst
  
  ${params.vcf2geno_exec} --inVcf ${bed.getBaseName()}_autosome.vcf.gz --peopleIncludeFile x00.lst --out ${bed.getBaseName()}_x00.lst
  rm x00.lst

  ## PCA projection

  ${params.trace_exec} -s ${bed.getBaseName()}_x00.lst.geno -g REF.geno -o 00 -k ${params.K} -K ${params.K} -nt ${params.threads_S5} ### produces RefPC.coord used later


    """
}

// Parallelize PCA projection by individuals
process S5_PCAprojection {
  label 'S5_PCAprojection'
  memory { 16.GB + 5.GB * task.attempt }
  errorStrategy 'retry'
  errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate' }
  maxRetries 5

  input: 
  path(list)
  each path(vcf_autosome)
  each path(ref_PCA)
  each path(REF)
  each path(REF_site)
  
  output:
  path("*.ProPC.coord"), emit : coord

  script:
  """
  ${params.vcf2geno_exec} --inVcf ${vcf_autosome} --peopleIncludeFile ${list} --out ${vcf_autosome.getBaseName()}_${list.getBaseName()} 

  ${params.trace_exec} -s ${vcf_autosome.getBaseName()}_${list.getBaseName()}.geno -c ${ref_PCA} -g ${REF} -o ${vcf_autosome.getBaseName()}_${list.getBaseName()} -k ${params.K} -K ${params.K} -nt ${params.threads_S5} ### produces RefPC.coord used later

  """
}

// Create unrelated individuals list
process unrelated {
  label 'unrelated'
  cache 'lenient'

  input :
  tuple path(bed), path(bim), path(fam)

  output :
  tuple path("*unrelated.vcf.gz"), path("*unrelated.vcf.gz.csi"), emit : unrelated
  path("*_tmp.king.cutoff.in.id"), emit: kingcut
  """
  plink2 --bfile ${bed.getBaseName()} --alt1-allele ${bim} 6 2 --output-chr chrMT --export vcf bgz id-paste=iid --out ${bed.getBaseName()} --threads 12

  plink2 --vcf ${bed.getBaseName()}.vcf.gz --king-cutoff ${params.King_cutOff} --output-chr chrMT --make-pgen --out ${bed.getBaseName()}_tmp --threads 12

  plink2 --pfile ${bed.getBaseName()}_tmp --alt1-allele ${bim} 6 2 --output-chr chrMT --keep ${bed.getBaseName()}_tmp.king.cutoff.in.id --export vcf bgz id-paste=iid --out ${bed.getBaseName()}_unrelated --threads 12

  bcftools index ${bed.getBaseName()}_unrelated.vcf.gz

  """
}

// LRT of Array 
process S5_PCAtest {
  label 'S5_PCAtest'
  cache 'lenient'

  input:
  tuple val(chrom), path(bed), path(bim), path(fam), path(vcf), path(index)
  path(kingcut)
  path(CaG_array)
  path(coord)
  path(coord1)


  output:
  path("AF_ftest_${params.K}PC*"), emit : Tests

  publishDir "${params.OutDir}/S2", pattern: "AF_ftest_${params.K}PC", mode: "copy"
  script:
  """
  files=(\$(ls ${coord}))
  base_file="\${files[0]%.*}"
  unset files[0] # pop out the first file

  cat ${coord1} > ${params.Merged_name}.ProPC.coord
  for file in ${coord} ; do sed 1d \$file >> ${params.Merged_name}.ProPC.coord ; done
  ## Create files of unrelated CaG individuals

  ## Created Label files

  for FAM in ${CaG_array}; do 
    cut -f 1 -d ' ' \${FAM} | sed "s/\$/\t\${FAM%.*}/g" >> label;
  done

  cut -f 1,5 ${bed.getBaseName()}.fam | grep 2\$ | cut -f 1 > female.filter
  cut -f 1,5 ${bed.getBaseName()}.fam | grep 1\$ | cut -f 1 > male.filter
  cut -f 1  ${kingcut} >  unrelated_king.txt

  grep -f unrelated_king.txt label > label_unrelated
  grep -f male.filter label_unrelated > label_unrelated_male
  grep -f female.filter label_unrelated > label_unrelated_female

  sed "s/indivID/IID/g" ${params.Merged_name}.ProPC.coord > ${params.Merged_name}.PCA

  PAR1='chrX:10001-2781479'       # PAR1 hg38
  PAR2='chrX:155701383-156030895' # PAR2 hg38
  nPAR='chrX:2781479-155701383'   # non-par hg38

  if [ "${chrom}" -eq 23 ]; then
  compare_ancestry_adjusted_af.py -v ${vcf} -l label_unrelated -r \${PAR1} -p  ${params.Merged_name}.PCA -k ${params.K} -o AF_ftest_${params.K}PC_PAR1
  compare_ancestry_adjusted_af.py -v ${vcf} -l label_unrelated -r \${PAR2} -p  ${params.Merged_name}.PCA -k ${params.K} -o AF_ftest_${params.K}PC_PAR2
  compare_ancestry_adjusted_af.py -v ${vcf} -l label_unrelated_male -r \${nPAR}  -p  ${params.Merged_name}.PCA -k ${params.K} -o AF_ftest_${params.K}PC_Xmale
  compare_ancestry_adjusted_af.py -v ${vcf} -l label_unrelated_female -r \${nPAR}  -p  ${params.Merged_name}.PCA -k ${params.K} -o AF_ftest_${params.K}PC_Xfemale
  else
  compare_ancestry_adjusted_af.py -v ${vcf} -l label_unrelated -r chr${chrom} -p  ${params.Merged_name}.PCA -k ${params.K} -o AF_ftest_${params.K}PC_chr${chrom}
  fi

  """
}

// Produce Harmonious vcf
process S5_Filtering {
  label 'S5_Filtering'
  cache 'lenient'
  
  input:
  tuple path(bed), path(bim), path(fam)
  path(AF_Test)


  output:
   tuple path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit : vcf
   path("Variants_failing_ancestry_LRTtest.txt")

  publishDir "${params.OutDir}/S5/VCF", pattern: "*.vcf.gz", mode: "copy"
  publishDir "${params.OutDir}/S5/Index", pattern: "*.tbi", mode: "copy"
  publishDir "${params.OutDir}/S5/Unharmonious_variant", pattern: "Variants_failing_ancestry_LRTtest.txt", mode: "copy"

  """
    ## concat results
  cat AF_ftest_${params.K}PC_chr1 > AF_ftest_${params.K}PC
  for i in chr{2..22} PAR1 PAR2 Xfemale Xmale
  do
  sed 1d AF_ftest_${params.K}PC_\${i} >> AF_ftest_${params.K}PC
  done

  S5_scrpit.py -i AF_ftest_${params.K}PC -p ${params.threshold}
  plink -bfile ${bed.getBaseName()} --a1-allele ${bed.getBaseName()}.bim 6 2  --output-chr chrMT --exclude Variants_failing_ancestry_LRTtest.txt --recode vcf-iid bgz --out temp_${params.Merged_name}_Harmonized
  bcftools +fill-tags temp_${params.Merged_name}_Harmonized.vcf.gz -Oz -o ${params.Merged_name}_Harmonized.vcf.gz
  bcftools index -t ${params.Merged_name}_Harmonized.vcf.gz
  rm temp*
  """
  }


workflow {
  //Input files
     G_Arrays = Channel.fromPath(params.beds_initial).map(f -> [f, file("${f.getParent()}/${f.getBaseName()}.bim"), file("${f.getParent()}/${f.getBaseName()}.fam")])
     Consent = Channel.fromPath(params.Consent)
     REF_Projection = Channel.fromPath(params.REF_projection).map(f -> [f, file("${f}.csi")])
     chromosome = channel.from( 1..23 )
  // Ressource file
     Chain = Channel.fromPath(params.chain)
     Genome = Channel.fromPath(params.genome)
     REF_ID=Channel.fromPath(params.REF_ID)
     // LongLD=Channel.fromPath(params.LongLD)
  // Steps
    //Step 0
        S0 = S0_preprocess(G_Arrays, Consent)
    //Step 1
        S1 = S1_LiftOver(S0.Preprocessed, Chain)
    // Step 2-3
        S2 = S2_Alignement_S3_HHmissing(S1.Lifted, Genome)
    // Step 4
        S4 = S4_Merge(S2.Aligned.collect())
    // Step 5.1
        S5_0 = REF_projection_prep(REF_Projection.combine(S4.unharmonized), REF_ID)
        S5_1 = S5_PCAprojection_prep(S4.unharmonized, S5_0.REF_VCFs.collect(sort:true), S5_0.REF_index.collect())
    // Step 5.2 
    	
	S5_2 = S5_PCAprojection(S5_1.lists.flatten(),S5_1.vcf_autosome, S5_1.REF_PCA, S5_1.REF, S5_1.REF_site)
    // Step 5.3
    	unrelated = unrelated(S4.unharmonized)
        ch1=S4.unharmonized.concat(unrelated.unrelated).flatten().collect()
	chromosome.combine(ch1).flatten().collate( 6 ).set { result }
	S5_3 = S5_PCAtest(result, unrelated.kingcut,S4.Family_Arrays.collect(),S5_2.coord.collect(), S5_1.coord)
    // Step 5.4
        S5_4 = S5_Filtering(S4.unharmonized , S5_3.Tests.collect())

}
