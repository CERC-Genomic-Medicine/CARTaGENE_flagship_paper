#!/usr/bin/env nextflow

/*
* AUTHOR: Daniel Taliun, PhD <daniel.taliun@mcgill.ca>
* VERSION: 1.0
* YEAR: 2024
*/


process run_stepwise_conditional_analysis {
   cache "lenient"
   //scratch false

   cpus 1
   memory "8GB"
   time "8h"

   container "${params.regenie_container}"

   input:
   tuple val(chrom), val(start), val(end), val(phenotype_name), path(gwas_results), path(loco_pred), path(gwas_genotypes_file), path(sample_file), path(variants_file)
   each path(phenotypes_file)
   each path(covariates_file)
  
   output:       
   tuple val(phenotype_name), val(chrom), val(start), val(end), path("conditioning_variants.txt"), emit: out
   tuple path("*_step*"), path("index_variants.txt"), emit: aux

   publishDir "${params.output_dir}/${phenotype_name}/${chrom}_${start}_${end}/stepwise_conditional/", pattern: "*_step*", mode: "copy"
   publishDir "${params.output_dir}/${phenotype_name}/${chrom}_${start}_${end}/stepwise_conditional/", pattern: "index_variants.txt", mode: "copy"

   """
   # Initialize the common Regenie imput options
   if [ ${gwas_genotypes_file.getExtension()} = "pgen" ]; then
      input_genotypes="--pgen ${gwas_genotypes_file.getBaseName()}"
   else
      input_genotypes="--bgen ${gwas_genotypes_file} --sample ${sample_file}"
   fi

   if [ -z "${params.categorical_covariates}" ]; then
      categorical_covariates=""
   else
      categorical_covariates="--catCovarList ${params.categorical_covariates}"
   fi

   # Create file to map phenotype name to the corresponding genomic predictions:
   echo "${phenotype_name} ${loco_pred}" > loco_pred_list.txt

   # Get the top significant hit using the LOG10P column from the Regenie output
   signif_threshold=7.301029995663981
   gzip -dc ${gwas_results} | tail -n+2 | awk -v t=\${signif_threshold} '{if (\$13 > t) print \$3, \$13}' | sort -k2,2 -gr | head -n1 | cut -f1 -d" " > conditioning_variant.txt

   i=1
   gzip -dc ${gwas_results} | head -n1 | sed "s/\$/\tPHENOTYPE\tANALYSIS/" > index_variants.txt
   while [ -s conditioning_variant.txt ]; do
      cat conditioning_variant.txt >> conditioning_variants.txt
      
      if [ \${i} -eq 1 ]; then
         gzip -dc ${gwas_results} | grep -wFf conditioning_variant.txt | sed "s/\$/\t${phenotype_name}\tPRIMARY/" >> index_variants.txt
      else
         gzip -dc ${gwas_results} | grep -wFf conditioning_variant.txt | sed "s/\$/\t${phenotype_name}\tSECONDARY/" >> index_variants.txt
      fi

      regenie \
    		--step 2 \
		--gz \
		--loocv \
		--bsize ${params.block_size} \
		--phenoFile ${phenotypes_file} --strict \
		--covarFile ${covariates_file} \${categorical_covariates} \
		--condition-list conditioning_variants.txt \
		\${input_genotypes} ${params.apply_rint} \
		--range ${chrom.split('_')[0]}:${start}-${end} \
		--out "${chrom}_${start}_${end}_step\${i}" \
		--pred loco_pred_list.txt \
		--threads 1 \
		--minMAC ${params.min_mac} \
		--minINFO ${params.min_info} \
		--lowmem

      signif_threshold=5
      gzip -dc ${chrom}_${start}_${end}_step\${i}_${phenotype_name}.regenie.gz | tail -n+2 | awk -v t=\${signif_threshold} '{if (\$13 > t) print \$3, \$13}' | sort -k2,2 -gr | head -n1 | cut -f1 -d" " > conditioning_variant.txt
      i=\$((i+1))
   done
   """
}


process run_gwas_by_index_variant {
   cache "lenient"
   //scratch false

   cpus 1
   memory "8GB"
   time "8h"

   container "${params.regenie_container}"

   input:
   tuple val(chrom), val(start), val(end), val(phenotype_name), path(index_variants), path(loco_pred), path(gwas_genotypes_file), path(sample_file), path(variants_file)
   each path(phenotypes_file)
   each path(covariates_file)

   output:
   tuple val(chrom), val(start), val(end), val(phenotype_name), path("*_index_*.regenie.gz"), emit: out
   path "*_index_*", emit: intermediate

   publishDir "${params.output_dir}/${phenotype_name}/${chrom}_${start}_${end}/gwas_by_index_variant/", pattern: "*_index_*.*", mode: "copy"
   
   """
   # Initialize the common Regenie imput options
   if [ ${gwas_genotypes_file.getExtension()} = "pgen" ]; then
      input_genotypes="--pgen ${gwas_genotypes_file.getBaseName()}"
   else
      input_genotypes="--bgen ${gwas_genotypes_file} --sample ${sample_file}"
   fi

   if [ -z "${params.categorical_covariates}" ]; then
      categorical_covariates=""
   else
      categorical_covariates="--catCovarList ${params.categorical_covariates}"
   fi

   # Create file to map phenotype name to the corresponding genomic predictions:
   echo "${phenotype_name} ${loco_pred}" > loco_pred_list.txt

   while read -r variant; do
      grep -vF \${variant} conditioning_variants.txt > conditioning_variants_index_\${variant//:/_}.txt
      
      if [ -s conditioning_variants_index_\${variant//:/_}.txt ]; then
         conditioning="--condition-list conditioning_variants_index_\${variant//:/_}.txt"
      else
         conditioning=""
      fi

      regenie \
                --step 2 \
                --gz \
                --loocv \
                --bsize ${params.block_size} \
                --phenoFile ${phenotypes_file} --strict \
                --covarFile ${covariates_file} \${categorical_covariates} \${conditioning} \
		\${input_genotypes} ${params.apply_rint} \
                --range ${chrom.split('_')[0]}:${start}-${end} \
                --out "${chrom}_${start}_${end}_index_\${variant//:/_}" \
                --pred loco_pred_list.txt \
                --threads 1 \
		--minMAC ${params.min_mac} \
		--minINFO ${params.min_info} \
                --lowmem
   done < ${index_variants} || true # We do here this trick because `read` returns exit code 1 if it encounters EOF.
   """
}


process compute_credible_sets {
   cache "lenient"
   //scratch false

   cpus 1
   memory "4GB"
   time "1h"

   input:
   tuple val(chrom), val(start), val(end), val(phenotype_name), path(gwas_results)

   output:
   tuple val(chrom), val(start), val(end), val(phenotype_name), path("*.set99.txt"), emit: out

   publishDir "${params.output_dir}/${phenotype_name}/${chrom}_${start}_${end}/credible_set_by_index_variant/", pattern: "*.set99.txt", mode: "copy"

   """
   get_credible_set.py -g ${gwas_results} -s 0.99 -o ${gwas_results.getBaseName()}.set99.txt
   sed -i "s/PIP_CUMSUM\$/PIP_CUMSUM\tINDEX_ID\tPHENOTYPE/" ${gwas_results.getBaseName()}.set99.txt # append column name
   
   index_variant=${gwas_results.getSimpleName()}
   index_variant=`echo \${index_variant} | sed 's/.*_index_//' | sed 's/_${phenotype_name}//' | sed 's/_/:/g'`

   sed -i "s/[0-9]\$/\t\${index_variant}\t${phenotype_name}/" ${gwas_results.getBaseName()}.set99.txt # appned column values
   """
}


workflow {
   // Load phenotypes
   phenotypes = Channel.fromPath(params.phenotypes_file)
   //phenotypes.view()

   // Load covariates
   covariates = Channel.fromPath(params.covariates_file)
   //covariates.view()

   // Load loci: phenotype, chrom, start, end, gwas_results.
   loci = Channel.from(file(params.loci_file).readLines()).map { line -> def fields = line.split(); [ fields[0], fields[1], fields[2], fields[3], file(fields[4]) ] }
   //loci.view()

   // Load genomic predictions
   // 1) Read all *_pred.list and extract phenotype names and names of *.loco.gz files.
   loco2pheno = Channel.fromPath(params.genomic_predictions_file).map( it -> it.readLines()).flatten().map { line -> def fields = line.split(); [fields[1], fields[0]]}
   //loco2pheno.view()

   // 2) Get full paths of all *.loco.gz files
   loco2path = Channel.fromPath(params.genomic_predictions_file).map(it -> files(it.toString().replaceAll(/_pred.list$/, "_*.loco.gz"), checkIfExists: true)).flatten().map { it -> [it.getName(), it] }
   //loco2path.view()

   // 3) Merge
   genomic_predictions = loco2pheno.combine(loco2path, by: [0]).map(it -> [it[1], it[2]])
   //genomic_predictions.view()

   // Load genotype files. Here we assume that the filenames start with the chromosome name e.g. "chr1.study.pgen" or "1.study.pgen".
   gwas_genotypes = Channel.fromPath(params.gwas_genotypes_file).map(f -> f.getExtension() == "pgen" ? [(f.getBaseName() =~ /^(?:chr)(X_par1|X_par2|X_nonpar|X|[1-9][0-9]?)/)[0][1], f, file("${f.getParent()}/${f.getBaseName()}.psam"), file("${f.getParent()}/${f.getBaseName()}.pvar")] : [f, file("${f.getParent()}/${f.getBaseName()}.sample"), f + ".bgi"])
   //gwas_genotypes.view()

   index_variants = run_stepwise_conditional_analysis(loci.combine(genomic_predictions, by: [0]).map {it -> [it[1], it[2], it[3], it[0], it[4], it[5]]}.combine(gwas_genotypes, by: [0]), phenotypes, covariates).out
   //index_variants.view()

   gwas_by_index_variant = run_gwas_by_index_variant(index_variants.combine(genomic_predictions, by: [0]).map {it -> [it[1], it[2], it[3], it[0], it[4], it[5]]}.combine(gwas_genotypes, by: [0]), phenotypes, covariates).out
   //gwas_by_index_variant.view()

   compute_credible_sets(gwas_by_index_variant.transpose())
}
