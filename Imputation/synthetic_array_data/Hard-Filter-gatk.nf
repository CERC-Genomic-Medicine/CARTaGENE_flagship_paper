process HardFilter {
   errorStrategy 'finish'
   maxRetries 1
   memory "8 GB"
   cache "lenient"
   cpus 1
   time "3h"

   container "${params.gatkContainer}"

   publishDir "VCF_gatk_hard_filter/", pattern: "${name}.gatk_hard_filter.vcf.gz*", mode: "copy"

   input:
   tuple val(name), path(vcf), path(vcf_index)

   output:
   //tuple path("${name}.gatk_hard_filter.vcf.gz"), path("${name}.gatk_hard_filter.vcf.gz.tbi")
   path("${name}.gatk_hard_filter.vcf.gz"), emit : filtered_vcf
   path("${name}.gatk_hard_filter.vcf.gz.tbi"), emit : filtered_tbi

   """
   gatk --java-options -Xmx4G SelectVariants -V ${vcf} -select-type SNP -O snps.vcf.gz

   gatk --java-options -Xmx4G VariantFiltration -V snps.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ${name}.gatk_hard_filter.vcf.gz

   """ 
}


workflow {
   vcfs = Channel.fromPath(params.inputFiles, checkIfExists:true) \
    | map { file -> [ file.getName().toString().tokenize('.').get(0), file, file + ".tbi"] } 
   
   HardFilter(vcfs)
}
