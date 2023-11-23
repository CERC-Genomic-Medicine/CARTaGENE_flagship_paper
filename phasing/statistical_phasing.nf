#!/usr/bin/env nextflow
/*
* AUTHOR: Mohadese Sayahian Dehkordi <mohadese.sayahiandehkordi@mail.mcgill.ca>, Daniel Taliun <daniel.taliun@mcgill.ca>
* VERSION: 2.0
* YEAR: 2023
*/

// How to run:
// nextflow run statistical_phasing.nf --vcf_path /path/to/chr2.vcf.gz --eagle /path/to//Eagle_v2.4.1/eagle --genetic_map /path/to/genetic_map_hg38_withX.txt.gz -with-report report.html

//params.vcf_path = "/path/to/data/*.bcf" // Absolute path to the Input VCF/BCF files split by chromosome
//params.eagle = "/path/to/executable/eagle" // Absolute path to the Eagle executable
//params.genetic_map = "/path/to/genetic_map/genetic_map_hg38_withX.txt.gz" // Absolute path to the genetic map from Eagle software

process get_chr_name {
    //errorStrategy 'retry'
    //maxRetries 3
    cache "lenient"
    cpus 1
    memory "4GB"
    time "00:15:00"
    input:
    tuple path(snv_vcf), path(snv_index)

    output:
    tuple stdout, path(snv_vcf), path(snv_index)

    """
    chrom=`bcftools index -s ${snv_vcf} | cut -f1`
	printf "\${chrom}"
    """
}

process beagle_statistical_phasing {
    //errorStrategy 'retry'
    //maxRetries 3
    cache "lenient"
    
    //executor 'slurm'
    //clusterOptions '--account=rrg-vmooser'

    cpus 8
    memory "32GB"
    time "10h"
    //scratch "$SLURM_TMPDIR"

    input:
    tuple val(chromosome), path(vcf), path(vcf_index), path(genetic_map)
    
    output:
    path("*beagle_phased*")
    
    publishDir "phased/", pattern: "*beagle_phased*", mode: "copy"
    
    """
    java -jar -Xmx32g ${params.beagle} window=25.0 overlap=2.5 nthreads=8 gt=${vcf} map=${genetic_map} out=${vcf.getBaseName()}.beagle_phased
    """
}

process recalculate_AF_phased {
    errorStrategy 'retry'
    maxRetries 3
    cache "lenient"
    cpus 1
    memory "8GB"
    time "4h"
    input:
    tuple path(snv_vcf), path(snv_index)
    
    output:
    tuple path("*.AF_calculated.vcf.gz"), path("*.AF_calculated.vcf.gz.tbi")
    
    publishDir "phased/recalculated_AF_vcfs/", pattern: "*.vcf.gz", mode: "copy"   
    publishDir "phased/recalculated_AF_vcfs/", pattern: "*.vcf.gz.tbi", mode: "copy"    

    """
    bcftools +fill-tags $snv_vcf -Oz -o ${snv_vcf.getSimpleName()}.AF_calculated.vcf.gz -- -t AN,AC,AF,NS
    bcftools tabix --tbi ${snv_vcf.getSimpleName()}.AF_calculated.vcf.gz
    """
}
process remove_singletons {
    cache "lenient"
    cpus 1
    memory "16GB"
    time "5h"
    input:
    tuple path(chr), path(log)
    
    output:
    tuple path("*.vcf.gz"), path("*.vcf.gz.tbi")
    
    publishDir "phased/ref_withoutsingletons_vcfs/", pattern: "*.vcf.gz", mode: "copy"
    publishDir "phased/ref_withoutsingletons_vcfs/", pattern: "*.vcf.gz.tbi", mode: "copy"

    """
    bcftools view $chr -c 2 -Oz -o  ${chr.getSimpleName()}.with_out_singletons.vcf.gz
    bcftools index --tbi ${chr.getSimpleName()}.with_out_singletons.vcf.gz
    """
}



workflow {
    genetic_map_ch = Channel.fromPath(params.genetic_map_path).map { file -> [ file.name.toString().tokenize('_').get(1), file] }
    vcf_ch = Channel.fromPath(params.vcf_path).map{ vcf -> [ vcf, vcf + ".tbi" ] }
    vcf_with_chr_name = get_chr_name(vcf_ch)
    stat_phasing_ch = vcf_with_chr_name.join(genetic_map_ch)
	beagle_statistical_phasing(stat_phasing_ch)
}