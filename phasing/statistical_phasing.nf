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
    tuple path("*beagle_phased.reheader.vcf.gz"), path("*beagle_phased.reheader.vcf.gz.tbi")
    
    publishDir "phased/", pattern: "*beagle_phased.reheader.vcf.gz*", mode: "copy"
    
    """
    # we need to modify the header after phasing since the beagle tool removes the original header
    bcftools view -h ${vcf} > original_header.txt
    head -n -1 original_header.txt > original_header_without_sample_names.txt
    java -jar -Xmx32g ${params.beagle} window=25.0 overlap=2.5 nthreads=8 gt=${vcf} map=${genetic_map} out=${vcf.getBaseName()}.beagle_phased
    bcftools view -h ${vcf.getBaseName()}.beagle_phased.vcf.gz > new_header.txt
    cat original_header_without_sample_names.txt new_header.txt > combined_header.txt
    bcftools reheader -h combined_header.txt -o ${vcf.getBaseName()}.beagle_phased.reheader.vcf.gz ${vcf.getBaseName()}.beagle_phased.vcf.gz
    bcftools index --tbi ${vcf.getBaseName()}.beagle_phased.reheader.vcf.gz
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
    tuple path(vcf), path(vcf_index)
    
    output:
    tuple path("*.AF_calculated.vcf.gz"), path("*.AF_calculated.vcf.gz.tbi")
    
    publishDir "phased/recalculated_AF_vcfs/", pattern: "*.vcf.gz*", mode: "copy"   

    """
    bcftools +fill-tags $vcf -Oz -o ${vcf.getSimpleName()}.AF_calculated.vcf.gz -- -t AN,AC,AF,NS
    bcftools tabix --tbi ${vcf.getSimpleName()}.AF_calculated.vcf.gz
    """
}
process remove_singletons {
    cache "lenient"
    cpus 1
    memory "16GB"
    time "5h"
    input:
    tuple path(vcf), path(vcf_index)
    
    output:
    tuple path("*.vcf.gz"), path("*.vcf.gz.tbi")
    
    publishDir "phased/ref_withoutsingletons_vcfs/", pattern: "*.vcf.gz*", mode: "copy"

    """
    bcftools view $vcf -c 2 -Oz -o  ${vcf.getSimpleName()}.with_out_singletons.vcf.gz
    bcftools index --tbi ${vcf.getSimpleName()}.with_out_singletons.vcf.gz
    """
}



workflow {
    genetic_map_ch = Channel.fromPath(params.genetic_map_path).map { file -> [ file.name.toString().tokenize('_').get(1), file] }
    vcf_ch = Channel.fromPath(params.vcf_path).map{ vcf -> [ vcf, vcf + ".tbi" ] }
    vcf_with_chr_name = get_chr_name(vcf_ch)
    vcf_with_chr_name.view{"Result: ${it}"}
    genetic_map_ch.view{"Result: ${it}"}
    //stat_phasing_ch = vcf_with_chr_name.join(genetic_map_ch)
	//phased_ch = beagle_statistical_phasing(stat_phasing_ch)
    //recalculated_AF_ch = recalculate_AF_phased(phased_ch)
    //remove_singletons(recalculated_AF_ch)
}