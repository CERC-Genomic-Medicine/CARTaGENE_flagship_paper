#!/usr/bin/env nextflow
/*
* AUTHOR: Mohadese Sayahian Dehkordi <mohadese.sayahiandehkordi@mail.mcgill.ca>, Daniel Taliun <daniel.taliun@mcgill.ca>
* VERSION: 2.0
* YEAR: 2023
*/

// How to run:
// nextflow run statistical_phasing.nf --vcf_path /path/to/chr*.vcf.gz --beagle /path/to/BEAGLE5/beagle --genetic_map_path /path/to/plink.chr*.map -with-report report.html

//params.vcf_path = "/path/to/chr*.vcf.gz" // Absolute path to the Input VCF/BCF files split by chromosome (chromosome X splitted to PAR1, PAR2, non_PAR)
//params.beagle = "/path/to/executable/beagle" // Absolute path to the Beagle executable
//params.genetic_map_path = "/path/to/plink.chr*.map" // Absolute path to the genetic map for Beagle split by chromosome (chromosome X splitted to PAR1, PAR2, non_PAR)
//params.iteration = 12 // number of iteration for the BEAGLE statistical phasing, default value is 12.
//params.seed = 3849 // seed for the Hidden Markov Model of BEAGLE, it can be any random number.

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
    # Beagle removes all the original headers and adds new header lines. 
    # This will be problematic when we have SVs in our VCF file because after phasing the end tag for VCF won't be recognizable.
    # First, we need to extract the header from the original VCF input.
    bcftools view -h ${vcf} > original_header.txt
    # Second, we need to remove the sample names and column names from the header, i.e. last line of the VCF header since we want to merge this header to the top of phased VCF headers.
    head -n -1 original_header.txt > original_header_without_sample_names.txt
    java -jar -Xmx32g ${params.beagle} window=25.0 overlap=2.5 nthreads=8 gt=${vcf} map=${genetic_map} seed=${params.seed} iterations=${params.iterations} out=${vcf.getBaseName()}.beagle_phased
    # Third, we need to extract the headers from the phased VCFs.
    bcftools view -h ${vcf.getBaseName()}.beagle_phased.vcf.gz > new_header.txt
    # Fourth, we need to combine the original header with the new header.
    cat original_header_without_sample_names.txt new_header.txt > combined_header.txt
    # Fifth, we need to change the header of phased VCFs using our new combined VCFs.
    bcftools reheader -h combined_header.txt -o ${vcf.getBaseName()}.beagle_phased.reheader.vcf.gz ${vcf.getBaseName()}.beagle_phased.vcf.gz
    bcftools index --tbi ${vcf.getBaseName()}.beagle_phased.reheader.vcf.gz
    """
}

process recalculate_AF_phased {
    errorStrategy 'retry'
    maxRetries 3
    cache "lenient"

    //executor 'slurm'
    //clusterOptions '--account=rrg-vmooser'

    cpus 1
    memory "8GB"
    time "4h"

    //scratch "$SLURM_TMPDIR"

    input:
    tuple path(vcf), path(vcf_index)
    
    output:
    tuple path("*.AF_calculated.vcf.gz"), path("*.AF_calculated.vcf.gz.tbi")
    
    publishDir "phased/recalculated_AF_vcfs/", pattern: "*.vcf.gz*", mode: "copy"   

    """
    # Beagle removes the allele frequency columns. Therefore if we want to have them in our final phased VCF files, we need to recalculate those tags.
    bcftools +fill-tags $vcf -Oz -o ${vcf.getSimpleName()}.AF_calculated.vcf.gz -- -t AN,AC,AF,NS
    bcftools tabix --tbi ${vcf.getSimpleName()}.AF_calculated.vcf.gz
    """
}
process remove_singletons {
    cache "lenient"
    cpus 1

    //executor 'slurm'
    //clusterOptions '--account=rrg-vmooser'

    memory "16GB"
    time "5h"

    //scratch "$SLURM_TMPDIR"
    input:
    tuple path(vcf), path(vcf_index)
    
    output:
    tuple path("*.vcf.gz"), path("*.vcf.gz.tbi")
    
    publishDir "phased/ref_withoutsingletons_vcfs/", pattern: "*.vcf.gz*", mode: "copy"

    """
    # After phasing, we need to remove singletons as part of imputation reference panel construction process.
    bcftools view $vcf -c 2 -Oz -o  ${vcf.getSimpleName()}.with_out_singletons.vcf.gz
    bcftools index --tbi ${vcf.getSimpleName()}.with_out_singletons.vcf.gz
    """
}



workflow {
    genetic_map_ch = Channel.fromPath(params.genetic_map_path).map { file -> [ file.name.toString().tokenize('.').get(1), file] }
    vcf_ch = Channel.fromPath(params.vcf_path).map{ vcf -> [vcf.getSimpleName(), vcf, vcf + ".tbi" ] }
    stat_phasing_ch = vcf_ch.join(genetic_map_ch)
	phased_ch = beagle_statistical_phasing(stat_phasing_ch)
    recalculated_AF_ch = recalculate_AF_phased(phased_ch)
    remove_singletons(recalculated_AF_ch)
}