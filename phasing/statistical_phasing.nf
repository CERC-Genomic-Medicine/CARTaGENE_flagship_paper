#!/usr/bin/env nextflow
/*
* AUTHOR: Mohadese Sayahian Dehkordi <mohadese.sayahiandehkordi@mail.mcgill.ca>, Daniel Taliun <daniel.taliun@mcgill.ca>
* VERSION: 2.0
* YEAR: 2023
*/

// How to run:
// nextflow run statistical_phasing.nf --vcf_path /path/to/chr2.vcf.gz --eagle /path/to//Eagle_v2.4.1/eagle --genetic_map /path/to/genetic_map_hg38_withX.txt.gz -with-report report.html

params.vcf_path = "/path/to/data/*.bcf" // Absolute path to the Input VCF/BCF files split by chromosome
params.eagle = "/path/to/executable/eagle" // Absolute path to the Eagle executable
params.genetic_map = "/path/to/genetic_map/genetic_map_hg38_withX.txt.gz" // Absolute path to the genetic map from Eagle software

process eagle_statistical_phasing {
    errorStrategy 'retry'
    maxRetries 3
    cache 'lenient'

    executor 'slurm'
    clusterOptions '--account=rrg-vmooser'

    cpus 4
    memory "16GB"
    time "6h"
    scratch '$SLURM_TMPDIR'

    input:
    path vcf

    output:
    path("*eagle_phased*")

    publishDir "phased/", pattern: "*eagle_phased*", mode: "copy"

    """
    ${params.eagle} --numThreads 4 --geneticMapFile ${params.genetic_map} --vcf ${vcf} --vcfOutFormat z --outPrefix ${vcf.getBaseName()}.eagle_phased
    """
}


process beagle_statistical_phasing {
    errorStrategy 'retry'
    maxRetries 3
    cache "lenient"
    
    executor 'slurm'
    clusterOptions '--account=rrg-vmooser'

    cpus 8
    memory "32GB"
    time "10h"
    scratch "$SLURM_TMPDIR"

    input:
    path vcf
    
    output:
    path("*beagle_phased*")
    
    publishDir "phased/", pattern: "*beagle_phased*", mode: "copy"
    
    """
    java -jar -Xmx32g ${params.beagle} window=25.0 overlap=2.5 nthreads=8 gt=${vcf} out=${vcf.getBaseName()}.beagle_phased
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

	eagle_statistical_phasing(Channel.fromPath(params.vcf_path))
	beagle_statistical_phasing(Channel.fromPath(params.vcf_path))

        //concat_ch = Channel.fromPath(params.snv_vcf_path).map{ vcf -> [vcf, vcf + ".tbi" ] }
        //sv_ch = Channel.fromPath(params.sv_vcf_path).map{ vcf -> [vcf, vcf + ".tbi" ] }

        //setGT_snv_ch = setGT_non_PASS_GT_SNVs(snv_ch)
        //recalculated_ch = recalculate_AF_SNVs(setGT_snv_ch)
        //left_aligned_ch = left_align_SNVs(recalculated_ch)
        //rm_dup_ch = remove_duplicates_SNVs(left_aligned_ch)
        //prep_snv_ch = filter_based_on_AC_SNVs(rm_dup_ch)
        //prep_sv_ch = prep_SVs(sv_ch)
        //filled_ref_ch = fill_REF_SVs(prep_sv_ch)

        //snv_with_chr_name_ch = get_chr_name_SNVs(snv_ch)
        //sv_with_chr_name_ch = get_chr_name_SVs(filled_ref_ch)

        //stat_phasing_ch = snv_with_chr_name_ch.join(sv_with_chr_name_ch)
        //stat_phasing_ch_combine = concat_vcfs(stat_phasing_ch)

        //phased_vcfs = beagle_statistical_phasing(concat_ch)
        //recal_phased = recalculate_AF_phased(phased_vcfs)

        //remove_singletons(recal_phased)
}
