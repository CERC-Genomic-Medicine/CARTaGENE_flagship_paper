#!/usr/bin/env nextflow
/*
* AUTHOR: Mohadese Sayahian Dehkordi, <mohadese.sayahiandehkordi@mail.mcgill.ca>
* VERSION: 1.0
* YEAR: 2023
*/

process beagle_statistical_phasing {
    errorStrategy 'retry'
    maxRetries 3
    cache "lenient"
    cpus 8
    memory "32GB"
    time "10h"
    input:
    tuple path(chr), path(index)
    
    output:
    tuple path("*.ref.vcf.gz"), path("*.ref.vcf.gz.tbi"), path("*.ref.log")
    
    publishDir "phased/ref_vcfs/", pattern: "*.vcf.gz", mode: "copy"
    publishDir "phased/ref_vcfs/", pattern: "*.vcf.gz.tbi", mode: "copy"
    publishDir "phased/ref_logs/", pattern: "*.ref.log", mode: "copy"
    
    """
    java -jar -Djava.io.tmpdir=./temp/ -Xmx32g ${params.beagle} window=25.0 overlap=2.5 nthreads=8 gt=$chr out=${chr.getSimpleName()}.ref 
    bcftools index --tbi ${chr.getSimpleName()}.ref.vcf.gz
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
        concat_ch = Channel.fromPath(params.snv_vcf_path).map{ vcf -> [vcf, vcf + ".tbi" ] }
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

        phased_vcfs = beagle_statistical_phasing(concat_ch)
        recal_phased = recalculate_AF_phased(phased_vcfs)

        remove_singletons(recal_phased)
}