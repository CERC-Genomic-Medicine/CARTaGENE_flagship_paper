#!/usr/bin/env nextflow
/*
* AUTHOR: Mohadese Sayahian Dehkordi, <mohadese.sayahiandehkordi@mail.mcgill.ca>, Rose Laflamme, <rose.laflamme@umontreal.ca>
* VERSION: 1.0
* YEAR: 2023
*/
process fill_REF_SVs {
    errorStrategy 'retry'
    maxRetries 3
    cache "lenient"
    cpus 1
    memory "16GB"
    time "7h"
    input:
    tuple path(sv_vcf), path(sv_index)
    
    output:
    tuple path("*.filled_REF.vcf.gz"), path("*.filled_REF.vcf.gz.tbi")
    
    publishDir "Filled_REF_SV_vcfs/", pattern: "*.vcf.gz", mode: "copy"   
    publishDir "Filled_REF_SV_vcfs/", pattern: "*.vcf.gz.tbi", mode: "copy"    

    """
    bcftools +fill-from-fasta $sv_vcf -Oz -o ${sv_vcf.getBaseName()}.filled_REF.vcf.gz -- -c REF -f ${params.ref} 
    bcftools tabix --tbi ${sv_vcf.getBaseName()}.filled_REF.vcf.gz
    """
}
process handle_duplicates_SVs {
    //errorStrategy 'retry'
    //maxRetries 3
    cache "lenient"
    cpus 1
    memory "16GB"
    time "1h"
    input:
    tuple path(sv_vcf), path(sv_index)
    
    output:
    tuple path("*.handle_dup.vcf.gz"), path("*.handle_dup.vcf.gz.tbi")
    
    publishDir "handle_dp_vcfs/", pattern: "*.vcf.gz", mode: "copy"   
    publishDir "handle_dp_vcfs/", pattern: "*.vcf.gz.tbi", mode: "copy"    

    """
    handle_dp_SVs.py -i $sv_vcf -o ${sv_vcf.getBaseName()}.handle_dp.vcf.gz
    bcftools view ${sv_vcf.getBaseName()}.handle_dp.vcf.gz -Oz -o ${sv_vcf.getBaseName()}.handle_dup.vcf.gz
    bcftools index --tbi ${sv_vcf.getBaseName()}.handle_dup.vcf.gz
    """
}

process get_chr_name_SNVs {
    //errorStrategy 'retry'
    //maxRetries 3
    cache "lenient"
    cpus 1
    memory "4GB"
    time "1h"
    input:
    tuple path(snv_vcf), path(snv_index)

    output:
    tuple stdout, path(snv_vcf), path(snv_index)

    """
    tabix -l ${snv_vcf} 
    """
}

process get_chr_name_SVs {
    //errorStrategy 'retry'
    //maxRetries 3
    cache "lenient"
    cpus 1
    memory "4GB"
    time "1h"
    input:
    tuple path(sv_vcf), path(sv_index)

    output:
    tuple stdout, path(sv_vcf), path(sv_index)

    """
    tabix -l ${sv_vcf} 
    """

}
process concat_vcfs {
    //errorStrategy 'retry'
    //maxRetries 3
    cache "lenient"
    cpus 1
    memory "64GB"
    time "8h"
    input:
    tuple val(chr_name), path(snv_vcf), path(snv_index), path(sv_vcf), path(sv_index)

    output:
    tuple path("*.combined.vcf.gz"), path("*.combined.vcf.gz.tbi")

    publishDir "combined_vcfs/", pattern: "*.vcf.gz", mode: "copy"   
    publishDir "combined_vcfs/", pattern: "*.vcf.gz.tbi", mode: "copy"    

    """
    bcftools concat $snv_vcf $sv_vcf | bcftools sort -m 64G -Oz -o ${snv_vcf.getSimpleName()}.combined.vcf.gz
    bcftools index --tbi ${snv_vcf.getSimpleName()}.combined.vcf.gz
    """
}
workflow {
        snv_ch = Channel.fromPath(params.snv_vcf_path).map{ vcf -> [vcf, vcf + ".tbi" ] }
        sv_ch = Channel.fromPath(params.sv_vcf_path).map{ vcf -> [vcf, vcf + ".tbi" ] }

        //setGT_snv_ch = setGT_non_PASS_GT_SNVs(snv_ch)
        //recalculated_ch = recalculate_AF_SNVs(setGT_snv_ch)
        //left_aligned_ch = left_align_SNVs(recalculated_ch)
        //rm_dup_ch = remove_duplicates_SNVs(left_aligned_ch)
        //prep_snv_ch = filter_based_on_AC_SNVs(rm_dup_ch)
        //prep_sv_ch = prep_SVs(sv_ch)
        //recalculated_sv_ch = recalculate_AF_SVs(prep_sv_ch)
        filled_ref_ch = fill_REF_SVs(sv_ch)
        handle_dp_ch = handle_duplicates_SVs(filled_ref_ch)

        snv_with_chr_name_ch = get_chr_name_SNVs(snv_ch)
        sv_with_chr_name_ch = get_chr_name_SVs(handle_dp_ch)
        stat_phasing_ch = snv_with_chr_name_ch.join(sv_rename_ch)
        stat_phasing_ch_combine = concat_vcfs(stat_phasing_ch)

        //phased_vcfs = beagle_statistical_phasing(stat_phasing_ch_combine)
        //recal_phased = recalculate_AF_phased(phased_vcfs)

        //remove_singletons(recal_phased)
}
