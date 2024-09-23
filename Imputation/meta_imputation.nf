/*
* AUTHOR: Mohadese Sayahian Dehkordi <mohadese.sayahiandehkordi@mail.mcgill.ca>, Daniel Taliun <daniel.taliun@mcgill.ca>
* VERSION: 2.0
* YEAR: 2024
*/


// Example:
// nextflow meta_imputation.nf --panel1_name="Study_1" --imputed_vcfs_panel1="path/to/imputed_panel1/chr*.{dose,empiricalDose}.vcf.gz" --panel2_name="Study_2" --imputed_vcfs_panel1="path/to/imputed_panel2/chr*.{dose,empiricalDose}.vcf.gz" --metaminimac2="/path/to/executable/MetaMinimac2"

params.panel1_name = "Study_1" // Name of the first imputation panel. This name will be used in the final VCF header. No whitespaces or special characters allowed. 
params.imputed_vcfs_panel1 = "/path/to/imputed_panel1/chr*.{dose,empiricalDose}.vcf.gz" // VCFs/BCFs with genotypes imputed using Panel 1. One file per chromosome. Index (tbi/csi) is not required.
params.panel2_name = "Study_2" // Name of the second imputation panel. This name will be used in the final VCF header. No whitespaces or special characters  allowed.
params.imputed_vcfs_panel2 = "/path/to/imputed_panel2/chr*.{dose,empiricalDose}.vcf.gz" // VCFs/BCFs with genotypes imputed using Panel 2. One file per chromosome. Index (tbi/csi) is not required.
params.metaminimac2 = "/path/to/executable/MetaMinimac2" // Absolute path to the MetaMinimac2 executable


process meta_impute {
    errorStrategy 'finish'
    cache "lenient"

    executor 'slurm'
    clusterOptions '--account=rrg-vmooser'
    
    cpus 1
    memory "8GB"
    time "7d"
    scratch "$SLURM_TMPDIR"

    input:
    tuple val(chr), path(dose_vcf_panel1, stageAs: "${params.panel1_name}.dose.vcf.gz"), path(empirical_dose_vcf_panel1, stageAs: "${params.panel1_name}.empiricalDose.vcf.gz"), path(dose_vcf_panel2, stageAs: "${params.panel2_name}.dose.vcf.gz"), path(empirical_dose_vcf_panel2, stageAs: "${params.panel2_name}.empiricalDose.vcf.gz")
     
    output:
    tuple path("*.metaDose.vcf.gz"), path("*.metaWeights.gz")

    publishDir "meta_imputed_vcfs/", pattern: "*.gz"
    
    """
    ${params.metaminimac2} -i ${params.panel1_name}:${params.panel2_name} --weight -o ${chr}
    """
}



workflow {
    // Chromosome X must be split into files starring with chrX_nonpar, chrX_par1, and chrX_par2.
    // Output tuple: [chr, dose VCF/BCF,  empiricalDose VCF/BCF]. 
    imputed_vcfs_panel1 = Channel.fromFilePairs("${params.imputed_vcfs_panel1}", size: -1) { file -> file.getSimpleName() }.map { it -> [(it[0] =~ /chrX_nonpar|chrX_par1|chrX_par2|chr[1-9][0-9]?/)[0], it[1][0], it[1][1]] }
    imputed_vcfs_panel2 = Channel.fromFilePairs("${params.imputed_vcfs_panel2}", size: -1) { file -> file.getSimpleName() }.map { it -> [(it[0] =~ /chrX_nonpar|chrX_par1|chrX_par2|chr[1-9][0-9]?/)[0], it[1][0], it[1][1]] }
    //imputed_vcfs_panel1.view()
    //imputed_vcfs_panel2.view()

    joined_by_chr = imputed_vcfs_panel1.join(imputed_vcfs_panel2)
    //joined_by_chr.view()

    meta_impute(joined_by_chr)
}
