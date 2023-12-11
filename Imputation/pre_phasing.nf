#!/usr/bin/env nextflow
/*
* AUTHOR: Mohadese Sayahian Dehkordi <mohadese.sayahiandehkordi@mail.mcgill.ca>, Daniel Taliun <daniel.taliun@mcgill.ca>
* VERSION: 4.0
* YEAR: 2023
*/

// This pipeline performs reference-based phasing using EAGLE.

// How to run:
// nextflow run pre_phasing.nf --ref_vcf_path /path/to/*.vcf.gz  --study_vcf_path /path/to/*.vcf.gz --eagle /path/to/Eagle_v2.4.1/eagle --eagle_genetic_map /path/to/genetic_map_hg38_withX.txt.gz  -with-report report.html

params.ref_vcf_path = "/path/to/data/*.vcf.gz" // Absolute path to the haplotype reference panel VCF/BCF file (WGS) split by the chromosome
params.study_vcf_path = "/path/to/data/*.vcf.gz" // Absolute path to the GWAS study VCF/BCF files (genotyping array) split by the chromosome

params.eagle = "/path/to/executable/eagle" // Absolute path to the Eagle executable
params.eagle_genetic_map = "/path/to/genetic_map/genetic_map_hg38_withX.txt.gz" // Absolute path to the genetic map from Eagle software

process eagle_phased{
    errorStrategy 'retry'
    maxRetries 3
    cache "lenient"

    executor 'slurm'
    clusterOptions '--account=rrg-vmooser'

    cpus 8
    memory "128GB"
    time "48:00:00"
    scratch "$SLURM_TMPDIR"

    input:
    tuple val(chromosome), path(ref_vcf), path(ref_vcf_index), path(study_vcf), path(study_vcf_index) 
        
    output:
    path("*.phased.final.vcf.gz*")
    publishDir "phased/", pattern: "*.vcf.gz*", mode: "copy"

    script:
    """
    ${params.eagle} --vcfRef ${ref_vcf}  --vcfTarget ${study_vcf} --numThreads 8 --geneticMapFile ${params.eagle_genetic_map} --outPrefix ${study_vcf.getBaseName()}.phased.final --allowRefAltSwap --vcfOutFormat z
    bcftools index --tbi ${study_vcf.getBaseName()}.phased.final.vcf.gz
    """
}

workflow {
        ref_ch = Channel.fromPath(params.ref_vcf_path).map{ vcf -> [ vcf.name.toString().tokenize('.')[0], vcf, vcf + ".tbi" ] }
        study_ch = Channel.fromPath(params.study_vcf_path).map{ vcf -> [ vcf.name.toString().tokenize('.')[0], vcf, vcf + ".tbi" ] }
   
        phased = eagle_phased(ref_ch.combine(study_ch, by:[0]))
}