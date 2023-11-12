#!/usr/bin/env nextflow
/*
* AUTHOR: Mohadese Sayahian Dehkordi <mohadese.sayahiandehkordi@mail.mcgill.ca>, Daniel Taliun <daniel.taliun@mcgill.ca>
* VERSION: 2.0
* YEAR: 2023
*/

// Pre-requisites/assumptions:
// 1. Script was designed for autosomal chromosomes (i.e. 1-22) and chromosome X
// 2. Upstream variant caller represented males as diploids (i.e. 0/0 or 1/1) in non-PAR regions of chromosome X
// 3. Multi-allelic variants were split into multiple VCF entries (one entry per alternate allele)
// 4. Make sure that bcftools is installed/loaded on your system
// 5. Compile/install vt (https://github.com/atks/vt)

// How to execute:
// nextflow run prepare_SNVs_indels.nf --vcf_path "/home/user/data/chr*.vcf.gz" --ref "/home/user/reference/GRCh38.fa" --vt "/home/user/tools/vt/vt"

params.vcf_path = "/path/to/data/*.vcf.gz" // Absolute path to the Input VCF/BCF files split by chromosome
params.ref = "/path/to/reference/Homo_sapiens.GRCh38.fa" // Absolute path to the FASTA file with the human genome reference
params.vt = "/path/to/executable/vt" // Absolute path to the vt executable

process clean_VCF {
    errorStrategy 'retry'
    maxRetries 3
    cache "lenient"
    
    executor "slurm"
    clusterOptions "--account=rrg-vmooser"

    cpus 1
    memory "16GB"
    time "6h"
    scratch '$SLURM_TMPDIR'

    input:
    path vcf
    
    output:
    tuple path("*.norm.dedup.vcf.gz"), path("*.norm.dedup.vcf.gz.tbi")
    
    publishDir "SNVs_indels/", pattern: "*.vcf.gz*", mode: "copy"

    """
    # a. set GT to missing if it didn't pass QC
    # b. change `.` to `./.` in autosomal chromosomes (`.` is an artifact from variant calling pipeline). When no options are given to bcftools fixploidy, then all samples are treated as females and, thus, all samples will be set to diploids on chromosome X.
    # c. remove all INFO and FORMAT fields which will be not needed in downstream analyses
    # d. recalculate allele counts after setting some GT fields to mssing
    # e. remove any monomorphic variants which were introduced in previous steps
    # f. left-align indels
    # g. remove any duplicated records after left-alignment (see vt documentation)
    # h. update ID field to reflect new reference and alternate alleles after left-alignment
    
    bcftools +setGT -Ou ${vcf}  -- -t q -n . -i 'FT!="PASS"' | bcftools +fixploidy -Ou | bcftools annotate -x INFO,^FORMAT/GT -Ou | bcftools +fill-tags -Ou -- -t AN,AC,AF,NS | bcftools view -c 1 -Ou | ${params.vt} normalize - -r ${params.ref} -o + | ${params.vt} uniq + -o + | bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' -Oz -o ${vcf.getSimpleName()}.norm.dedup.vcf.gz
    bcftools tabix --tbi ${vcf.getSimpleName()}.norm.dedup.vcf.gz
    """
}

workflow {
	clean_VCF(Channel.fromPath(params.vcf_path))
}
