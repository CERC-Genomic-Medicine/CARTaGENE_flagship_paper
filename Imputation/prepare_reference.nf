#!/usr/bin/env nextflow
/*
* AUTHOR: Mohadese Sayahian Dehkordi, <mohadese.sayahiandehkordi@mail.mcgill.ca>, Daniel Taliun, PhD <daniel.taliun@mcgill.ca>
* VERSION: 5.0
* YEAR: 2024
*/

// This pipeline builds reference panel for mputation with Minimac 4 from phased genotypes in VCF/BCF format.

// How to run:
// nextflow run prepare_reference.nf --phased_vcf_path "/path/to/chr*.vcf.gz*"  --study_vcf_path /path/to/*.vcf.gz --minimac4 /path/to/minimac4

params.phased_vcf_path = "/path/to/chr*.vcf.gz*" // Absolute path to the phased VCF/BCF file(s) with the corresponding TBI/CSI index file(s). One file per chromosome is expected. Chromosome X splitted into PAR1, PAR2, nonPAR.
params.minimac4 = "/path/to/minimac4" // Path to Minimac 4 executable

process filter_and_convert {
    errorStrategy 'finish'
    cache "lenient"

    executor 'slurm'
    clusterOptions '--account=def-vmooser'
    //scratch '$SLURM_TMPDIR'

    cpus 1
    memory "8GB"
    time "4h"


   input:
   tuple val(chrom), path(vcf), path(vcf_index)

   output:
   path("${chrom}.reference.*")

   publishDir "reference", mode: "copy"

   """
   # Remove singletons (recommended), and remove unnecessary FORMAT fields
   bcftools view -c2 ${vcf} -Ou | bcftools annotate -x FORMAT/PP -Oz -o ${chrom}.reference.vcf.gz
   bcftools index -t ${chrom}.reference.vcf.gz

   # Convert to Minimac 4 format
   ${params.minimac4} --compress-reference ${chrom}.reference.vcf.gz > ${chrom}.reference.msav
   """
}


workflow {
   // Here, we assume that the chromosome name is encoded in the file prefix e.g. chr1.study_name.vcf.gz, chr2.study_name.vcf.gz, ..., chrX_nonPAR.study_name.vcf.gz, chrX_PAR1.study_name.vcf.gz, chrX_PAR2.study_name.vcf.gz
   study_vcfs = Channel.fromFilePairs("${params.phased_vcf_path}", size: -1) { file -> file.getName().replaceAll("(.tbi|.csi)\$", "") }.map {it -> [ (it[0] =~ /chrX{1}[.]{1}|chrX_nonpar|chrX_par1|chrX_par2|chr[1-9][0-9]?/)[0], it[1][0], it[1][1]] }
   
   filter_and_convert(study_vcfs)
}
