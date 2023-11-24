#!/usr/bin/env nextflow
/*
* AUTHOR: Mohadese Sayahian Dehkordi <mohadese.sayahiandehkordi@mail.mcgill.ca>, Rose Laflamme <rose.laflamme@umontreal.ca>, Daniel Taliun <daniel.taliun@mcgill.ca>
* VERSION: 2.0
* YEAR: 2023
*/

// This pipeline merges together VCFs with SNVs and VCFs with SVs.
// Pre-requisites/assumptions:
// 1) Script was designed for autosomal chromosomes, since at the time there were no SVs data available for chromosome X.
// 2) The SNV file names are prefixed with "chr", e.g. chr1.snvs.vcf.gz, chr2.snvs.vcf.gz, ... .  
// 3) Only samples that are present in both VCFs will be kept.
// 4) Input VCFs must be split by chromosome and indexed.
// 5) Uses bcftools

// How to run:
// nextflow run prepare_and_merge_SVs.nf --snv_vcf_path="/path/to/snv/*.vcf.gz" --sv_vcf_path="/path/to/sv/*.vcf.gz" --ref /path/to/Homo_sapiens.GRCh38.fa

params.snv_vcf_path = "/path/to/data/*.vcf.gz" // Absolute path to the Input VCF/BCF files split by chromosome. The index files must be located in the same directory.
params.sv_vcf_path = "/path/to/data/*.vcf.gz" // Absolute path to the Input VCF/BCF files split by chromosome. The index files must be located in the same directory.
params.ref = "/path/to/reference/Homo_sapiens.GRCh38.fa" // Absolute path to the FASTA file with the human genome reference


process clean_SV_VCF {
    cache "lenient"
    
    executor "slurm"
    clusterOptions "--account=rrg-vmooser"

    cpus 1
    memory "8GB"
    time "1h"
    scratch '$SLURM_TMPDIR'

    input:
    tuple path(sv_vcf), path(sv_index)
    
    output:
    tuple stdout, path("*.SVs.cleaned.vcf.gz"), path("*.SVs.cleaned.vcf.gz.tbi")
    
    """
    # Get the chromosome name for this VCF/BCF
    chr=`bcftools index -s ${sv_vcf} | cut -f1`
    
    # a) Compute missigness for each variant
    # b) Remove variants with missigness >0.1
    # c) Fill the missing REF alleles with the allele from the reference genome (upstream SV calling tools set all REF alleles to `.` which is not in line with the VCF specs)
    bcftools +fill-tags ${sv_vcf} -Ou -- -t F_MISSING | bcftools view -e 'F_MISSING>0.1' -Ou | bcftools +fill-from-fasta -Oz -o \${chr}.SVs.cleaned.vcf.gz -- -c REF -f ${params.ref} 
    bcftools index -t \${chr}.SVs.cleaned.vcf.gz
    echo -n "\${chr}"
    """
}


process concat_VCF {
    cache "lenient"
  
    executor "slurm"
    clusterOptions "--account=rrg-vmooser"

    cpus 1
    memory "8GB"
    time "2h"
    scratch '$SLURM_TMPDIR'

    input:
    tuple val(chr_name), path(snv_vcf), path(snv_index), path(sv_vcf), path(sv_index)

    output:
    path("*.SNVs_Indels_SVs.vcf.gz*")

    publishDir "combined_vcfs/", pattern: "*.SNVs_Indels_SVs.vcf.gz*", mode: "copy"

    """
    # First, we need to make sure that the sample names and samples' order in both files match, otherwise the bcftools concat will fail
    bcftools query -l ${snv_vcf} | sort > samples1.txt
    bcftools query -l ${sv_vcf} | sort > samples2.txt
    comm --check-order -12 samples1.txt samples2.txt > samples_union.txt
    
    # Then, subset and re-order samples in both files
    bcftools view -S samples_union.txt ${snv_vcf} -Ob -o temp1.bcf
    bcftools index temp1.bcf
    bcftools view -S samples_union.txt ${sv_vcf} -Ob -o temp2.bcf
    bcftools index temp2.bcf

    # Now, we can concatenate files. The --allow-overlaps option will ensure that the resulting concatenated VCF is sorted by position.
    bcftools concat --allow-overlaps temp1.bcf temp2.bcf -Oz -o ${chr_name}.SNVs_Indels_SVs.vcf.gz
    bcftools index --tbi ${chr_name}.SNVs_Indels_SVs.vcf.gz
    """
}


workflow {
    snv_ch = Channel.fromPath(params.snv_vcf_path).map{ vcf -> [vcf.getSimpleName(), vcf, vcf + ".tbi" ] }
    sv_ch = Channel.fromPath(params.sv_vcf_path).map{ vcf -> [vcf, vcf + ".tbi" ] }
    concat_VCF(snv_ch.combine(clean_SV_VCF(sv_ch), by: 0))
}
