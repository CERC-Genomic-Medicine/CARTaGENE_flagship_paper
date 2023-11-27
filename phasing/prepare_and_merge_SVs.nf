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

//params.snv_vcf_path = "/path/to/data/*.vcf.gz" // Absolute path to the Input VCF/BCF files split by chromosome. The index files must be located in the same directory.
//params.sv_vcf_path = "/path/to/data/*.vcf.gz" // Absolute path to the Input VCF/BCF files split by chromosome. The index files must be located in the same directory.
//params.ref = "/path/to/reference/Homo_sapiens.GRCh38.fa" // Absolute path to the FASTA file with the human genome reference


process clean_SV_VCF {
    cache "lenient"
    
    //executor "slurm"
    //clusterOptions "--account=rrg-vmooser"

    cpus 1
    memory "8GB"
    time "1h"
    //scratch '$SLURM_TMPDIR'

    input:
    tuple path(sv_vcf), path(sv_index)
    
    output:
    tuple stdout, path("*.SVs.cleaned.vcf.gz"), path("*.SVs.cleaned.vcf.gz.tbi")
    
    publishDir "SV_vcfs/", pattern: "*.SVs.cleaned.vcf.gz*", mode: "copy"
    """
    chr=`bcftools index -s ${sv_vcf} | cut -f1`
    if [ "\$chr" = "chrX" ]; then
    # Process the PAR regions (where males and females are diploid) in a similar way as autosomal chromosomes below.
    bcftools view ${sv_vcf} ${params.par1_region} -Ou | bcftools +fill-tags -Ou -- -t F_MISSING | bcftools view -e 'F_MISSING>0.1' -Ou | bcftools +fixploidy -Ou | bcftools +fill-from-fasta -Oz -o chrX_PAR1.SVs.cleaned.vcf.gz -- -c REF -f ${params.ref} 
    bcftools index -t chrX_PAR1.SVs.cleaned.vcf.gz

    bcftools view ${sv_vcf} ${params.par2_region} -Ou | bcftools +fill-tags -Ou -- -t F_MISSING | bcftools view -e 'F_MISSING>0.1' -Ou | bcftools +fixploidy -Ou | bcftools +fill-from-fasta -Oz -o chrX_PAR2.SVs.cleaned.vcf.gz -- -c REF -f ${params.ref} 
    bcftools index -t chrX_PAR2.SVs.cleaned.vcf.gz

    # Process non-PAR part of chromosome X. Sample must be either haploid or diploid across all non-PAR variants.
    # First, remove PAR regions:
    bcftools view -t ^${params.par1_region},${params.par2_region} ${sv_vcf} -Ob -o temp.bcf
    bcftools index temp.bcf
    # Second, check that ploidy is the same across all variants for each individual
    # The bcftools +check-ploidy will output multiple rows per individual if the ploidy was different. If ploidy is the same, then it will output only one row per individual.
    bcftools +check-ploidy temp.bcf | tail -n+2 | cut -f1 > ploidy.txt
    bcftools query -l temp.bcf > all_sample_names.txt
    if ! cmp -s ploidy.txt all_sample_names.txt; then
            exit 1 # there are individuals with different ploidy at different variants
    fi
    # Now we can apply filters similar to the autosomal chromosomes
    bcftools +fill-tags temp.bcf -Ou -- -t F_MISSING | bcftools view -e 'F_MISSING>0.1' -Ou | bcftools +fixploidy -Ou | bcftools +fill-from-fasta -Oz -o chrX_nonPAR.SVs.cleaned.vcf.gz -- -c REF -f ${params.ref} 
    bcftools index -t chrX_nonPAR.SVs.cleaned.vcf.gz
    echo -n "\${chr}"
    else
    # a) Compute missigness for each variant
    # b) Remove variants with missigness >0.1
    # c) Fill the missing REF alleles with the allele from the reference genome (upstream SV calling tools set all REF alleles to `.` which is not in line with the VCF specs)
    bcftools +fill-tags ${sv_vcf} -Ou -- -t F_MISSING | bcftools view -e 'F_MISSING>0.1' -Ou | bcftools +fixploidy -Ou | bcftools +fill-from-fasta -Oz -o \${chr}.SVs.cleaned.vcf.gz -- -c REF -f ${params.ref} 
    bcftools index -t \${chr}.SVs.cleaned.vcf.gz
    echo -n "\${chr}"
    fi
    """
}


process concat_VCF {
    cache "lenient"
  
    //executor "slurm"
    //clusterOptions "--account=rrg-vmooser"

    cpus 1
    memory "8GB"
    time "2h"
    //scratch '$SLURM_TMPDIR'

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
    cleaned_sv_ch = clean_SV_VCF(sv_ch)

    sv_chrX_vcfs = cleaned_sv_ch.filter { it[0].contains('chrX') }
    sv_chrX_ch = sv_chrX_vcfs.flatMap { chr, vcfs, vcf_indices ->
         vcfs.collect { vcf -> [vcf.getSimpleName(), vcf, vcf.toString() + ".tbi"] }
         }
    sv_autosomal_vcfs = cleaned_sv_ch.filter { !it[0].contains('chrX')}

    snv_chrX_vcfs = snv_ch.filter { it[0].contains('chrX') }
    snv_autosomal_vcfs = snv_ch.filter { !it[0].contains('chrX')}


    combined_autosomal = snv_autosomal_vcfs.combine(sv_autosomal_vcfs, by: 0)
    combined_chrX = snv_chrX_vcfs.combine(sv_chrX_ch, by: 0)
    all_vcfs = combined_autosomal.concat(combined_chrX)

    concat_VCF(all_vcfs)


}
