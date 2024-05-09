#!/usr/bin/env nextflow
/*
* AUTHOR: Mohadese Sayahian Dehkordi <mohadese.sayahiandehkordi@mail.mcgill.ca>, Daniel Taliun <daniel.taliun@mcgill.ca>
* VERSION: 5.0
* YEAR: 2024
*/

// This pipeline performs imputation using Minimac 4. If input data is not phased, then the reference-based phasing using Eagle 4.2.1 will be performed before imputation.
// Note: this pipeline is intended for benchmarking imputation in small sample (N<1000) using local reference panel (N<10000) and, thus, the parallelization is done per chromosome rather than per chromosomal region (i.e. by chunking).


// How to run:
// nextflow run benchmark_imputation.nf --study_vcf_path "/path/to/*.vcf.gz*" --imp_ref_msav_path "/path/to/*.msav*" --minimac4 /path/to/minimac4 --phasing_ref_vcf_path "/path/to/*.vcf.gz*" --eagle /path/to/Eagle_v2.4.1/eagle --eagle_genetic_map /path/to/genetic_map_hg38_withX.txt.gz  -with-report report.html

params.study_vcf_path = "/path/to/data/*.vcf.gz*" // Absolute path to the input study array-based VCF/BCF file(s) with the corresponding TBI/CSI index file(s). One file per chromosome is expected.
params.imp_ref_msav_path = "/path/to/data/*.msav*" // Absolute path to the reference files *.msav, *.msav.INFO.tsv.gz, *.msav.INFO.txt for Minimac 4. One triplet (i.e. *.msav, *.msav.INFO.tsv.gz, *.msav.INFO.txt) per chromosome is expected; chrX is split into chrX_nonpar.[suffix], chrX_par1.[suffix] and chrX_par2.[suffix] files. These reference files will be used for the imputation.
params.minimac4 = "/path/to/minimac4"

// If your data is already phased, then set this to empty string ""
params.phasing_ref_vcf_path = "/path/to/data/*.vcf.gz*" // Absolute path to the reference VCF/BCF file(s) with the TBI/CSI index file(s). One file per chromosome is expected; chrX is split into chrX_nonpar.[suffix], chrX_par1.[suffix] and chrX_par2.[suffix] files. This reference files will be used for the reference-based phasing.

// If your data is already phased, then you don't need to use these arguments.
params.eagle = "/path/to/executable/eagle" // Absolute path to the Eagle 2.4.1 executable
params.eagle_genetic_map = "/path/to/genetic_map/genetic_map_hg38_withX.txt.gz" // Absolute path to the genetic map from Eagle software

// PAR region coordinates for GRCh38. Change only when working with different human genome reference build.
params.par1_region = "chrX:10001-2781479"
params.par2_region = "chrX:155701383-156030895"


process phase_with_ref {
    errorStrategy 'finish'
    cache "lenient"

    executor 'slurm'
    clusterOptions '--account=def-vmooser'

    cpus 8
    memory "32GB"
    time "4h"
    //scratch "$SLURM_TMPDIR"

    input:
    tuple val(chromosome), path(study_vcf), path(study_vcf_index), val(ref_prefix), path(ref_vcf), path(ref_vcf_index)
        
    output:
    tuple val(ref_prefix), path("*.phased.vcf.gz"), path("*.phased.vcf.gz.tbi")
    
    script:
    """
    # For the reference-based phasing we don't need to split chromosome X PAR and non-PAR regions as long as the reference panel is split itself.
    ${params.eagle} --vcfRef ${ref_vcf}  --vcfTarget ${study_vcf} --numThreads 8 --geneticMapFile ${params.eagle_genetic_map} --outPrefix ${ref_prefix}.phased --allowRefAltSwap --vcfOutFormat z
    bcftools index --tbi ${ref_prefix}.phased.vcf.gz
    """
}


process impute {
    errorStrategy 'finish'
    cache "lenient"

    executor 'slurm'
    clusterOptions '--account=def-vmooser'

    cpus 8
    memory "32GB"
    time "2h"
    //scratch "$SLURM_TMPDIR"

    input:
    tuple val(ref_prefix), path(study_vcf), path(study_vcf_index), path(ref_msav), path(ref_msav_info), path(ref_msav_info_index), path(ref_msav_info_header)

    output:
    tuple val(ref_prefix), path("*.imputed.vcf.gz"), path("*.imputed.vcf.gz.tbi")

    publishDir "imputed", mode: "copy"

    script:
    """
    # Impute
    if [ "${ref_prefix}" = "chrX_par1" ]; then
       bcftools view ${study_vcf} ${params.par1_region} -Ob -o temp.bcf
       bcftools index temp.bcf
       ${params.minimac4} ${ref_msav} temp.bcf -t 8 -f GT,DS,HDS -o imputed_noinfo.bcf
    elif [ "${ref_prefix}" = "chrX_par2" ]; then
       bcftools view ${study_vcf} ${params.par2_region} -Ob -o temp.bcf
       bcftools index temp.bcf
       ${params.minimac4} ${ref_msav} temp.bcf -t 8 -f GT,DS,HDS -o imputed_noinfo.bcf
    elif [ "${ref_prefix}" = "chrX_nonpar" ]; then
       bcftools view -t ^${params.par1_region},${params.par2_region} ${study_vcf} -Ob -o temp.bcf
       bcftools index temp.bcf
       ${params.minimac4} ${ref_msav} temp.bcf -t 8 -f GT,DS,HDS -o imputed_noinfo.bcf
    else
       ${params.minimac4} ${ref_msav} ${study_vcf} -t 8 -f GT,DS,HDS -o imputed_noinfo.bcf
    fi
    
    # Index
    bcftools index imputed_noinfo.bcf

    # Add back INFO fields and relevant meta-lines for the SVs
    bcftools annotate -a ${ref_msav_info} -h ${ref_msav_info_header} -c CHROM,POS,~ID,REF,ALT,INFO/END,INFO/SVTYPE,INFO/SVLEN imputed_noinfo.bcf -Oz -o ${ref_prefix}.imputed.vcf.gz
    bcftools index -t ${ref_prefix}.imputed.vcf.gz
    """
}


workflow {
   study_vcfs = Channel.fromFilePairs("${params.study_vcf_path}", size: -1) { file -> file.getName().replaceAll("(.tbi|.csi)\$", "") }.map {it -> [ (it[0] =~ /chrX{1}|chr[1-9][0-9]?/)[0], it[1][0], it[1][1]] }
   //study_vcfs.view()

   // Input tuple: [msav, msav.INFO.tsv.gz, msav.INFO.tsv.gz.tbi, msav.INFO.txt]
   imp_ref_msav = Channel.fromFilePairs("${params.imp_ref_msav_path}", size: -1) { file -> file.getName().replaceAll("\\.msav.*\$", "") }.map { it -> [ (it[0] =~ /chrX{1}|chr[1-9][0-9]?/)[0], (it[0] =~ /chrX{1}[.]{1}|chrX_nonpar|chrX_par1|chrX_par2|chr[1-9][0-9]?/)[0], it[1][0], it[1][1], it[1][2], it[1][3] ] }
   //imp_ref_msav.view()

   if (params.phasing_ref_vcf_path != "") {
      phasing_ref_vcfs = Channel.fromFilePairs("${params.phasing_ref_vcf_path}", size: -1) { file -> file.getName().replaceAll("(.tbi|.csi)\$", "") }.map {it -> [ (it[0] =~ /chrX{1}|chr[1-9][0-9]?/)[0], (it[0] =~ /chrX{1}[.]{1}|chrX_nonpar|chrX_par1|chrX_par2|chr[1-9][0-9]?/)[0], it[1][0], it[1][1]]}
      phased_with_ref = phase_with_ref(study_vcfs.combine(phasing_ref_vcfs, by: 0))
      
      imputed = impute(phased_with_ref.combine(imp_ref_msav.map {it -> [it[1], it[2], it[3], it[4], it[5]]}, by: 0))
   } else {
      phased_with_ref = study_vcfs.combine(imp_ref_msav, by: 0).map { it -> [it[3], it[1], it[2], it[4], it[5], it[6], it[7]] }
      imputed = impute(phased_with_ref)
   }
}
