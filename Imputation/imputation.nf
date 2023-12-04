#!/usr/bin/env nextflow
/*
* AUTHOR: Mohadese Sayahian Dehkordi, <mohadese.sayahiandehkordi@mail.mcgill.ca>, Daniel Taliun, PhD <daniel.taliun@mcgill.ca>
* Adapted and modified get_chunk and concat from https://github.com/CERC-Genomic-Medicine/shapeit4_pipeline/blob/main/Phasing.nf by Daniel Taliun, PhD <daniel.taliun@mcgill.ca>
* VERSION: 4.0
* YEAR: 2023
*/

// This pipeline performs imputation using minimac 4.

// How to run:
// nextflow run imputation.nf --ref_vcf_path /path/to/*.vcf.gz  --study_vcf_path /path/to/*.vcf.gz --minimac3 /path/to/minimac3 --minimac4 /path/to/minimac4  --window window_size --flank flank_size -with-report report.html

params.ref_vcf_path = "/path/to/data/*.vcf.gz" // Absolute path to the Input VCF/BCF file split by the chromosome
params.study_vcf_path = "/path/to/data/*.vcf.gz" // Absolute path to the Input VCF/BCF files split by the chromosome

params.minimac3 = "/path/to/executable/minimac3" // Absolute path to the minimac3 executable
params.minimac4 = "/path/to/executable/minimac4" // Absolute path to the minimac4 executable 

params.window = 2000000 // Size of sliding window that we want to perform imputation in.
params.flank = 0 //Size of the overlap of sliding windows.

process convert_ref_vcfs{
    errorStrategy 'retry'
    maxRetries 3
    cache "lenient"
    cpus 4

    executor 'slurm'
    clusterOptions '--account=rrg-vmooser'

    memory "64GB"
    time "5h"
    scratch "$SLURM_TMPDIR"

    input:
    tuple val(chr_name), path(ref_vcf), path(ref_vcf_index)
    
    output:
    tuple val(chr_name), path( "*.m3vcf.gz"), path("*.vcf.gz"), path("*.vcf.gz.tbi")
    
    publishDir "minimac_m3vcfs/", pattern: "*.m3vcf.gz", mode: "copy"

    script:
    """
    # First, Minimac requires us to remove "chr" from the name of the chromosome in the VCF file. i.e if the chromosome name in the VCF looks like chr9 we need to convert it to 9 for chromosome 9.
    
    echo "${chr_name}" > chroms1.txt
    chr=${chr_name}
    echo "\${chr:3}" > chroms2.txt

    paste chroms1.txt chroms2.txt > chr_name_conv.txt  

    bcftools annotate --rename-chrs chr_name_conv.txt ${ref_vcf} -Oz -o ${ref_vcf.getBaseName()}.vcf.gz
    bcftools index --tbi ${ref_vcf.getBaseName()}.vcf.gz

    # Second, Minimac require us to convert the reference haplotypes to m3vcf format which can be accuired by using Minimac 3 software. This format allow faster analysis.
    ${params.minimac3} --refHaps ${ref_vcf.getBaseName()}.vcf.gz --cpus 4 --processReference --rsid --prefix ${ref_vcf.getBaseName()}
    """
}


process rm_chr_name_study{
    errorStrategy 'retry'
    maxRetries 3
    cache "lenient"
    cpus 1

    executor 'slurm'
    clusterOptions '--account=rrg-vmooser'


    memory "4GB"
    time "00:30:00"
    scratch "$SLURM_TMPDIR"

    input:
    tuple val(chr_name), path(vcf), path(vcf_index)

    output:
    tuple val(chr_name), path("*.vcf.gz"), path("*.vcf.gz.tbi")
    
    script:
    """
    # Minimac requires us to remove "chr" from the name of the chromosome in the VCF file. i.e if the chromosome name in the VCF looks like chr9 we need to convert it to 9 for chromosome 9.

    echo "${chr_name}" > chroms1.txt
    chr=${chr_name}
    echo "\${chr:3}" > chroms2.txt

    paste chroms1.txt chroms2.txt > chr_name_conv.txt  

    bcftools annotate --rename-chrs chr_name_conv.txt $vcf -Oz -o ${vcf.getBaseName()}.vcf.gz
    bcftools index --tbi ${vcf.getBaseName()}.vcf.gz
    """
}

process get_chunks {
    errorStrategy 'retry'
    maxRetries 3
    cache "lenient"
    cpus 1

    executor 'slurm'
    clusterOptions '--account=rrg-vmooser'

    memory "4GB"
    time "00:30:00"
    scratch "$SLURM_TMPDIR"

    input:
    tuple val(chr_name), path(study_vcf), path(study_vcf_index), path(ref_m3vcf), path(ref_vcf), path(ref_vcf_index)

    output:
    tuple path("*.chunk"), val(chr_name), path(study_vcf), path(study_vcf_index), path(ref_m3vcf)

    """
    # To be able to perform imputation for large scale genotyping array data, we would need to split each chromosome to smaller chunks and
    # perform the imputation on each of the chunks separately and combine the imputed chunks at the end.

    # We extract the chromosome number, start and end position of the chromosome in the reference panel.
    chrom=`bcftools index -s ${ref_vcf} | cut -f1`
    start_bp=`bcftools view -HG ${ref_vcf} | head -n1 | cut -f2`
    stop_bp=`bcftools index -s ${ref_vcf} | cut -f2`

    # We identify the start and end position of each chunk by iterating over the start-end interval. And we save the information of each non-empty chunk in a txt file.
    extend=0
    for i in `seq \${start_bp} ${params.window} \${stop_bp}`; do
        if [ \${extend} -eq 0 ]; then
            chunk_start=\$((${params.flank} > i ? 1 : i - ${params.flank}))
        fi
        chunk_stop=\$((i + ${params.window} + ${params.flank}))
        n=`bcftools view -HG ${study_vcf} \${chrom}:\${chunk_start}-\${chunk_stop} | wc -l`
        if [ \${n} -gt 0 ]; then
            printf "\${chrom}\t\${chunk_start}\t\${chunk_stop}\t\${n}\n" > \${chrom}_\${chunk_start}_\${chunk_stop}.chunk
            extend=0
        else
            extend=1
        fi
    done
    if [ \${extend} -eq 1 ]; then
        printf "\${chrom}\t\${chunk_start}\t\${chunk_stop}\t\${n}\n" > \${chrom}_\${chunk_start}_\${chunk_stop}.chunk
    fi
    """
}



process minimac_imputation{
    errorStrategy 'retry'
    maxRetries 3
    cache "lenient"
    cpus 8

    executor 'slurm'
    clusterOptions '--account=rrg-vmooser'

    memory "128GB"
    time "36h"
    scratch "$SLURM_TMPDIR"

    input:
    tuple val(chr_name), path(chunk), file(study_vcf), file(study_vcf_index), file(ref_vcf)
     
    output:
    tuple val(chr_name), path("*.imp.dose.vcf.gz"), path("*.imp.dose.vcf.gz.tbi"), path("*.imp.empiricalDose.vcf.gz"), path("*.imp.empiricalDose.vcf.gz.tbi"), path("*.info")

    publishDir "imputed_dose_vcfs/", pattern: "*.imp.dose.vcf.gz*", mode: "copy"
    publishDir "empirical_dose_vcfs/", pattern: "*.imp.empiricalDose.vcf.gz*", mode: "copy"
    publishDir "imputed_info/", pattern: "*.info", mode: "copy"

    script:
    """
    # We extract the start and end position of each chunk from the chunk info file.
    chrom=`head -n1 ${chunk} | cut -f1`
    start_bp=`head -n1 ${chunk} | cut -f2`
    stop_bp=`head -n1 ${chunk} | cut -f3`

    # We extract the chunk interval from the input genotyping array VCF.
    bcftools view -r \${chrom}:\${start_bp}-\${stop_bp} ${study_vcf} -Oz -o study.\${chrom}_\${start_bp}_\${stop_bp}.vcf.gz 
    bcftools index --tbi study.\${chrom}_\${start_bp}_\${stop_bp}.vcf.gz

    # Minimac allow us to input the whole vcf for the reference panel and specify the start and end of the chunk we want to impute. This way we don't need to split the ref VCF.
    ${params.minimac4} --refHaps ${ref_vcf} --rsid --haps study.\${chrom}_\${start_bp}_\${stop_bp}.vcf.gz  --chr \${chrom} --start \${start_bp} --end \${stop_bp} --minRatio 0.00001 --prefix study.\${chrom}_\${start_bp}_\${stop_bp} --cpus 8  --meta --ignoreDuplicates

    # After performing imputation we need to reannotate the name of the chromosome in the final imputed VCF. i.e. from 9 to chr9.
    echo "${chr_name}" > chroms1.txt
    chr=${chr_name}
    echo "\${chr:3}" > chroms2.txt

    paste chroms2.txt chroms1.txt > chr_name_conv.txt  

    bcftools annotate --rename-chrs chr_name_conv.txt study.\${chrom}_\${start_bp}_\${stop_bp}.dose.vcf.gz -Oz -o study.\${chrom}_\${start_bp}_\${stop_bp}.imp.dose.vcf.gz
    bcftools index --tbi study.\${chrom}_\${start_bp}_\${stop_bp}.imp.dose.vcf.gz

    bcftools annotate --rename-chrs chr_name_conv.txt study.\${chrom}_\${start_bp}_\${stop_bp}.empiricalDose.vcf.gz -Oz -o study.\${chrom}_\${start_bp}_\${stop_bp}.imp.empiricalDose.vcf.gz
    bcftools index --tbi study.\${chrom}_\${start_bp}_\${stop_bp}.imp.empiricalDose.vcf.gz
    """
}

process concat_vcfs {
    errorStrategy 'retry'
    maxRetries 3
    cache "lenient"
    scratch true
    cpus 1
    executor 'slurm'
    clusterOptions '--account=rrg-vmooser'

    memory "512G"
    time "24h"
    scratch "$SLURM_TMPDIR"

    input:
    tuple val(chromosome), path(imputed_vcfs), path(imputed_vcfs_index), path(imputed_emp_vcfs), path(imputed_emp_vcfs_index), path(info)

    output:
    tuple path("*.imputed.dose.vcf.gz"), path("*.imputed.dose.vcf.gz.tbi"), path("*.imputed.empiricalDose.vcf.gz"), path("*.imputed.empiricalDose.vcf.gz.tbi")

    publishDir "final_imputed_vcfs/", pattern: "*.vcf.gz*", mode: "copy"

    script:
    """
    # We need to concat the imputed and meta-imputed chunks for each chromosome.
    chrom=\$(echo -n "${chromosome}" | tr -d '\n')
    for f in ${imputed_vcfs}; do echo \${f}; done | sort -V > files_list.txt
    bcftools concat -f files_list.txt  -Oz -o \${chrom}.imputed.dose.vcf.gz
    bcftools index --tbi \${chrom}.imputed.dose.vcf.gz
    for f in ${imputed_emp_vcfs}; do echo \${f}; done | sort -V > files_list.txt
    bcftools concat -f files_list.txt -Oz -o \${chrom}.imputed.empiricalDose.vcf.gz
    bcftools index --tbi \${chrom}.imputed.empiricalDose.vcf.gz
    """
}


workflow {

        ref_ch = Channel.fromPath(params.ref_vcf_path).map{ vcf -> [ vcf.name.toString().tokenize('.').get(0), vcf, vcf + ".tbi" ] }
        study_ch = Channel.fromPath(params.study_vcf_path).map{ vcf -> [ vcf.name.toString().tokenize('.').get(0), vcf, vcf + ".tbi" ] }

        ref_cnv_vcfs = convert_ref_vcfs(ref_ch)
        study_rm_chr_vcfs = rm_chr_name_study(study_ch)
        imputation_ch = study_rm_chr_vcfs.combine(ref_cnv_vcfs, by:[0])
        chunks = get_chunks(imputation_ch)
         
        all_chunks = chunks.flatMap { chunks, chromosome, study_vcf, study_vcf_index, ref_m3vcf ->
        chunks.collect { chunk -> [chromosome, chunk, study_vcf, study_vcf_index, ref_m3vcf] }
        }

        imputed_chunks = minimac_imputation(all_chunks)
        concat_vcfs(imputed_chunks.groupTuple(by:[0]))
}