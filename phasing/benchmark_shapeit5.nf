#!/usr/bin/env nextflow
/*
* AUTHOR: Mohadese Sayahian Dehkordi <mohadese.sayahiandehkordi@mail.mcgill.ca>, Daniel Taliun <daniel.taliun@mcgill.ca>
* VERSION: 3.0
* YEAR: 2023
*/

// This pipeline benchmarks phasing quality with ShapeIt and external reference panel, using samples with phase information inferred from parents.
// For each sample, the pipeline merges it's unphased genotypes with the unphased study VCF. Then it performes phasing of the merged VCFs.
// Pre-requisites/assumptions:
// 1. Only autosomal chromosomes (sex chromosomes were not tested) were used to perform the benchmarking.

// How to run:
// nextflow run benchmark_shapeit5.nf --vcf_path "/path/to/chr[1-9]*_*.vcf.gz*" --ref_vcf_path "/path/to/chr[1-9]*_*.vcf.gz*" --ref /path/to/reference.fa --vt /path/to/executable/vt --shapeit5_common /path/to/phase_common_static --shapeit5_rare /path/tophase_rare_static --shapeit5_genetic_map /path/to/chr2.b38.gmap.gz --shapeit5_chunks /path/to/chunks/b38/4cM/chunks_chr2.txt --shapeit5_switch /path/to/executable/switch_static -with-report report.html

params.vcf_path = "/path/to/chr[1-9]*_*.vcf.gz*" // Absolute path to the Input VCF/BCF file with the TBI/CSI indices for the benchmarking
params.ref_vcf_path = "/path/to/*chr[1-9]*_.vcf.gz*" // Absolute path to the Reference VCF/BCF file with the TBI/CSI indices
params.ref = "/path/to/reference/Homo_sapiens.GRCh38.fa" // Absolute path to the FASTA file with the human genome reference
params.vt = "/path/to/executable/vt" // Absolute path to the vt executable
params.shapeit5_common = "/path/to/executable/phase_common_static" // Absolute path to the Shapeit5 executable for common variants phasing
params.shapeit5_rare = "/path/to/executable/phase_rare_static" // Absolute path to the Shapeit5 executable for rare variants phasing
params.shapeit5_genetic_map = "/path/to/genetic_map/chr*.b38.gmap.gz" // Absolute path to the genetic map from Shapeit5 software
params.shapeit5_chunks = "/path/to/chunks/b38/4cM/chunks_chr2.txt" // Absolute path to the chunks file from Shapeit5 software
params.shapeit5_switch = "/path/to/executable/switch_static" // Absolute path the "switch" program executable from the Shapeit5


process prepare_target_vcf {
   cache 'lenient'

   executor 'slurm'
   clusterOptions '--account=def-vmooser'
   //scratch '$SLURM_TMPDIR'

   cpus 1
   memory "4GB"
   time "2h"

   input:
   tuple path(vcf), path(vcf_inddex)

   output:
   path("*.target_only.norm.bcf*")

   """
   # We expect one chromosome per file
   chr=`bcftools index --stats ${vcf} | cut -f1`

   # Note: we do not remove monomorphic variants, because they will be helpful when merging with benchmark samples to increase overlap
   bcftools +setGT -Ou ${vcf} -- -t q -n . -i 'FT!="PASS"' | bcftools +fixploidy -Ou | bcftools annotate -x INFO,^FORMAT/GT -Ou | ${params.vt} normalize - -r ${params.ref} -o + | ${params.vt} uniq + -o + | bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' -Ob -o \${chr}.target_only.norm.bcf
   bcftools index \${chr}.target_only.norm.bcf
   """
}


process merge_sample {
   cache 'lenient'

   executor 'slurm'
   clusterOptions '--account=def-vmooser'
   //scratch '$SLURM_TMPDIR'

   cpus 1
   memory "4GB"
   time "2h"

   input:
   tuple val(chr), path(vcf), path(vcf_index), path(benchmark_unphased_sample_vcf), path(benchmark_unphased_sample_vcf_index)

   output:
   tuple val(chr), val("${benchmark_unphased_sample_vcf.getSimpleName()}"), path("${chr}.*.with_target.unphased.bcf"), path("${chr}.*.with_target.unphased.bcf.csi")

   publishDir "benchmark_shapeit5/merged_with_target/", pattern: "*.with_target.unphased.bcf*", mode: "copy"

   """
   bcftools merge ${benchmark_unphased_sample_vcf} ${vcf} -Ou | bcftools annotate -x INFO -Ou | bcftools view -m2 -M2 -Ou | bcftools +fill-tags -Ou -- -t AN,AC,AF,NS,F_MISSING | bcftools view -e "F_MISSING > 0.1" -Ou | bcftools view -i "GT[0]!='./.'" -Ob -o ${chr}.${benchmark_unphased_sample_vcf.getSimpleName()}.with_target.unphased.bcf

   bcftools index ${chr}.${benchmark_unphased_sample_vcf.getSimpleName()}.with_target.unphased.bcf
   """
}


process shapeit5_phasing_with_reference {
    errorStrategy 'finish'
    cache "lenient"

    executor 'slurm'
    clusterOptions '--account=def-vmooser'
    scratch '$SLURM_TMPDIR'

    cpus 8
    memory "32GB"
    time "12h"

    input:
    tuple val(chr), val(benchmark_sample_name), path(target_vcf), path(target_vcf_index), path(ref_vcf), path(ref_vcf_index)
    path genetic_map_files
    path chunks_files

    output:
    tuple val(chr), val(benchmark_sample_name), path("${chr}.${benchmark_sample_name}.shapeit5_phased.bcf"), path("${chr}.${benchmark_sample_name}.shapeit5_phased.bcf.csi")

    """
    ${params.shapeit5_common} --thread 8 --input ${target_vcf} --reference ${ref_vcf} --map ${chr}.b38.gmap.gz --region ${chr} --filter-maf 0.001 --output common.scaffold_with_ref.bcf
    ${params.shapeit5_common} --thread 8 --input ${target_vcf} --scaffold common.scaffold_with_ref.bcf --map ${chr}.b38.gmap.gz --region ${chr} --filter-maf 0.001 --output common.phased.bcf

    while read LINE; do
        ID=\$(echo \$LINE | awk '{ print \$1; }')
        SRG=\$(echo \$LINE | awk '{ print \$3; }')
        IRG=\$(echo \$LINE | awk '{ print \$4; }')
        ${params.shapeit5_rare} --thread 8 --input ${target_vcf} --scaffold common.phased.bcf --map ${chr}.b38.gmap.gz --input-region \${IRG} --scaffold-region \${SRG} --output chunk\${ID}.phased.bcf
    done < chunks_${chr}.txt

    ls -1v chunk*.phased.bcf > files.txt

    bcftools concat --naive -f files.txt -o ${chr}.${benchmark_sample_name}.with_target.shapeit5_phased.bcf --threads 8
    bcftools index ${chr}.${benchmark_sample_name}.with_target.shapeit5_phased.bcf
    bcftools view -s ${benchmark_sample_name} ${chr}.${benchmark_sample_name}.with_target.shapeit5_phased.bcf -Ob -o ${chr}.${benchmark_sample_name}.shapeit5_phased.bcf
    bcftools index ${chr}.${benchmark_sample_name}.shapeit5_phased.bcf
    """
}


process count_errors {
    errorStrategy 'finish'
    cache "lenient"

    cpus 1
    memory "4GB"
    time "1h"

    input:
    tuple val(benchmark_sample_name), val(chr), path(unphased_with_target_vcf), path(unphased_with_target_index), path(phased_vcf), path(phased_index), path(benchmark_vcf), path(benchmark_index)

    output:
    path "*.variant.switch.txt.gz"
    
    publishDir "benchmark_shapeit5/with_ref/", pattern: "*.variant.switch.txt.gz", mode: "copy"

    """
    ${params.shapeit5_switch} -F ${unphased_with_target_vcf} -E ${phased_vcf} -V ${benchmark_vcf} --region ${chr} --output ${benchmark_sample_name}.${chr}
    
    rm *.block.switch.txt.gz
    rm *.calibration.switch.txt.gz
    rm *.flipsAndSwitches.txt.gz
    rm *.frequency.switch.txt.gz
    rm *.sample.switch.txt.gz
    rm *.sample.typing.txt.gz
    rm *.type.switch.txt.gz
    rm *.variant.typing.txt.gz
    """
}


workflow {
    // tuple: [vcf, tbi/csi]
    study_vcfs = Channel.fromFilePairs("${params.vcf_path}", size: -1) { file -> file.getName().replaceAll("(.tbi|.csi)\$", "") }.map {it -> it[1] }
    prepared_study_vcfs = prepare_target_vcf(study_vcfs)

    benchmark_unphased_samples_vcfs = Channel.fromPath("${workflow.projectDir}/Trios/*.unphased.bcf").map{ vcf -> [vcf, vcf + ".csi"]}
    merged_vcfs = merge_sample(prepared_study_vcfs.map { it -> [it[0].getSimpleName(), it[0], it[1]]}.combine(benchmark_unphased_samples_vcfs))

    ref_vcfs = Channel.fromFilePairs("${params.ref_vcf_path}", size: -1) { file -> file.getName().replaceAll("(.tbi|.csi)\$", "") }.map { it -> [ (it[0] =~ /chr[1-9][0-9]?/)[0], it[1][0], it[1][1] ] }

    shapeit5_phased = shapeit5_phasing_with_reference(merged_vcfs.combine(ref_vcfs, by: 0), Channel.fromPath(params.shapeit5_genetic_map).collect(), Channel.fromPath(params.shapeit5_chunks).collect())


    benchmark_phased_samples_vcfs = Channel.fromPath("${workflow.projectDir}/Trios/*.phased.bcf").map{ vcf -> [vcf.getSimpleName(), vcf, vcf + ".csi"]}
    count_errors(merged_vcfs.combine(shapeit5_phased, by: [0, 1]).map { it -> [it[1], it[0], it[2], it[3], it[4], it[5]] }.combine(benchmark_phased_samples_vcfs, by: 0))
}
