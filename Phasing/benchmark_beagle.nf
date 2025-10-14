#!/usr/bin/env nextflow
/*
* AUTHOR: Mohadese Sayahian Dehkordi <mohadese.sayahiandehkordi@mail.mcgill.ca>, Daniel Taliun <daniel.taliun@mcgill.ca>
* VERSION: 2.0
* YEAR: 2023
*/

// This pipeline benchmarks phasing quality with Beagle at different "iteratation" values using samples with phase information inferred from parents.
// For each sample, the pipeline merges it's unphased genotypes with the unphased study VCF. Then, it performes phasing of the merged VCF with Beagle with different "iteration" values.
// Pre-requisites/assumptions:
// 1. Only autosomal chromosomes (sex chromosomes were not tested) were used to perform the benchmarking.

// How to run:
// nextflow run benchmark_beagle.nf --vcf_path "/path/to/chr[1-9]*_*.vcf.gz*" --ref /home/user/reference/GRCh38.fa --vt /home/user/tools/vt/vt --beagle /path/to/beagle.jar --beagle_genetic_map /path/to/plink.chr*.GRCh38.map --shapeit5_switch /path/to/executable/switch_static -with-report report.html

params.vcf_path = "/path/to/chr[1-9]*_*.vcf.gz*" // Absolute path to the Input VCF/BCF file with the TBI/CSI indices for the benchmarking

params.ref = "/path/to/reference/Homo_sapiens.GRCh38.fa" // Absolute path to the FASTA file with the human genome reference
params.vt = "/path/to/executable/vt" // Absolute path to the vt executable

params.beagle = "/path/to/executable/beagle.jar" // Absolute path to the Beagle JAR executable
params.beagle_genetic_map = "/path/to/genetic_map/plink.chr*.GRCh38.map" // Absolute path to the genetic map from Beagle software. If Beagle GRCh38 genetic map files do not have "chr" prefix, then run "for f in *.map; do sed -i 's/^/chr/' ${f}; done"
params.shapeit5_switch = "/path/to/executable/switch_static" // Absolute path the "switch" program executable from the Shapeit5


params.beagle_iterations = [8, 12, 16, 20, 24, 28, 32]


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

   """
   bcftools merge ${benchmark_unphased_sample_vcf} ${vcf} -Ou | bcftools annotate -x INFO -Ou | bcftools view -m2 -M2 -Ou | bcftools +fill-tags -Ou -- -t AN,AC,AF,NS,F_MISSING | bcftools view -e "F_MISSING > 0.1" -Ou | bcftools view -i "GT[0]!='./.'" -Ob -o ${chr}.${benchmark_unphased_sample_vcf.getSimpleName()}.with_target.unphased.bcf

   bcftools index ${chr}.${benchmark_unphased_sample_vcf.getSimpleName()}.with_target.unphased.bcf
   """
}


process beagle_statistical_phasing {
    errorStrategy 'finish'
    cache "lenient"

    executor 'slurm'
    clusterOptions '--account=def-vmooser'
    scratch '$SLURM_TMPDIR'

    cpus 8
    memory "32GB"
    time "24h"

    input:
    tuple val(chr), val(benchmark_sample_name), path(vcf), path(index)
    path genetic_map_files
    each n_iterations

    output:
    tuple val(chr), val(benchmark_sample_name), val("${n_iterations}"), path("${chr}.${benchmark_sample_name}.${n_iterations}.beagle_phased.bcf"), path("${chr}.${benchmark_sample_name}.${n_iterations}.beagle_phased.bcf.csi")

    """
    # Fix the random seed and nthreads to ensure that any observed differences in phasing were only because of different iterations
    bcftools view ${vcf} | java -jar -Xmx32g ${params.beagle} seed=1985 window=25.0 overlap=2.5 iterations=${n_iterations} nthreads=8 gt=/dev/stdin map=plink.${chr}.GRCh38.map out=${chr}.${benchmark_sample_name}.${n_iterations}.with_target.beagle_phased
    bcftools index ${chr}.${benchmark_sample_name}.${n_iterations}.with_target.beagle_phased.vcf.gz
    bcftools view -s ${benchmark_sample_name} ${chr}.${benchmark_sample_name}.${n_iterations}.with_target.beagle_phased.vcf.gz -Ob -o ${chr}.${benchmark_sample_name}.${n_iterations}.beagle_phased.bcf
    bcftools index ${chr}.${benchmark_sample_name}.${n_iterations}.beagle_phased.bcf
    """
}


process count_errors {
    errorStrategy 'finish'
    cache "lenient"

    cpus 1
    memory "4GB"
    time "1h"

    input:
    tuple val(benchmark_sample_name), val(chr), path(unphased_with_target_vcf), path(unphased_with_target_index), val(n_iterations), path(phased_vcf), path(phased_index), path(benchmark_vcf), path(benchmark_index)
    
    output:
    tuple path("*.sample.switch.txt.gz"), path("*.frequency.switch.txt.gz")

    publishDir "benchmark_beagle/${n_iterations}/", pattern: "*.sample.switch.txt.gz", mode: "copy"
    publishDir "benchmark_beagle/${n_iterations}/", pattern: "*.frequency.switch.txt.gz", mode: "copy"

    """
    ${params.shapeit5_switch} -F ${unphased_with_target_vcf} -E ${phased_vcf} -V ${benchmark_vcf} --region ${chr} --output ${benchmark_sample_name}.${chr}.${n_iterations}
    """
}


workflow {
    // tuple: [vcf, tbi/csi]
    study_vcfs = Channel.fromFilePairs("${params.vcf_path}", size: -1) { file -> file.getName().replaceAll("(.tbi|.csi)\$", "") }.map {it -> it[1] }
    prepared_study_vcfs = prepare_target_vcf(study_vcfs)

    benchmark_unphased_samples_vcfs = Channel.fromPath("${workflow.projectDir}/Trios/*.unphased.bcf").map{ vcf -> [vcf, vcf + ".csi"]}
    merged_vcfs = merge_sample(prepared_study_vcfs.map { it -> [it[0].getSimpleName(), it[0], it[1]]}.combine(benchmark_unphased_samples_vcfs))

    beagle_phased = beagle_statistical_phasing(merged_vcfs, Channel.fromPath(params.beagle_genetic_map).collect(), Channel.from(params.beagle_iterations))


    benchmark_phased_samples_vcfs = Channel.fromPath("${workflow.projectDir}/Trios/*.phased.bcf").map{ vcf -> [vcf.getSimpleName(), vcf, vcf + ".csi"]}

    count_errors(merged_vcfs.combine(beagle_phased, by: [0, 1]).map { it -> [it[1], it[0], it[2], it[3], it[4], it[5], it[6]] }.combine(benchmark_phased_samples_vcfs, by: 0))
}
