#!/usr/bin/env nextflow
/*
* AUTHOR: Mohadese Sayahian Dehkordi <mohadese.sayahiandehkordi@mail.mcgill.ca>, Daniel Taliun <daniel.taliun@mcgill.ca>
* VERSION: 2.0
* YEAR: 2023
*/

// This pipeline benchmarks phasing quality with Beagle at different "iteratation" values using samples with curated phase information from GIAB.
// For each GIAB sample, the pipeline merges it's unphased genotypes with the unphased study VCF. Then, it performes phasing of the merged VCF with Beagle with different "iteration" values.
// Pre-requisites/assumptions:
// 1. Only one autosomal chromosome (sex chromosomes were not tested) is used to perform the benchmarking.

// How to run:
// nextflow run benchmark_beagle.nf --vcf_path /path/to/chr2.vcf.gz --beagle /path/to/beagle.jar --beagle_genetic_map /path/to/plink.chr2.GRCh38.map --shapeit5_switch /path/to/executable/switch_static -with-report report.html

params.vcf_path = "/path/to/data/*.vcf.gz" // Absolute path to the Input VCF/BCF file with the chromosome of choice for the benchmarking

params.beagle = "/path/to/executable/beagle.jar" // Absolute path to the Beagle JAR executable
params.beagle_genetic_map = "/path/to/genetic_map/plink.chr2.GRCh38.map" // Absolute path to the genetic map from Beagle software. If Beagle GRCh38 genetic map files do not have "chr" prefix, then run "for f in *.map; do sed -i 's/^/chr/' ${f}; done"
params.shapeit5_switch = "/path/to/executable/switch_static" // Absolute path the "switch" program executable from the Shapeit5


params.beagle_iterations = [8, 12, 16, 20, 24, 28, 32]


process benchmark_sample {
   cache 'lenient'

   executor 'slurm'
   clusterOptions '--account=rrg-vmooser'

   cpus 1
   memory "4GB"
   time "4h"
   scratch '$SLURM_TMPDIR'

   input:
   tuple path(benchmark_unphased_sample_vcf), path(benchmark_unphased_sample_vcf_tbi)
   tuple path(vcf), path(vcf_tbi)

   output:
   path("*.unphased.merged.vcf.gz*")

   """
   bcftools merge -Ou ${benchmark_unphased_sample_vcf} ${vcf} | bcftools annotate -x INFO -Ou | bcftools view -m2 -M2 -Ou | bcftools +fill-tags -Ou -- -t AN,AC,AF,NS | bcftools view -e "NS < 2" -Ou | bcftools view -i "GT[0]!='./.'" -Oz -o ${benchmark_unphased_sample_vcf.getSimpleName()}.unphased.merged.vcf.gz
   bcftools index -t ${benchmark_unphased_sample_vcf.getSimpleName()}.unphased.merged.vcf.gz
   """
}


process beagle_statistical_phasing {
    errorStrategy 'finish'
    cache "lenient"
    
    executor 'slurm'
    clusterOptions '--account=rrg-vmooser'

    cpus 8
    memory "32GB"
    time "24h"
    scratch '$SLURM_TMPDIR'

    input:
    tuple path(vcf), path(tbi)
    each n_iterations

    output:
    tuple val("${vcf.getSimpleName()}"), val("${n_iterations}"), path("*beagle_phased.merged.vcf.gz")
    
    """
    # Fix the random seed and nthreads to ensure that any observed differences in phasing were only because of different iterations
    java -jar -Xmx32g ${params.beagle} seed=1985 window=25.0 overlap=2.5 iterations=${n_iterations} nthreads=8 gt=${vcf} map=${params.beagle_genetic_map} out=${vcf.getSimpleName()}.${n_iterations}.beagle_phased.merged
    """
}


process count_errors {
    errorStrategy 'finish'
    cache "lenient"

    cpus 1
    memory "4GB"
    time "1h"

    input:
    tuple val(sample_name), val(n_iterations), path(phased_vcf), path(true_phases_vcf)

    output:
    path("*.sample.switch.txt.gz")

    publishDir "benchmark_beagle/", pattern: "*.sample.switch.txt.gz", mode: "copy"

    """
    bcftools index ${phased_vcf}
    bcftools index ${true_phases_vcf}
    chr=`bcftools index -s ${phased_vcf} | cut -f1`
    ${params.shapeit5_switch} -V ${true_phases_vcf} -E ${phased_vcf} --region \${chr} --output ${sample_name}.${n_iterations}
    """
}


workflow {
    benchmark_unphased_samples_vcfs = Channel.fromPath("${workflow.projectDir}/GIAB/HG*.unphased.vcf.gz").map{ vcf -> [vcf, vcf + ".tbi"]}
    benchmark_phased_samples_vcfs = Channel.fromPath("${workflow.projectDir}/GIAB/HG*.phased.vcf.gz").map{ vcf -> [vcf.getSimpleName(), vcf]}
    study_vcf = Channel.fromPath("${params.vcf_path}").map{ vcf -> [vcf, vcf + ".tbi"]}.first() // we benchmark only one chromosome, thus only one file in the Channel is allowed

    merged_vcfs = benchmark_sample(benchmark_unphased_samples_vcfs, study_vcf)

    beagle_phased = beagle_statistical_phasing(merged_vcfs, Channel.from(params.beagle_iterations))
    count_errors(beagle_phased.combine(benchmark_phased_samples_vcfs, by: 0))
}
