#!/usr/bin/env nextflow
/*
* AUTHOR: Mohadese Sayahian Dehkordi <mohadese.sayahiandehkordi@mail.mcgill.ca>, Daniel Taliun <daniel.taliun@mcgill.ca>
* VERSION: 2.0
* YEAR: 2023
*/

// This pipeline benchmarks phasing quality with Beagle, Eagle, and ShapeIt using samples with curated phase information from GIAB.
// For each GIAB sample, the pipeline merges it's unphased genotypes with the unphased study VCF. Then it performes phasing of the merged VCF with each of the three phasing tools.
// Pre-requisites/assumptions:
// 1. Only one autosomal chromosome (sex chromosomes were not tested) is used to perform the benchmarking.

// How to run:
// nextflow run statistical_phasing.nf --vcf_path /path/to/chr2.vcf.gz --eagle /path/to/Eagle_v2.4.1/eagle --eagle_genetic_map /path/to/genetic_map_hg38_withX.txt.gz --beagle /path/to/beagle.jar --beagle_genetic_map /path/to/plink.chr2.GRCh38.map --shapeit5_common /path/to/phase_common_static --shapeit5_rare /path/tophase_rare_static --shapeit5_genetic_map /path/to/chr2.b38.gmap.gz --shapeit5_chunks /path/to/chunks/b38/4cM/chunks_chr2.txt --shapeit5_switch /path/to/executable/switch_static -with-report report.html

params.vcf_path = "/path/to/data/*.vcf.gz" // Absolute path to the Input VCF/BCF file with the chromosome of choice for the benchmarking

params.eagle = "/path/to/executable/eagle" // Absolute path to the Eagle executable
params.eagle_genetic_map = "/path/to/genetic_map/genetic_map_hg38_withX.txt.gz" // Absolute path to the genetic map from Eagle software
params.beagle = "/path/to/executable/beagle.jar" // Absolute path to the Beagle JAR executable
params.beagle_genetic_map = "/path/to/genetic_map/plink.chr2.GRCh38.map" // Absolute path to the genetic map from Beagle software. If Beagle GRCh38 genetic map files do not have "chr" prefix, then run "for f in *.map; do sed -i 's/^/chr/' ${f}; done"
params.shapeit5_common = "/path/to/executable/phase_common_static" // Absolute path to the Shapeit5 executable for common variants phasing
params.shapeit5_rare = "/path/to/executable/phase_rare_static" // Absolute path to the Shapeit5 executable for rare variants phasing
params.shapeit5_genetic_map = "/path/to/genetic_map/chr2.b38.gmap.gz" // Absolute path to the genetic map from Shapeit5 software
params.shapeit5_chunks = "/path/to/chunks/b38/4cM/chunks_chr2.txt" // Absolute path to the chunks file from Shapeit5 software
params.shapeit5_switch = "/path/to/executable/switch_static" // Absolute path the "switch" program executable from the Shapeit5


process benchmark_sample {
   cache 'lenient'

   executor 'slurm'
   clusterOptions '--account=rrg-vmooser'

   cpus 1
   memory "4GB"
   time "2h"
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


process eagle_statistical_phasing {
    errorStrategy 'finish'
    cache 'lenient'

    executor 'slurm'
    clusterOptions '--account=rrg-vmooser'

    cpus 8
    memory "32GB"
    time "12h"
    scratch '$SLURM_TMPDIR'

    input:
    tuple path(vcf), path(tbi)

    output:
    path("*eagle_phased.merged.vcf.gz")

    """
    ${params.eagle} --numThreads 8 --geneticMapFile ${params.eagle_genetic_map} --vcf ${vcf} --vcfOutFormat z --outPrefix ${vcf.getSimpleName()}.eagle_phased.merged
    """
}


process beagle_statistical_phasing {
    errorStrategy 'finish'
    cache "lenient"
    
    executor 'slurm'
    clusterOptions '--account=rrg-vmooser'

    cpus 8
    memory "32GB"
    time "12h"
    scratch '$SLURM_TMPDIR'

    input:
    tuple path(vcf), path(tbi)
    
    output:
    path("*beagle_phased.merged.vcf.gz")
    
    """
    java -jar -Xmx32g ${params.beagle} window=25.0 overlap=2.5 nthreads=8 gt=${vcf} map=${params.beagle_genetic_map} out=${vcf.getSimpleName()}.beagle_phased.merged
    """
}


process shapeit5_statistical_phasing {
    errorStrategy 'finish'
    cache "lenient"
    
    executor 'slurm'
    clusterOptions '--account=rrg-vmooser'

    cpus 8
    memory "32GB"
    time "12h"
    scratch '$SLURM_TMPDIR'

    input:
    tuple path(vcf), path(tbi)
    
    output:
    path("*shapeit5_phased.merged.bcf")
    
    """
    chr=`bcftools index -s ${vcf} | cut -f1`
    ${params.shapeit5_common} --thread 8 --input ${vcf} --map ${params.shapeit5_genetic_map} --region \${chr} --filter-maf 0.001 --output common.phased.bcf

    while read LINE; do
        ID=\$(echo \$LINE | awk '{ print \$1; }')
        SRG=\$(echo \$LINE | awk '{ print \$3; }')
        IRG=\$(echo \$LINE | awk '{ print \$4; }')
        ${params.shapeit5_rare} --thread 8 --input ${vcf} --scaffold common.phased.bcf --map ${params.shapeit5_genetic_map} --input-region \${IRG} --scaffold-region \${SRG} --output chunk\${ID}.phased.bcf
    done < ${params.shapeit5_chunks}

    ls -1v chunk*.phased.bcf > files.txt

    bcftools concat --naive -f files.txt -o ${vcf.getSimpleName()}.shapeit5_phased.merged.bcf --threads 8
    """
}


process count_errors {
    errorStrategy 'finish'
    cache "lenient"

    cpus 1
    memory "4GB"
    time "1h"

    input:
    tuple val(sample_name), val(tool_name), path(phased_vcf), path(true_phases_vcf)

    output:
    path("*.sample.switch.txt.gz")

    publishDir "benchmark_phasing/", pattern: "*.sample.switch.txt.gz", mode: "copy"

    """
    bcftools index ${phased_vcf}
    bcftools index ${true_phases_vcf}
    chr=`bcftools index -s ${phased_vcf} | cut -f1`
    ${params.shapeit5_switch} -V ${true_phases_vcf} -E ${phased_vcf} --region \${chr} --output ${tool_name}.${sample_name}
    """
}


workflow {
    benchmark_unphased_samples_vcfs = Channel.fromPath("${workflow.projectDir}/GIAB/HG*.unphased.vcf.gz").map{ vcf -> [vcf, vcf + ".tbi"]}
    benchmark_phased_samples_vcfs = Channel.fromPath("${workflow.projectDir}/GIAB/HG*.phased.vcf.gz").map{ vcf -> [vcf.getSimpleName(), vcf]}
    study_vcf = Channel.fromPath("${params.vcf_path}").map{ vcf -> [vcf, vcf + ".tbi"]}.first() // we benchmark only one chromosome, thus only one file in the Channel is allowed

    merged_vcfs = benchmark_sample(benchmark_unphased_samples_vcfs, study_vcf) 

    shapeit5_phased = shapeit5_statistical_phasing(merged_vcfs)
    beagle_phased = beagle_statistical_phasing(merged_vcfs)
    eagle_phased = eagle_statistical_phasing(merged_vcfs)

    all_phased = beagle_phased.map{ vcf -> [vcf.getSimpleName(), 'beagle', vcf]}.combine(benchmark_phased_samples_vcfs, by: 0).concat(
        eagle_phased.map{ vcf -> [vcf.getSimpleName(), 'eagle', vcf]}.combine(benchmark_phased_samples_vcfs, by: 0),
	shapeit5_phased.map{ vcf -> [vcf.getSimpleName(), 'shapeit5', vcf]}.combine(benchmark_phased_samples_vcfs, by: 0),
    )

    count_errors(all_phased)
}
