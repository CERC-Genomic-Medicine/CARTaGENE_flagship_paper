#!/usr/bin/env nextflow
/*
* AUTHOR: Mohadese Sayahian Dehkordi <mohadese.sayahiandehkordi@mail.mcgill.ca>, Daniel Taliun <daniel.taliun@mcgill.ca>
* VERSION: 3.0
* YEAR: 2023
*/

// This pipeline benchmarks phasing quality with Beagle, Eagle, and ShapeIt using samples with phase information inferred from parents.
// For each sample, the pipeline merges it's unphased genotypes with the unphased study VCF. Then it performes phasing of the merged VCF with each of the three phasing tools.
// Pre-requisites/assumptions:
// 1. Only one autosomal chromosome (sex chromosomes were not tested) is used to perform the benchmarking.

// How to run:
// nextflow run benchmark_phasing.nf --vcf_path "/path/to/chr[1-9]*_*.vcf.gz*" --eagle /path/to/executable/eagle --eagle_genetic_map /path/to/genetic_map/genetic_map_hg38_withX.txt.gz --beagle /path/to/jar/beagle.jar --beagle_genetic_map /path/to/genetic_map/plink.chr*.GRCh38.map --shapeit5_common /path/to/executable/phase_common_static --shapeit5_rare /path/to/executable/phase_rare_static --shapeit5_genetic_map /path/to/genetic_map//chr*.b38.gmap.gz --shapeit5_chunks /path/to/chromosome_chunks/b38/4cM/chunks_chr*.txt --shapeit5_switch /path/to/executable/switch_static  --ref /path/to/reference.fa --vt /path/to/executable/vt -with-report report.html

params.vcf_path = "/path/to/data/*.vcf.gz*" // Absolute path to the input VCF/BCF and TBI/CSI files with the chromosome(s) of choice for the benchmarking

params.ref = "/path/to/reference/Homo_sapiens.GRCh38.fa" // Absolute path to the FASTA file with the human genome reference
params.vt = "/path/to/executable/vt" // Absolute path to the vt executable

// If you set all the parameters below to an empty string (i.e. ""), then eagle will not be run.
params.eagle = "/path/to/executable/eagle" // Absolute path to the Eagle executable
params.eagle_genetic_map = "/path/to/genetic_map/genetic_map_hg38_withX.txt.gz" // Absolute path to the genetic map from Eagle software

// If you set all the parameters below to an empty string (i.e. ""), then beagle will not be run. 
params.beagle = "/path/to/executable/beagle.jar" // Absolute path to the Beagle JAR executable
params.beagle_genetic_map = "/path/to/genetic_map/plink.chr*.GRCh38.map" // Absolute path to the genetic map from Beagle software. If Beagle GRCh38 genetic map files do not have "chr" prefix, then run "for f in *.map; do sed -i 's/^/chr/' ${f}; done"

// If you set all the parameters below to an empty string (i.e. ""), then shapeit5 will not be run.
params.shapeit5_common = "/path/to/executable/phase_common_static" // Absolute path to the Shapeit5 executable for common variants phasing
params.shapeit5_rare = "/path/to/executable/phase_rare_static" // Absolute path to the Shapeit5 executable for rare variants phasing
params.shapeit5_genetic_map = "/path/to/genetic_map/chr*.b38.gmap.gz" // Absolute path to the genetic map from Shapeit5 software
params.shapeit5_chunks = "/path/to/chunks/b38/4cM/chunks_chr*.txt" // Absolute path to the chunks file from Shapeit5 software

// The switch tool from the shapit5 toolset is used for counting switch errors.
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

   publishDir "benchmark_phasing/target_only/", pattern: "*.target_only.norm.bcf*", mode: "copy"

   """
   # We expect one chromosome per file
   chr=`bcftools index --stats ${vcf} | cut -f1`

   # Note: we do not remove monomorphic variants, because they will be helpful when merging with benchmark samples to increase overlap
   bcftools +fixploidy -Ou ${vcf} | bcftools annotate -x INFO,^FORMAT/GT -Ou | ${params.vt} normalize - -r ${params.ref} -o + | ${params.vt} uniq + -o + | bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' -Ob -o \${chr}.target_only.norm.bcf

   # For DRAGEN output
   #bcftools +setGT -Ou ${vcf} -- -t q -n . -i 'FT!="PASS"' | bcftools +fixploidy -Ou | bcftools annotate -x INFO,^FORMAT/GT -Ou | ${params.vt} normalize - -r ${params.ref} -o + | ${params.vt} uniq + -o + | bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' -Ob -o \${chr}.target_only.norm.bcf
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

   publishDir "benchmark_phasing/merged_with_target/", pattern: "*.with_target.unphased.bcf*", mode: "copy"

   """
   bcftools merge ${benchmark_unphased_sample_vcf} ${vcf} -Ou | bcftools annotate -x INFO -Ou | bcftools view -m2 -M2 -Ou | bcftools +fill-tags -Ou -- -t AN,AC,AF,NS,F_MISSING | bcftools view -e "F_MISSING > 0.1" -Ou | bcftools view -i "GT[0]!='./.'" -Ob -o ${chr}.${benchmark_unphased_sample_vcf.getSimpleName()}.with_target.unphased.bcf
   
   bcftools index ${chr}.${benchmark_unphased_sample_vcf.getSimpleName()}.with_target.unphased.bcf
   """
}



process eagle_statistical_phasing {
    errorStrategy 'retry'
    cache 'lenient'

    executor 'slurm'
    clusterOptions '--account=def-vmooser'
    scratch '$SLURM_TMPDIR'

    cpus 8
    memory "32GB"
    time "12h"
    
    input:
    tuple val(chr), val(benchmark_sample_name), path(vcf), path(index)
    path genetic_map_file

    output:
    tuple val(chr), val(benchmark_sample_name), path("${chr}.${benchmark_sample_name}.eagle_phased.bcf"), path("${chr}.${benchmark_sample_name}.eagle_phased.bcf.csi")

    publishDir "benchmark_phasing/eagle_phased/", pattern: "${chr}.${benchmark_sample_name}.eagle_phased.bcf*", mode: "copy"

    """
    ${params.eagle} --numThreads 8 --geneticMapFile ${genetic_map_file} --vcf ${vcf} --vcfOutFormat b --outPrefix ${chr}.${benchmark_sample_name}.with_target.eagle_phased
    bcftools index ${chr}.${benchmark_sample_name}.with_target.eagle_phased.bcf
    bcftools view -s ${benchmark_sample_name} ${chr}.${benchmark_sample_name}.with_target.eagle_phased.bcf -Ob -o ${chr}.${benchmark_sample_name}.eagle_phased.bcf
    bcftools index ${chr}.${benchmark_sample_name}.eagle_phased.bcf
    """
}


process beagle_statistical_phasing {
    errorStrategy 'retry'
    cache "lenient"
    
    executor 'slurm'
    clusterOptions '--account=def-vmooser'
    scratch '$SLURM_TMPDIR'

    cpus 8
    memory "32GB"
    time "12h"
    
    input:
    tuple val(chr), val(benchmark_sample_name), path(vcf), path(index)
    path genetic_map_files
    
    output:
    tuple val(chr), val(benchmark_sample_name), path("${chr}.${benchmark_sample_name}.beagle_phased.bcf"), path("${chr}.${benchmark_sample_name}.beagle_phased.bcf.csi")
    
    publishDir "benchmark_phasing/beagle_phased/", pattern: "${chr}.${benchmark_sample_name}.beagle_phased.bcf*", mode: "copy"

    """
    bcftools view ${vcf} | java -jar -Xmx32g ${params.beagle} window=25.0 overlap=2.5 nthreads=8 gt=/dev/stdin map=plink.${chr}.GRCh38.map out=${chr}.${benchmark_sample_name}.with_target.beagle_phased
    bcftools index ${chr}.${benchmark_sample_name}.with_target.beagle_phased.vcf.gz
    bcftools view -s ${benchmark_sample_name} ${chr}.${benchmark_sample_name}.with_target.beagle_phased.vcf.gz -Ob -o ${chr}.${benchmark_sample_name}.beagle_phased.bcf
    bcftools index ${chr}.${benchmark_sample_name}.beagle_phased.bcf
    """
}


process shapeit5_statistical_phasing {
    errorStrategy 'retry'
    cache "lenient"
    
    executor 'slurm'
    clusterOptions '--account=def-vmooser'
    scratch '$SLURM_TMPDIR'

    cpus 8
    memory "32GB"
    time "12h"
    
    input:
    tuple val(chr), val(benchmark_sample_name), path(vcf), path(index)
    path genetic_map_files
    path chunks_files
    
    output:
    tuple val(chr), val(benchmark_sample_name), path("${chr}.${benchmark_sample_name}.shapeit5_phased.bcf"), path("${chr}.${benchmark_sample_name}.shapeit5_phased.bcf.csi")
   
    publishDir "benchmark_phasing/shapeit5_phased/", pattern: "${chr}.${benchmark_sample_name}.shapeit5_phased.bcf*", mode: "copy"

    """
    ${params.shapeit5_common} --thread 8 --input ${vcf} --map ${chr}.b38.gmap.gz --region ${chr} --filter-maf 0.001 --output common.phased.bcf

    while read LINE; do
        ID=\$(echo \$LINE | awk '{ print \$1; }')
        SRG=\$(echo \$LINE | awk '{ print \$3; }')
        IRG=\$(echo \$LINE | awk '{ print \$4; }')
        ${params.shapeit5_rare} --thread 8 --input ${vcf} --scaffold common.phased.bcf --map ${chr}.b38.gmap.gz --input-region \${IRG} --scaffold-region \${SRG} --output chunk\${ID}.phased.bcf
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
    tuple val(benchmark_sample_name), val(tool_name), val(chr), path(unphased_with_target_vcf), path(unphased_with_target_index), path(phased_vcf), path(phased_index), path(benchmark_vcf), path(benchmark_index)

    output:
    path "*.variant.switch.txt.gz"

    publishDir "benchmark_phasing/${tool_name}/", pattern: "*.variant.switch.txt.gz", mode: "copy"
    
    """
    ${params.shapeit5_switch} -F ${unphased_with_target_vcf} -E ${phased_vcf} -V ${benchmark_vcf} --region ${chr} --output ${tool_name}.${benchmark_sample_name}.${chr}

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
   
    if ((params.beagle != "") && (params.beagle_genetic_map != "")) {
    	beagle_phased = beagle_statistical_phasing(merged_vcfs, Channel.fromPath(params.beagle_genetic_map).collect())
    } else {
    	beagle_phased = Channel.empty()
    }
    if ((params.eagle != "") && (params.eagle_genetic_map != "")) {
	eagle_phased = eagle_statistical_phasing(merged_vcfs, Channel.fromPath(params.eagle_genetic_map).collect())
    } else {
    	eagle_phased = Channel.empty()
    }
    if ((params.shapeit5_common != "") && (params.shapeit5_rare != "") && (params.shapeit5_genetic_map != "") && (params.shapeit5_chunks != "")) {
    	shapeit5_phased = shapeit5_statistical_phasing(merged_vcfs, Channel.fromPath(params.shapeit5_genetic_map).collect(), Channel.fromPath(params.shapeit5_chunks).collect())
    } else {
    	shapeit5_phased = Channel.empty()
    }

    benchmark_phased_samples_vcfs = Channel.fromPath("${workflow.projectDir}/Trios/*.phased.bcf").map{ vcf -> [vcf.getSimpleName(), vcf, vcf + ".csi"]}

    // sample, tool_name, chr, unphased_with_target_vcf, unphased_with_target_index, phased_vcf, phased_index, benchmark_vcf, benchmark_index
    all_phased = merged_vcfs.combine(beagle_phased, by: [0, 1]).map { it -> [it[1], 'beagle', it[0], it[2], it[3], it[4], it[5]] }.combine(benchmark_phased_samples_vcfs, by: 0).concat(
       merged_vcfs.combine(eagle_phased, by: [0, 1]).map { it -> [it[1], 'eagle', it[0], it[2], it[3], it[4], it[5]] }.combine(benchmark_phased_samples_vcfs, by: 0),
       merged_vcfs.combine(shapeit5_phased, by: [0, 1]).map { it -> [it[1], 'shapeit5', it[0], it[2], it[3], it[4], it[5]] }.combine(benchmark_phased_samples_vcfs, by: 0)
    )
   
    count_errors(all_phased)
}
