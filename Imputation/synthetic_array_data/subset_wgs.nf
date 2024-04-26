#!/usr/bin/env nextflow
/*
* AUTHOR: Peyton McClelland <peyton.mcclelland@mail.mcgill.ca>, Daniel Taliun <daniel.taliun@mcgill.ca>
* VERSION: 2.0
* YEAR: 2024
*/

// How to run:
// nextflow run subset_wgs.nf --sample_vcf_path "/path/to/*.vcf.gz*" --array_positions_file /path/to/array_positions --haploids_file /path/to/haploids.txt

params.sample_vcf_path = "/path/to/chr*.vcf.gz*" // One VCF/BCF file per individuals. Single VCF/BCF file must contain all chromosomes and be indexed with TBI/CSI.
params.array_positions_file = "/path/to/array_positions.txt" // File with array positions. Has two tab-delimited columns: chromosome, position. No header.
params.haploids_file = "/path/to/haploids.txt" // File with individuals who are haploids in chrX non-PAR region. One sample name per line. No header.

params.par1_region = "chrX:10001-2781479"
params.par2_region = "chrX:155701383-156030895"
params.fixploidy_region1 = "chrX 1 10000 M 1"
params.fixploidy_region2 = "chrX 2781480 155701382 M 1"

process subset_sample_wgs {
   errorStrategy 'finish'
   cache "lenient"

   executor 'slurm'
   clusterOptions '--account=def-vmooser'

   cpus 1
   memory "4GB"
   time "1h"

   input:
   tuple val(sample), path(sample_vcf), path(sample_vcf_index)
   path array_positions

   output:
   path("${sample}.array.bcf*")

   //publishDir "array", mode: "copy"

   """
   # Count number of SNPs with missing FILTER
   n=`bcftools view -m2 -M2 -v snps -HG -f. ${sample_vcf} | wc -l`
   
   if [[ \${n} -eq 0 ]]; then # FILTER is set for all SNPs
       # Subset to bi-allelic array positions; then, set GT to missing for SNPs with non-PASS filter; then remove unnecessary annotations.
       bcftools view -m2 -M2 -v snps -R ${array_positions} ${sample_vcf} -Ou | bcftools +setGT -Ou -- -t q -n . -e 'FILTER="PASS"' | bcftools annotate -xFILTER,INFO,^FORMAT/GT -Ob -o ${sample}.array.bcf
   else # FILTER not set
       # Here, we assume that the data comes from GATK DRAGEN and we apply DRAGEN's hard filters
       bcftools view -m2 -M2 -v snps -R ${array_positions} ${sample_vcf} -Ou | bcftools +setGT -Ou -- -t q -n . -i 'QUAL<10.4139 || INFO/DP<=1 || FORMAT/GQ==0 || LOD<6.3' | bcftools annotate -xINFO,^FORMAT/GT -Ob -o ${sample}.array.bcf
   fi
   bcftools index ${sample}.array.bcf
   """
}


process merge_samples {
   cache "lenient"

   cpus 1
   memory "4GB"
   time "1h"

   input:
   path samples

   output:
   path("chr*.merged.bcf*")

   //publishDir "array", mode: "copy"

   """
   find . -name "*.array.bcf" | sort > bcf_list.txt

   # Merge samples together:
   # 1) For variants that are not present in a sample, genotype is set to 0/0 in that sample
   # 2) All multi-allelic entries are removed after merging
   # 3) QUAL field is set to missing, because it has no use at this point
   bcftools merge -l bcf_list.txt -0 -Ou | bcftools view -m2 -M2 -Ou | bcftools annotate -xQUAL -Ob -o merged.bcf 
   bcftools index merged.bcf

   # Split by chromosome. We are interested in autosomals and X.
   for i in {1..22} X; do
      if `bcftools index -s merged.bcf | grep -q "chr\${i}"`; then
         bcftools view merged.bcf chr\${i} -Ob -o chr\${i}.merged.bcf
         bcftools index chr\${i}.merged.bcf
      fi
   done
   """
}


process finalize_chromosome {
   cache "lenient"

   cpus 1
   memory "4GB"
   time "1h"

   input:
   tuple val(chrom), path(vcf), path(vcf_index)
   path haploids

   output:
   path("${chrom}.array.qc.vcf.gz*")

   publishDir "array", mode: "copy"

   """
   if [ "${chrom}" = "chrX" ]; then
      # Process X pseudo-autosomal regions in the same way as autosomal chromosomes
      bcftools view ${vcf} ${params.par1_region} -Ou | bcftools +fixploidy -Ou | bcftools +fill-tags -Ou -- -t AN,AC,AF,NS,F_MISSING | bcftools view -e 'F_MISSING>0.1' -Ou | bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' -Ou | bcftools filter -s PASS -Oz -o ${chrom}_par1.array.qc.vcf.gz
      bcftools index ${chrom}_par1.array.qc.vcf.gz
      bcftools view ${vcf} ${params.par2_region} -Ou | bcftools +fixploidy -Ou | bcftools +fill-tags -Ou -- -t AN,AC,AF,NS,F_MISSING | bcftools view -e 'F_MISSING>0.1' -Ou | bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' -Ou | bcftools filter -s PASS -Oz -o ${chrom}_par2.array.qc.vcf.gz
      bcftools index ${chrom}_par2.array.qc.vcf.gz

      # Process non-PAR part of chromosome X. Sample must be either haploid or diploid across all non-PAR variants.
      # First, remove PAR regions:
      bcftools view -t ^${params.par1_region},${params.par2_region} chrX.merged.bcf -Ob -o temp_nonpar.bcf
      bcftools index temp_nonpar.bcf

      # Second, prepare files with diploids and input files for bcftools +fixploidy
      bcftools query -l temp_nonpar.bcf | grep -vFf ${haploids} > diploids.txt
      cat <(sed s"/\$/\tM/"g ${haploids}) <(sed s"/\$/\tF/"g diploids.txt) > sex.txt
      
      echo "${params.fixploidy_region1}" > fixploidy_regions.txt
      echo "${params.fixploidy_region2}" >> fixploidy_regions.txt

      # Third, resolve umbiguous genotypes:
      # 1) Set haploid genotypes (GT=0 or GT=1) in diploids (i.e. females) to missing.
      # 2) Set heterozygous genotypes (GT=0/1) in haploids (i.e. males) to missing.
      bcftools +setGT temp_nonpar.bcf -Ou -- -i "GT[@diploids.txt]='1'" -t q -n . | bcftools +setGT -Ou -- -i "GT[@diploids.txt]='0'" -t q -n c:0/0 |  bcftools +setGT -Ou -- -i "GT[@${haploids}]='het'" -t q -n . | bcftools +fixploidy -Ou -- -p fixploidy_regions.txt -s sex.txt | bcftools +fill-tags -Ou -- -t AN,AC,AF,NS,F_MISSING | bcftools view -e 'F_MISSING>0.1 || AC < 1' -Ou | bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' -Ou | bcftools filter -s PASS -Oz -o ${chrom}_nonpar.array.qc.vcf.gz
      bcftools index ${chrom}_nonpar.array.qc.vcf.gz

      # Merge nonPAR and PAR regions
      bcftools concat -a -Oz -o chrX.array.qc.vcf.gz chrX_par1.array.qc.vcf.gz chrX_nonpar.array.qc.vcf.gz chrX_par2.array.qc.vcf.gz
      bcftools index chrX.array.qc.vcf.gz
   else
      # Force ploidy to 2; then, compute AN, AC, AF, and other stats; then, remove variants with high missigness; then, set variant ID; then, set all variants to PASS.
      bcftools +fixploidy ${vcf} -Ou | bcftools +fill-tags -Ou -- -t AN,AC,AF,NS,F_MISSING | bcftools view -e 'F_MISSING>0.1 || AC < 1' -Ou | bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' -Ou | bcftools filter -s PASS -Oz -o ${chrom}.array.qc.vcf.gz
      bcftools index ${chrom}.array.qc.vcf.gz
   fi
   """
}


workflow {
   study_vcfs = Channel.fromFilePairs("${params.sample_vcf_path}", size: -1) { file -> file.getName().replaceAll("(.tbi|.csi)\$", "") }.map {it -> [it[1][0].getSimpleName(), it[1][0], it[1][1]] }
   subset_per_sample = subset_sample_wgs(study_vcfs, Channel.fromPath(params.array_positions_file).collect())
   merged_per_chromosome = merge_samples(subset_per_sample.collect())
   finalize_chromosome(merged_per_chromosome.flatten().map{ it -> [it.getSimpleName(), it]}.groupTuple().map{it -> [it[0], it[1][0], it[1][1]]}, Channel.fromPath(params.haploids_file).collect())
}
