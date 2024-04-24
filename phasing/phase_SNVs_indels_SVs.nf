#!/usr/bin/env nextflow
/*
* AUTHOR: Mohadese Sayahian Dehkordi <mohadese.sayahiandehkordi@mail.mcgill.ca>, Daniel Taliun <daniel.taliun@mcgill.ca>
* VERSION: 1.0
* YEAR: 2024
*/

// This pipeline performs reference-based phasing using shapeit5.1.1.
// The pipeline was designed for small sequencing studies (i.e. study N<10,000, reference N<10,000) and, thus, chromosomes are phased in parallel but without additional chunking (with exception to rare variants).
// How to run:
// nextflow run phase_SNVs_indels_SVs.nf --vcf_path "/path/to/chr*.vcf.gz*" --ref_vcf_path "/path/to/chr*_.vcf.gz*" --males_list "/path/to/males.txt"  --shapeit5_common /path/to/shapeit/phase_common_static --shapeit5_rare /path/to/shapeit/phase_rare_static --shapeit5_genetic_map "/path/to/shapeit/genetic_map/chr*.b38.gmap.gz" --shapeit5_chunks "/path/to/chunks/b38/4cM/chunks_chr*.txt"  -with-report report.html

params.vcf_path = "/path/to/chr*.vcf.gz*" // Absolute path to the input study VCF/BCF file(s) with the corresponding TBI/CSI index file(s). One file per chromosome is expected. Chromosome X splitted into PAR1, PAR2, nonPAR.
params.ref_vcf_path = "/path/to/chr*_.vcf.gz*" // Absolute path to the reference VCF/BCF file(s) with the TBI/CSI index file(s). One file per chromosome is expected.
params.males_list = "/path/to/file.txt" // Absolute path to the file listing all males (i.e. haploid samples for chromosome X non-PAR region). One sample ID per line.

params.shapeit5_common = "/path/to/executable/phase_common_static" // Absolute path to the Shapeit5 executable for common variants phasing
params.shapeit5_rare = "/path/to/executable/phase_rare_static" // Absolute path to the Shapeit5 executable for rare variants phasing
params.shapeit5_genetic_map = "/path/to/genetic_map/chr*.b38.gmap.gz" // Absolute path to the genetic map from Shapeit5 software
params.shapeit5_chunks = "/path/to/chunks/b38/4cM/chunks_chr*.txt" // Absolute path to the chunks file from Shapeit5 software


process shapeit5_phasing_with_reference {
    errorStrategy 'finish'
    cache "lenient"

    executor 'slurm'
    clusterOptions '--account=def-vmooser'
    scratch '$SLURM_TMPDIR'

    cpus 8
    memory "64GB"
    time "48h"


    input:
    tuple val(chr), val(prefix), path(study_vcf), path(study_vcf_index), path(ref_vcf), path(ref_vcf_index)
    path genetic_map_files
    path chunks_files

    output:
    path("*.phased_using_shapeit5_with_1KGP_HGDP.bcf*")

    publishDir "phased", mode: "copy"

    """
    # The INFO/END field for SVs will be dropped by the Shapeit5.1.1 from the output bcf. Without the INFO/END field, the next step, phase_rare, will fail because underlying code will not be able to determine which SVs were already phased.
    # Thus, we will subset the SVs' INFO annotations for later use.
    bcftools query -i 'INFO/END!="."' -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%END\t%SVTYPE\t%SVLEN" ${study_vcf} | bgzip -c > SVs_annotations.tsv.gz
    tabix -s1 -b2 -e2 SVs_annotations.tsv.gz
    bcftools view -h ${study_vcf} | grep -E "ID=(END,|SVTYPE,|SVLEN,)|##ALT=" > SVs_annotations_header.txt

    ${params.shapeit5_common} --thread 8 --input ${study_vcf} --reference ${ref_vcf} --map ${prefix}.b38.gmap.gz --region ${chr} --filter-maf 0.001 --output common.scaffold_with_ref.bcf
    ${params.shapeit5_common} --thread 8 --input ${study_vcf} --scaffold common.scaffold_with_ref.bcf --map ${prefix}.b38.gmap.gz --region ${chr} --filter-maf 0.001 --output common.phased.bcf

    # Adding INFO/END and other related INFO fields back.    
    bcftools annotate -a SVs_annotations.tsv.gz -h SVs_annotations_header.txt -c CHROM,POS,~ID,REF,ALT,+INFO/END,+INFO/SVTYPE,+INFO/SVLEN common.phased.bcf -Ob -o common.phased_with_SVs_annotations.bcf
    bcftools index common.phased_with_SVs_annotations.bcf

    while read LINE; do
        ID=\$(echo \$LINE | awk '{ print \$1; }')
        SRG=\$(echo \$LINE | awk '{ print \$3; }')
        IRG=\$(echo \$LINE | awk '{ print \$4; }')
        
	n=`bcftools view -HG  common.phased_with_SVs_annotations.bcf \${IRG} | wc -l`
	if [[ \${n} -eq 0 ]]; then
	    echo "skipping chunk \${ID}, \${IRG}"
	    continue # nothing to phase. This should mainly happen when working with chromosome X PAR regions, for chunks which are outside them.
	fi

	${params.shapeit5_rare} --thread 8 --input ${study_vcf} --scaffold common.phased_with_SVs_annotations.bcf --map ${chr}.b38.gmap.gz --input-region \${IRG} --scaffold-region \${SRG} --output chunk\${ID}.phased.bcf
    done < chunks_${chr}.txt

    ls -1v chunk*.phased.bcf > files.txt

    bcftools concat --naive -f files.txt -o all_chunks.phased.bcf --threads 8
    bcftools index all_chunks.phased.bcf

    # Adding INFO/END and other SVs related INFO fields back.
    bcftools annotate -a SVs_annotations.tsv.gz -h SVs_annotations_header.txt -c CHROM,POS,~ID,REF,ALT,+INFO/END,+INFO/SVTYPE,+INFO/SVLEN all_chunks.phased.bcf -Ob -o ${prefix}.phased_using_shapeit5_with_1KGP_HGDP.bcf
    bcftools index ${prefix}.phased_using_shapeit5_with_1KGP_HGDP.bcf
    """
}


process shapeit5_Xnonpar_phasing_with_reference {
    errorStrategy 'finish'
    cache "lenient"

    executor 'slurm'
    clusterOptions '--account=def-vmooser'
    scratch '$SLURM_TMPDIR'

    cpus 8
    memory "64GB"
    time "48h"

    input:
    tuple val(chr), val(prefix), path(study_vcf), path(study_vcf_index), path(ref_vcf), path(ref_vcf_index)
    path males_ids_file
    path genetic_map_files
    path chunks_files

    output:
    path("*.phased_using_shapeit5_with_1KGP_HGDP.bcf*")

    publishDir "phased", mode: "copy"

    """
    # Shapeit5.1.1 didn’t work properly when males were encoded as haploids (i.e. 0 and 1). Let’s recode them to diploids (i.e. 0/0 and 1/1):
    # 1. For all males change GT=1 to GT=1/1
    # 2. For all males change GT=0 to GT=0/0
    # 3. After first 2 steps, only missing '.' are left and we can safely run +fixploidy command forcing diploid GT.
    # 4. After these transformation, the AC, AN, and AF values may change and must be recomputed.
    bcftools +setGT ${study_vcf} -Ou -- -i 'GT[@${males_ids_file}]="1"' -t q -n c:"1/1" | bcftools +setGT -Ou -- -i 'GT[@${males_ids_file}]="0"' -t q -n c:"0/0" | bcftools +fixploidy -Ou | bcftools annotate -x ^INFO/END,INFO/SVTYPE,INFO/SVLEN -Ou | bcftools +fill-tags -Ob -o unphased.males2diploids.bcf -- -t AN,AC,AF,NS,F_MISSING

    # The INFO/END field for SVs will be dropped by the Shapeit5.1.1 from the output bcf. Without the INFO/END field, the next step, phase_rare, will fail because underlying code will not be able to determine which SVs were already phased.
    # Thus, we will subset the SVs' INFO annotations for later use.
    bcftools query -i 'INFO/END!="."' -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%END\t%SVTYPE\t%SVLEN" unphased.males2diploids.bcf | bgzip -c > SVs_annotations.tsv.gz
    tabix -s1 -b2 -e2 SVs_annotations.tsv.gz
    bcftools view -h unphased.males2diploids.bcf | grep -E "ID=(END,|SVTYPE,|SVLEN,)|##ALT=" > SVs_annotations_header.txt

    ${params.shapeit5_common} --thread 8 --input unphased.males2diploids.bcf --haploids ${males_ids_file} --reference ${ref_vcf} --map ${chr}.b38.gmap.gz --region ${chr} --filter-maf 0.001 --output common.scaffold_with_ref.bcf
    ${params.shapeit5_common} --thread 8 --input unphased.males2diploids.bcf --haploids ${males_ids_file} --scaffold common.scaffold_with_ref.bcf --map ${chr}.b38.gmap.gz --region ${chr} --filter-maf 0.001 --output common.phased.bcf

    # Adding INFO/END and other related INFO fields back.
    bcftools annotate -a SVs_annotations.tsv.gz -h SVs_annotations_header.txt -c CHROM,POS,~ID,REF,ALT,+INFO/END,+INFO/SVTYPE,+INFO/SVLEN common.phased.bcf -Ob -o common.phased_with_SVs_annotations.bcf
    bcftools index common.phased_with_SVs_annotations.bcf

    while read LINE; do
        ID=\$(echo \$LINE | awk '{ print \$1; }')
        SRG=\$(echo \$LINE | awk '{ print \$3; }')
        IRG=\$(echo \$LINE | awk '{ print \$4; }')

        n=`bcftools view -HG  common.phased_with_SVs_annotations.bcf \${IRG} | wc -l`
        if [[ \${n} -eq 0 ]]; then
            echo "skipping chunk \${ID}, \${IRG}"
            continue # nothing to phase. This should mainly happen if chunk completely falls into the chromosome X PAR regions.
        fi

        ${params.shapeit5_rare} --thread 8 --input unphased.males2diploids.bcf --haploids ${males_ids_file} --scaffold common.phased_with_SVs_annotations.bcf --map ${chr}.b38.gmap.gz --input-region \${IRG} --scaffold-region \${SRG} --output chunk\${ID}.phased.bcf
    done < chunks_${chr}.txt

    ls -1v chunk*.phased.bcf > files.txt

    bcftools concat --naive -f files.txt -o all_chunks.phased.bcf --threads 8
    bcftools index all_chunks.phased.bcf

    # Check if there are any phased heterozygous genotypes for males.
    bcftools query -S ${males_ids_file} -i 'GT="het"' -f '[%CHROM:%POS %SAMPLE %GT\n]' all_chunks.phased.bcf > males_hets.txt
    n=`wc -l males_hets.txt`
    if [[ \${n} -ne 0 ]]; then
        echo "There are \${n} phased heterozygous genotypes for males. Check your input data!"
        exit 1 # there are phased heterozygous genotypes for males
    fi

    # Recode males from 0/0s and 1/1s back to 0s and 1s to comply with the VCF spec.
    bcftools +setGT all_chunks.phased.bcf -Ou -- -i 'GT[@${males_ids_file}]="1|1"' -t q  -n c:"1" | bcftools +setGT -Ou -- -i 'GT[@${males_ids_file}]="0|0"' -t q  -n c:"0" | bcftools +fill-tags -Ob -o all_chunks.phased.males_recoded.bcf  -- -t AN,AC
    bcftools index all_chunks.phased.males_recoded.bcf

    # Adding INFO/END and other SVs related fields back.
    bcftools annotate -a SVs_annotations.tsv.gz -h SVs_annotations_header.txt -c CHROM,POS,~ID,REF,ALT,+INFO/END,+INFO/SVTYPE,+INFO/SVLEN all_chunks.phased.males_recoded.bcf -Ob -o ${prefix}.phased_using_shapeit5_with_1KGP_HGDP.bcf
    bcftools index ${prefix}.phased_using_shapeit5_with_1KGP_HGDP.bcf
    """
}


workflow {
    // Here, we assume that the chromosome name is encoded in the file prefix e.g. chr1.study_name.vcf.gz, chr2.study_name.vcf.gz, ..., chrX_nonPAR.study_name.vcf.gz, chrX_PAR1.study_name.vcf.gz, chrX_PAR2.study_name.vcf.gz
    study_vcfs = Channel.fromFilePairs("${params.vcf_path}", size: -1) { file -> file.getName().replaceAll("(.tbi|.csi)\$", "") }.map {it -> [ (it[0] =~ /chrX{1}|chr[1-9][0-9]?/)[0], (it[0] =~ /chrX{1}[.]{1}|chrX_nonPAR|chrX_PAR1|chrX_PAR2|chr[1-9][0-9]?/)[0].replace("PAR", "par"), it[1][0], it[1][1]] }
    //study_vcfs.view()

    // Here, we assume that the chromosome name is encoded in the file in the form of chr1, chr2, ..., chrX.
    ref_vcfs = Channel.fromFilePairs("${params.ref_vcf_path}", size: -1) { file -> file.getName().replaceAll("(.tbi|.csi)\$", "") }.map { it -> [ (it[0] =~ /chrX{1}|chr[1-9][0-9]?/)[0], it[1][0], it[1][1] ] }
    //ref_vcfs.view()


    //study_vcfs.combine(ref_vcfs, by: 0).view()

    input_branched = study_vcfs.combine(ref_vcfs, by: 0).branch { 
    	nonpar: it[1] == "chrX_nonpar"
    	auto: true
    }

    // Phase autosomal chromosomex and chromosome X PAR1 and PAR2 pseudo-autosomal regions. 
    shapeit5_phasing_with_reference(input_branched.auto, Channel.fromPath(params.shapeit5_genetic_map).collect(), Channel.fromPath(params.shapeit5_chunks).collect())

    // Phase X non-PAR region
    shapeit5_Xnonpar_phasing_with_reference(input_branched.nonpar, Channel.fromPath(params.males_list).collect(), Channel.fromPath(params.shapeit5_genetic_map).collect(), Channel.fromPath(params.shapeit5_chunks).collect())
}
