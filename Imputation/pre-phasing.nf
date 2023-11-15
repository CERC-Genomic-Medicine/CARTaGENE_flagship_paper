#!/usr/bin/env nextflow
/*
* AUTHOR: Mohadese Sayahian Dehkordi, <mohadese.sayahiandehkordi@mail.mcgill.ca>
* VERSION: 3.0
* YEAR: 2023
*/


process get_ref_chr_names{
    cache "lenient"
    cpus 1
    memory "4GB"
    time "00:30:00"
    scratch true
    input:
    tuple val(sex_id), path(vcf), path(vcf_index)

    output:
    tuple stdout, val(sex_id), path(vcf), path(vcf_index)

    """
    chrom=`bcftools index -s ${vcf} | cut -f1`
	printf "\${chrom}"
    """
}


process get_study_chr_names{
    cache "lenient"
    cpus 1
    memory "4GB"
    time "00:30:00"
    scratch true
    input:
    tuple val(sex_id), path(vcf), path(vcf_index)

    output:
    tuple stdout, val(sex_id), path(vcf), path(vcf_index)

    """
    chrom=`bcftools index -s ${vcf} | cut -f1`
	printf "\${chrom}"
    """
}

process get_chunks {
	//executor "local"
	cache "lenient"
	cpus 1
	memory "4GB"
    time "00:30:00"

	input:
	tuple val(chromosome), val(sex_id), path(vcf), path(vcf_index)

    output:
    tuple path("*.chunk"), val(chromosome), val(sex_id), path(vcf), path(vcf_index)

	"""
    chrom=`bcftools index -s ${vcf} | cut -f1`
    start_bp=`bcftools view -HG ${vcf} | head -n1 | cut -f2`
	stop_bp=`bcftools index -s ${vcf} | cut -f2`
       
    extend=0
	for i in `seq \${start_bp} ${params.window} \${stop_bp}`; do
		if [ \${extend} -eq 0 ]; then
			chunk_start=\$((${params.flank} > i ? 0 : i - ${params.flank}))
		fi
		chunk_stop=\$((i + ${params.window} + ${params.flank}))
		n=`bcftools view -HG ${vcf} \${chrom}:\${chunk_start}-\${chunk_stop} | wc -l`
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

process eagle_phased{
    //errorStrategy 'retry'
    //maxRetries 3
    cache "lenient"
    cpus 8
    memory "128GB"
    time "48:00:00"
    scratch true
    input:
    tuple val(chromosome), val(sex_id), path(ref_vcf), path(ref_vcf_index), path(study_vcf), path(study_vcf_index) 
        
    output:
    tuple val(chromosome), val(sex_id), path("*.phased.final.vcf.gz"), path("*.phased.final.vcf.gz.tbi")
    publishDir "phased/", pattern: "*.vcf.gz*", mode: "copy"

    script:
    script:
    if(params.chromosomeX == true){
    """
    echo "${chromosome}" > chroms1.txt
    echo "chrX" > chroms2.txt

    paste chroms1.txt chroms2.txt > chr_name_conv.txt   
    
    bcftools annotate --rename-chrs chr_name_conv.txt $study_vcf -Oz -o ${study_vcf.getBaseName()}.vcf.gz
    bcftools index --tbi ${study_vcf.getBaseName()}.vcf.gz

    bcftools annotate --rename-chrs chr_name_conv.txt $ref_vcf -Oz -o ${ref_vcf.getBaseName()}.vcf.gz
    bcftools index --tbi ${ref_vcf.getBaseName()}.vcf.gz

    ${params.eagle} --vcfRef ${ref_vcf.getBaseName()}.vcf.gz  --vcfTarget ${study_vcf.getBaseName()}.vcf.gz --numThreads 4 --geneticMapFile ${params.genetic_map} --outPrefix ${study_vcf.getBaseName()}.phased --allowRefAltSwap --vcfOutFormat z
    bcftools index --tbi ${study_vcf.getBaseName()}.phased.vcf.gz
    paste chroms2.txt chroms1.txt > chr_name_conv.txt   

    bcftools annotate --rename-chrs chr_name_conv.txt ${study_vcf.getBaseName()}.phased.vcf.gz -Oz -o ${study_vcf.getBaseName()}.phased.final.vcf.gz
    bcftools index --tbi ${study_vcf.getBaseName()}.phased.final.vcf.gz

    """
    } else {
    """
    ${params.eagle} --vcfRef ${ref_vcf}  --vcfTarget ${study_vcf} --numThreads 8 --geneticMapFile ${params.genetic_map} --outPrefix ${study_vcf.getBaseName()}.phased.final --allowRefAltSwap --vcfOutFormat z
    bcftools index --tbi ${study_vcf.getBaseName()}.phased.final.vcf.gz
    """
    }
}

workflow {
        ref_ch = Channel.fromPath(params.ref_vcf_path).map{ vcf -> [ vcf.name.toString().tokenize('.').contains('female'), vcf, vcf + ".tbi" ] }
        study_ch = Channel.fromPath(params.study_vcf_path).map{ vcf -> [ vcf.name.toString().tokenize('.').contains('female'), vcf, vcf + ".tbi" ] }

        ref_vcfs = get_ref_chr_names(ref_ch)
        study_vcfs = get_study_chr_names(study_ch)
   
        phasing_ch = ref_vcfs.combine(study_vcfs, by:[0, 1])
        phased_chunks = eagle_phased(phasing_ch)
}