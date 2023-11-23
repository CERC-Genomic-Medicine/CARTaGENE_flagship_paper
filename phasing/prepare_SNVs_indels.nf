#!/usr/bin/env nextflow
/*
* AUTHOR: Mohadese Sayahian Dehkordi <mohadese.sayahiandehkordi@mail.mcgill.ca>, Daniel Taliun <daniel.taliun@mcgill.ca>
* VERSION: 2.0
* YEAR: 2023
*/

// Pre-requisites/assumptions:
// 1. Script was designed for autosomal chromosomes (i.e. 1-22) and chromosome X
// 2. Upstream variant caller represented males as haploids (i.e. 0 or 1) in non-PAR regions of chromosome X
// 3. Multi-allelic variants were split into multiple VCF entries (one entry per alternate allele)
// 4. Make sure that bcftools is installed/loaded on your system
// 5. Compile/install vt (https://github.com/atks/vt)

// How to execute:
// nextflow run prepare_SNVs_indels.nf --vcf_path "/home/user/data/chr*.vcf.gz" --ref "/home/user/reference/GRCh38.fa" --vt "/home/user/tools/vt/vt"

params.vcf_path = "/path/to/data/*.vcf.gz" // Absolute path to the Input VCF/BCF files split by chromosome
params.ref = "/path/to/reference/Homo_sapiens.GRCh38.fa" // Absolute path to the FASTA file with the human genome reference
params.vt = "/path/to/executable/vt" // Absolute path to the vt executable

// PAR region coordinates for GRCh38. Change only when working with different human genome reference build.
params.par1_region = "chrX:10001-2781479"
params.par2_region = "chrX:155701383-156030895"

process clean_VCF {
    errorStrategy 'retry'
    maxRetries 3
    cache "lenient"
    
    executor "slurm"
    clusterOptions "--account=rrg-vmooser"

    cpus 1
    memory "16GB"
    time "8h"
    scratch '$SLURM_TMPDIR'

    input:
    tuple path(vcf), path(vcf_index)
    
    output:
    tuple path("*.norm.dedup.vcf.gz"), path("*.norm.dedup.vcf.gz.tbi")
    
    publishDir "SNVs_indels/", pattern: "*.vcf.gz*", mode: "copy"

    """
    chr=`bcftools index -s ${vcf} | cut -f1`
    if [ "\$chr" = "chrX" ]; then
        # Process the PAR regions (where males and females are diploid) in a similar way as autosomal chromosomes below.
        bcftools view ${vcf} ${params.par1_region} -Ou | bcftools +setGT -Ou -- -t q -n . -i 'FT!="PASS"' | bcftools +fixploidy -Ou | bcftools annotate -x INFO,^FORMAT/GT -Ou | bcftools +fill-tags -Ou -- -t AN,AC,AF,NS,F_MISSING | bcftools view -e 'F_MISSING>0.1 || AC < 1' -Ou | ${params.vt} normalize - -r ${params.ref} -o + | ${params.vt} uniq + -o + | bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' -Oz -o chrX_PAR1.norm.dedup.vcf.gz
	bcftools index -t chrX_PAR1.norm.dedup.vcf.gz
        
	bcftools view ${vcf} ${params.par2_region} -Ou | bcftools +setGT -Ou -- -t q -n . -i 'FT!="PASS"' | bcftools +fixploidy -Ou | bcftools annotate -x INFO,^FORMAT/GT -Ou | bcftools +fill-tags -Ou -- -t AN,AC,AF,NS,F_MISSING | bcftools view -e 'F_MISSING>0.1 || AC < 1' -Ou | ${params.vt} normalize - -r ${params.ref} -o + | ${params.vt} uniq + -o + | bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' -Oz -o chrX_PAR2.norm.dedup.vcf.gz
	bcftools index -t chrX_PAR2.norm.dedup.vcf.gz

        # Process non-PAR part of chromosome X. Sample must be either haploid or diploid across all non-PAR variants.
	# First, remove PAR regions:
	bcftools view -t ^${params.par1_region},${params.par2_region} ${vcf} -Ob -o temp.bcf
	bcftools index temp.bcf
	# Second, check that ploidy is the same across all variants for each individual
	# The bcftools +check-ploidy will output multiple rows per individual if the ploidy was different. If ploidy is the same, then it will output only one row per individual.
	bcftools +check-ploidy temp.bcf | tail -n+2 | cut -f1 > ploidy.txt
	bcftools query -l temp.bcf > all_sample_names.txt
	if ! cmp -s ploidy.txt all_sample_names.txt; then
            exit 1 # there are individuals with different ploidy at different variants
	fi
	# Now we can apply filters similar to the autosomal chromosomes
	bcftools +setGT -Ou temp.bcf -- -t q -n . -i 'FT!="PASS"' | bcftools annotate -x INFO,^FORMAT/GT -Ou | bcftools +fill-tags -Ou -- -t AN,AC,AF,NS,F_MISSING | bcftools view -e 'F_MISSING>0.1 || AC < 1' -Ou | ${params.vt} normalize - -r ${params.ref} -o + | ${params.vt} uniq + -o + | bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' -Oz -o chrX_nonPAR.norm.dedup.vcf.gz
	bcftools index -t chrX_nonPAR.norm.dedup.vcf.gz
    else
        # a. set GT to missing if it didn't pass QC
        # b. change `.` to `./.` in autosomal chromosomes (`.` is an artifact from variant calling pipeline). When no options are given to bcftools fixploidy, then all samples are treated as females and, thus, all samples will be set to diploids on chromosome X.
        # c. remove all INFO and FORMAT fields which will be not needed in downstream analyses
        # d. recalculate allele counts after setting some GT fields to missing
        # e. remove any monomorphic variants which were introduced in previous steps or variants with missingness >0.1
        # f. left-align indels
        # g. remove any duplicated records after left-alignment (see vt documentation)
        # h. update ID field to reflect new reference and alternate alleles after left-alignment
    
        bcftools +setGT -Ou ${vcf}  -- -t q -n . -i 'FT!="PASS"' | bcftools +fixploidy -Ou | bcftools annotate -x INFO,^FORMAT/GT -Ou | bcftools +fill-tags -Ou -- -t AN,AC,AF,NS,F_MISSING | bcftools view -e 'F_MISSING>0.1 || AC < 1' -Ou | ${params.vt} normalize - -r ${params.ref} -o + | ${params.vt} uniq + -o + | bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' -Oz -o \${chr}.norm.dedup.vcf.gz
        bcftools tabix --tbi \${chr}.norm.dedup.vcf.gz
    fi
    """
}

workflow {
	clean_VCF(Channel.fromPath(params.vcf_path).map { vcf -> [vcf, vcf + ".tbi" ] })
}
