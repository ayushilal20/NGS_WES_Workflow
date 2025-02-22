process SAMTOOLS_MPILEUP {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    
    output:
    tuple val(sample_id), path("${sample_id}.mpileup"), emit: mpileup
    
    script:
    """
    # Create reference index if it doesn't exist
    if [ ! -f ${reference}.fai ]; then
        samtools faidx ${reference}
    fi
    
    # Generate mpileup
    samtools mpileup -f ${reference} -o ${sample_id}.mpileup ${bam}
    """
}

process VARSCAN_CALL {
    tag "${sample_id}"
    publishDir "${params.outdir}/variants", mode: 'copy'
    
    input:
    tuple val(sample_id), path(mpileup)
    
    output:
    tuple val(sample_id), path("${sample_id}.vcf"), emit: vcf
    
    script:
    """
    # Call variants with VarScan
    varscan mpileup2snp ${mpileup} \
        --min-coverage ${params.varscan_min_coverage} \
        --min-var-freq ${params.varscan_min_var_freq} \
        --p-value ${params.varscan_p_value} \
        --output-vcf 1 > ${sample_id}.vcf
    """
}
