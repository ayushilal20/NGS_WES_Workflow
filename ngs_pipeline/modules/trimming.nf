process FASTP {
    tag "${sample_id}"
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads1), path(reads2)
    
    output:
    tuple val(sample_id), path("${params.outdir}/trimmed/${sample_id}_1_trimmed.fastq.gz"), path("${params.outdir}/trimmed/${sample_id}_2_trimmed.fastq.gz"), emit: trimmed_reads
    
    script:
    """
    fastp \
        -i ${reads1} \
        -I ${reads2} \
        -o ${params.outdir}/trimmed/${sample_id}_1_trimmed.fastq.gz \
        -O ${params.outdir}/trimmed/${sample_id}_2_trimmed.fastq.gz \
        -f 13 \
        -t 10
    """
}

