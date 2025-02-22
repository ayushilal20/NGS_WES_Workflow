process FASTQ_DUMP {
    tag "${sra_file.baseName}"
    publishDir "${params.outdir}/fastq", mode: 'copy'
    
    input:
    path sra_file
    
    output:
    tuple val("${sra_file.simpleName}"), path("*_R{1,2}.fastq.gz"), emit: reads
    
    script:
    """
    fastq-dump --split-files --gzip ${sra_file}
    mv ${sra_file.simpleName}_1.fastq.gz ${sra_file.simpleName}_R1.fastq.gz
    mv ${sra_file.simpleName}_2.fastq.gz ${sra_file.simpleName}_R2.fastq.gz
    """
}
