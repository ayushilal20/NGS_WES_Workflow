process FASTP {
    tag "${sample_id}"
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_trimmed_*.fastq.gz"), emit: trimmed_reads
    path "${sample_id}_fastp.html", emit: html_report
    path "${sample_id}_fastp.json", emit: json_report
    
    script:
    // Handle both paired-end and single-end data
    if (reads.size() == 2) {
        """
        fastp \
            -i ${reads[0]} \
            -I ${reads[1]} \
            -o ${sample_id}_trimmed_1.fastq.gz \
            -O ${sample_id}_trimmed_2.fastq.gz \
            -f ${params.fastp_trim_front1} \
            -t ${params.fastp_trim_tail1} \
            -h ${sample_id}_fastp.html \
            -j ${sample_id}_fastp.json \
            --thread ${task.cpus}
        """
    } else {
        """
        fastp \
            -i ${reads[0]} \
            -o ${sample_id}_trimmed_1.fastq.gz \
            -f ${params.fastp_trim_front1} \
            -t ${params.fastp_trim_tail1} \
            -h ${sample_id}_fastp.html \
            -j ${sample_id}_fastp.json \
            --thread ${task.cpus}
        """
    }
}
