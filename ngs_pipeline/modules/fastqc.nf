process FASTQC {
    tag "${sample_id}"
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "*.html", emit: html_reports
    path "*.zip", emit: zip_reports
    
    script:
    """
    fastqc -t ${task.cpus} ${reads}
    """
}
