process MINIMAP2_INDEX {
    tag "${reference.baseName}"
    publishDir "${params.outdir}/reference", mode: 'copy'
    
    input:
    path reference
    
    output:
    tuple path(reference), path("${reference}.mmi"), emit: minimap_index
    
    script:
    """
    minimap2 -d ${reference}.mmi ${reference}
    """
}

process MINIMAP2_ALIGN {
    tag "${sample_id}"
    publishDir "${params.outdir}/aligned", mode: 'copy', pattern: "*.sam"
    
    input:
    tuple val(sample_id), path(reads)
    tuple path(reference), path(index)
    
    output:
    tuple val(sample_id), path("${sample_id}.sam"), emit: sam_file
    
    script:
    if (reads.size() == 2) {
        """
        # Align paired-end reads
        minimap2 -ax sr -t ${task.cpus} ${index} ${reads[0]} ${reads[1]} > ${sample_id}.sam
        """
    } else {
        """
        # Align single-end reads
        minimap2 -ax sr -t ${task.cpus} ${index} ${reads[0]} > ${sample_id}.sam
        """
    }
}

process SAM_TO_BAM {
    tag "${sample_id}"
    publishDir "${params.outdir}/aligned", mode: 'copy'
    
    input:
    tuple val(sample_id), path(sam_file)
    
    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: sorted_bam
    
    script:
    """
    samtools view -bS ${sam_file} | samtools sort -o ${sample_id}.sorted.bam
    samtools index ${sample_id}.sorted.bam
    """
}
