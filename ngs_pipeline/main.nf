#!/usr/bin/env nextflow
// Import modules
include { FASTQC } from './modules/fastqc'
include { FASTP } from './modules/fastp'
include { MINIMAP2_INDEX; MINIMAP2_ALIGN; SAM_TO_BAM } from './modules/alignment'
include { SAMTOOLS_MPILEUP; VARSCAN_CALL } from './modules/variant_calling'

params.help = false

if (params.help) {
    log.info"""
    =========================================
    NGS Analysis Pipeline
    =========================================
    Usage:
    nextflow run main.nf --reads '/path/to/reads/*_{1,2}.fastq.gz' --reference reference.fasta
    
    Mandatory arguments:
      --reads           Path to FASTQ files (must be quoted and support glob patterns)
                        Example: '/path/to/reads/*_{1,2}.fastq.gz'
      --reference       Path to reference genome FASTA file
    
    Optional arguments:
      --outdir          Output directory (default: results)
      --trim_minlen     Minimum read length after trimming (default: 50)
      --trim_quality    Minimum base quality during trimming (default: 20)
    """.stripIndent()
    exit 0
}

// Parameter validation
if (!params.reads) {
    error "Please provide path to FASTQ files with --reads"
}
if (!params.reference) {
    error "Please provide a reference genome with --reference"
}

// Define input channels with better handling for mixed file types
reads_ch = Channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .map { sample_id, files -> 
        println "Found sample: ${sample_id} with files: ${files}"
        return tuple(sample_id, files)
    }
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }

reference_ch = Channel
    .fromPath(params.reference, checkIfExists: true)
    .ifEmpty { error "Cannot find reference file: ${params.reference}" }

// Define workflow
workflow {
    // QC Steps
    FASTQC(reads_ch)
    FASTP(reads_ch)
    
    // Reference preparation and alignment
    MINIMAP2_INDEX(reference_ch)
    MINIMAP2_ALIGN(FASTP.out.trimmed_reads, MINIMAP2_INDEX.out.minimap_index)
    SAM_TO_BAM(MINIMAP2_ALIGN.out.sam_file)
    
    // Variant calling with samtools/varscan
    SAMTOOLS_MPILEUP(SAM_TO_BAM.out.sorted_bam, reference_ch)
    VARSCAN_CALL(SAMTOOLS_MPILEUP.out.mpileup)
}

// Completion handler
workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Time taken: $workflow.duration"
    log.info "Execution status: ${ workflow.success ? 'SUCCESS' : 'FAILED' }"
}
