# Nextflow Pipeline for Whole Exome Sequencing (WES) Analysis

## Overview
This Nextflow pipeline performs NGS data analysis including quality control, read trimming, alignment, and variant calling. The pipeline uses Docker containers for reproducibility and supports both single-end and paired-end sequencing data.

## Prerequisites
- Nextflow (version 22.10.0 or later)
- Docker
- At least 16GB RAM (recommended)
- Multi-core processor (recommended: 4+ cores)
- Sufficient storage space for NGS data processing

## Installation

1. Clone the repository:
```bashOptional arguments:

    --outdir: Output directory (default: results)
    --trim_minlen: Minimum read length after trimming (default: 50)
    --trim_quality: Minimum base quality during trimming (default: 20
git clone <repository-url>
cd ngs-pipeline
```

2. Ensure Docker is installed and running on your system

## Pipeline Steps

1. Quality Control
   - Runs FastQC for initial read quality assessment
   - Uses FastP for read trimming and filtering

2. Read Alignment
   - Indexes the reference genome with Minimap2
   - Aligns reads to the reference genome

3. BAM Processing
   - Converts SAM to BAM and sorts the alignments
   - Indexes BAM files for variant calling

4. Variant Calling
   - Uses Samtools to generate mpileup files
   - Calls variants using VarScan

## Resource Configuration

The pipeline resources can be configured in `nextflow.config`. Adjust these settings based on your system specifications:

### Basic Resource Configuration
```nextflow
process {
    cpus = 4          // Number of CPU cores per process
    memory = '8 GB'   // Memory per process
    time = '2h'       // Time limit per process
}
```

### Process-Specific Resources
For optimal performance, configure resources for specific processes:

```nextflow
process {
    withName: 'MINIMAP2_ALIGN' {
        cpus = 8
        memory = '16 GB'
    }
    withName: 'VARSCAN_CALL' {
        cpus = 4
        memory = '8 GB'
    }
}
```

## Usage

### Basic Command
```bash
nextflow run main.nf \
    --reads "/path/to/reads/*_{1,2}.fastq.gz" \
    --reference /path/to/reference.fa \
    --outdir results
```

### Available Parameters
- `--reads`: Path to input FASTQ files (required)
- `--reference`: Path to reference genome (required)
- `--outdir`: Output directory (default: results)
- `--trim_minlen`: Minimum read length after trimming (default: 50)
- `--trim_quality`: Minimum base quality during trimming (default: 20

### Execution Profiles
- Default: Standard execution
- Dryrun: Validate pipeline without execution

Example with profile:
```bash
nextflow run main.nf \
    --reads "/path/to/reads/*_{1,2}.fastq.gz" \
    --reference /path/to/reference.fa \
    -profile dryrun
```

## Output Structure
```
results/
├── fastqc/           # Quality control reports
├── trimmed/          # Trimmed reads
├── aligned/          # Alignment files
│   ├── *.sam
│   ├── *.bam
│   └── *.bai
├── variants/         # Called variants
│   └── *.vcf
└── reference/        # Indexed reference
```

## Performance Optimization

### Memory Usage
- For whole genome sequencing: Recommend 32GB+ RAM
- For exome/targeted sequencing: Minimum 16GB RAM
- Adjust process-specific memory in `nextflow.config`

### Storage Requirements
- Raw data: ~10GB per sample
- Intermediate files: 2-3x raw data size
- Final results: ~5GB per sample

### CPU Utilization
- Recommended: 8+ cores for optimal performance
- Minimum: 4 cores
- Configure process-specific CPU usage in config file

## Troubleshooting

### Common Issues
1. Insufficient memory
   - Solution: Adjust memory settings in config
   - Example: `memory = '32 GB'`

2. Docker permission errors
   - Solution: Add user to docker group
   - Command: `sudo usermod -aG docker $USER`

3. Pipeline execution errors
   - Check error logs in `.nextflow.log`
   - Use `-resume` flag to restart from last successful step

## Support
For issues and questions, please open an issue in the repository's issue tracker.
