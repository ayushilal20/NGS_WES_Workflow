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
git clone <repository-url>
cd ngs_pipeline
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
Modify the nextflow.config to include resource specifications:

```nextflow
// Existing configuration
params {
    reads = null          
    reference = null
    outdir = "results"
    help = false
    
    // fastp parameters
    fastp_trim_front1 = 13    
    fastp_trim_tail1 = 10     
    fastp_trim_front2 = 13    
    fastp_trim_tail2 = 10
    
    // VarScan parameters
    varscan_min_coverage = 8
    varscan_min_var_freq = 0.1
    varscan_p_value = 0.05
}

// Docker configuration
docker {
    enabled = true
    runOptions = '-u $(id -u):$(id -g)'
}

// Default process configuration
process {
    executor = 'local'
    errorStrategy = { task.exitStatus in [143,137,104,134,139] ? 'retry' : 'finish' }
    maxRetries = 3

    // Default resources for all processes
    cpus = 2
    memory = 4.GB
    time = '2h'

    // Process-specific resources
    withName: 'FASTQC' {
        cpus = 2
        memory = 4.GB
        time = '1h'
    }

    withName: 'FASTP' {
        cpus = 4
        memory = 8.GB
        time = '2h'
    }

    withName: 'MINIMAP2_INDEX' {
        cpus = 2
        memory = 8.GB
        time = '2h'
    }

    withName: 'MINIMAP2_ALIGN' {
        cpus = 8
        memory = 16.GB
        time = '4h'
    }

    withName: 'SAM_TO_BAM' {
        cpus = 4
        memory = 8.GB
        time = '2h'
    }

    withName: 'SAMTOOLS_MPILEUP' {
        cpus = 4
        memory = 8.GB
        time = '2h'
    }

    withName: 'VARSCAN_CALL' {
        cpus = 4
        memory = 8.GB
        time = '2h'
    }
}

// Including process-specific configurations
includeConfig 'process.config'

// Execution profiles
profiles {
    // Standard profile
    standard {
        process {
            cpus = 2
            memory = 4.GB
        }
    }

    // High-performance profile
    highperf {
        process {
            cpus = 8
            memory = 32.GB
        }
    }

    // Low-resource profile
    lowmem {
        process {
            cpus = 1
            memory = 2.GB
        }
    }

    test {
        params.reads = "$baseDir/test/*_{1,2}.fastq.gz"
        params.reference = "$baseDir/test/reference.fa"
    }
    
    dryrun {
        process.executor = 'local'
        executor.queueSize = 1
        executor.submitRateLimit = '1/s'
        dag.enabled = true
        executor.pollInterval = '10sec'
        dryRun = true
    }
}
```

Key points about this configuration:

1. **Default Resources**: Set baseline resources for all processes
   ```nextflow
   process {
       cpus = 2
       memory = 4.GB
       time = '2h'
   }
   ```

2. **Process-Specific Resources**: Each process can have its own resource requirements
   ```nextflow
   withName: 'MINIMAP2_ALIGN' {
       cpus = 8
       memory = 16.GB
       time = '4h'
   }
   ```

3. **Resource Profiles**: Different profiles for various computing environments
   - standard: Default resources
   - highperf: High-performance computing
   - lowmem: Low-resource environments

To use these configurations, you can run the pipeline with:

```bash
# Use standard profile
nextflow run main.nf --reads "..." --reference "..."

# Use high-performance profile
nextflow run main.nf --reads "..." --reference "..." -profile highperf

# Use low-memory profile
nextflow run main.nf --reads "..." --reference "..." -profile lowmem
```

The resources are specified using these units:
- Memory: Can use GB, MB, KB (e.g., 4.GB, 512.MB)
- Time: Can use h, m, s (e.g., '2h', '30m')
- CPUs: Whole numbers

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
