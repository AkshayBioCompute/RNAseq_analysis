# WDL file: pipeline.wdl

workflow RNASeqPipeline {
    input {
        File reference_genome_gz      # Reference genome file in gzip format
        File annotation_gtf           # Annotation file in GTF format
        Array[File] fastq_files       # Array of input FASTQ files
        Int threads = 8               # Number of threads for parallel execution
    }

    # Unzip the reference genome if needed
    call UnzipReference {
        input:
            reference_genome_gz = reference_genome_gz
    }

    # Build Hisat2 index
    call Hisat2Build {
        input:
            reference_genome = UnzipReference.unzipped_reference
    }

    # Generate MD5 checksums for all FASTQ files
    call GenerateMD5 {
        input:
            fastq_files = fastq_files
    }

    # Run FastQC and MultiQC on all FASTQ files
    scatter (fastq in fastq_files) {
        call FastQC {
            input:
                fastq_file = fastq
        }
    }

    call MultiQC {
        input:
            fastqc_reports = FastQC.fastqc_outputs
    }

    # Trim paired-end reads with Fastp
    scatter (fastq_pair in paired_fastqs(fastq_files)) {
        call Fastp {
            input:
                fastq_r1 = fastq_pair.r1,
                fastq_r2 = fastq_pair.r2,
                threads = threads
        }
    }

    # Align reads with Hisat2
    scatter (trimmed_pair in Fastp.trimmed_pairs) {
        call Hisat2Align {
            input:
                r1 = trimmed_pair.r1,
                r2 = trimmed_pair.r2,
                index_prefix = Hisat2Build.index_prefix
        }
    }

    # Convert SAM to BAM, sort, and generate statistics
    scatter (sam_file in Hisat2Align.sam_outputs) {
        call SAMtoSortedBAM {
            input:
                sam_file = sam_file,
                threads = threads
        }
    }

    # Run featureCounts on all sorted BAM files
    call FeatureCounts {
        input:
            sorted_bams = SAMtoSortedBAM.sorted_bams,
            annotation_gtf = annotation_gtf,
            threads = threads
    }

    output {
        Array[File] md5_checksums = GenerateMD5.md5_outputs
        Array[File] fastqc_reports = FastQC.fastqc_outputs
        File multiqc_report = MultiQC.multiqc_report
        Array[File] trimmed_fastqs = Fastp.trimmed_fastqs
        Array[File] sorted_bams = SAMtoSortedBAM.sorted_bams
        File gene_counts = FeatureCounts.gene_counts
    }
}

# Task to unzip the reference genome
task UnzipReference {
    input {
        File reference_genome_gz
    }

    command {
        gunzip -c ~{reference_genome_gz} > reference_genome.fa
    }

    output {
        File unzipped_reference = "reference_genome.fa"
    }
}

# Task to build Hisat2 index
task Hisat2Build {
    input {
        File reference_genome
    }

    command {
        hisat2-build ~{reference_genome} index_prefix
    }

    output {
        Array[File] index_prefix = glob("index_prefix.*.ht2")
    }
}

# Task to generate MD5 checksums
task GenerateMD5 {
    input {
        Array[File] fastq_files
    }

    command {
        md5sum ~{sep=' ' fastq_files} > md5_checksums.txt
    }

    output {
        File md5_outputs = "md5_checksums.txt"
    }
}

# Task to run FastQC
task FastQC {
    input {
        File fastq_file
    }

    command {
        fastqc -o . ~{fastq_file}
    }

    output {
        Array[File] fastqc_outputs = glob("*_fastqc*")
    }
}

# Task to run MultiQC
task MultiQC {
    input {
        Array[File] fastqc_reports
    }

    command {
        multiqc . -o .
    }

    output {
        File multiqc_report = "multiqc_report.html"
    }
}

# Task to trim paired-end reads using Fastp
task Fastp {
    input {
        File fastq_r1
        File fastq_r2
        Int threads
    }

    command {
        fastp -i ~{fastq_r1} -I ~{fastq_r2} -o trimmed_~{basename(fastq_r1)} -O trimmed_~{basename(fastq_r2)} --thread ~{threads} -h fastp.html -j fastp.json
    }

    output {
        File r1_trimmed = "trimmed_~{basename(fastq_r1)}"
        File r2_trimmed = "trimmed_~{basename(fastq_r2)}"
        File html_report = "fastp.html"
        File json_report = "fastp.json"
    }
}

# Task to align reads with Hisat2
task Hisat2Align {
    input {
        File r1
        File r2
        Array[File] index_prefix
    }

    command {
        hisat2 -x index_prefix -1 ~{r1} -2 ~{r2} -S aligned.sam
    }

    output {
        File sam_output = "aligned.sam"
    }
}

# Task to convert SAM to sorted BAM and generate stats
task SAMtoSortedBAM {
    input {
        File sam_file
        Int threads
    }

    command {
        samtools view -bS ~{sam_file} | samtools sort -o sorted.bam
        samtools flagstat sorted.bam > alignment_statistics.txt
    }

    output {
        File sorted_bam = "sorted.bam"
        File alignment_stats = "alignment_statistics.txt"
    }
}

# Task to count features with featureCounts
task FeatureCounts {
    input {
        Array[File] sorted_bams
        File annotation_gtf
        Int threads
    }

    command {
        featureCounts -T ~{threads} -a ~{annotation_gtf} -o gene_counts.txt -g gene_id -t exon -s 1 -p ~{sep=' ' sorted_bams}
    }

    output {
        File gene_counts = "gene_counts.txt"
    }
}
