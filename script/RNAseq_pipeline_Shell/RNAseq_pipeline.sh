#!/bin/bash

# Navigate to the directory containing the fastq.gz files
cd /path/to/your/working/directory  # Change to your working directory

# Generate MD5 checksums for all FASTQ files
md5sum *.fastq.gz > md5_checksums.txt

# Unzip all FASTQ files
gunzip *.fastq.gz

# Create directories for outputs
mkdir -p QC_reports

# Run FastQC on all FASTQ files
fastqc -o QC_reports/ *.fastq

# Generate MultiQC summary report
multiqc -o QC_reports/ QC_reports/

# Trim paired-end reads using Fastp (loops over all *_R1.fastq files)
for r1 in *_R1.fastq; do
    r2=${r1/_R1.fastq/_R2.fastq}  # Find the matching R2 file
    fastp -i "$r1" -I "$r2" -o "trimmed_${r1}" -O "trimmed_${r2}" \
    -h "QC_reports/fastp_${r1%%_R1.fastq}.html" -j "QC_reports/fastp_${r1%%_R1.fastq}.json" --thread 8
done

# Unzip the reference genome if not already unzipped
gunzip reference_genome.fa.gz

# Build the Hisat2 index
hisat2-build reference_genome.fa index_prefix

# Align reads with Hisat2 for all paired-end fastq files
for r1 in trimmed_*_R1.fastq; do
    r2=${r1/_R1.fastq/_R2.fastq}  # Find the matching R2 file
    output_sam=${r1%%_R1.fastq}.sam
    hisat2 -x index_prefix -1 "$r1" -2 "$r2" -S "$output_sam"
done

# Convert SAM to BAM, sort, and generate statistics for all samples
for sam_file in *.sam; do
    bam_file=${sam_file%.sam}_sorted.bam
    samtools view -bS "$sam_file" | samtools sort -o "$bam_file"
    samtools flagstat "$bam_file" > "${bam_file%.bam}_alignment_statistics.txt"
done

# Run featureCounts on all sorted BAM files
featureCounts -T 8 -a annotation.gtf -o gene_counts.txt \
  -g gene_id -t exon -s 1 -p *.bam
