import os
import subprocess
import glob

# Set the working directory
working_dir = "/path/to/your/working/directory"
os.chdir(working_dir)

# Generate MD5 checksums for all FASTQ files
with open("md5_checksums.txt", "w") as f:
    subprocess.run(["md5sum"] + glob.glob("*.fastq.gz"), stdout=f)

# Unzip all FASTQ files
subprocess.run(["gunzip", "*.fastq.gz"])

# Create directories for outputs
os.makedirs("QC_reports", exist_ok=True)

# Run FastQC on all FASTQ files
subprocess.run(["fastqc", "-o", "QC_reports/"] + glob.glob("*.fastq"))

# Generate MultiQC summary report
subprocess.run(["multiqc", "-o", "QC_reports/", "QC_reports/"])

# Trim paired-end reads using Fastp
for r1 in glob.glob("*_R1.fastq"):
    r2 = r1.replace("_R1.fastq", "_R2.fastq")  # Find matching R2 file
    output_r1 = f"trimmed_{r1}"
    output_r2 = f"trimmed_{r2}"
    html_report = f"QC_reports/fastp_{r1.split('_R1.fastq')[0]}.html"
    json_report = f"QC_reports/fastp_{r1.split('_R1.fastq')[0]}.json"
    subprocess.run(["fastp", "-i", r1, "-I", r2, "-o", output_r1, "-O", output_r2,
                    "-h", html_report, "-j", json_report, "--thread", "8"])

# Unzip the reference genome if not already unzipped
if os.path.exists("reference_genome.fa.gz"):
    subprocess.run(["gunzip", "reference_genome.fa.gz"])

# Build the Hisat2 index
subprocess.run(["hisat2-build", "reference_genome.fa", "index_prefix"])

# Align reads with Hisat2 for all paired-end fastq files
for r1 in glob.glob("trimmed_*_R1.fastq"):
    r2 = r1.replace("_R1.fastq", "_R2.fastq")  # Find matching R2 file
    output_sam = r1.replace("_R1.fastq", ".sam")
    subprocess.run(["hisat2", "-x", "index_prefix", "-1", r1, "-2", r2, "-S", output_sam])

# Convert SAM to BAM, sort, and generate statistics for all samples
for sam_file in glob.glob("*.sam"):
    bam_file = sam_file.replace(".sam", "_sorted.bam")
    with open(f"{bam_file.replace('_sorted.bam', '_alignment_statistics.txt')}", "w") as f:
        subprocess.run(["samtools", "view", "-bS", sam_file], stdout=subprocess.PIPE)
        subprocess.run(["samtools", "sort", "-o", bam_file], stdout=subprocess.PIPE)
        subprocess.run(["samtools", "flagstat", bam_file], stdout=f)

# Run featureCounts on all sorted BAM files
sorted_bam_files = glob.glob("*.bam")
subprocess.run(["featureCounts", "-T", "8", "-a", "annotation.gtf", "-o", "gene_counts.txt",
                "-g", "gene_id", "-t", "exon", "-s", "1", "-p"] + sorted_bam_files)
