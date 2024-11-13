#!/bin/bash
# Converts SAM to sorted BAM and generates alignment statistics

samtools view -bS "$1" | samtools sort -o sorted.bam
samtools flagstat sorted.bam > alignment_statistics.txt
