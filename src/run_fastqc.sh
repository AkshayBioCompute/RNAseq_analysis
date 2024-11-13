#!/bin/bash
# Runs FastQC on each FASTQ file

fastqc -o . "$1"
