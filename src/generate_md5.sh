#!/bin/bash
# Generates MD5 checksums for the input FASTQ files

md5sum "$@" > md5_checksums.txt
