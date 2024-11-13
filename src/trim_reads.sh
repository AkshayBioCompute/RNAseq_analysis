#!/bin/bash
# Trims paired-end reads using Fastp

fastp -i "$1" -I "$2" -o trimmed_"$(basename "$1")" -O trimmed_"$(basename "$2")" --thread "$3" -h fastp.html -j fastp.json
