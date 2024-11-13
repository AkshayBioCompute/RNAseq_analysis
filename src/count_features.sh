#!/bin/bash
# Counts features using featureCounts

featureCounts -T "$3" -a "$2" -o gene_counts.txt -g gene_id -t exon -s 1 -p "$@"
