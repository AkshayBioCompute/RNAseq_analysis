#!/bin/bash
# Aligns reads with HISAT2

hisat2 -x index_prefix -1 "$1" -2 "$2" -S aligned.sam
