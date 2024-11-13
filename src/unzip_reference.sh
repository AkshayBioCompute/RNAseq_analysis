#!/bin/bash
# Unzips the reference genome file

gzip -d -c "$1" > reference_genome.fa
