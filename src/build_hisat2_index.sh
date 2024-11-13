#!/bin/bash
# Builds a HISAT2 index for the reference genome

hisat2-build "$1" index_prefix
