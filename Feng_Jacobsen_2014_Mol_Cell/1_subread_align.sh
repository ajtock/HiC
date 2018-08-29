#!/bin/bash

# Subread v1.6.2

# Use Subread to align reads in fastq format to the indexed TAIR10 assembly
# Supply the "--trim3 50" option where reads are longer than 50 bp

i=$1

echo "$i"
subread-align -i /projects/ajt200/TAIR10/TAIR10_chr_all_GI.nix -r "${i}.fastq" -t 1 \
              -o "${i}.bam" -P 3 -M 0 -I 0 -T 24
