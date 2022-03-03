#!/bin/bash

cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --minimum-length 25 \
    -q 22 \ 
    -o trimmed.R1.fastq.gz -p trimmed.R2.fastq.gz \
    reads.R1.fastq.gz reads.R2.fastq.gz


