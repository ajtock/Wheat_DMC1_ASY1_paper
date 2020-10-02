#!/bin/bash

# Build index for wheat reference genome (wheat_v1.0.fa)
# using bowtie2 version 2.3.4.3

# Usage:
# ./bowtie2_index.sh wheat_v1.0.fa wheat_v1.0

genome=$1
idxBaseName=$2

/home/ajt200/anaconda3/envs/ChIPseq_mapping/bin/bowtie2-build --threads 48 $genome $idxBaseName
