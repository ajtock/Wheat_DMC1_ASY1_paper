#!/bin/bash

# Extract nucleotide sequence at given coordinates (specified in a BED file)
# in the genome in fasta format

# Note that BED files use 0-based half-open coordinates;
# start coordinates are 0-based and end coordinates are 1-based,
# such that a feature START coordinate that cooresponds to the first base in
# a chromosome is numbered 0, while a feature END coordinate that corresponds
# to the second base in a chromosome is numbered 2 in the BED file:
# see https://genome.ucsc.edu/FAQ/FAQformat.html#format1

# Usage:
# ./bedtools_getfasta.sh wheat_v1.0.fa DMC1_Rep1_ChIP_peaks_minuslog10Qsorted_in_Agenome_euchromatin_summits200bp 

genome=$1
prefix=$2

/home/ajt200/anaconda3/bin/bedtools getfasta -fi ${genome} \
                                             -bed ${prefix}.bed \
                                             -fo ${prefix}.fa \
                                             -name
