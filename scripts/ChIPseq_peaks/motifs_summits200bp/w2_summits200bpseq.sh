#!/bin/bash

# Evaluate the top 10,000 compartmentalized DMC1 ChIP-seq peaks
# in each subgenome (ordered by decreasing -log10(ranger-assigned FDR))
# for over-representation of DNA sequence motifs using Weeder v2.0

# Usage:
# ./w2_summits200bpseq.sh DMC1_Rep1_ChIP_peaks_minuslog10Qsorted_in_Agenome_euchromatin_summits200bp

prefix=$1
toolDir="/home/ajt200/tools/Weeder2.0"

$toolDir/weeder2 -f ${prefix}.fa \
                 -O TA -chipseq -top 10000
