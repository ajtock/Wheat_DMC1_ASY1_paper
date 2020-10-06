#!/bin/bash

# Generate background oligonucleotide frequency files for Weeder v2.0 analyses,
# based on the sequences of randomly positioned loci within the same
# subgenomic compartments as the peak sets to be analyzed

# Usage:
# ./w2frequency_maker_ranLoc_200bpseq.sh ASY1_CS_Rep1_ChIP_peaks_minuslog10Qsorted_in_Agenome_euchromatin_summits200bp

prefix=$1
toolDir="/home/ajt200/tools/Weeder2.0"

[ -d FreqFiles ] || mkdir FreqFiles

$toolDir/w2frequency_maker ${prefix}_randomLoci.fa TA ds

mv *.freq FreqFiles
