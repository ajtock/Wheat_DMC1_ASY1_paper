#!/bin/bash

# Call peaks in ChIP-seq signal using the ranger tool from PeakRanger v1.18
# ( http://ranger.sourceforge.net/manual1.18.html ),
# providing alignments for an ChIP input library as a background control

# Usage:
# ./1_identify_ChIPseq_peaks.sh DMC1 DMC1_Rep1_ChIP input input_SRR6350669 0.001 0.01 200 48

ChIP=$1
ChIPLibName=$2
control=$3
controlLibName=$4
pval=$5
qval=$6
fragLen=$7
threads=$8

ChIPDir="/home/ajt200/analysis/wheat/"$ChIP"/snakemake_ChIPseq/mapped/both"
controlDir="/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/"$control"/snakemake_ChIPseq/mapped/both"
rangerDir="/home/ajt200/tools/PeakRanger-1.18/bin"

$rangerDir/peakranger ranger -d $ChIPDir/$ChIPLibName"_MappedOn_wheat_v1.0_lowXM_both_sort.bam" \
                             -c $controlDir/$controlLibName"_MappedOn_wheat_v1.0_lowXM_both_sort.bam" \
                             --format bam -o $ChIPLibName"_peaks_peakranger_ranger_p"$pval"_q"$qval \
                             -p $pval -q $qval -l $fragLen --pad -t $threads --verbose
