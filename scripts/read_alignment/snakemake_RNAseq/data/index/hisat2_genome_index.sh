#!/bin/bash

# Generate wheat_v1.0.fa reference genome index
# with IWGSC_v1.1_HC_20170706.gtf annotation
# using hisat2 version 2.1.0

THREADS=48
GENOME_PREFIX=wheat_v1.0
GTF_PREFIX=IWGSC_v1.1_HC_20170706

HISAT2_DIR=/home/ajt200/anaconda3/envs/RNAseq_mapping/bin

HISAT2_BUILD_EXE=$HISAT2_DIR/hisat2-build
if [ ! -x "$HISAT2_BUILD_EXE" ] ; then
  echo "Could not find this hisat2-build installation"
  exit 1
fi

HISAT2_SS_SCRIPT=$HISAT2_DIR/hisat2_extract_splice_sites.py
if [ ! -x "$HISAT2_SS_SCRIPT" ] ; then
  echo "Could not find this hisat2_extract_splice_sites.py installation"
  exit 1
fi

HISAT2_EXON_SCRIPT=$HISAT2_DIR/hisat2_extract_exons.py
if [ ! -x "$HISAT2_EXON_SCRIPT" ] ; then
  echo "Could not find this hisat2_extract_exons.py installation"
  exit 1
fi

if [ -f ${GTF_PREFIX}.gtf ] ; then
  ${HISAT2_SS_SCRIPT} ${GTF_PREFIX}.gtf > ${GTF_PREFIX}.ss
  ${HISAT2_EXON_SCRIPT} ${GTF_PREFIX}.gtf > ${GTF_PREFIX}.exon
else
  echo "Could not find the GTF-format annotation file"
  exit 1
fi

if [ ! -f ${GENOME_PREFIX}.fa ] ; then
  echo "Could not find the reference genome file"
  exit 1
fi

# hisat2 USAGE:
# hisat2-build [options]* <reference_in> <ht2_base>
CMD="${HISAT2_BUILD_EXE} -p ${THREADS} ${GENOME_PREFIX}.fa ${GENOME_PREFIX}"
#CMD="${HISAT2_BUILD_EXE} -p ${THREADS} --ss ${GTF_PREFIX}.ss --exon ${GTF_PREFIX}.exon ${GENOME_PREFIX}.fa ${GENOME_PREFIX}"
echo Running $CMD
if $CMD ; then
  echo "Index built successfully"
else
  echo "Index building failed; see error message"
fi
