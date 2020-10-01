#!/bin/bash

conda activate BSseq_mapping
snakemake -p --cores 48
conda deactivate
