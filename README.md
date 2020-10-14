# Wheat DMC1 and ASY1 ChIP-seq custom scripts

This repository contains custom scripts used in the preparation of the paper entitled "Wheat meiotic recombination hotspots are marked by DMC1 recombinase, the chromosome axis, facultative heterochromatin and signatures of adaptation" by Andrew J. Tock, Daniel M. Holland, Wei Jiang, Kim Osman, Eugenio Sanchez-Moran, James D. Higgins, Keith J. Edwards, Cristobal Uauy, F.C.H. Franklin, and Ian R. Henderson.

These files can be downloaded together by cloning the repository:

```
git clone https://github.com/ajtock/Wheat_DMC1_ASY1_paper/
```

Alternatively, individual files (e.g., `Snakefile`) can be downloaded using `wget`:

```
wget https://raw.githubusercontent.com/ajtock/Wheat_DMC1_ASY1_paper/scripts/read_alignment/snakemake_ChIPseq_MNaseseq/Snakefile
```

## Alignment workflows

Workflows for processing and aligning next-generation sequencing (NGS) reads were developed using [Snakemake](https://snakemake.readthedocs.io/en/stable/) v3.13.3 in conjunction with [conda](https://conda.io/en/latest/) v4.6.9 package and environment manager.
These workflows are located in [scripts/read_alignment/](https://github.com/ajtock/Wheat_DMC1_ASY1_paper/tree/master/scripts/read_alignment) in this repository.

## Feature analyses

Scripts for analyzing ChIP-seq peaks and genes are located in [scripts/ChIPseq_peaks/](https://github.com/ajtock/Wheat_DMC1_ASY1_paper/tree/master/scripts/ChIPseq_peaks) and [scripts/genes/](https://github.com/ajtock/Wheat_DMC1_ASY1_paper/tree/master/scripts/genes).

## Chromosome compartments

Previously defined coordinates delimiting compartments of the chromosomes of hexaploid wheat landrace Chinese Spring are detailed in [chromosome_compartments/](https://github.com/ajtock/Wheat_DMC1_ASY1_paper/tree/master/chromosome_compartments).

