# snakemake\_RNAseq

## Automated workflow for paired-end RNA-seq data processing and alignment

This is a Snakemake workflow for automated processing and alignment of paired-end RNA-seq data.

### Requirements

- Installation of [Snakemake](https://snakemake.readthedocs.io/en/stable/) and optionally [conda](https://conda.io/docs/)
- Demultiplexed paired-end reads in gzipped FASTQ format located in the `data/` directory. These should be named according to the following naming convention: `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz`
- A samtools-indexed reference genome in FASTA format and a chromosome sizes file (e.g., `wheat_v1.0.fa`, `wheat_v1.0.fa.fai`, and `wheat_v1.0.fa.sizes`, the latter two of which generated with `samtools faidx wheat_v1.0.fa; cut -f1,2 wheat_v1.0.fa.fai > wheat_v1.0.fa.sizes`), each located in `data/index/`
- Gene annotations in GTF format (e.g., `IWGSC_v1.1_HC_20170706.gtf`) and a reference genome index with exon and splice site annotations for hisat2, each located in `data/index/` (see `data/index/hisat2_genome_index.sh` for an example of how to generate these)
- Lists of potential contaminant sequences to be removed, located in `contaminants/`, provided or described in this repository
- `Snakefile` in this repository. This contains "rules" that each execute a step in the workflow
- `config.yaml` in this repository. This contains customizable parameters including `reference_prefix`, which should be the reference genome file name without the `.fa` extension (e.g., `wheat_v1.0`)
- Optional: `environment.yaml` in this repository, used to create the software environment if conda is used
- If conda is not used, the tools listed in environment.yaml must be specified in the PATH variable

This repository can be downloaded with:

```
git clone https://github.com/ajtock/Wheat_DMC1_ASY1_paper/scripts/read_alignment/snakemake_RNAseq/
```

Alternatively, individual files (e.g., `Snakefile`) can be downloaded using `wget`:

```
wget https://raw.githubusercontent.com/ajtock/Wheat_DMC1_ASY1_paper/scripts/read_alignment/snakemake_RNAseq/Snakefile
```

### Creating the conda environment

```
conda env create --file environment.yaml --name RNAseq_mapping
```

### Usage

In a Unix shell, navigate to the base directory containing `Snakefile`, `config.yaml`, `environment.yaml`, and the `data\` subdirectory, which should have a directory tree structure like this:

```
.
├── config.yaml
├── contaminants
│   ├── cat_all_and_TruSeq_Single_Indexes.fa
│   ├── contaminants_list_fastqc.txt
│   ├── ribokmers.fa.gz
│   └── ribokmers_README.txt
├── data
│   ├── index
│   │   ├── hisat2_genome_index.sh
│   │   ├── IWGSC_v1.1_HC_20170706.exon
│   │   ├── IWGSC_v1.1_HC_20170706.gtf
│   │   ├── IWGSC_v1.1_HC_20170706.ss
│   │   ├── samtools_faidx_chr_sizes.sh
│   │   ├── wheat_v1.0.1.ht2l
│   │   ├── wheat_v1.0.2.ht2l
│   │   ├── wheat_v1.0.3.ht2l
│   │   ├── wheat_v1.0.4.ht2l
│   │   ├── wheat_v1.0.5.ht2l
│   │   ├── wheat_v1.0.6.ht2l
│   │   ├── wheat_v1.0.7.ht2l
│   │   ├── wheat_v1.0.8.ht2l
│   │   ├── wheat_v1.0.fa
│   │   ├── wheat_v1.0.fa.fai
│   │   └── wheat_v1.0.fa.sizes
│   ├── WT_RNAseq_Rep1_ERR2402974_R1.fastq.gz
│   ├── WT_RNAseq_Rep1_ERR2402974_R2.fastq.gz
│   ├── WT_RNAseq_Rep2_ERR2402973_R1.fastq.gz
│   ├── WT_RNAseq_Rep2_ERR2402973_R2.fastq.gz
│   ├── WT_RNAseq_Rep3_ERR2402972_R1.fastq.gz
│   └── WT_RNAseq_Rep3_ERR2402972_R2.fastq.gz
├── environment.yaml
├── README.md
├── scripts
│   └── keepPaired.py
└── Snakefile
```

Then run the following commands in the base directory (`--cores` should match the `THREADS` parameter in `config.yaml`):

```
conda activate RNAseq_mapping
snakemake -p --cores 48
conda deactivate
```

### Useful Snakemake parameters

- `--cores` specifies the maximum number of threads
- `-n` performs a dry run
- `-p` prints commands
- `--use-conda`
- `--conda-prefix ~/.myconda`
- `--forcerun calc_coverage` forces rerun of a given rule (e.g., `calc_coverage`)

### Outputs

Below is the directory tree structure including files generated once the Snakemake workflow has run to completion.

```
.
├── config.yaml
├── contaminants
│   ├── cat_all_and_TruSeq_Single_Indexes.fa
│   ├── contaminants_list_fastqc.txt
│   ├── ribokmers.fa.gz
│   └── ribokmers_README.txt
├── data
│   ├── index
│   │   ├── IWGSC_v1.1_HC_20170706.exon
│   │   ├── IWGSC_v1.1_HC_20170706.gtf
│   │   ├── IWGSC_v1.1_HC_20170706.nss
│   │   ├── IWGSC_v1.1_HC_20170706.ss
│   │   ├── samtools_faidx_chr_sizes.sh
│   │   ├── wheat_v1.0.1.ht2l
│   │   ├── wheat_v1.0.2.ht2l
│   │   ├── wheat_v1.0.3.ht2l
│   │   ├── wheat_v1.0.4.ht2l
│   │   ├── wheat_v1.0.5.ht2l
│   │   ├── wheat_v1.0.6.ht2l
│   │   ├── wheat_v1.0.7.ht2l
│   │   ├── wheat_v1.0.8.ht2l
│   │   ├── wheat_v1.0.fa
│   │   ├── wheat_v1.0.fa.fai
│   │   └── wheat_v1.0.fa.sizes
│   ├── trimmed
│   │   ├── WT_RNAseq_Rep1_ERR2402974_R1_rRNAremoved_trimmed.fastq.gz
│   │   ├── WT_RNAseq_Rep1_ERR2402974_R1_rRNAremoved_trimmed_unpaired.fastq.gz
│   │   ├── WT_RNAseq_Rep1_ERR2402974_R2_rRNAremoved_trimmed.fastq.gz
│   │   ├── WT_RNAseq_Rep1_ERR2402974_R2_rRNAremoved_trimmed_unpaired.fastq.gz
│   │   ├── WT_RNAseq_Rep2_ERR2402973_R1_rRNAremoved_trimmed.fastq.gz
│   │   ├── WT_RNAseq_Rep2_ERR2402973_R1_rRNAremoved_trimmed_unpaired.fastq.gz
│   │   ├── WT_RNAseq_Rep2_ERR2402973_R2_rRNAremoved_trimmed.fastq.gz
│   │   ├── WT_RNAseq_Rep2_ERR2402973_R2_rRNAremoved_trimmed_unpaired.fastq.gz
│   │   ├── WT_RNAseq_Rep3_ERR2402972_R1_rRNAremoved_trimmed.fastq.gz
│   │   ├── WT_RNAseq_Rep3_ERR2402972_R1_rRNAremoved_trimmed_unpaired.fastq.gz
│   │   ├── WT_RNAseq_Rep3_ERR2402972_R2_rRNAremoved_trimmed.fastq.gz
│   │   └── WT_RNAseq_Rep3_ERR2402972_R2_rRNAremoved_trimmed_unpaired.fastq.gz
│   ├── WT_RNAseq_Rep1_ERR2402974_R1.fastq.gz
│   ├── WT_RNAseq_Rep1_ERR2402974_R1_rRNA.fastq.gz
│   ├── WT_RNAseq_Rep1_ERR2402974_R1_rRNAremoved.fastq.gz
│   ├── WT_RNAseq_Rep1_ERR2402974_R2.fastq.gz
│   ├── WT_RNAseq_Rep1_ERR2402974_R2_rRNA.fastq.gz
│   ├── WT_RNAseq_Rep1_ERR2402974_R2_rRNAremoved.fastq.gz
│   ├── WT_RNAseq_Rep2_ERR2402973_R1.fastq.gz
│   ├── WT_RNAseq_Rep2_ERR2402973_R1_rRNA.fastq.gz
│   ├── WT_RNAseq_Rep2_ERR2402973_R1_rRNAremoved.fastq.gz
│   ├── WT_RNAseq_Rep2_ERR2402973_R2.fastq.gz
│   ├── WT_RNAseq_Rep2_ERR2402973_R2_rRNA.fastq.gz
│   ├── WT_RNAseq_Rep2_ERR2402973_R2_rRNAremoved.fastq.gz
│   ├── WT_RNAseq_Rep3_ERR2402972_R1.fastq.gz
│   ├── WT_RNAseq_Rep3_ERR2402972_R1_rRNA.fastq.gz
│   ├── WT_RNAseq_Rep3_ERR2402972_R1_rRNAremoved.fastq.gz
│   ├── WT_RNAseq_Rep3_ERR2402972_R2.fastq.gz
│   ├── WT_RNAseq_Rep3_ERR2402972_R2_rRNA.fastq.gz
│   └── WT_RNAseq_Rep3_ERR2402972_R2_rRNAremoved.fastq.gz
├── environment.yaml
├── logs
│   ├── bamCoverage
│   │   ├── WT_RNAseq_Rep1_ERR2402974_MappedOn_wheat_v1.0_both_sort_norm_binSize10kb.log
│   │   ├── WT_RNAseq_Rep1_ERR2402974_MappedOn_wheat_v1.0_both_sort_norm_binSize1Mb.log
│   │   ├── WT_RNAseq_Rep1_ERR2402974_MappedOn_wheat_v1.0_both_sort_norm.log
│   │   ├── WT_RNAseq_Rep1_ERR2402974_MappedOn_wheat_v1.0_unique_sort_norm_binSize10kb.log
│   │   ├── WT_RNAseq_Rep1_ERR2402974_MappedOn_wheat_v1.0_unique_sort_norm_binSize1Mb.log
│   │   ├── WT_RNAseq_Rep1_ERR2402974_MappedOn_wheat_v1.0_unique_sort_norm.log
│   │   ├── WT_RNAseq_Rep2_ERR2402973_MappedOn_wheat_v1.0_both_sort_norm_binSize10kb.log
│   │   ├── WT_RNAseq_Rep2_ERR2402973_MappedOn_wheat_v1.0_both_sort_norm_binSize1Mb.log
│   │   ├── WT_RNAseq_Rep2_ERR2402973_MappedOn_wheat_v1.0_both_sort_norm.log
│   │   ├── WT_RNAseq_Rep2_ERR2402973_MappedOn_wheat_v1.0_unique_sort_norm_binSize10kb.log
│   │   ├── WT_RNAseq_Rep2_ERR2402973_MappedOn_wheat_v1.0_unique_sort_norm_binSize1Mb.log
│   │   ├── WT_RNAseq_Rep2_ERR2402973_MappedOn_wheat_v1.0_unique_sort_norm.log
│   │   ├── WT_RNAseq_Rep3_ERR2402972_MappedOn_wheat_v1.0_both_sort_norm_binSize10kb.log
│   │   ├── WT_RNAseq_Rep3_ERR2402972_MappedOn_wheat_v1.0_both_sort_norm_binSize1Mb.log
│   │   ├── WT_RNAseq_Rep3_ERR2402972_MappedOn_wheat_v1.0_both_sort_norm.log
│   │   ├── WT_RNAseq_Rep3_ERR2402972_MappedOn_wheat_v1.0_unique_sort_norm_binSize10kb.log
│   │   ├── WT_RNAseq_Rep3_ERR2402972_MappedOn_wheat_v1.0_unique_sort_norm_binSize1Mb.log
│   │   └── WT_RNAseq_Rep3_ERR2402972_MappedOn_wheat_v1.0_unique_sort_norm.log
│   ├── bbduk
│   │   ├── WT_RNAseq_Rep1_ERR2402974_rRNAremoved.log
│   │   ├── WT_RNAseq_Rep2_ERR2402973_rRNAremoved.log
│   │   └── WT_RNAseq_Rep3_ERR2402972_rRNAremoved.log
│   ├── fastqc
│   │   ├── raw
│   │   │   ├── WT_RNAseq_Rep1_ERR2402974_R1_fastqc.html
│   │   │   ├── WT_RNAseq_Rep1_ERR2402974_R1_fastqc.zip
│   │   │   ├── WT_RNAseq_Rep1_ERR2402974_R1.log
│   │   │   ├── WT_RNAseq_Rep1_ERR2402974_R2_fastqc.html
│   │   │   ├── WT_RNAseq_Rep1_ERR2402974_R2_fastqc.zip
│   │   │   ├── WT_RNAseq_Rep1_ERR2402974_R2.log
│   │   │   ├── WT_RNAseq_Rep2_ERR2402973_R1_fastqc.html
│   │   │   ├── WT_RNAseq_Rep2_ERR2402973_R1_fastqc.zip
│   │   │   ├── WT_RNAseq_Rep2_ERR2402973_R1.log
│   │   │   ├── WT_RNAseq_Rep2_ERR2402973_R2_fastqc.html
│   │   │   ├── WT_RNAseq_Rep2_ERR2402973_R2_fastqc.zip
│   │   │   ├── WT_RNAseq_Rep2_ERR2402973_R2.log
│   │   │   ├── WT_RNAseq_Rep3_ERR2402972_R1_fastqc.html
│   │   │   ├── WT_RNAseq_Rep3_ERR2402972_R1_fastqc.zip
│   │   │   ├── WT_RNAseq_Rep3_ERR2402972_R1.log
│   │   │   ├── WT_RNAseq_Rep3_ERR2402972_R2_fastqc.html
│   │   │   ├── WT_RNAseq_Rep3_ERR2402972_R2_fastqc.zip
│   │   │   └── WT_RNAseq_Rep3_ERR2402972_R2.log
│   │   └── trimmed
│   │       ├── WT_RNAseq_Rep1_ERR2402974_R1_rRNAremoved_trimmed_fastqc.html
│   │       ├── WT_RNAseq_Rep1_ERR2402974_R1_rRNAremoved_trimmed_fastqc.zip
│   │       ├── WT_RNAseq_Rep1_ERR2402974_R1_rRNAremoved_trimmed.log
│   │       ├── WT_RNAseq_Rep1_ERR2402974_R2_rRNAremoved_trimmed_fastqc.html
│   │       ├── WT_RNAseq_Rep1_ERR2402974_R2_rRNAremoved_trimmed_fastqc.zip
│   │       ├── WT_RNAseq_Rep1_ERR2402974_R2_rRNAremoved_trimmed.log
│   │       ├── WT_RNAseq_Rep2_ERR2402973_R1_rRNAremoved_trimmed_fastqc.html
│   │       ├── WT_RNAseq_Rep2_ERR2402973_R1_rRNAremoved_trimmed_fastqc.zip
│   │       ├── WT_RNAseq_Rep2_ERR2402973_R1_rRNAremoved_trimmed.log
│   │       ├── WT_RNAseq_Rep2_ERR2402973_R2_rRNAremoved_trimmed_fastqc.html
│   │       ├── WT_RNAseq_Rep2_ERR2402973_R2_rRNAremoved_trimmed_fastqc.zip
│   │       ├── WT_RNAseq_Rep2_ERR2402973_R2_rRNAremoved_trimmed.log
│   │       ├── WT_RNAseq_Rep3_ERR2402972_R1_rRNAremoved_trimmed_fastqc.html
│   │       ├── WT_RNAseq_Rep3_ERR2402972_R1_rRNAremoved_trimmed_fastqc.zip
│   │       ├── WT_RNAseq_Rep3_ERR2402972_R1_rRNAremoved_trimmed.log
│   │       ├── WT_RNAseq_Rep3_ERR2402972_R2_rRNAremoved_trimmed_fastqc.html
│   │       ├── WT_RNAseq_Rep3_ERR2402972_R2_rRNAremoved_trimmed_fastqc.zip
│   │       └── WT_RNAseq_Rep3_ERR2402972_R2_rRNAremoved_trimmed.log
│   ├── hisat2
│   │   ├── WT_RNAseq_Rep1_ERR2402974_MappedOn_wheat_v1.0.log
│   │   ├── WT_RNAseq_Rep2_ERR2402973_MappedOn_wheat_v1.0.log
│   │   └── WT_RNAseq_Rep3_ERR2402972_MappedOn_wheat_v1.0.log
│   ├── samtools
│   │   ├── stats
│   │   │   ├── WT_RNAseq_Rep1_ERR2402974_MappedOn_wheat_v1.0_both_sort_flagstat.log
│   │   │   ├── WT_RNAseq_Rep1_ERR2402974_MappedOn_wheat_v1.0_both_sort_idxstats.log
│   │   │   ├── WT_RNAseq_Rep1_ERR2402974_MappedOn_wheat_v1.0_unique_sort_flagstat.log
│   │   │   ├── WT_RNAseq_Rep1_ERR2402974_MappedOn_wheat_v1.0_unique_sort_idxstats.log
│   │   │   ├── WT_RNAseq_Rep2_ERR2402973_MappedOn_wheat_v1.0_both_sort_flagstat.log
│   │   │   ├── WT_RNAseq_Rep2_ERR2402973_MappedOn_wheat_v1.0_both_sort_idxstats.log
│   │   │   ├── WT_RNAseq_Rep2_ERR2402973_MappedOn_wheat_v1.0_unique_sort_flagstat.log
│   │   │   ├── WT_RNAseq_Rep2_ERR2402973_MappedOn_wheat_v1.0_unique_sort_idxstats.log
│   │   │   ├── WT_RNAseq_Rep3_ERR2402972_MappedOn_wheat_v1.0_both_sort_flagstat.log
│   │   │   ├── WT_RNAseq_Rep3_ERR2402972_MappedOn_wheat_v1.0_both_sort_idxstats.log
│   │   │   ├── WT_RNAseq_Rep3_ERR2402972_MappedOn_wheat_v1.0_unique_sort_flagstat.log
│   │   │   └── WT_RNAseq_Rep3_ERR2402972_MappedOn_wheat_v1.0_unique_sort_idxstats.log
│   │   ├── WT_RNAseq_Rep1_ERR2402974_MappedOn_wheat_v1.0_both_sort.log
│   │   ├── WT_RNAseq_Rep1_ERR2402974_MappedOn_wheat_v1.0_unique_sort.log
│   │   ├── WT_RNAseq_Rep2_ERR2402973_MappedOn_wheat_v1.0_both_sort.log
│   │   ├── WT_RNAseq_Rep2_ERR2402973_MappedOn_wheat_v1.0_unique_sort.log
│   │   ├── WT_RNAseq_Rep3_ERR2402972_MappedOn_wheat_v1.0_both_sort.log
│   │   └── WT_RNAseq_Rep3_ERR2402972_MappedOn_wheat_v1.0_unique_sort.log
│   └── trimmomatic
│       ├── WT_RNAseq_Rep1_ERR2402974_rRNAremoved_trimmed.log
│       ├── WT_RNAseq_Rep2_ERR2402973_rRNAremoved_trimmed.log
│       └── WT_RNAseq_Rep3_ERR2402972_rRNAremoved_trimmed.log
├── mapped
│   ├── both
│   │   ├── bg
│   │   │   ├── WT_RNAseq_Rep1_ERR2402974_MappedOn_wheat_v1.0_both_sort_norm.bedgraph
│   │   │   ├── WT_RNAseq_Rep1_ERR2402974_MappedOn_wheat_v1.0_both_sort_norm_binSize10kb.bedgraph
│   │   │   ├── WT_RNAseq_Rep1_ERR2402974_MappedOn_wheat_v1.0_both_sort_norm_binSize1Mb.bedgraph
│   │   │   ├── WT_RNAseq_Rep2_ERR2402973_MappedOn_wheat_v1.0_both_sort_norm.bedgraph
│   │   │   ├── WT_RNAseq_Rep2_ERR2402973_MappedOn_wheat_v1.0_both_sort_norm_binSize10kb.bedgraph
│   │   │   ├── WT_RNAseq_Rep2_ERR2402973_MappedOn_wheat_v1.0_both_sort_norm_binSize1Mb.bedgraph
│   │   │   ├── WT_RNAseq_Rep3_ERR2402972_MappedOn_wheat_v1.0_both_sort_norm.bedgraph
│   │   │   ├── WT_RNAseq_Rep3_ERR2402972_MappedOn_wheat_v1.0_both_sort_norm_binSize10kb.bedgraph
│   │   │   └── WT_RNAseq_Rep3_ERR2402972_MappedOn_wheat_v1.0_both_sort_norm_binSize1Mb.bedgraph
│   │   ├── bw
│   │   │   ├── WT_RNAseq_Rep1_ERR2402974_MappedOn_wheat_v1.0_both_sort_norm.bw
│   │   │   ├── WT_RNAseq_Rep2_ERR2402973_MappedOn_wheat_v1.0_both_sort_norm.bw
│   │   │   └── WT_RNAseq_Rep3_ERR2402972_MappedOn_wheat_v1.0_both_sort_norm.bw
│   │   ├── WT_RNAseq_Rep1_ERR2402974_MappedOn_wheat_v1.0_both_sort.bam
│   │   ├── WT_RNAseq_Rep1_ERR2402974_MappedOn_wheat_v1.0_both_sort.bam.csi
│   │   ├── WT_RNAseq_Rep2_ERR2402973_MappedOn_wheat_v1.0_both_sort.bam
│   │   ├── WT_RNAseq_Rep2_ERR2402973_MappedOn_wheat_v1.0_both_sort.bam.csi
│   │   ├── WT_RNAseq_Rep3_ERR2402972_MappedOn_wheat_v1.0_both_sort.bam
│   │   └── WT_RNAseq_Rep3_ERR2402972_MappedOn_wheat_v1.0_both_sort.bam.csi
│   ├── unique
│   │   ├── bg
│   │   │   ├── WT_RNAseq_Rep1_ERR2402974_MappedOn_wheat_v1.0_unique_sort_norm.bedgraph
│   │   │   ├── WT_RNAseq_Rep1_ERR2402974_MappedOn_wheat_v1.0_unique_sort_norm_binSize10kb.bedgraph
│   │   │   ├── WT_RNAseq_Rep1_ERR2402974_MappedOn_wheat_v1.0_unique_sort_norm_binSize1Mb.bedgraph
│   │   │   ├── WT_RNAseq_Rep2_ERR2402973_MappedOn_wheat_v1.0_unique_sort_norm.bedgraph
│   │   │   ├── WT_RNAseq_Rep2_ERR2402973_MappedOn_wheat_v1.0_unique_sort_norm_binSize10kb.bedgraph
│   │   │   ├── WT_RNAseq_Rep2_ERR2402973_MappedOn_wheat_v1.0_unique_sort_norm_binSize1Mb.bedgraph
│   │   │   ├── WT_RNAseq_Rep3_ERR2402972_MappedOn_wheat_v1.0_unique_sort_norm.bedgraph
│   │   │   ├── WT_RNAseq_Rep3_ERR2402972_MappedOn_wheat_v1.0_unique_sort_norm_binSize10kb.bedgraph
│   │   │   └── WT_RNAseq_Rep3_ERR2402972_MappedOn_wheat_v1.0_unique_sort_norm_binSize1Mb.bedgraph
│   │   ├── bw
│   │   │   ├── WT_RNAseq_Rep1_ERR2402974_MappedOn_wheat_v1.0_unique_sort_norm.bw
│   │   │   ├── WT_RNAseq_Rep2_ERR2402973_MappedOn_wheat_v1.0_unique_sort_norm.bw
│   │   │   └── WT_RNAseq_Rep3_ERR2402972_MappedOn_wheat_v1.0_unique_sort_norm.bw
│   │   ├── WT_RNAseq_Rep1_ERR2402974_MappedOn_wheat_v1.0_unique_sort.bam
│   │   ├── WT_RNAseq_Rep1_ERR2402974_MappedOn_wheat_v1.0_unique_sort.bam.csi
│   │   ├── WT_RNAseq_Rep2_ERR2402973_MappedOn_wheat_v1.0_unique_sort.bam
│   │   ├── WT_RNAseq_Rep2_ERR2402973_MappedOn_wheat_v1.0_unique_sort.bam.csi
│   │   ├── WT_RNAseq_Rep3_ERR2402972_MappedOn_wheat_v1.0_unique_sort.bam
│   │   └── WT_RNAseq_Rep3_ERR2402972_MappedOn_wheat_v1.0_unique_sort.bam.csi
│   ├── WT_RNAseq_Rep1_ERR2402974_MappedOn_wheat_v1.0.bam
│   ├── WT_RNAseq_Rep2_ERR2402973_MappedOn_wheat_v1.0.bam
│   └── WT_RNAseq_Rep3_ERR2402972_MappedOn_wheat_v1.0.bam
├── scripts
│   └── keepPaired.py
├── README.md
└── Snakefile
```

### Updating the conda environment

```
conda env update --file environment.yaml --name RNAseq_mapping
```
