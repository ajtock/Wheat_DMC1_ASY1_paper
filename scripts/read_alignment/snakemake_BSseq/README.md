# snakemake\_BSseq

## Automated workflow for whole-genome BS-seq (Illumina) and EM-Seq (NEB) data

This is a Snakemake workflow for automated processing of DNA methylation data derived from whole-genome bisulfite sequencing (BS-seq or WGBS) or NEBNext Enzymatic Methyl-seq (EM-seq):
- [**BS-seq**](https://en.wikipedia.org/wiki/Bisulfite_sequencing)
- [**EM-seq**](https://international.neb.com/about-neb/news-and-press-releases/new-england-biolabs-to-present-latest-innovations-for-ngs-sample-preparation-at-agbt-2017)

Both BS-seq and EM-seq produce the same kind of data, so no adjustment for post-processing is needed.

***IMPORTANT***: This Snakemake pipeline should not be run with Bismark version 0.21.0 or later due to the addition of HISAT2 support, which requires Python 2, which conflicts with the Python 3 requirements of other parts of this pipeline (e.g., the pigz part of the `trim_galore` rule).

### Requirements

- Installation of [snakemake](https://snakemake.readthedocs.io/en/stable/) and optionally [conda](https://conda.io/docs/)
- Demultiplexed paired-end reads in gzipped FASTQ format located in the `data/` directory. These should be named according to the following naming convention: `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz`
- A reference genome in FASTA format (e.g., `wheat_v1.0_incl_organelles_controls.fa`) and a chromosome sizes file (e.g., `wheat_v1.0_incl_organelles_controls.fa.sizes`, generated with `samtools faidx wheat_v1.0_incl_organelles_controls.fa; cut -f1,2 wheat_v1.0_incl_organelles_controls.fa.fai > wheat_v1.0_incl_organelles_controls.fa.sizes`), both located in `data/index/`. No other FASTA format files should be located in this directory
- A bisulfite-converted reference genome index for bowtie2, located in `data/index/` (see `data/index/bismark_genome_preparation.sh` for an example of how to generate `data/index/Bisulfite_Genome/`)
- `Snakefile` in this repository. This contains "rules" that each execute a step in the workflow
- `config.yaml` in this repository. This contains customizable parameters including `reference_prefix`, which should be the reference genome file name without the `.fa` extension (e.g., `wheat_v1.0_incl_organelles_controls`)
- Optional: `environment.yaml` in this repository, used to create the software environment if conda is used
- If conda is not used, `bismark`, `bowtie2`, `fastqc`, `trim_galore`, `samtools`, `deeptools`, `ucsc-bedgraphtobigwig` and `python3` must be specified in the PATH variable

This repository can be downloaded with:

```
git clone https://github.com/ajtock/Wheat_DMC1_ASY1_paper/scripts/read_alignment/snakemake_BSseq/
```

Alternatively, individual files (e.g., `Snakefile`) can be downloaded using `wget`:

```
wget https://raw.githubusercontent.com/ajtock/Wheat_DMC1_ASY1_paper/scripts/read_alignment/snakemake_BSseq/Snakefile
```

### Creating the conda environment

```
conda env create --file environment.yaml --name BSseq_mapping
```

### Usage

In a Unix shell, navigate to the base directory containing `Snakefile`, `config.yaml`, `environment.yaml`, and the `data\` subdirectory, which should have a directory tree structure like this:

```
.
├── config.yaml
├── data
│   ├── BSseq_Rep8a_SRR6792678_R1.fastq.gz
│   ├── BSseq_Rep8a_SRR6792678_R2.fastq.gz
│   └── index
│       ├── bismark_genome_preparation.sh
│       ├── Bisulfite_Genome
│       │   ├── CT_conversion
│       │   │   ├── BS_CT.1.bt2l
│       │   │   ├── BS_CT.2.bt2l
│       │   │   ├── BS_CT.3.bt2l
│       │   │   ├── BS_CT.4.bt2l
│       │   │   ├── BS_CT.rev.1.bt2l
│       │   │   ├── BS_CT.rev.2.bt2l
│       │   │   └── genome_mfa.CT_conversion.fa
│       │   └── GA_conversion
│       │       ├── BS_GA.1.bt2l
│       │       ├── BS_GA.2.bt2l
│       │       ├── BS_GA.3.bt2l
│       │       ├── BS_GA.4.bt2l
│       │       ├── BS_GA.rev.1.bt2l
│       │       ├── BS_GA.rev.2.bt2l
│       │       └── genome_mfa.GA_conversion.fa
│       ├── cat_wheat_nuclear_chloroplast_mitochondrial_genomes_and_control_sequences.sh
│       ├── genomic_nucleotide_frequencies.txt
│       ├── README_chloroplast_Lambda_pUC19_control_sequences.txt
│       ├── samtools_faidx_chr_sizes.sh
│       ├── separate_nuclear_organelles_controls
│       │   ├── lambda_NEB.fa
│       │   ├── pUC19_NEB.fa
│       │   ├── wheat_CS_chloroplast_genome.fa
│       │   ├── wheat_CS_mitochondrial_genome.fa
│       │   └── wheat_v1.0.fa
│       ├── wheat_v1.0_incl_organelles_controls.fa
│       ├── wheat_v1.0_incl_organelles_controls.fa.fai
│       └── wheat_v1.0_incl_organelles_controls.fa.sizes
├── environment.yaml
├── README.md
└── Snakefile
```

Then run the following commands in the base directory (`--cores` should match the `THREADS` parameter in `config.yaml`):

```
conda activate BSseq_mapping
snakemake -p --cores 48
conda deactivate
```

### Useful Snakemake parameters

- `--cores` specifies the maximum number of threads
- `-n` performs a dry run
- `-p` prints commands
- `--use-conda`
- `--conda-prefix ~/.myconda`
- `--forcerun bismark2bedGraph` forces rerun of a given rule (e.g., `bismark2bedGraph`)

### Outputs

Below is the directory tree structure including files generated once the Snakemake workflow has run to completion.

```
.
├── config.yaml
├── coverage
│   ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_CHG
│   ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_CHG.gz
│   ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_CHG.gz.bismark.cov
│   ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_CHG.gz.bismark.cov.gz
│   ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_CHH
│   ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_CHH.gz
│   ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_CHH.gz.bismark.cov
│   ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_CHH.gz.bismark.cov.gz
│   ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_CpG
│   ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_CpG.gz
│   ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_CpG.gz.bismark.cov
│   ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_CpG.gz.bismark.cov.gz
│   ├── bw
│   │   ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_CHG.bw
│   │   ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_CHH.bw
│   │   └── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_CpG.bw
│   └── report
│       ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_CHG.CX_report.txt.gz
│       ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_CHH.CX_report.txt.gz
│       └── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_CpG.CX_report.txt.gz
├── data
│   ├── BSseq_Rep8a_SRR6792678_R1.fastq.gz
│   ├── BSseq_Rep8a_SRR6792678_R2.fastq.gz
│   └── index
│       ├── bismark_genome_preparation.sh
│       ├── Bisulfite_Genome
│       │   ├── CT_conversion
│       │   │   ├── BS_CT.1.bt2l
│       │   │   ├── BS_CT.2.bt2l
│       │   │   ├── BS_CT.3.bt2l
│       │   │   ├── BS_CT.4.bt2l
│       │   │   ├── BS_CT.rev.1.bt2l
│       │   │   ├── BS_CT.rev.2.bt2l
│       │   │   └── genome_mfa.CT_conversion.fa
│       │   └── GA_conversion
│       │       ├── BS_GA.1.bt2l
│       │       ├── BS_GA.2.bt2l
│       │       ├── BS_GA.3.bt2l
│       │       ├── BS_GA.4.bt2l
│       │       ├── BS_GA.rev.1.bt2l
│       │       ├── BS_GA.rev.2.bt2l
│       │       └── genome_mfa.GA_conversion.fa
│       ├── cat_wheat_nuclear_chloroplast_mitochondrial_genomes_and_control_sequences.sh
│       ├── genomic_nucleotide_frequencies.txt
│       ├── README_chloroplast_Lambda_pUC19_control_sequences.txt
│       ├── samtools_faidx_chr_sizes.sh
│       ├── separate_nuclear_organelles_controls
│       │   ├── lambda_NEB.fa
│       │   ├── pUC19_NEB.fa
│       │   ├── wheat_CS_chloroplast_genome.fa
│       │   ├── wheat_CS_mitochondrial_genome.fa
│       │   └── wheat_v1.0.fa
│       ├── wheat_v1.0_incl_organelles_controls.fa
│       ├── wheat_v1.0_incl_organelles_controls.fa.fai
│       └── wheat_v1.0_incl_organelles_controls.fa.sizes
├── environment.yaml
├── logs
│   ├── bamCoverage
│   │   └── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_sort_bw.log
│   ├── bedGraphToBigWig
│   │   ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_CHG_bedGraphToBigWig.log
│   │   ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_CHH_bedGraphToBigWig.log
│   │   └── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_CpG_bedGraphToBigWig.log
│   ├── bismark
│   │   ├── BSseq_Rep2a_SRR6792682_MappedOn_wheat_v1.0_incl_organelles_controls.log
│   │   ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls.log
│   │   ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls.nucleotide_stats.txt
│   │   └── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_report.txt
│   ├── bismark2bedGraph
│   │   ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_CHG_bismark2bedGraph.log
│   │   ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_CHH_bismark2bedGraph.log
│   │   └── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_CpG_bismark2bedGraph.log
│   ├── coverage2cytosine
│   │   ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_CHG_coverage2cytosine.log
│   │   ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_CHH_coverage2cytosine.log
│   │   └── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_CpG_coverage2cytosine.log
│   ├── deduplicate_bismark
│   │   ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup.log
│   │   └── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_report.txt
│   ├── fastqc
│   │   ├── raw
│   │   │   ├── BSseq_Rep2a_SRR6792682_R1_fastqc.html
│   │   │   ├── BSseq_Rep2a_SRR6792682_R1_fastqc.zip
│   │   │   ├── BSseq_Rep2a_SRR6792682_R1.log
│   │   │   ├── BSseq_Rep2a_SRR6792682_R2_fastqc.html
│   │   │   ├── BSseq_Rep2a_SRR6792682_R2_fastqc.zip
│   │   │   ├── BSseq_Rep2a_SRR6792682_R2.log
│   │   │   ├── BSseq_Rep8a_SRR6792678_R1_fastqc.html
│   │   │   ├── BSseq_Rep8a_SRR6792678_R1_fastqc.zip
│   │   │   ├── BSseq_Rep8a_SRR6792678_R1.log
│   │   │   ├── BSseq_Rep8a_SRR6792678_R2_fastqc.html
│   │   │   ├── BSseq_Rep8a_SRR6792678_R2_fastqc.zip
│   │   │   └── BSseq_Rep8a_SRR6792678_R2.log
│   │   └── trimmed
│   │       ├── BSseq_Rep2a_SRR6792682_R1_trimmed_fastqc.html
│   │       ├── BSseq_Rep2a_SRR6792682_R1_trimmed_fastqc.zip
│   │       ├── BSseq_Rep2a_SRR6792682_R1_trimmed.log
│   │       ├── BSseq_Rep2a_SRR6792682_R2_trimmed_fastqc.html
│   │       ├── BSseq_Rep2a_SRR6792682_R2_trimmed_fastqc.zip
│   │       ├── BSseq_Rep2a_SRR6792682_R2_trimmed.log
│   │       ├── BSseq_Rep8a_SRR6792678_R1_trimmed_fastqc.html
│   │       ├── BSseq_Rep8a_SRR6792678_R1_trimmed_fastqc.zip
│   │       ├── BSseq_Rep8a_SRR6792678_R1_trimmed.log
│   │       ├── BSseq_Rep8a_SRR6792678_R2_trimmed_fastqc.html
│   │       ├── BSseq_Rep8a_SRR6792678_R2_trimmed_fastqc.zip
│   │       └── BSseq_Rep8a_SRR6792678_R2_trimmed.log
│   ├── methylation_extractor
│   │   └── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_methylation_extractor.log
│   └── trim_galore
│       ├── BSseq_Rep2a_SRR6792682_R1.fastq.gz_trimming_report.txt
│       ├── BSseq_Rep2a_SRR6792682_R2.fastq.gz_trimming_report.txt
│       ├── BSseq_Rep2a_SRR6792682_trimmed.log
│       ├── BSseq_Rep8a_SRR6792678_R1.fastq.gz_trimming_report.txt
│       ├── BSseq_Rep8a_SRR6792678_R2.fastq.gz_trimming_report.txt
│       └── BSseq_Rep8a_SRR6792678_trimmed.log
├── mapped
│   ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls.bam
│   └── dedup
│       ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup.bam
│       ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_sort.bam
│       ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_sort.bam.csi
│       └── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_sort.bw
├── methylation_extracted
│   ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup.M-bias.txt
│   ├── BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_splitting_report.txt
│   ├── CHG_context_BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup.txt.gz
│   ├── CHH_context_BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup.txt.gz
│   └── CpG_context_BSseq_Rep8a_SRR6792678_MappedOn_wheat_v1.0_incl_organelles_controls_dedup.txt.gz
├── methylation_extractor_help.txt
├── README.md
├── Snakefile
└── trimmed
    ├── BSseq_Rep2a_SRR6792682_R1_trimmed.fastq.gz
    ├── BSseq_Rep2a_SRR6792682_R2_trimmed.fastq.gz
    ├── BSseq_Rep8a_SRR6792678_R1_trimmed.fastq.gz
    └── BSseq_Rep8a_SRR6792678_R2_trimmed.fastq.gz
```

### Updating the conda environment

```
conda env update --file environment.yaml --name BSseq_mapping
```
