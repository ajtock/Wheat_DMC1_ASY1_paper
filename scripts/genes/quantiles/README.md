# Gene quantiles

This subdirectory contains the R script `group_genes_into_quantiles_and_popgenetics.R` for dividing genes into groups corresponding to those in given percentile ranges (e.g., 100th–75th (Quantile 1), 75th–50th (Quantile 2), 50th–25th (Quantile 3) and 25th–0th (Quantile 4)) with regard to various ordering factors, such as mean crossover recombination rate (cM/Mb) values, derived from a Chinese Spring × Renan genetic map, ChIP-seq signal, or population genetics statistics.
This script extends gene boundaries by 1 kb on each side, and computes mean cM/Mb values within these intervals using previously calculated mean recombination rates in 10-Mb sliding windows with a 1-Mb step (`iwgsc_refseqv1.0_recombination_rate.txt`, available as part of the [IWGSC RefSeq v1.0 annotation](https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.0/)).
`group_genes_into_quantiles_and_popgenetics.R` requires coordinates for representative gene models from the [IWGSC RefSeq v1.1 annotation](https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.1/) and, separately, random loci in BED6 format: column 1 = chromosome ID; column 2 = 0-based start coordinates; column 3 = 1-based end coordinates; column 4 = sequential or otherwise unique numbers; column 5 = fill with NA; column 6 = strand.
To compute population genetics statistics for each gene, this script also requires a variant call format (VCF) file containing ~3 million exome sequencing-derived SNP sites ([all.GP08_mm75_het3_publication01142019.vcf.gz](http://wheatgenomics.plantpath.ksu.edu/1000EC/)), with geographical information about each accession obtained from [He et al. (2019) *Nat. Genet.* **51**. DOI: 10.1038/s41588-019-0382-2](https://www.nature.com/articles/s41588-019-0382-2).

`metaprofiles/gene_quantile_metaprofiles.R` calculates and plots metaprofiles of ChIP-seq, MNase-seq and RNA-seq signals, DNA methylation proportions, and SNP and transposon frequencies (gene windowed means and 95% confidence intervals, CIs) for each group of genes, defined either by decreasing recombination rate, for example, or randomly.

`stat_means_95CIs/gene_quantile_cMMb_density_mean_95CI_plot.R` plots density and means with 95% CIs of recombination rate for each group of genes.

`stat_means_95CIs/gene_quantile_popgenetics_stats_density_mean_95CI_plot.R` plots density and means with 95% CIs of a given population genetics statistic for each group of genes.

`GO_term_enrichment/topGO_gene_quantiles.R` and `GO_term_enrichment/topGOslim_gene_quantiles.R` use the Bioconductor package [topGO](https://bioconductor.org/packages/release/bioc/html/topGO.html) to evaluate groups of genes ("quantiles", defined by decreasing crossover recombination rate or ChIP-seq signal) for over-representation gene ontology (GO) and high-level GO (GO slim) terms, respectively.

`hypergeometric_tests/proportion_*_in_gene_quantiles_hypergeometricTest.R` scripts evaluate different gene categories for over- and under-representation in each group of genes by applying hypergeometric tests.
These scripts also sample from the hypergeometric distribution 100,000 times to obtain a probability distribution for plotting the log<sub>2</sub>(observed/expected) ratio and significance threshold.
