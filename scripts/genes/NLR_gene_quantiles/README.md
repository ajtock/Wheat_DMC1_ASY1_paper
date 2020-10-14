# Gene quantiles

This subdirectory contains the R script `group_NLR_genes_into_quantiles.R` for dividing NLR-encoding genes into groups corresponding to those in given percentile ranges (e.g., 100th–75th (Quantile 1), 75th–50th (Quantile 2), 50th–25th (Quantile 3) and 25th–0th (Quantile 4)) with regard to different ordering factors, such as mean crossover recombination rate (cM/Mb) values, derived from a Chinese Spring × Renan genetic map, ChIP-seq signal, or physical cluster size.

`metaprofiles/gene_quantile_metaprofiles.R` calculates and plots metaprofiles of ChIP-seq, MNase-seq and RNA-seq signals, DNA methylation proportions, and SNP and transposon frequencies (gene windowed means and 95% confidence intervals, CIs) for each group of NLR genes, defined either by decreasing recombination rate, for example, or randomly.

`stat_means_95CIs/gene_quantile_cMMb_density_mean_95CI_plot.R` plots density and means with 95% CIs of recombination rate for each group of genes.

`stat_means_95CIs/gene_quantile_popgenetics_stats_density_mean_95CI_plot.R` plots density and means with 95% CIs of a given population genetics statistic for each group of genes.

`GO_term_enrichment/topGO_gene_quantiles.R` and `GO_term_enrichment/topGOslim_gene_quantiles.R` use the Bioconductor package [topGO](https://bioconductor.org/packages/release/bioc/html/topGO.html) to evaluate groups of genes ("quantiles", defined by decreasing crossover recombination rate or ChIP-seq signal) for over-representation gene ontology (GO) and high-level GO (GO slim) terms, respectively.

`hypergeometric_tests/proportion_*_in_gene_quantiles_hypergeometricTest.R` scripts evaluate different gene categories for over- and under-representation in each group of genes by applying hypergeometric tests.
These scripts also sample from the hypergeometric distribution 100,000 times to obtain a probability distribution for plotting the log<sub>2</sub>(observed/expected) ratio and significance threshold.
