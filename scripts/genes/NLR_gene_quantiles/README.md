# NLR-encoding gene quantiles

This subdirectory contains the R script `group_NLR_genes_into_quantiles.R` for dividing NLR-encoding genes into groups corresponding to those in given percentile ranges (e.g., 100th–75th (Quantile 1), 75th–50th (Quantile 2), 50th–25th (Quantile 3) and 25th–0th (Quantile 4)) with regard to different ordering factors, such as mean crossover recombination rate (cM/Mb) values, derived from a Chinese Spring × Renan genetic map, ChIP-seq signal, or physical cluster size.

`metaprofiles/NLR_gene_quantile_metaprofiles.R` calculates and plots metaprofiles of ChIP-seq, MNase-seq and RNA-seq signals, DNA methylation proportions, and SNP and transposon frequencies (gene windowed means and 95% confidence intervals, CIs) for each group of NLR genes, defined either by decreasing recombination rate, for example, or randomly.

`stat_means_95CIs/NLR_gene_quantile_cMMb_density_mean_95CI_plot.R` and `stat_means_95CIs/NLR_gene_quantile_cluster_size_density_mean_95CI_plot.R` plot density and means with 95% CIs of recombination rate and physical cluster size, respectively, for each group of NLR genes. 

`hypergeometric_tests/proportion_*_in_NLR_quantiles_hypergeometricTest.R` scripts evaluate different NLR gene categories for over- and under-representation in each group of NLR genes by applying hypergeometric tests.
These scripts also sample from the hypergeometric distribution 100,000 times to obtain a probability distribution for plotting the log<sub>2</sub>(observed/expected) ratio and significance threshold.
