#!/applications/R/R-3.3.2/bin/Rscript

# Convert narrowPeak into GRanges object (with overlapping peaks merged)
# to be provided as an input file for genomic feature overlap analyses
# with the Bioconductor package regioneR

# Usage:
# ./4_convert_narrowPeakToGRanges_for_regioneR.R DMC1_Rep1_ChIP 0.001 0.01

args <- commandArgs(trailingOnly = TRUE)
ChIPLibName <- args[1]
pval <- as.character(args[2])
qval <- as.character(args[3])

library(GenomicRanges)

rangerPeaks <- read.table(paste0(ChIPLibName,
                                 "_peaks_peakranger_ranger_p", pval,
                                 "_q", qval, "_TreadsNormCreads.narrowPeak"))
rangerPeaks <- cbind(rangerPeaks[,1:3],
                     rangerPeaks[,7:10])
colnames(rangerPeaks) <- c("chr", "start0based", "end",
                           "sigval", "pval", "qval", "summit0based")
rangerPeaks <- data.frame(chr          = as.character(rangerPeaks$chr),
                          start        = as.integer(rangerPeaks$start0based+1),
                          end          = as.integer(rangerPeaks$end),
                          sigval       = as.numeric(rangerPeaks$sigval),
                          pval         = as.numeric(rangerPeaks$pval),
                          qval         = as.numeric(rangerPeaks$qval),
                          summit0based = as.integer(rangerPeaks$summit0based))

# Create GRanges objects sorted by decreasing -log10(qval)
rangerPeaksGR <- sort(GRanges(seqnames = rangerPeaks$chr,
                              ranges = IRanges(start = rangerPeaks$start,
                                               end = rangerPeaks$end),
                              strand = "*",
                              sigval = rangerPeaks$sigval,
                              pval = rangerPeaks$pval,
                              qval = rangerPeaks$qval,
                              summit0based = rangerPeaks$summit0based),
                      by = ~ qval, decreasing = T)
# Merge overlapping peaks
rangerPeaksGRmergedOverlaps <- reduce(rangerPeaksGR)
save(rangerPeaksGRmergedOverlaps,
     file = paste0(ChIPLibName,
                   "_rangerPeaksGRmergedOverlaps_minuslog10_p", pval, "_q", qval, "_noMinWidth.RData"))
