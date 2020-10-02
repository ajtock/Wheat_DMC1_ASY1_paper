#!/applications/R/R-3.3.2/bin/Rscript

# -log10 transform p-values and q-values for consistency with narrowPeak format conventions,
# and to be used for sorting by decreasing values for downstream motif enrichment analyses
# Define narrowPeak "signalValue" (column 7) as:
# signalValue = log2((region ChIP reads+1)/(region control reads+1))
# However, refrain from using this signalValue as it is unknown whether these read counts
# are normalized by library size by the PeakRanger ranger tool

# Usage:
# ./3_narrowPeak_minuslog10PQ.R DMC1_Rep1_ChIP 0.001 0.01

args <- commandArgs(trailingOnly = TRUE)
ChIPLibName <- args[1]
pval <- as.character(args[2])
qval <- as.character(args[3])

peaks <- read.table(paste0(ChIPLibName, "_peaks_peakranger_ranger_p",
                           pval, "_q", qval,
                           "_Treads_Creads.narrowPeak.UntransformedPQ"))
colnames(peaks) <- c("chr", "start0based", "end",
                     "name", "score", "strand",
                     "signalVal", "pValUntrans", "qValUntrans",
                     "summit0based", "treads", "creads")
peaks <- cbind(peaks[,1:6], log2((peaks[,11]+1)/(peaks[,12]+1)),
               -log10(peaks[,8]), -log10(peaks[,9]), peaks[,10])
colnames(peaks) <- c("chr", "start0based", "end",
                     "name", "score", "strand",
                     "TreadsNormCreads", "pVal", "qVal",
                     "summit0based")
peaks$pVal[which(!is.finite(peaks$pVal))] <- 323
peaks$qVal[which(!is.finite(peaks$qVal))] <- 323
head(peaks)
write.table(peaks, file = paste0(ChIPLibName, "_peaks_peakranger_ranger_p",
                                 pval, "_q", qval,
                                 "_log2TreadsNormCreads.narrowPeak"),
            col.names = F, row.names = F, sep = "\t", quote = F)
