#!/applications/R/R-3.5.0/bin/Rscript

# Divide ChIP-seq peaks into four groups corresponding to those within the
# 100th--75th (Quantile 1), 75th--50th (Quantile 2),
# 50th--25th (Quantile 3) and 25th--0th (Quantile 4) percentiles
# with regard to their mean crossover recombination rate (cM/Mb) values,
# derived from a Chinese Spring x Renan genetic map.
#
# This script extends ChIP-seq peak boundaries by 1 kb on each side,
# and computes mean cM/Mb values within these intervals using previously
# calculated mean recombination rates in 10-Mb sliding windows with a 1-Mb step
# (iwgsc_refseqv1.0_recombination_rate.txt), available as part of the IWGSC RefSeq v1.0 annotation:
# https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.0/
#
# It requires peak coordinates and, separately, random locus coordinates in BED6 format:
# column 1 = chromosome ID; column 2 = 0-based start coordinate; column 3 = 1-based end coordinates;
# column 4 = sequential or otherwise unique numbers (this will speed up computation; see comment from dpryan79 on 13/09/2018 under GitHub issue [computeMatrix has problem with multi processors #760](https://github.com/deeptools/deepTools/issues/760)]);
# column 5 = fill with NA; column 6 = fill with \*.

# Extract and save feature IDs for each quantile for further analyses
# (e.g., mean and 95% CI metaprofile calculation and plotting).

# Usage:
# /applications/R/R-3.5.0/bin/Rscript group_peaks_into_cMMb_quantiles.R DMC1_Rep1_ChIP DMC1 'Agenome_euchromatin,Bgenome_euchromatin,Dgenome_euchromatin' 4

#libName <- "DMC1_Rep1_ChIP"
#dirName <- "DMC1"
#featureName <- unlist(strsplit("Agenome_euchromatin,Bgenome_euchromatin,Dgenome_euchromatin",
#                               split = ","))
#quantiles <- 4

args <- commandArgs(trailingOnly = T)
libName <- args[1]
dirName <- args[2]
featureName <- unlist(strsplit(args[3],
                               split = ","))
quantiles <- as.numeric(args[4])

library(GenomicRanges)
library(dplyr)
library(parallel)

# Load features in BED format and convert into GRanges
# Note addition of 1 to 0-based BED start coordinates
features <- lapply(seq_along(featureName), function(y) {
  tmp <- read.table(paste0("/home/ajt200/analysis/wheat/", dirName,
                           "/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.001_q0.01/",
                           libName,
                           "_rangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth_in_",
                           featureName[y], ".bed"),
                    header = F)
  data.frame(tmp,
             V7 = paste0(featureName[y], "_", tmp$V4),
             stringsAsFactors = F) 
})
# If features from multiple subgenomes and/or compartments are to be analysed,
# concatenate the corresponding feature data.frames
if(length(featureName) > 1) {
  features <- do.call(rbind, features)
} else {
  features <- features[[1]]
}
colnames(features) <- c("chr", "start", "end", "name", "score", "strand", "featureID")
featuresGR <- GRanges(seqnames = features$chr,
                      ranges = IRanges(start = features$start+1,
                                       end = features$end),
                      strand = features$strand,
                      featureID = features$featureID)
# Extend feature boundaries to include 1000 bp upstream
# and downstream for calculation of mean cM/Mb around loci
featuresGR_ext <- GRanges(seqnames = seqnames(featuresGR),
                          ranges = IRanges(start = start(featuresGR)-1000,
                                           end = end(featuresGR)+1000),
                          strand = strand(featuresGR),
                          featureID = featuresGR$featureID)
print(featuresGR_ext)

# Load ranLocs in BED format and convert into GRanges
# Note addition of 1 to 0-based BED start coordinates
ranLocs <- lapply(seq_along(featureName), function(y) {
  tmp <- read.table(paste0("/home/ajt200/analysis/wheat/", dirName,
                           "/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.001_q0.01/",
                           libName,
                           "_rangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth_in_",
                           featureName[y], "_randomLoci.bed"),
                    header = F)
  data.frame(tmp,
             V7 = paste0(featureName[y], "_", tmp$V4),
             stringsAsFactors = F) 
})
# If ranLocs from multiple subgenomes and/or compartments are to be analysed,
# concatenate the corresponding ranLoc data.frames
if(length(featureName) > 1) {
  ranLocs <- do.call(rbind, ranLocs)
} else {
  ranLocs <- ranLocs[[1]]
}
colnames(ranLocs) <- c("chr", "start", "end", "name", "score", "strand", "ranLocID")
ranLocsGR <- GRanges(seqnames = ranLocs$chr,
                     ranges = IRanges(start = ranLocs$start+1,
                                      end = ranLocs$end),
                     strand = ranLocs$strand,
                     ranLocID = ranLocs$ranLocID)
# Extend ranLoc boundaries to include 1000 bp upstream
# and downstream for calculation of mean cM/Mb around loci
ranLocsGR_ext <- GRanges(seqnames = seqnames(ranLocsGR),
                         ranges = IRanges(start = start(ranLocsGR)-1000,
                                          end = end(ranLocsGR)+1000),
                         strand = strand(ranLocsGR),
                         ranLocID = ranLocsGR$ranLocID)
print(ranLocsGR_ext)

# Convert windowed recombination rate into GRanges
cMMb <- read.table("iwgsc_refseqv1.0_recombination_rate.txt",
                   header = T)
cMMbGR <- GRanges(seqnames = cMMb$chromosome,
                  ranges = IRanges(start = cMMb$intervalStart,
                                   end = cMMb$intervalEnd),
                  strand = "*",
                  cMMb = cMMb$recombinationRate)

# Obtain cMMb values for each feature extended by 1000 bp on each side
# Where features overlap more than one winName window, calculate mean cMMb
feature_cMMb_overlaps <- findOverlaps(query = featuresGR_ext,
                                      subject = cMMbGR,
                                      type = "any",
                                      select = "all",
                                      ignore.strand = TRUE)
feature_cMMb_overlapsList <- lapply(seq_along(featuresGR_ext), function(x) {
  subjectHits(feature_cMMb_overlaps)[queryHits(feature_cMMb_overlaps) == x]
})
feature_cMMb <- sapply(feature_cMMb_overlapsList,
                       function(x) mean(cMMbGR$cMMb[x], na.rm = TRUE))
featuresGR <- GRanges(featuresGR,
                      featureID = featuresGR$featureID,
                      cMMb = feature_cMMb)

# Obtain cMMb values for each ranLoc extended by 1000 bp on each side
# Where ranLocs overlap more than one winName window, calculate mean cMMb
ranLoc_cMMb_overlaps <- findOverlaps(query = ranLocsGR_ext,
                                     subject = cMMbGR,
                                     type = "any",
                                     select = "all",
                                     ignore.strand = TRUE)
ranLoc_cMMb_overlapsList <- lapply(seq_along(ranLocsGR_ext), function(x) {
  subjectHits(ranLoc_cMMb_overlaps)[queryHits(ranLoc_cMMb_overlaps) == x]
})
ranLoc_cMMb <- sapply(ranLoc_cMMb_overlapsList,
                      function(x) mean(cMMbGR$cMMb[x], na.rm = TRUE))
ranLocsGR <- GRanges(ranLocsGR,
                     ranLocID = ranLocsGR$ranLocID,
                     cMMb = ranLoc_cMMb)


# Define orderingFactor to be used for grouping features into quantiles
orderingFactor <- "cMMb"
outDir <- paste0("quantiles_by_", orderingFactor, "/")
plotDir_list <- lapply(seq_along(outDir), function(w) {
  paste0(outDir[w], "plots/")
})

sapply(seq_along(outDir), function(w) {
 system(paste0("[ -d ", outDir[w], " ] || mkdir ", outDir[w]))
})
sapply(seq_along(outDir), function(w) {
  system(paste0("[ -d ", plotDir_list[[w]], " ] || mkdir ", plotDir_list[[w]]))
})

# For each population, divide features into quantiles based on decreasing orderingFactor
# all subgenomes
featuresDF <- data.frame(featuresGR,
                         quantile = as.character(""),
                         stringsAsFactors = F)
mclapply(seq_along(orderingFactor), function(w) {
  print(orderingFactor[w])
  # Assign 0s to NA values only for coverage data
  if(grepl("_in_", orderingFactor[w])) {
    featuresDF[,which(colnames(featuresDF) == orderingFactor[w])][which(is.na(featuresDF[,which(colnames(featuresDF) == orderingFactor[w])]))] <- 0
  }
  if(grepl("HudsonRM", orderingFactor[w])) {
    quantiles <- 2
  }
  quantilesStats <- data.frame()
  for(k in 1:quantiles) {
    # First quantile should span 1 to greater than, e.g., 0.75 proportions of features
    if(k < quantiles) {
      featuresDF[ !is.na(featuresDF[,which(colnames(featuresDF) == orderingFactor[w])]) &
                  percent_rank(featuresDF[,which(colnames(featuresDF) == orderingFactor[w])]) <= 1-((k-1)/quantiles) &
                  percent_rank(featuresDF[,which(colnames(featuresDF) == orderingFactor[w])]) >  1-(k/quantiles), ]$quantile <- paste0("Quantile ", k)
    } else {
    # Final quantile should span 0 to, e.g., 0.25 proportions of features
      featuresDF[ !is.na(featuresDF[,which(colnames(featuresDF) == orderingFactor[w])]) &
                  percent_rank(featuresDF[,which(colnames(featuresDF) == orderingFactor[w])]) <= 1-((k-1)/quantiles) &
                  percent_rank(featuresDF[,which(colnames(featuresDF) == orderingFactor[w])]) >= 1-(k/quantiles), ]$quantile <- paste0("Quantile ", k)
    }
    write.table(featuresDF[featuresDF$quantile == paste0("Quantile ", k),],
                file = paste0(outDir[w],
                              "quantile", k, "_of_", quantiles,
                              "_by_", orderingFactor[w],
                              "_of_", libName, "_peaks_in_",
                              paste0(featureName,
                                     collapse = "_"), ".txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)
    stats <- data.frame(quantile = as.integer(k),
                        n = as.integer(dim(featuresDF[featuresDF$quantile == paste0("Quantile ", k),])[1]),
                        mean_width = as.integer(round(mean(featuresDF[featuresDF$quantile == paste0("Quantile ", k),]$width, na.rm = T))),
                        total_width = as.integer(sum(featuresDF[featuresDF$quantile == paste0("Quantile ", k),]$width, na.rm = T)),
                        mean_orderingFactor = as.numeric(mean(featuresDF[featuresDF$quantile == paste0("Quantile ", k),][,which(colnames(featuresDF) == orderingFactor[w])], na.rm = T)))
    quantilesStats <- rbind(quantilesStats, stats)
  }
  write.table(quantilesStats,
              file = paste0(outDir[w],
                            "summary_", quantiles, "quantiles",
                             "_by_", orderingFactor[w],
                             "_of_", libName, "_peaks_in_",
                             paste0(featureName,
                                    collapse = "_"), ".txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(featuresDF,
              file = paste0(outDir[w],
                            "features_", quantiles, "quantiles",
                             "_by_", orderingFactor[w],
                             "_of_", libName, "_peaks_in_",
                             paste0(featureName,
                                    collapse = "_"), ".txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)

  # Divide ranLocs into quantiles based on feature quantile indices
  ranLocsDF <- data.frame(ranLocsGR,
                          random = as.character(""),
                          stringsAsFactors = F)
  # Get row indices for each feature quantile
  quantileIndices <- lapply(1:quantiles, function(k) {
    which(featuresDF$quantile == paste0("Quantile ", k))
  })
  for(k in 1:quantiles) {
    ranLocsDF[quantileIndices[[k]],]$random <- paste0("Random ", k)
  }
  write.table(ranLocsDF,
              file = paste0(outDir[w],
                            "features_", quantiles, "quantiles",
                             "_by_", orderingFactor[w],
                             "_of_", libName, "_peaks_in_",
                             paste0(featureName,
                                    collapse = "_"), "_ranLocs.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)
}, mc.cores = length(orderingFactor), mc.preschedule = F)
