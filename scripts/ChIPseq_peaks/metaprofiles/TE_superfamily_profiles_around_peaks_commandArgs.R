#!/applications/R/R-3.4.0/bin/Rscript

# Profile TE frequency around compartmentalised peaks and random loci

# Wheat subgenome compartments (a.k.a. partitions):
# 1. R1 and R3 (distal or "euchromatin")
# 2. R2a and R2b (interstitial)
# 3. C (proximal)
# 4. heterochromatin (interstitial and proximal)
# 5. centromeres (defined by IWGSC (2018) Science 361 using CENH3 ChIP-seq data from Guo et al. (2016) PLOS Genet. 12)

# Usage:
# /applications/R/R-3.4.0/bin/Rscript ./TE_superfamily_profiles_around_peaks_commandArgs.R DMC1_Rep1_ChIP euchromatin A 400 2000 2kb 20

#libName <- "DMC1_Rep1_ChIP"
#region <- "euchromatin"
#genomeName <- "A"
#bodyLength <- 400
#flankSize <- 2000
#flankName <- "2kb"
#winSize <- 20

args <- commandArgs(trailingOnly = T)
libName <- args[1]
region <- args[2]
genomeName <- args[3]
bodyLength <- as.numeric(args[4])
flankSize <- as.numeric(args[5])
flankName <- as.character(args[6])
winSize <- as.numeric(args[7])

library(EnrichedHeatmap)
library(parallel)

matDir <- paste0("matrices/")
system(paste0("[ -d ", matDir, " ] || mkdir ", matDir))

# Genomic definitions
chrs <- as.vector(read.table("wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrStart <- c(rep(1, times = length(chrs)))
chrLens <- as.vector(read.table("wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table(paste0("chromosome_compartments/",
                                               "centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt"))[,2])
centromereEnd <- as.vector(read.table(paste0("chromosome_compartments/",
                                             "centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt"))[,3])
chrPartitions <- read.table(paste0("chromosome_compartments/",
                                   "chromosome_partitions_IWGSC_2018_Science_Table_S29.txt"),
                            header = TRUE)
genomeGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = chrStart,
                                     end = chrLens),
                    strand = "*")
genomeGR <- genomeGR[grep(genomeName,
                          seqnames(genomeGR))@values]

# Define region to be analysed
if(region == "euchromatin") {
  regionGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                      ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                                 chrPartitions$R2b_R3),
                                       end = c(chrPartitions$R1_R2a,
                                               chrLens)),
                      strand = "*")
  regionGR <- regionGR[grep(genomeName,
                            seqnames(regionGR))@values]
} else if(region == "interstitial") {
  regionGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                      ranges = IRanges(start = c(chrPartitions$R1_R2a+1,
                                                 chrPartitions$C_R2b),
                                       end = c(chrPartitions$R2a_C,
                                               chrPartitions$R2b_R3-1)),
                      strand = "*")
  regionGR <- regionGR[grep(genomeName,
                            seqnames(regionGR))@values]
} else if(region == "proximal") {
  regionGR <- GRanges(seqnames = chrPartitions$chrom,
                      ranges = IRanges(start = chrPartitions$R2a_C+1,
                                       end = chrPartitions$C_R2b-1),
                      strand = "*")
  regionGR <- regionGR[grep(genomeName,
                            seqnames(regionGR))@values]
} else if(region == "heterochromatin") {
  regionGR <- GRanges(seqnames = chrPartitions$chrom,
                      ranges = IRanges(start = chrPartitions$R1_R2a+1,
                                       end = chrPartitions$R2b_R3-1),
                      strand = "*")
  regionGR <- regionGR[grep(genomeName,
                            seqnames(regionGR))@values]
} else if(region == "centromeres") {
  regionGR <- GRanges(seqnames = chrs,
                      ranges = IRanges(start = centromereStart,
                                       end = centromereEnd),
                      strand = "*")
  regionGR <- regionGR[grep(genomeName,
                            seqnames(regionGR))@values]
} else if(region == "genomewide") {
  regionGR <- GRanges(seqnames = chrs,
                      ranges = IRanges(start = rep(1, length(chrs)),
                                       end = chrLens),
                      strand = "*")
  regionGR <- regionGR[grep(genomeName,
                            seqnames(regionGR))@values]
} else {
  stop("region is not euchromatin, interstitial, proximal, heterochromatin, centromeres, or genomewide")
}

# Define region to be masked out of analysis
if(region == "euchromatin") {
  maskGR <- GRanges(seqnames = chrPartitions$chrom,
                    ranges = IRanges(start = chrPartitions$R1_R2a+1,
                                     end = chrPartitions$R2b_R3-1),
                    strand = "*")
  maskGR <- maskGR[grep(genomeName,
                        seqnames(maskGR))@values]
} else if(region == "interstitial") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 3),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$R2a_C+1,
                                               chrPartitions$R2b_R3),
                                     end = c(chrPartitions$R1_R2a,
                                             chrPartitions$C_R2b-1,
                                             chrLens)),
                    strand = "*")
  maskGR <- maskGR[grep(genomeName,
                        seqnames(maskGR))@values]
} else if(region == "proximal") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$C_R2b),
                                     end = c(chrPartitions$R2a_C,
                                             chrLens)),
                    strand = "*")
  maskGR <- maskGR[grep(genomeName,
                        seqnames(maskGR))@values]
} else if(region == "heterochromatin") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$R2b_R3),
                                     end = c(chrPartitions$R1_R2a,
                                             chrLens)),
                    strand = "*")
  maskGR <- maskGR[grep(genomeName,
                        seqnames(maskGR))@values]
} else if(region == "centromeres") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               centromereEnd+1),
                                     end = c(centromereStart-1,
                                             chrLens)),
                    strand = "*")
  maskGR <- maskGR[grep(genomeName,
                        seqnames(maskGR))@values]
} else if(region == "genomewide") {
  maskGR <- GRanges()
  maskGR <- maskGR[grep(genomeName,
                        seqnames(maskGR))@values]
} else {
  stop("region is not euchromatin, interstitial, proximal, heterochromatin, centromeres, or genomewide")
}

# Load peaks in BED format and convert into GRanges
# Note addition of 1 to 0-based BED start coordinates
peaks <- read.table(paste0("/home/ajt200/analysis/wheat/DMC1/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.001_q0.01/",
                           libName,
                           "_rangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth_in_",
                           genomeName, "genome_", region, ".bed"),
                    header = F)
colnames(peaks) <- c("chr", "start", "end", "name", "score", "strand")
peaksGR <- GRanges(seqnames = peaks$chr,
                   ranges = IRanges(start = peaks$start+1,
                                    end = peaks$end),
                   strand = peaks$strand,
                   number = peaks$name)
peaksGR <- peaksGR[seqnames(peaksGR) != "chrUn"]
# Subset to include only those not overlapping masked region
mask_peaks_overlap <- findOverlaps(query = maskGR,
                                   subject = peaksGR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
if(length(mask_peaks_overlap) > 0) {
  peaksGR <- peaksGR[-subjectHits(mask_peaks_overlap)]
}
# Load ranLoc in BED format and convert into GRanges
# Note addition of 1 to 0-based BED start coordinates
ranLoc <- read.table(paste0("/home/ajt200/analysis/wheat/DMC1/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.001_q0.01/",
                            libName,
                            "_rangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth_in_",
                            genomeName, "genome_", region, "_randomLoci.bed"),
                     header = F)
colnames(ranLoc) <- c("chr", "start", "end", "name", "score", "strand")
ranLocGR <- GRanges(seqnames = ranLoc$chr,
                    ranges = IRanges(start = ranLoc$start+1,
                                     end = ranLoc$end),
                    strand = ranLoc$strand,
                    number = ranLoc$name)
ranLocGR <- ranLocGR[seqnames(ranLocGR) != "chrUn"]
# Subset to include only those not overlapping masked region
mask_ranLoc_overlap <- findOverlaps(query = maskGR,
                                    subject = ranLocGR,
                                    type = "any",
                                    select = "all",
                                    ignore.strand = TRUE)
if(length(mask_ranLoc_overlap) > 0) {
  ranLocGR <- ranLocGR[-subjectHits(mask_ranLoc_overlap)]
}

# Load TE superfamily BED files
superfamCode <- c("RLG",
                  "RLC",
                  "RLX",
                  "RIX",
                  "SIX",
                  "DTC",
                  "DTM",
                  "DTX",
                  "DTH",
                  "DMI",
                  "DTT",
                  "DXX",
                  "DTA",
                  "DHH",
                  "XXX")
superfamName <- c("Gypsy_LTR",
                  "Copia_LTR",
                  "Unclassified_LTR",
                  "LINE",
                  "SINE",
                  "CACTA",
                  "Mutator",
                  "Unclassified_with_TIRs",
                  "Harbinger",
                  "MITE",
                  "Mariner",
                  "Unclassified_class_2",
                  "hAT",
                  "Helitrons",
                  "Unclassified_repeats")

superfamListGR <- mclapply(seq_along(superfamName), function(x) {
  superfam <- read.table(paste0("iwgsc_refseqv1.0_TransposableElements_2017Mar13_superfamily_",
                                superfamName[x], "_", superfamCode[x], ".bed"),
                         header = F)
  colnames(superfam) <- c("chr", "start", "end", "name", "score", "strand")
  superfamGR <- GRanges(seqnames = superfam$chr,
                        ranges = IRanges(start = superfam$start+1,
                                         end = superfam$end),
                        strand = "*",
                        number = superfam$name,
                        coverage = rep(1, dim(superfam)[1]))
  superfamGR <- superfamGR[seqnames(superfamGR) != "chrUn"]
  # Subset to include only those not overlapping masked region
  mask_superfam_overlap <- findOverlaps(query = maskGR,
                                        subject = superfamGR,
                                        type = "any",
                                        select = "all",
                                        ignore.strand = TRUE)
  if(length(mask_superfam_overlap) > 0) {
    superfamGR <- superfamGR[-subjectHits(mask_superfam_overlap)]
  }
  superfamGR
}, mc.cores = length(superfamName))

# Define matrix and column mean outfiles
outDF <- lapply(seq_along(superfamName), function(x) {
  list(paste0(matDir, superfamName[x], "_", superfamCode[x],
              "_around_DMC1_peaks_in_", genomeName, "genome_", region,
              "_matrix_bin", winSize, "bp_flank", flankName, ".tab"),
       paste0(matDir, superfamName[x], "_", superfamCode[x],
              "_around_DMC1_peaks_in_", genomeName, "genome_", region,
              "_ranLoc_matrix_bin", winSize, "bp_flank", flankName, ".tab"))
})
outDFcolMeans <- lapply(seq_along(superfamName), function(x) {
  list(paste0(matDir, superfamName[x], "_", superfamCode[x],
              "_around_DMC1_peaks_in_", genomeName, "genome_", region,
              "_matrix_bin", winSize, "bp_flank", flankName, "_colMeans.tab"),
       paste0(matDir, superfamName[x], "_", superfamCode[x],
              "_around_DMC1_peaks_in_", genomeName, "genome_", region,
              "_ranLoc_matrix_bin", winSize, "bp_flank", flankName, "_colMeans.tab"))
})

# Function to create TE frequency matrices for
# feature loci and random loci (incl. flanking regions)
# and to calculate mean profiles across all feature loci and random loci
covMatrix <- function(signal,
                      feature,
                      ranLoc,
                      featureSize,
                      flankSize,
                      winSize,
                      outDF,
                      outDFcolMeans) {
  # feature loci
  set.seed(2840)
  feature_smoothed <- normalizeToMatrix(signal = signal,
                                        target = feature,
                                        value_column = "coverage",
                                        extend = flankSize,
                                        mean_mode = "w0",
                                        w = winSize,
                                        background = 0,
                                        smooth = TRUE,
                                        include_target = TRUE,
                                        target_ratio = featureSize/(featureSize+(flankSize*2)))
  print("feature_smoothed")
  print(feature_smoothed)
  print("feature_smoothed rows = ")
  print(length(feature_smoothed)/round((featureSize/winSize)+((flankSize*2)/winSize)))
  feature_smoothed_DF <- data.frame(feature_smoothed)
  feature_smoothed_DF_colMeans <- as.vector(colMeans(feature_smoothed_DF,
                                                     na.rm = T))
  write.table(feature_smoothed_DF,
              file = outDF[[1]],
              quote = F, sep = "\t", row.names = F, col.names = T)
  write.table(feature_smoothed_DF_colMeans,
              file = outDFcolMeans[[1]],
              quote = F, sep = "\t", row.names = F, col.names = T)

  # random loci
  set.seed(8472)
  ranLoc_smoothed <- normalizeToMatrix(signal = signal,
                                       target = ranLoc,
                                       value_column = "coverage",
                                       extend = flankSize,
                                       mean_mode = "w0",
                                       w = winSize,
                                       background = 0,
                                       smooth = TRUE,
                                       include_target = TRUE,
                                       target_ratio = featureSize/(featureSize+(flankSize*2)))
  print("ranLoc_smoothed")
  print(ranLoc_smoothed)
  print("ranLoc_smoothed rows = ")
  print(length(ranLoc_smoothed)/round((featureSize/winSize)+((flankSize*2)/winSize)))
  ranLoc_smoothed_DF <- data.frame(ranLoc_smoothed)
  ranLoc_smoothed_DF_colMeans <- as.vector(colMeans(ranLoc_smoothed_DF,
                                                    na.rm = T))
  write.table(ranLoc_smoothed_DF,
              file = outDF[[2]],
              quote = F, sep = "\t", row.names = F, col.names = T)
  write.table(ranLoc_smoothed_DF_colMeans,
              file = outDFcolMeans[[2]],
              quote = F, sep = "\t", row.names = F, col.names = T)
}

# Run covMatrix() function on each feature GRanges object to obtain matrices
# containing normalised feature frequency values around target and random loci
mclapply(seq_along(superfamName), function(x) {
  covMatrix(signal = superfamListGR[[x]],
            feature = peaksGR,
            ranLoc = ranLocGR,
            featureSize = bodyLength,
            flankSize = flankSize,
            winSize = winSize,
            outDF = outDF[[x]],
            outDFcolMeans = outDFcolMeans[[x]])
  print(paste0(superfamName[x], "_", superfamCode[x],
               "_around_DMC1_peaks_in_", genomeName, "genome_", region,
               " profile calculation complete"))
}, mc.cores = length(superfamName))
