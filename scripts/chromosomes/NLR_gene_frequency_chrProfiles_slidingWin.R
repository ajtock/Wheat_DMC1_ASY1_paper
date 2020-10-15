#!/applications/R/R-3.5.0/bin/Rscript

# Calculate NLR-encoding gene frequency in
# sliding windows along each wheat chromosome

# Usage:
# ./NLR_gene_frequency_chrProfiles_slidingWin.R 10Mb 10000000 1Mb 1000000

#winName <- "10Mb"
#winSize <- 10000000
#stepName <- "1Mb"
#stepSize <- 1000000

args <- commandArgs(trailingOnly = T)
winName <- args[1]
winSize <- as.numeric(args[2])
stepName <- args[3]
stepSize <- as.numeric(args[4])

library(rtracklayer)
library(parallel)
library(GenomicRanges)

inDir <- "/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/"
outDir <- "/home/ajt200/analysis/wheat/chromosomeProfiles/genes/"
genesGFF <- readGFF(paste0(inDir, "NLRs_Steuernagel_Wulff_2020_Plant_Physiol/NLR_genes_complete_representative_mRNA.gff3"))

# Genomic definitions
chrs <- as.vector(read.table("wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrLens <- as.vector(read.table("wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]

# Define windows as GRanges object
windowsGR <- GRanges()
for(i in 1:length(chrs)) {
  # Define sliding windows of width winSize nt,
  # with a step of stepSize nt
  winStarts <- seq(from = 1,
                   to = chrLens[i]-winSize,
                   by = stepSize)
  if(chrLens[i] - winStarts[length(winStarts)] >= winSize) {
    winStarts <- c(winStarts,
                   winStarts[length(winStarts)]+stepSize)
  }
  winEnds <- seq(from = winStarts[1]+winSize-1,
                 to = chrLens[i],
                 by = stepSize)
  winEnds <- c(winEnds,
               rep(chrLens[i], times = length(winStarts)-length(winEnds)))
  stopifnot(length(winStarts) == length(winEnds))

  chrWindowsGR <- GRanges(seqnames = chrs[i],
                          ranges = IRanges(start = winStarts,
                                           end = winEnds),
                          strand = "*")
  print(chrWindowsGR)
  try( if( end(chrWindowsGR)[length(chrWindowsGR)] != chrLens[i] )
       {
         stop("End coordinate of last window != chromosome length")
       }
     )
  windowsGR <- append(windowsGR, chrWindowsGR)
}

featureProfile <- NULL
for(i in 1:length(chrs)) {
  # Count genes within windows
  chrgenes <- genesGFF[genesGFF$seqid == chrs[i],]
  chrgenesGR <- GRanges(seqnames = chrs[i],
                        ranges = IRanges(start = chrgenes$start,
                                         end = chrgenes$end),
                        strand = chrgenes$strand)
  chrWindowsGR <- windowsGR[seqnames(windowsGR) == chrs[i]]
  wingenes <- countOverlaps(chrWindowsGR,
                            chrgenesGR,
                            ignore.strand = T)
  chrProfile <- data.frame(chr = as.character(chrs[i]),
                           window = as.integer(start(chrWindowsGR)),
                           features = as.integer(wingenes))
  featureProfile <- rbind(featureProfile, chrProfile)
}
write.table(featureProfile,
            file = paste0(outDir,
                          "NLR_gene_frequency_per_", winName,
                          "_step", stepName, ".txt"),
            sep = "\t", quote = F, row.names = F, col.names = T)
