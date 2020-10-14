#!/applications/R/R-3.5.0/bin/Rscript

# For a given wheat subgenome or across all subgenomes,
# analyse genes in each grouping (quantile) for over-representation of
# high-level gene ontology (GO slim) terms relative to their representation among
# all genes in the given subgenome(s).
# Use the default algorithm implemented in topGO ("weight01") to identify
# over-represented GO slim terms, with Fisher's exact test P-values <=sigLevel

# This script is based on a post by Avril Coghlan:
# http://avrilomics.blogspot.co.uk/2015/07/using-topgo-to-test-for-go-term.html

# Usage:
# /applications/R/R-3.5.0/bin/Rscript ./topGOslim_gene_quantiles.R BP 0.05 'DMC1_in_genes' 1 4 'genomewide' 'Agenome_Bgenome_Dgenome'

#ont <- "BP"
#sigLevel <- 0.05
#quantileBy <- "DMC1_in_genes"
#quantileFirst <- 1
#quantileLast <- 4
#region <- "genomewide"
#genomeName <- "Agenome_Bgenome_Dgenome"

args <- commandArgs(trailingOnly = TRUE)
ont <- args[1]
sigLevel <- args[2]
quantileBy <- args[3]
quantileFirst <- as.integer(args[4])
quantileLast <- as.integer(args[5])
region <- args[6]
genomeName <- args[7]

suppressMessages(library(topGO))

targets <- lapply(seq_along(quantileFirst:quantileLast), function(x) {
  paste0("quantiles_by_", quantileBy, "/",
         "featureIDs_quantile", as.character(x),
         "_of_", as.character(quantileLast),
         "_by_", quantileBy, "_of_genes_in_",
         genomeName, "_", region, ".txt")
})

# Read in GO slim annotations for genes to define "gene universe".
# GO and GO slim annotations were obtained from "OntologiesForGenes.rds"
# (available at https://opendata.earlham.ac.uk/wheat/under_license/toronto/Ramirez-Gonzalez_etal_2018-06025-Transcriptome-Landscape/data/TablesForExploration/),
# after retaining only rows containing "IWGSC+Stress" or "slim_IWGSC+Stress",
# respectively, in the "ontology" column of this table
geneID2GO <- readMappings(file = paste0("/home/ajt200/analysis/wheat/annotation/RamirezGonzalez_2018_Science_GO_anno/RamirezGonzalez_2018_iwgsc_refseqv1.0_OntologiesForGenes_FunctionalAnnotation_HCgenes_in_",
                                        genomeName, "_", region,
                                        "_GOslim_IDs_no_chrUn.tsv"))
geneUniverse <- names(geneID2GO)

genesetGO <- function(target) {
  # Define list of genes of interest; file should contain a single column of gene identifiers
  genesOfInterest <- as.character(read.table(target)$V1)
  
  # Specify where genes of interest appear in the gene universe vector
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
  
  # Build GOdata object in topGO
  capture.output(GOdata <- new("topGOdata", description = "Grouped genes", ontology = ont,
                               allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO),
                 file="/dev/null")
  
  ## Access list of genes of interest
  #sg <- sigGenes(GOdata)
  #print(str(sg))
  #print(numSigGenes(GOdata))
  
  # Run Fisher's exact tests to determine GOslim term enrichment
  capture.output(resultClassic <- runTest(GOdata, algorithm = "classic", statistic = "fisher"),
                 file="/dev/null")
  capture.output(resultElim <- runTest(GOdata, algorithm = "elim", statistic = "fisher"),
                 file="/dev/null")
  capture.output(resultWeight <- runTest(GOdata, algorithm = "weight", statistic = "fisher"),
                 file="/dev/null")
  capture.output(resultTopGO <- runTest(GOdata, algorithm = "weight01", statistic = "fisher"),
                 file="/dev/null")
  
  # Count number of results where weight01 gives a P-value <= sigLevel
  mySummary <- summary(attributes(resultTopGO)$score <= as.numeric(sigLevel))
  if(length(mySummary) > 2) {
    numSignif <- as.integer(mySummary[[3]])
  } else {
    numSignif <- as.integer(0)
  }
  if(numSignif < 2) {
    numSignif <- as.integer(2)
  }

  # List significant results and write to file
  capture.output(enrichRes <- GenTable(GOdata,
                                       classicFisher = resultClassic,
                                       elimFisher = resultElim,
                                       weightFisher = resultWeight,
                                       topGOFisher = resultTopGO,
                                       orderBy = "topGOFisher",
                                       ranksOf = "weightFisher",
                                       topNodes = numSignif),
                 file="/dev/null")
  
  # WARNING: ASSUMES INPUT FILE HAS A 3-LETTER EXTENSION
  basename <- basename(target)
  len <- nchar(basename)
  basename <- substr(basename, 1, len-4)
  
  out_name <- paste(basename, "GOslim", ont, "enrichment.tsv", sep="_")
  folder <- paste0(dirname(target), "/GOslim")
  system(paste0("[ -d ", folder, " ] || mkdir ", folder))
  folder2 <- paste0(folder, "/", basename, "_GOslim_", ont)
  system(paste0("[ -d ", folder2, " ] || mkdir ", folder2))
  
  capture.output(write.table(enrichRes, file = file.path(folder, out_name), sep = "\t",
                             row.names = FALSE, col.names = TRUE, quote = FALSE),
                 file="/dev/null")
  
  ## Visualise the positions of the top 5 statistically significant GOslim terms in the GO hierarchy
  #out_name2 <- paste(basename, "GOslim", ont, "enrichment", sep="_")
  #printGraph(GOdata, resultTopGO, firstSigNodes = 5,
  #           fn.prefix = file.path(folder, out_name2), useInfo = "all", pdfSW = TRUE)
  
  # Extract gene IDs annotated with significantly enriched GOslim terms
  myTerms <- enrichRes$GO.ID
  myGenes <- genesInTerm(GOdata, myTerms)
  for(i in 1:length(myTerms)) {
    myTerm <- myTerms[i]
    myGenesForTerm <- myGenes[myTerm][[1]]
    myFactor <- myGenesForTerm %in% genesOfInterest
    myGenesForTermT <- myGenesForTerm[myFactor == TRUE]
    myGenesForTermT <- paste(myGenesForTermT, collapse = ",")
    myGenesForTermT <- paste(myTerm, myGenesForTermT, sep = "\t")
    out_name3 <- paste0(basename, "_GOslim_", ont, "_enrichment_", myTerm, ".txt")  
    write(myGenesForTermT, file = file.path(folder2, out_name3))
  }
}

# Apply genesetGO() function to each target group of genes
lapply(seq_along(targets), function(x) {
  genesetGO(target = targets[[x]])
})
