# This R pipeline was used to call tag clusters (TCs) in TSS-Seq data using CAGEfightR
# Differentially expressed TCs were found using DESeq2 and then exported to BED files

# Install CAGEfightR from GitHub:
library(devtools)
install_github("MalteThodberg/CAGEfightR")
packageVersion("CAGEfightR")		# 0.99.0

# Load the required libraries:
library(CAGEfightR)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
library(BiocParallel)
register(MulticoreParam(4), default=TRUE) # On Windows use snowParam()
library(BSgenome.Athaliana.TAIR.TAIR9)
library(tibble)
library(dplyr)
library(edgeR)
library(DESeq2)
source("assignTxType_custom.R")

# Provide CAGEfightR with BigWig files from 5Cap-Seq experiment:
setwd("path/to/expanded/bigwig/files")
bw_plus_filenames <- list.files(".", pattern="fw_cov2.bw$")
bw_minus_filenames <- list.files(".", pattern="rev_cov2.bw$")
bw_plus <- BigWigFileList(bw_plus_filenames)
bw_minus <- BigWigFileList(bw_minus_filenames)
sample_names <- sub('_fw_cov2.bw', '', bw_plus_filenames)
names(bw_plus) <- sample_names
names(bw_minus) <- sample_names

# Make the appropriate Seqinfo object and the design matrix:
tair_seqinfo <- seqinfo(Athaliana)
seqlevels(tair_seqinfo) <- c("1", "2", "3", "4", "5", "Mt", "Pt")
design_matrix <- data.frame("Name"=sample_names, "BigWigPlus"=bw_plus_filenames, "BigWigMinus"=bw_minus_filenames, 
row.names=sample_names, genotype=factor(rep(c("wt", "hen2"), 4), levels=c("wt", "hen2")), 
temp=factor(rep(rep(c("rt", "cold"), each=2), 2), levels=c("rt", "cold")))

# Quantify all tag clusters (TCs):
ctss <- quantifyCTSSs(plusStrand=bw_plus, minusStrand=bw_minus, design=design_matrix, genome=tair_seqinfo)

# Call candidate TSS:
tss <- quickTSSs(ctss)

# Call candidate enhancers:
enhancers <- quickEnhancers(ctss)

# Annotate TSS and enhancers by genomic features (observe that a custom assignTxType() function was used):
rowRanges(tss)$txType <- suppressWarnings(assignTxType_custom(rowRanges(tss), txdb=txdb, asFactor=TRUE))
rowRanges(enhancers)$txType <- suppressWarnings(assignTxType_custom(rowRanges(enhancers), txdb=txdb, asFactor=TRUE))

# Remove enhancers overlapping known transcripts:
enhancers <- subset(enhancers, txType %in% c("intergenic", "intron"))

# Combine candidate TSS and enhancers into a single RangedSummarizedExperiment object:
rowRanges(tss)$clusterType <- "TSS"
rowRanges(enhancers)$clusterType <- "enhancer"
rse <- combineClusters(tss, enhancers, removeIfOverlapping="object1")

# Remove low expressed TCs:
rse <- subsetBySupport(rse, inputAssay = "counts", outputColumn = "support", unexpressed = 0, minSamples = 1)

# Annotate TCs by gene IDs:
rse <- suppressWarnings(assignGeneID(rse, geneModels=txdb))

# Annotate TCs by gene names:
tair_ann <- import.gff3('Arabidopsis_thaliana.TAIR10.26.gff3')
mapping <- data.frame('geneID' = sub('gene:', '', tair_ann$ID), 'name'=tair_ann$external_name)
tmp <- left_join(x=tibble(geneID=rowRanges(rse)$geneID), y=mapping, by="geneID", na_matches = "never")
mcols(rse) <- DataFrame(mcols(rse), tmp[,-1])
rm(tmp)

# For DE calling consider only strong peaks (TPM >= 1 in at least 2 samples):
rse2 <- rse[rowSums(cpm(assay(rse)) >= 1) >= 2]

# Use DESeq2:
dds <- DESeqDataSet(rse2, design = ~ genotype + temp + genotype:temp)
dds <- DESeq(dds)
resultsNames(dds)	# c("Intercept", "genotype_hen2_vs_wt", "temp_cold_vs_rt", "genotypehen2.tempcold")

# Extract DE results:
my_names <- c("hen2", "cold", "interaction")
res <- vector("list", 3)
for (i in 1:3) {
  contr <- rep(0, 4)
  contr[[i+1]] <- 1
  de <- results(dds, contrast=contr)
  df <- as.data.frame(de)[,c(2,6)]
  colnames(df) <- paste0(colnames(df), "_", my_names[[i]])
  de_lfc1 <- results(dds, contrast=contr, lfcThreshold=1)
  df_lfc1 <- as.data.frame(de_lfc1)[,c(2,6)]
  colnames(df_lfc1) <- paste0(colnames(df_lfc1), "_", my_names[[i]], "_lfc")
  out <- cbind(df, df_lfc1)
  out <- out[,c(1,2,4)]
  colnames(out) <- gsub("FoldChange", "FC", colnames(out))
  res[[i]] <- out
}
m <- cbind(res[[1]], res[[2]], res[[3]])
m[is.na(m)] <- 1 # replace padj=NA with padj=1

# Expand DE results to the original TC number:
m$key <- rownames(m)
orig <- data.frame("key"=names(rowRanges(rse)))
n <- left_join(orig, m, by=c("key"), all.x=TRUE)

# Add DE results to the RSE object:
mcols(rse) <- cbind(mcols(rse), n[-1])
saveRDS(rse, "rse.RData")

# Generate BED files (to highlight differentially expressed tag clusters in genome browsers):
for (i in 1:3) {
  name <- my_names[[i]]
  colname <- paste0("padj_", name)
  colname_lfc <- paste0("padj_", name, "_lfc")
  pvals <- mcols(rse)[, colname]
  pvals_lfc <- mcols(rse)[, colname_lfc]
  gr <- rowRanges(rse)[!is.na(pvals) & pvals < 0.05]
  gr_lfc <- rowRanges(rse)[!is.na(pvals_lfc) & pvals_lfc < 0.05]
  rtracklayer::export(granges(gr), paste0(name, "_DE.bed"))
  rtracklayer::export(granges(gr_lfc), paste0(name, "_DE_lfc1.bed"))
}


