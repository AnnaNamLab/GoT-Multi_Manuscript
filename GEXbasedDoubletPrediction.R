##################################################################
## Identify doublets based on aberrant gene expression profiles ##
##################################################################

# This code aims at identifying doublets based on aberrant gene expression profiles using two independent methods (DoubletFinder & scDblFinder) and retaining only the cells called as "doublets" by both methods as potential doublets.

##################################################################
##       Step 1. Identifying doublets using DoubletFinder       ##
##################################################################

# Define variables
Variables <- list(
  # Doublet formation rate from 10X user guide
  pN = 0.075, 
  
  # number of PCs to keep for clustering analysis 
  ## choose final after seeing the elbow plot
  pca.dim = 20,
  
  # Path to seurat object after cell and gene filtering (see CellAndGeneFiltering.R code, UMI counts as input)
  InputPath = "SeuratObject_CellLine.rds",
  
  # prefix to add to output file names
  OutputName = "_CellLine"
)

# Load required packages
library(Seurat)
library(tidyverse)
library(ggplot2)
library(DoubletFinder)

# Load the Seurat object containing the GEX data (after filtering genes and cells)
res.data <- readRDS(Variables$InputPath)

# Pre-process using standard workflow
res.data <- NormalizeData(object = res.data) ## normalize data
res.data <- FindVariableFeatures(object = res.data) ## find most variable genes
res.data <- ScaleData(object = res.data) ## sacle data
res.data <-  RunPCA(object = res.data) # run linear dimension reduction

pdf(paste0(Sys.Date(), "_ElbowPlot.pdf"))
print(ElbowPlot(res.data))
dev.off()

res.data <-  FindNeighbors(object = res.data, dims = 1:Variables$pca.dim) ## Computes nearest neighbor graph and SNN
res.data <-  FindClusters(object = res.data) ## Run cell clustering analysis
res.data <-  RunUMAP(object = res.data, dims = 1:Variables$pca.dim) ## Run non-linear dimensional reduction

# pK Identification
sweep.res.list <- paramSweep_v3(res.data, PCs = 1:Variables$pca.dim, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
p <-   ggplot(bcmvn, aes(pK, BCmetric, group = 1)) + 
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(paste0(Sys.Date(), "_pKIdentificationPlot.pdf"))
print(p)
dev.off()

pK <- bcmvn %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

# Estimate expected number of doublet
## Exp = expected number of doublets. This info can be obtained from the 10X user guide.
## Homotypic Doublet Proportion Estimate 
annotations <- res.data@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)         
nExp_poi <- round(Variables$pN * nrow(res.data@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder with varying classification stringencies -
## pN = the percentage of artificial doublets. It doesn't affect performance of DoubletFinder. It is set to a default of 0.25 (25%).
seurat.obj <- doubletFinder_v3(res.data, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)

# Visualize doublets
seurat.obj@meta.data$DoubletStatus <- seurat.obj@meta.data[,ncol(seurat.obj@meta.data)]
pdf(paste0(Sys.Date(), "_UMAPPlot_HighlightingDoubletSingletIdentity.pdf"), width = 5, height = 5)
DimPlot(seurat.obj, reduction = "umap", group.by = "DoubletStatus")
dev.off()

pdf(file = paste0(Sys.Date(), "_Number of doublets.pdf"))
ggplot(seurat.obj@meta.data, aes(x=orig.ident, fill=DoubletStatus)) + geom_bar(position="dodge")
dev.off()

# save results
write.csv(table(seurat.obj@meta.data$DoubletStatus), paste0(Sys.Date(), "_NumberOfDoubletsAndSingletsEstimated.csv"))
write.csv(seurat.obj@meta.data, paste0(Sys.Date(), Variables$OutputName, "_metadataWithDoubletResult_DoubletFinder.csv"))

rm(list = ls()); gc()

##################################################################
##       Step 2. Identifying doublets using scDblFinder         ##
##################################################################

# Define variables
Variables <- list(
  # Doublet formation rate from 10X user guide
  pN = 0.075, 
  
  # number of PCs to keep for clustering analysis 
  ## choose final after seeing the elbow plot
  pca.dim = 20,
  
  # Path to seurat object after cell and gene filtering (see CellAndGeneFiltering.R code, UMI counts as input)
  InputPath = "SeuratObject_CellLine.rds",
  
  # prefix to add to output file names
  OutputName = "_CellLine"
)

# Load required packages
library(scDblFinder)
library(Seurat)
library(BiocParallel)
library(ggplot2)

# Import data
ExpData <- readRDS(Variables$InputPath)

# Convert Seurat object to a Single Cell Experiment object
ExpData.sce <- as.SingleCellExperiment(ExpData)

# Identify doublets per sample
doublet.res <- scDblFinder(ExpData.sce, samples="orig.ident", BPPARAM=MulticoreParam(4), dbr.sd=0.015)

# Extract results to evaluate doublet score and class for each cell
metadata <- as.data.frame(doublet.res@colData@listData)
rownames(metadata) <- doublet.res@colData@rownames
table(metadata$scDblFinder.class)

# Visualize results
pdf(file = paste0(Sys.Date(), "_Number of doublets.pdf"))
ggplot(metadata, aes(x=orig.ident, fill=scDblFinder.class)) + geom_bar(position="dodge")
dev.off()

# Export cells' metadata
write.csv(metadata, paste0(Sys.Date(), "_CellMetadata_WithDoubletIdentification_scDblFinder.csv"))

rm(list = ls()); gc()

##################################################################
##       Step 3. Identifying doublets common to both tools      ##
##################################################################

# Import results from DoubletFinder
temp <- list.files(pattern = "metadataWithDoubletResult_DoubletFinder")
res1 <- read.csv(temp, header = TRUE, row.names = 1,sep=",")

# Import results from scDblFinder
temp <- list.files(pattern = "WithDoubletIdentification_scDblFinder")
res2 <- read.csv(temp, header = TRUE, row.names = 1,sep=",")
res2 <- res2[rownames(res1),]

# Check number of doublets identified by each method
Doublet.res <- data.frame(cellID = rownames(res1), Doublet_DoubletFinder = res1$DoubletStatus, Doublet_scDblFinder = res2$scDblFinder.class)
table(Doublet.res$Doublet_DoubletFinder, Doublet.res$Doublet_scDblFinder)
table(Doublet.res$Doublet_DoubletFinder)
table(Doublet.res$Doublet_scDblFinder)

# Keep only cells predicted to be doublets by both methods as potential doublets
for(i in 1:nrow(Doublet.res)){
  if((Doublet.res[i, "Doublet_DoubletFinder"] == "Doublet") & (Doublet.res[i, "Doublet_scDblFinder"] == "doublet")){
    Doublet.res[i, "Doublet_BothMethods"] <- "doublet"
  } else {
    Doublet.res[i, "Doublet_BothMethods"] <- "singlet"
  }
}
rownames(Doublet.res) <- Doublet.res$cellID
Doublet.res <- Doublet.res[,-1]

# save results
write.csv(Doublet.res, paste0(Sys.Date(), "_DoubletIdentificationResults_BothMethods.csv"))
