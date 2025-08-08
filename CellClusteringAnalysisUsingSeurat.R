##################################
## Rscript for clustering cells ##
##################################

# This code aims at 
# A. Running cell clustering analysis without any prior sample integration, and visualising the results on a UMAP.
# B. Generating dot plots of expression levels of markers of each cell group,
# C. Assigning cell group identities, regenerating dot plot of marker gene expression and assessing number of cells per cell group, and
# D. Assessing cell proportion differences/similarities between samples
# The code has been used for generating Figure 1C, 4B-C, S1C-D and S3C-E.

# Define variables
Variables <- list(
  # Path to seurat object after cell and gene filtering (see CellAndGeneFiltering.R code for cell and gene filtering, UMI counts as input)
  InputPath = "SeuratObject_CellLine.rds",

  # Sample ID (# will be added to output file names)
  Sample = "CellLine",
  
  # number of PCs to keep for clustering analysis 
  ## choose final after seeing the heatmap of PC contribution and the elbow plot
  pca.dim = 20,
  
  # specify if HTO data is available
  HTO.available = FALSE,
  
  # colors to assign to each sample
  Sample.colors = NA, # NA if sample identity wasn't assigned based on HTO label, or provide a string of colors like c("lightblue", "lightpink")
  
  # Gene symbols of genotyping targets
  ## to assess expression level
  geno.targets = c("S100A10", "HIST1H1C"),
  
  # Path to directory containing doublet inference analysis (see GEXbasedDoubletPrediction.R code for doublet inference analysis)
  doublet.path = NA,
  
  # list of genes to visualize expression levels in a dot plot format
  ## example below is for known gene markers for the two cell lines profiled
  genes.interest = c(
    "CTAG2", "S100A16", "KRT19", "FXYD3", "SPINT2", "VAMP8", "SERF2", "MYC", "CD151", "ERBB2", "MED1", "STARD3", "GRB7", "ADIRF", "DHRS2", # SKBR3
    "TFF1", "SLC9A3R1", "STARD10", "CD47", "ALCAM", "ESR1", "CYBA", "GDF15", "HIST1H2BD", "HSPB1", "PPDPF", "GAL", "CYP1B1", "EPCAM", "GATA3"  # MCF7
  ),
  
  # list of genes to visualize expression levels on UMAP
  ## example below is for known gene markers for the two cell lines profiled
  genes.UMAP = c("ESR1", "ERBB2"),
  
  # list of clusters belonging to each group (e.g. clusters with SK-BR-3 markers, clusters with MCF-7 markers, ...)
  ## determined after looking at the marker genes of each cluster, known markers and QC metrics
  ## can be NULL or a list of named groups and the corresponding cluster numbers
  group.labels = list("SKBR3" = c(0, 4, 5, 6,7,9),
                   "MCF7" = c(1,2,3,8,10,12),
                   "Mixed" = c(11)),
  
  # group to exclude (e.g. cell cluster with aberrant expression profile + enriched in GEX-based predicted doublets)
  ## determined after looking at the marker genes of each cluster, known markers and QC metrics
  group.exclude = "Mixed"
)

# Load required packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(stringr)
library(RColorBrewer)
library(dplyr)

# Load the Seurat object containing the GEX data (after filtering genes and cells)
res.data <- readRDS(Variables$InputPath)

# Normalize data
## Default of NormalizeData function= LogNormalize --> Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor (default = 10000). This is then natural-log transformed using log1p.
res.data <- NormalizeData(res.data)

# Identify top 2000 variable features
## All parameters set to default
## vst method --> First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess). Then standardizes the feature values using the observed mean and expected variance (given by the fitted line). Feature variance is then calculated on the standardized values after clipping to a maximum (see clip.max parameter).
res.data <- FindVariableFeatures(res.data, selection.method = "vst", nfeatures = 2000)

# Identify the 25 most highly variable genes

# Plot the 25 most highly variable genes with and without labels
top25 <- head(VariableFeatures(res.data), 25)
pdf(paste0(Sys.Date(), "_VariableGenes_", Variables$Sample, ".pdf"), height = 5, width=6)
VariableFeaturePlot(res.data)
LabelPoints(plot = VariableFeaturePlot(res.data), points = top25, repel = TRUE)
dev.off()

# Run the standard workflow for clustering and visualization
## Scale data
res.data <- ScaleData(res.data, verbose = FALSE, features = rownames(res.data))

## Perform linear dimensional reduction (+Visualize gens defining the first 10 PCs)
res.data <- RunPCA(res.data, npcs = 30, verbose = FALSE, features = VariableFeatures(object = res.data))

pdf(paste0(Sys.Date(), "_GenesDefiningFirst10PCs_", Variables$Sample, ".pdf"), height = 5, width=6)
VizDimLoadings(res.data, dims = 1:2, reduction = "pca")
VizDimLoadings(res.data, dims = 3:4, reduction = "pca")
VizDimLoadings(res.data, dims = 5:6, reduction = "pca")
VizDimLoadings(res.data, dims = 7:8, reduction = "pca")
VizDimLoadings(res.data, dims = 9:10, reduction = "pca")
dev.off()

## PC contribution and elbow plot to decide number of meaningful PCs to keep for clustering analysis
pdf(paste0(Sys.Date(), "_HeatmapPCcontribution_", Variables$Sample, ".pdf"), height = 12, width=6)
DimHeatmap(res.data, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()

pdf(paste0(Sys.Date(), "_ElbowPlot_", Variables$Sample, ".pdf"), height = 12, width=6)
ElbowPlot(object = res.data, ndims = 30)
dev.off()

## Run non-linear dimensional reduction 
### Uses Annoy for neighbor search, n_neighbors = 30. Builds Annoy index with metric = cosine, n_trees = 50
set.seed(123)
res.data <- RunUMAP(res.data, reduction = "pca", dims = 1:Variables$pca.dim)

## Computes nearest neighbor graph and SNN
res.data <- FindNeighbors(res.data, reduction = "pca", dims = 1:Variables$pca.dim) 

## Run cell clustering analysis
### Runs original Louvain algorithm by default
res.data <- FindClusters(res.data, resolution = 0.5) 

## Define color palette
cluster.colors <- c(brewer.pal(n = 12, name = "Set3"), brewer.pal(n = 8, name = "Set1"))

## Visualizing results
pdf(file = paste0(Sys.Date(), "_UMAPplots_HighlightingEachClusterSeparately_", Variables$Sample, ".pdf"), width = (5*length(unique(res.data@meta.data$seurat_clusters))), height = 5)
  DimPlot(res.data, reduction = "umap", split.by = "seurat_clusters", cols=cluster.colors)  + NoLegend()
dev.off()

if(Variables$HTO.available ==TRUE){
  pdf(file = paste0(Sys.Date(), "_UMAPplots_HighlightingEachSampleSeparately_", Variables$Sample, ".pdf"), width = (5*length(unique(res.data@meta.data$HTO_SampleLabel))), height = 5)
  print(DimPlot(res.data, reduction = "umap", split.by = "HTO_SampleLabel", cols=cluster.colors)  + NoLegend())
  dev.off()
  
  pdf(file = paste0(Sys.Date(), "_UMAPplots_HighlightingSampleOriginAndClusterNumbers_", Variables$Sample, ".pdf"), width = 5, height = 5)
  DimPlot(res.data, reduction = "umap", label = TRUE, repel = TRUE, cols=cluster.colors)
  DimPlot(res.data, reduction = "umap", group.by = "HTO_SampleLabel", cols = Variables$Sample.colors)
  dev.off()

}else{
  pdf(file = paste0(Sys.Date(), "_UMAPplots_HighlightingClusterNumbers_", Variables$Sample, ".pdf"), width = 5, height = 5)
  DimPlot(res.data, reduction = "umap", label = TRUE, repel = TRUE, cols=cluster.colors)
  dev.off()
}

## Highlight QC metrics per cluster
pdf(file = paste0(Sys.Date(), "_ViolinPlots_HighlightQCmetrics_", Variables$Sample, ".pdf"), width = 5, height = 5)
VlnPlot(res.data, features = "nFeature_RNA", pt.size = 0.2, cols=cluster.colors) + NoLegend()
VlnPlot(res.data, features = "nCount_RNA", pt.size = 0.2, cols=cluster.colors) + NoLegend()
VlnPlot(res.data, features = "percent.mt", pt.size = 0.2, cols=cluster.colors) + NoLegend()
dev.off()

pdf(file = paste0(Sys.Date(), "_UMAPplots_HighlightQCmetrics_", Variables$Sample, ".pdf"), width = 5, height = 5)
FeaturePlot(res.data, features = "nCount_RNA")
FeaturePlot(res.data, features = "nFeature_RNA")
FeaturePlot(res.data, features = "percent.mt")
dev.off()

# Assess cell cycle state
## Plot expression of proliferation markers
ProlifMarkers <- c("MKI67", "PCNA", "CDK1", "CDC20", "TOP2A", "MCM6")

DefaultAssay(res.data) <- "RNA"
pdf(file = paste0(Sys.Date(), "_Plots_HighlightingOtherProliferationMarkers_", Variables$Sample, ".pdf"), width = 15, height = 10)
FeaturePlot(res.data, features = ProlifMarkers, max.cutoff = "q98", ncol = 3)
VlnPlot(res.data, features = ProlifMarkers, n=3, cols=cluster.colors)
dev.off()

## Assign Cell-Cycle Scores
### A list of cell cycle markers, from Tirosh et al 2015, is available in Seurat
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
res.data <- CellCycleScoring(res.data, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

pdf(file = paste0(Sys.Date(), "_HighlightingCellCycleStates_", Variables$Sample, ".pdf"), width = 5, height = 5)
DimPlot(res.data, group.by = "Phase")
FeaturePlot(res.data, features = "S.Score")
VlnPlot(res.data, features = "S.Score", cols=cluster.colors) + NoLegend()
FeaturePlot(res.data, features = "G2M.Score")
VlnPlot(res.data, features = "G2M.Score", cols=cluster.colors) + NoLegend()
dev.off()

Freq.CellCycleStage <- as.data.frame(table(res.data@meta.data$seurat_clusters, res.data@meta.data$Phase))
pdf(file = paste0(Sys.Date(), "_BarPlot_HighlightingCellProportionInEachCellCycleStage_", Variables$Sample, ".pdf"), width = 4, height = 3)
ggplot(Freq.CellCycleStage, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_bar(position="fill", stat="identity") +
  scale_y_continuous(labels = scales::percent_format())
dev.off()

# Exclude ribosomal genes from variable genes before running DE analysis
features.use <- VariableFeatures(res.data)
rpl_genes <- features.use[startsWith(features.use, "RPL")]
rps_genes <- features.use[startsWith(features.use, "RPS")]
rps.rpl_genes <- c(rpl_genes, rps_genes)
features.use <- features.use[!features.use %in% rps.rpl_genes]

# Find genes differentially expressed between every cluster and all remaining cells
res.data <- SetIdent(res.data, value = "seurat_clusters")
DefaultAssay(res.data) <- "RNA"
cluster.markers <- FindAllMarkers(res.data, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0, features = features.use)

# Adjust p-values using Benjamini-Hochberg method
cluster.markers$p_val_adj <- p.adjust(cluster.markers$p_val, method = "BH")
write.csv(cluster.markers, paste0(Sys.Date(), "_GenesDistinguishingCellsFromOneClusterFromAllOtherCells_", Variables$Sample, ".csv"))

# Identify top 10 marker genes
cluster.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

cluster.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

# Generates an expression heatmap for the top 10 markers for each cluster
pdf(paste0(Sys.Date(), "_Heatmap_Top10markersForEachCluster_", Variables$Sample, ".pdf"), width =30, height = 20)
DoHeatmap(res.data, features = top10$gene) + NoLegend()
dev.off()

# Assess expression levels of genotyping targets
pdf(file = paste0(Sys.Date(), "_ExpressionLevelsOfTargetedGenesForGenotyping_", Variables$Sample, ".pdf"), width = 5*length(Variables$geno.targets), height = 5)
## Normalized counts
FeaturePlot(res.data, features = Variables$geno.targets, ncol = length(Variables$geno.targets), max.cutoff = "q98")
VlnPlot(res.data, features = Variables$geno.targets, ncol = length(Variables$geno.targets), group.by ="seurat_clusters", cols= cluster.colors)

## raw UMI counts
FeaturePlot(res.data, features = Variables$geno.targets, ncol = length(Variables$geno.targets), max.cutoff = "q98", slot = "counts")
VlnPlot(res.data, features = Variables$geno.targets, ncol = length(Variables$geno.targets), group.by ="seurat_clusters", cols= cluster.colors, slot = "counts")
dev.off()

if(Variables$HTO.available ==TRUE){
pdf(paste0(Sys.Date(), "_ExpressionLevelsOfTargetedGenesForGenotypingPerSample_", Variables$Sample, ".pdf"), width = 5*length(Variables$geno.targets), height = 5)
print(VlnPlot(res.data, group.by = "HTO_SampleLabel", features = Variables$geno.targets, ncol = length(Variables$geno.targets)) +NoLegend())
dev.off()
}

# Bar plot to show proportion of cells in each cluster from each sample
if(Variables$HTO.available == TRUE){
  freq.table <- as.data.frame(table(res.data$HTO_SampleLabel, res.data$seurat_clusters))
  names(freq.table) <- c("HTO_SampleLabel", "seurat_clusters", "Freq")
  
  pdf(paste0(Sys.Date(), "_BarPlot_ProportionOfCellsFromEachClusterPerSample_", Variables$Sample, ".pdf"), width = 5, height = 5)
  print(ggplot(freq.table, aes(fill=HTO_SampleLabel, y=Freq, x=seurat_clusters)) + 
    geom_bar(position="stack", stat="identity") + scale_fill_manual(values= Variables$Sample.colors)  +
    ggtitle("Number of cells in each cluster from each sample") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  print(ggplot(freq.table, aes(fill=HTO_SampleLabel, y=Freq, x=seurat_clusters)) + 
    geom_bar(position="fill", stat="identity") + scale_fill_manual(values= Variables$Sample.colors) + ylab("Proportion")   +
    ggtitle("Proportion of cells in each cluster from each sample")  +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  
  pdf(paste0(Sys.Date(), "_BarPlot_ProportionOfCellsFromEachSamplePerCluster_", Variables$Sample, ".pdf"), width = 5, height = 5)
  print(ggplot(freq.table, aes(fill=seurat_clusters, y=Freq, x=HTO_SampleLabel)) + 
          geom_bar(position="stack", stat="identity") + scale_fill_manual(values= cluster.colors)  +
          ggtitle("Number of cells from each cluster in each sample") +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  print(ggplot(freq.table, aes(fill=seurat_clusters, y=Freq, x=HTO_SampleLabel)) + 
          geom_bar(position="fill", stat="identity") + scale_fill_manual(values= cluster.colors) + ylab("Proportion")   +
          ggtitle("Proportion of cells from each cluster in each sample")  +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
}

write.csv(table(res.data$seurat_clusters), paste0(Sys.Date(), "_Number of cells per cluster.csv"))


# Visualize GEX-based predicted doublets
if(is.na(Variables$doublet.path)==FALSE){
  ## Import GEX-based predicted doublet identity
  temp <- list.files(path = Variables$doublet.path, pattern = "DoubletIdentificationResults_BothMethods.csv")
  doublet.res <- read.csv(paste0(Variables$doublet.path, temp), header = TRUE, row.names = 1)
  
  ## Add doublets to cell metadata
  res.data <- AddMetaData(res.data, metadata = doublet.res)

  ## Visualize GEX-based predicted doublets
  Freq.Table <- as.data.frame(table(res.data@meta.data$seurat_clusters, res.data@meta.data$Doublet_BothMethods))
  names(Freq.Table) <- c("seurat_clusters", "DoubletStatus", "Freq")
  
  pdf(file = paste0(Sys.Date(), "_HighlightingDoublets.pdf"), width = 5, height = 5)
  DimPlot(res.data, group.by = "Doublet_BothMethods", cols = c("blue", "grey"))

  ggplot(Freq.Table, aes(fill=DoubletStatus, y=Freq, x=seurat_clusters)) + 
    geom_bar(position="fill", stat="identity") +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(values = c("blue", "grey"))
  
  dev.off()
  
}

# Visualize expression of known cell type or cell line markers
if(is.na(Variables$genes.interest)==FALSE){
  pdf(file = paste0(Sys.Date(), "_DotPlot_HighlightingKnownMarkers.pdf"), width = 9, height = 2.5)
  print(DotPlot(res.data,
                group.by = "seurat_clusters", 
                features = Variables$genes.interest,
                cols = c("yellow3", "purple4")) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          xlab("") +ylab("")
  )
  dev.off()
}

if(is.na(Variables$genes.UMAP)==FALSE){
pdf(file = paste0(Sys.Date(), "_UMAPplot_HighlightingMarkers.pdf"), width = 5, height = 5)
  for(i in 1:length(Variables$genes.UMAP)){
    print(FeaturePlot(res.data, features = Variables$genes.UMAP[i], pt.size = 2, cols = c("lightgrey", "darkgreen"), max.cutoff = "q98"))
  }
dev.off()
  
}

if(length(Variables$group.labels) != 0){
  
  # Assign cell group labels
  res.data@meta.data$CellGroupLabel <- NA
  for(i in 1:length(Variables$group.labels)){
    res.data@meta.data$CellGroupLabel[res.data@meta.data$seurat_clusters %in% Variables$group.labels[[i]]] <- names(Variables$group.labels)[i]
  }
  table(res.data@meta.data$CellGroupLabel)
  
  if(is.na(Variables$group.exclude) == FALSE){
    metadata <- res.data@meta.data[!(res.data@meta.data$CellGroupLabel %in% Variables$group.exclude), ]
    res.data <- subset(res.data, cells = rownames(metadata))
    print(table(res.data@meta.data$CellGroupLabel))
    rm(metadata); gc()
  }
  
  # Highlight labels on UMAP
  pdf(file = paste0(Sys.Date(), "_UMAPplot_HighlightingCellGroupLabels.pdf"), width = 5, height = 5)
  DimPlot(res.data, group.by = "CellGroupLabel", label =TRUE, label.box = TRUE)  + NoLegend()
  dev.off()
  
  # Bar plot showing number of cells per cell group
  cell.freq <- table(res.data@meta.data$CellGroupLabel) %>% as.data.frame()
  pdf(file = paste0(Sys.Date(), "_BarPlot_ShowingCellNumberPerLabelledCellGroup.pdf"), width = 3, height = 6)
  print(ggplot(data=cell.freq, aes(x=Var1, y=Freq)) +
          geom_bar(stat="identity", fill="grey94", color="grey") +
          geom_text(aes(label=Freq), hjust=2, colour = "black", position = "dodge") +
          scale_y_log10() +
          xlab("") + ylab("Cell Number") +
          coord_flip() +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  )
  dev.off()
  
  # Dot plot of marker genes
  res.data <- SetIdent(res.data, value = "CellGroupLabel")
  if(is.na(Variables$genes.interest)==FALSE){
    pdf(file = paste0(Sys.Date(), "_DotPlot_HighlightingMarkersPerLabelledCellGroup.pdf"), width = 9, height = 2.5)
    DefaultAssay(res.data) <- "RNA"
    print(DotPlot(res.data,
                  group.by = "CellGroupLabel", 
                  features = Variables$genes.interest,
                  cols = c("yellow3", "purple4")) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
            xlab("") +ylab("")
    )
    dev.off()
  }
  
  # Proportion of cells per cell group per sample/experiment
  if(Variables$HTO.available == TRUE){
    freq.table <- as.data.frame(table(res.data$HTO_SampleLabel, res.data$CellGroupLabel))
    names(freq.table) <- c("HTO_SampleLabel", "CellGroupLabel", "Freq")
    
    pdf(paste0(Sys.Date(), "_BarPlot_ProportionOfCellsFromEachLabelledCellGroupPerSample_", Variables$Sample, ".pdf"), width = 5, height = 5)
    print(ggplot(freq.table, aes(fill=HTO_SampleLabel, y=Freq, x=CellGroupLabel)) + 
            geom_bar(position="stack", stat="identity") + scale_fill_manual(values= Variables$Sample.colors)  +
            ggtitle("Number of cells in each group from each sample") +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
    print(ggplot(freq.table, aes(fill=HTO_SampleLabel, y=Freq, x=CellGroupLabel)) + 
            geom_bar(position="fill", stat="identity") + scale_fill_manual(values= Variables$Sample.colors) + ylab("Proportion")   +
            ggtitle("Proportion of cells in each group from each sample")  +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
    dev.off()
    
    pdf(paste0(Sys.Date(), "_BarPlot_ProportionOfCellsFromEachSamplePerLabelledCellGroup_", Variables$Sample, ".pdf"), width = 5, height = 5)
    print(ggplot(freq.table, aes(fill=CellGroupLabel, y=Freq, x=HTO_SampleLabel)) + 
            geom_bar(position="stack", stat="identity") + scale_fill_manual(values= cluster.colors)  +
            ggtitle("Number of cells from each group in each sample") +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
    print(ggplot(freq.table, aes(fill=CellGroupLabel, y=Freq, x=HTO_SampleLabel)) + 
            geom_bar(position="fill", stat="identity") + scale_fill_manual(values= cluster.colors) + ylab("Proportion")   +
            ggtitle("Proportion of cells from each group in each sample")  +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
    dev.off()
  }
  
}

# save Seurat object with clustering results
saveRDS(res.data, file = paste0(Sys.Date(), "_SeuratObject_CellClusteringResults_", Variables$Sample, ".rds"))

rm(list = ls()); gc(); gc()
