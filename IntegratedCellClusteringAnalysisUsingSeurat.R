################################################################################
## Rscript for clustering cells after integrating data from different samples ##
################################################################################

# This code aims at
# A. Running cell clustering analysis after integrating data from different samples, and visualising the results on a UMAP.
# B. Generating dot plots of expression levels of markers of each cell group,
# C. Assigning cell group identities, regenerating dot plot of marker gene expression and assessing number of cells per cell group,
# D. Assessing cell proportion differences/similarities between samples/phenotypes, and
# E. Evaluating signature scores per cell group (e.g. cell cycle, …)
# The code has been used for generating Figure 2A-D, S4C-H.

# Define variables
Variables <- list(
  # Path(s) to seurat objects after cell and gene filtering containing data from multiple samples (see CellAndGeneFiltering.R code for cell and gene filtering, UMI counts as input)
  InputPath1 = "SeuratObject_Primary.rds",
  InputPath2 = NA,
  
  # name of analysis (will also appear in output file name)
  AnalysisName = "IntegratedPrimarySamples",
  
  # Identity of any sample to exclude (e.g. if low quality, ...)
  sample.exclude = NA,

  # number of PCs to keep for clustering analysis 
  ## choose final after seeing the heatmap of PC contribution and the elbow plot
  dims = 1:20,
  
  # resolution to use for clustering
  resolution = 0.5,
  
  # seed to use for clustering
  seed.use = 42,
  
  # name of column containing sample identity (after HTO-based sample/phenotype demultiplexing)
  TumorCol.plot = "HTO_SampleLabel",
  
  # threshold for minimum number of cells in which a gene should be detected for it to be retained in the dataset
  GeneFilter = 3,
  
  # specify if HTO data is available
  HTO.available = TRUE,
  
  # Gene symbols of genotyping targets
  ## to assess expression level
  geno.targets = c("B2M", "TNFAIP3"),
  
  # Path to directory containing doublet inference analysis (see GEXbasedDoubletPrediction.R code for doublet inference analysis)
  doublet.path = "DoubletIdentificationResults_BothMethods.csv",
  
  # genes to visualize expression levels in a dot plot format
  genes.interest = c(
    "CD3D", "CD3E", "CD2", "CD7", "CD4", "CD8A", "CD8B",
    "CD14", "ITGAX", "APOE", "LYZ", "CD68",
    "IL3RA", "CLEC4C", "NRP1", "CLEC10A", "CLEC9A",
    "PRDM1", "JCHAIN", "XBP1", "CD38", "TNFRSF17", "CD55"
  ),
  
  # genes to visualize expression levels on UMAP
  genes.UMAP = c("MS4A1", "CD3D"),
  
  # list of clusters belonging to each group (e.g. clusters with SK-BR-3 markers, clusters with MCF-7 markers, ...)
  ## determined after looking at the marker genes of each cluster, known markers and QC metrics
  ## can be NULL or a list of named groups and the corresponding cluster numbers
  group.labels = list("Tcells" = c(5, 11, 13),
                      "MyeloidCells" = c(0),
                      "MalignantB" = c(1,2,3,4,6,7,8,10, 12),
                      "Mixed" = c(9)),

  # group to exclude (e.g. cell cluster with aberrant expression profile + enriched in GEX-based predicted doublets)
  ## determined after looking at the marker genes of each cluster, known markers and QC metrics
  group.exclude = "Mixed",
  
  # path to list of genes to evaluate from literature (+ groups to retain + sheet number containing genes)
  ## in this case, from Nadeu et al 2022
  signature.path = "NadeuNatMed2022_DetectionEarlyRichterTransformationInCLL/41591_2022_1927_MOESM3_ESM.xlsx",
  signature.sheets = 80:84,
  signature.group = c("CXCR4hiCD27lo", "CXCR4loCD27hi", "MIR155HGhi", "CCND2lo RT", "CCND2hi RT", "RT proliferative", "MZB1hiIGHMhiXBP1hi")
)

# Load required packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(stringr)
library(dplyr)
library(RColorBrewer)

# Load the Seurat object containing the GEX data for all samples (merged data after filtering genes and cells)
GEX.data1 <- readRDS(Variables$InputPath1)
if(is.na(Variables$InputPath2) ==FALSE){
  GEX.data2 <- readRDS(Variables$InputPath2)
  
  # merge data
  GEX.data <- merge(x=GEX.data1, y = c(GEX.data2), add.cell.ids = NULL, project = "AllSamples")
  rm(GEX.data1); rm(GEX.data2); gc()

}else{
  GEX.data <- GEX.data1
  rm(GEX.data1); gc()
}

GEX.data@meta.data$orig.ident <- gsub("Primary_", "", GEX.data@meta.data$orig.ident)
GEX.data@meta.data$HTO_SampleLabel <- paste0(GEX.data@meta.data$orig.ident, "_", GEX.data@meta.data$HTO_SampleLabel)
table(GEX.data@meta.data$HTO_SampleLabel)

# exclude samples
if(is.na(Variables$sample.exclude)==FALSE){
  metadata <- GEX.data@meta.data[!(GEX.data@meta.data$HTO_SampleLabel %in% Variables$sample.exclude),]
  GEX.data <- subset(GEX.data, cells = rownames(metadata))
  rm(metadata); gc()
  table(GEX.data@meta.data$HTO_SampleLabel)
}

# Filter out genes detected in < X cells
counts <- GetAssayData(GEX.data, slot="counts", assay="RNA")   
genes.detection <- rowSums(counts!= 0 )
genes.filter <- names(genes.detection[genes.detection >= Variables$GeneFilter])  #select genes detected in at least X cells
GEX.data <- subset(GEX.data, features = genes.filter)

# Split dataset into a list of seurat objects (one per tumor)
GEX.data.list <- SplitObject(GEX.data, split.by = Variables$TumorCol.plot)
rm(GEX.data); gc(); gc()

# Normalize data and identify variable features for each sample
GEX.data.list <- lapply(X = GEX.data.list, FUN = function(x) {
  
  ## NormalizeData function: default = LogNormalize --> Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor (default = 10000). This is then natural-log transformed using log1p.
  x <- NormalizeData(x)
  
  ## FindVariableFeatures function: all parameters set to default. selection.method (method used to choose top variable features) = vst --> First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess). Then standardizes the feature values using the observed mean and expected variance (given by the fitted line). Feature variance is then calculated on the standardized values after clipping to a maximum (see clip.max parameter).
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Select features that are repeatedly variable across patient tumors for integration
features <- SelectIntegrationFeatures(object.list = GEX.data.list)

# Identify a set of anchors that will be used to integrate the patient tumors
## 1. Perform dimensional reduction on the dataset pair using the CCA method
## 2. Identify anchors - pairs of cells from each dataset that are contained within each other's neighborhoods (also known as mutual nearest neighbors).
## 3. Filter low-confidence anchors to ensure anchors in the low dimension space are in broad agreement with the high dimensional measurements. This is done by looking at the neighbors of each query cell in the reference dataset using max.features to define this space. If the reference cell isn't found within the first k.filter neighbors, remove the anchor.
## 4. Assign each remaining anchor a score. For each anchor cell, determine the nearest k.score anchors within its own dataset and within its pair's dataset. Based on these neighborhoods, construct an overall neighbor graph and then compute the shared neighbor overlap between anchor and query cells (analogous to an SNN graph). We use the 0.01 and 0.90 quantiles on these scores to dampen outlier effects and rescale to range between 0-1.
integration.anchors <- FindIntegrationAnchors(object.list = GEX.data.list, anchor.features = features)

# Integrate the data form distinct patient tumors
integrated.data <- IntegrateData(anchorset = integration.anchors)

# Run the standard workflow for visualization and clustering on the integrated dataset
DefaultAssay(integrated.data) <- "integrated"
integrated.data <- ScaleData(integrated.data, verbose = FALSE)
integrated.data <- RunPCA(integrated.data, npcs = 30, verbose = FALSE)

pdf(paste0(Sys.Date(), "_HeatmapPCcontribution.pdf"), height = 12, width=6)
DimHeatmap(integrated.data, dims = 1:30, cells = 1000, balanced = TRUE)
dev.off()

pdf(paste0(Sys.Date(), "_ElbowPlot.pdf"), height = 6, width=6)
ElbowPlot(integrated.data)
dev.off()

set.seed(123)
integrated.data <- RunUMAP(integrated.data, reduction = "pca", dims = Variables$dims) ## Uses Annoy for neighbor search, n_neighbors = 30. Builds Annoy index with metric = cosine, n_trees = 50.
integrated.data <- FindNeighbors(integrated.data, reduction = "pca", dims = Variables$dims) ## Computes nearest neighbor graph and SNN
integrated.data <- FindClusters(integrated.data, resolution = Variables$resolution) ## Runs original Louvain algorithm (default)

# Define color palette
cluster.colors <- c(brewer.pal(n = 12, name = "Set3"), brewer.pal(n = 8, name = "Set1"))

# Visualizing results
pdf(file = paste0(Sys.Date(), "_UMAPplots_HighlightingSampleOriginAndClusterNumbers.pdf"), width = 5, height = 5)
DimPlot(integrated.data, reduction = "umap", label = TRUE, repel = TRUE, label.box = TRUE, cols= cluster.colors) + NoLegend()
DimPlot(integrated.data, reduction = "umap", label = FALSE, group.by = "orig.ident", cols= cluster.colors)
DimPlot(integrated.data, reduction = "umap", label = FALSE, group.by = Variables$TumorCol.plot, cols= cluster.colors)
if(Variables$HTO.available ==TRUE){
  DimPlot(integrated.data, reduction = "umap", label = FALSE, group.by = "HTO_Phenotype", cols= cluster.colors)
}
dev.off()

# Visualize patient tumors side-by-side
pdf(file = paste0(Sys.Date(), "_UMAPplots_HighlightingExperiments.pdf"), width = (5*length(unique(integrated.data@meta.data$orig.ident))), height = 5)
DimPlot(integrated.data, reduction = "umap", split.by = "orig.ident", cols= cluster.colors)
dev.off()

pdf(file = paste0(Sys.Date(), "_UMAPplots_HighlightingSampleOriginPerExperiment.pdf"), width = (5*length(unique(integrated.data@meta.data[,Variables$TumorCol.plot]))), height = 5)
DimPlot(integrated.data, reduction = "umap", split.by = Variables$TumorCol.plot, cols= cluster.colors)
dev.off()

pdf(file = paste0(Sys.Date(), "_UMAPplots_HighlightingClusterNumbers.pdf"), width = (5*length(unique(integrated.data@meta.data$seurat_clusters))), height = 5)
DimPlot(integrated.data, reduction = "umap", split.by = "seurat_clusters", cols= cluster.colors)  + NoLegend()
dev.off()

if(Variables$HTO.available ==TRUE){
  pdf(file = paste0(Sys.Date(), "_UMAPplots_HighlightingSampleOrigin.pdf"), width = (5*length(unique(integrated.data@meta.data$HTO_SampleLabel))), height = 5)
  DimPlot(integrated.data, reduction = "umap", split.by = "HTO_SampleLabel", group.by = "HTO_Phenotype", cols= cluster.colors)
  dev.off()
}

write.csv(integrated.data@meta.data, paste0(Sys.Date(), "_CellMetadataWithClusteringResults.csv"))

## Highlight QC metrics and CD3 expression
pdf(file = paste0(Sys.Date(), "_ViolinPlots_HighlightQCmetricsPerCluster.pdf"), width = 5, height = 5) ##
VlnPlot(integrated.data, features = "nFeature_RNA", pt.size = 0.2, cols= cluster.colors) + NoLegend()
VlnPlot(integrated.data, features = "nCount_RNA", pt.size = 0.2, cols= cluster.colors) + NoLegend()
VlnPlot(integrated.data, features = "percent.mt", pt.size = 0.2, cols= cluster.colors) + NoLegend()
dev.off()

pdf(file = paste0(Sys.Date(), "_UMAPplots_HighlightQCmetrics.pdf"), width = 5, height = 5)
FeaturePlot(integrated.data, features = "nCount_RNA")
FeaturePlot(integrated.data, features = "nFeature_RNA")
FeaturePlot(integrated.data, features = "percent.mt")
dev.off()

# Assess cell cycle state
## Plot expression of proliferation markers
ProlifMarkers <- c("MKI67", "PCNA", "CDK1", "CDC20", "TOP2A", "MCM6")

DefaultAssay(integrated.data) <- "RNA"
pdf(file = paste0(Sys.Date(), "_Plots_HighlightingOtherProliferationMarkers.pdf"), width = 15, height = 10)
FeaturePlot(integrated.data, features = ProlifMarkers, max.cutoff = "q98", ncol = 3)
VlnPlot(integrated.data, features = ProlifMarkers, n=3, cols= cluster.colors)
dev.off()

## Assign Cell-Cycle Scores
### A list of cell cycle markers, from Tirosh et al 2015, is available in Seurat
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
integrated.data <- CellCycleScoring(integrated.data, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

pdf(file = paste0(Sys.Date(), "_HighlightingCellCycleStates.pdf"), width = 5, height = 5)
DimPlot(integrated.data, group.by = "Phase")
FeaturePlot(integrated.data, features = "S.Score")
VlnPlot(integrated.data, features = "S.Score", cols= cluster.colors) + NoLegend()
FeaturePlot(integrated.data, features = "G2M.Score")
VlnPlot(integrated.data, features = "G2M.Score", cols= cluster.colors) + NoLegend()
dev.off()

Freq.CellCycleStage <- as.data.frame(table(integrated.data@meta.data$seurat_clusters, integrated.data@meta.data$Phase))
pdf(file = paste0(Sys.Date(), "_BarPlot_HighlightingCellProportionInEachCellCycleStage.pdf"), width = 4, height = 3)
ggplot(Freq.CellCycleStage, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_bar(position="fill", stat="identity") +
  scale_y_continuous(labels = scales::percent_format())
dev.off()

# Exclude ribosomal genes from variable genes before running DE analysis
DefaultAssay(integrated.data) <- "integrated"
features.use <- VariableFeatures(integrated.data)
rpl_genes <- features.use[startsWith(features.use, "RPL")]
rps_genes <- features.use[startsWith(features.use, "RPS")]
rps.rpl_genes <- c(rpl_genes, rps_genes)
features.use <- features.use[!features.use %in% rps.rpl_genes]

# Find genes differentially expressed between every cluster and all remaining cells
integrated.data <- SetIdent(integrated.data, value = "seurat_clusters")
DefaultAssay(integrated.data) <- "RNA"
cluster.markers <- FindAllMarkers(integrated.data, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0, features = features.use)
write.csv(cluster.markers, paste0(Sys.Date(), "_GenesDistinguishingCellsFromOneClusterFromAllOtherCells.csv"))

# Identify top 10 marker genes
cluster.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

cluster.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

# scale data for plotting heatmap
integrated.data <- ScaleData(integrated.data, verbose = FALSE)

# Generates an expression heatmap for the top 10 markers for each cluster
cluster.ids <- unique(cluster.markers$cluster)
pdf(paste0(Sys.Date(), "_Heatmap_Top10markersForEachCluster.pdf"), width =30, height = 20)
DoHeatmap(integrated.data, features = top10$gene, group.colors = cluster.colors) + NoLegend()
dev.off()

# Assess expression levels of genotyping targets
pdf(file = paste0(Sys.Date(), "_ExpressionLevelsOfTargetedGenesForGenotyping_", Variables$Sample, ".pdf"), width = 5*length(Variables$geno.targets), height = 5)
## Normalized counts
FeaturePlot(integrated.data, features = Variables$geno.targets, ncol = length(Variables$geno.targets), max.cutoff = "q98")
VlnPlot(integrated.data, features = Variables$geno.targets, ncol = length(Variables$geno.targets), group.by ="seurat_clusters", cols= cluster.colors)

## raw UMI counts
FeaturePlot(integrated.data, features = Variables$geno.targets, ncol = length(Variables$geno.targets), max.cutoff = "q98", slot = "counts")
VlnPlot(integrated.data, features = Variables$geno.targets, ncol = length(Variables$geno.targets), group.by ="seurat_clusters", cols= cluster.colors, slot = "counts")
dev.off()

if(Variables$HTO.available ==TRUE){
  pdf(paste0(Sys.Date(), "_ExpressionLevelsOfTargetedGenesForGenotypingPerSample.pdf"), width = 5*length(Variables$geno.targets), height = 5)
  print(VlnPlot(integrated.data, group.by = "HTO_SampleLabel", features = Variables$geno.targets, ncol = length(Variables$geno.targets)) +NoLegend())
  dev.off()
}

# Bar plot to show proportion of cells in each cluster from each sample
if(Variables$HTO.available == TRUE){
  freq.table <- as.data.frame(table(integrated.data$HTO_SampleLabel, integrated.data$seurat_clusters))
  names(freq.table) <- c("HTO_SampleLabel", "seurat_clusters", "Freq")
  
  pdf(paste0(Sys.Date(), "_BarPlot_ProportionOfCellsFromEachClusterPerSample.pdf"), width = 5, height = 5)
  print(ggplot(freq.table, aes(fill=HTO_SampleLabel, y=Freq, x=seurat_clusters)) + 
          geom_bar(position="stack", stat="identity") + scale_fill_manual(values= cluster.colors)  +
          ggtitle("Number of cells in each cluster from each sample") +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  print(ggplot(freq.table, aes(fill=HTO_SampleLabel, y=Freq, x=seurat_clusters)) + 
          geom_bar(position="fill", stat="identity") + scale_fill_manual(values=cluster.colors) + ylab("Proportion")   +
          ggtitle("Proportion of cells in each cluster from each sample")  +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  
  pdf(paste0(Sys.Date(), "_BarPlot_ProportionOfCellsFromEachSamplePerCluster.pdf"), width = 5, height = 5)
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

write.csv(table(integrated.data$seurat_clusters), paste0(Sys.Date(), "_Number of cells per cluster.csv"))

# Visualize GEX-based predicted doublets
if(is.na(Variables$doublet.path)==FALSE){
  ## Import GEX-based predicted doublet identity
  temp <- list.files(path = Variables$doublet.path, pattern = "DoubletIdentificationResults_BothMethods.csv")
  doublet.res <- read.csv(paste0(Variables$doublet.path, temp), header = TRUE, row.names = 1)
  
  ## Add doublets to cell metadata
  integrated.data <- AddMetaData(integrated.data, metadata = doublet.res)
  
  ## Visualize GEX-based predicted doublets
  Freq.Table <- as.data.frame(table(integrated.data@meta.data$seurat_clusters, integrated.data@meta.data$Doublet_BothMethods))
  names(Freq.Table) <- c("seurat_clusters", "DoubletStatus", "Freq")
  
  pdf(file = paste0(Sys.Date(), "_HighlightingDoublets.pdf"), width = 5, height = 5)
  DimPlot(integrated.data, group.by = "Doublet_BothMethods", cols = c("blue", "grey"))
  
  ggplot(Freq.Table, aes(fill=DoubletStatus, y=Freq, x=seurat_clusters)) + 
    geom_bar(position="fill", stat="identity") +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(values = c("blue", "grey"))
  
  dev.off()
  
}

# Visualize expression of known cell type markers
if(is.na(Variables$genes.interest)==FALSE){
  pdf(file = paste0(Sys.Date(), "_DotPlot_HighlightingKnownMarkers.pdf"), width = 9, height = 2.5)
  print(DotPlot(integrated.data,
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
    print(FeaturePlot(integrated.data, features = Variables$genes.UMAP[i], pt.size = 2, cols = c("lightgrey", "darkgreen"), max.cutoff = "q98"))
  }
  dev.off()
}

if(length(Variables$group.labels) != 0){
  
  # Assign cell group labels
  integrated.data@meta.data$CellGroupLabel <- NA
  for(i in 1:length(Variables$group.labels)){
    integrated.data@meta.data$CellGroupLabel[integrated.data@meta.data$seurat_clusters %in% Variables$group.labels[[i]]] <- names(Variables$group.labels)[i]
  }
  table(integrated.data@meta.data$CellGroupLabel)
  
  if(is.na(Variables$group.exclude) == FALSE){
    metadata <- integrated.data@meta.data[!(integrated.data@meta.data$CellGroupLabel %in% Variables$group.exclude), ]
    integrated.data <- subset(integrated.data, cells = rownames(metadata))
    print(table(integrated.data@meta.data$CellGroupLabel))
    rm(metadata); gc()
  }
  
  # Highlight labels on UMAP
  pdf(file = paste0(Sys.Date(), "_UMAPplot_HighlightingCellGroupLabels.pdf"), width = 5, height = 5)
  DimPlot(integrated.data, group.by = "CellGroupLabel", label =TRUE, label.box = TRUE)  + NoLegend()
  dev.off()
  
  # Bar plot showing number of cells per cell group
  cell.freq <- table(integrated.data@meta.data$CellGroupLabel) %>% as.data.frame()
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
  integrated.data <- SetIdent(integrated.data, value = "CellGroupLabel")
  if(is.na(Variables$genes.interest)==FALSE){
    pdf(file = paste0(Sys.Date(), "_DotPlot_HighlightingMarkersPerLabelledCellGroup.pdf"), width = 9, height = 2.5)
    DefaultAssay(integrated.data) <- "RNA"
    print(DotPlot(integrated.data,
                  group.by = "CellGroupLabel", 
                  features = Variables$genes.interest,
                  cols = c("yellow3", "purple4")) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
            xlab("") +ylab("")
    )
    dev.off()
  }
  
  # Dot plot of marker genes per labelled group per patient
  integrated.data@meta.data$CellGroupLabel_Patient <- paste0(integrated.data@meta.data$CellGroupLabel, "_", integrated.data@meta.data$HTO_SampleLabel)
  pdf(file = paste0(Sys.Date(), "_DotPlot_HighlightingMarkersPerLabelledCellGroupPerPatient.pdf"), width = 11, height = 11)
  print(DotPlot(integrated.data,
                group.by = "CellGroupLabel_Patient", 
                features = Variables$genes.interest,
                cols = c("yellow3", "purple4")) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          xlab("") +ylab("")
  )
  dev.off()
  
  # Bar plot showing number of cells per cell group per patient
  cell.freq <- table(integrated.data@meta.data$CellGroupLabel_Patient) %>% as.data.frame()
  pdf(file = paste0(Sys.Date(), "_BarPlot_ShowingCellNumberPerCellGroupPerPatient.pdf"), width = 3, height = 18)
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
  
  # Proportion of cells per cell group per sample/experiment
  if(Variables$HTO.available == TRUE){
    freq.table <- as.data.frame(table(integrated.data$HTO_SampleLabel, integrated.data$CellGroupLabel))
    names(freq.table) <- c("HTO_SampleLabel", "CellGroupLabel", "Freq")
    
    pdf(paste0(Sys.Date(), "_BarPlot_ProportionOfCellsFromEachLabelledCellGroupPerSample_", Variables$Sample, ".pdf"), width = 5, height = 5)
    print(ggplot(freq.table, aes(fill=HTO_SampleLabel, y=Freq, x=CellGroupLabel)) + 
            geom_bar(position="stack", stat="identity") + scale_fill_manual(values= cluster.colors)  +
            ggtitle("Number of cells in each group from each sample") +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
    print(ggplot(freq.table, aes(fill=HTO_SampleLabel, y=Freq, x=CellGroupLabel)) + 
            geom_bar(position="fill", stat="identity") + ylab("Proportion")   +scale_fill_manual(values= cluster.colors)  +
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
  
  # Proportion of cells from each cell group in CLL versus LBCL
  integrated.data_subset <- subset(integrated.data, subset = HTO_Phenotype %in% c("CLL", "LBCL"))
  cellNumbers <- as.data.frame(table(integrated.data_subset@meta.data$CellGroupLabel, integrated.data_subset@meta.data$HTO_Phenotype))
  names(cellNumbers) <- c("CellGroupLabel", "DiseasePhenotype", "freq")
  for(i in c("CLL", "LBCL")){
    cellNumbers_subset <- cellNumbers[cellNumbers$DiseasePhenotype == i,]
    sum.cells <- sum(cellNumbers_subset$freq)
    
    for(j in 1:nrow(cellNumbers)){
      if(cellNumbers[j, "DiseasePhenotype"] == i){
        cellNumbers[j, "Proportion_All"] <- (cellNumbers[j, "freq"] / sum.cells)
      }
    }
  }
  
  patient.ids <- unique(integrated.data_subset$HTO_SampleLabel)
  for(k in 1:length(patient.ids)){
    integrated.data_subset2 <- subset(integrated.data_subset, subset = HTO_SampleLabel == patient.ids[k])
    
    cellNumbers2 <- as.data.frame(table(integrated.data_subset2@meta.data$CellGroupLabel, integrated.data_subset2@meta.data$HTO_Phenotype))
    names(cellNumbers2) <- c("CellGroupLabel", "DiseasePhenotype", "freq")
    for(i in c("CLL", "LBCL")){
      cellNumbers_subset <- cellNumbers2[cellNumbers2$DiseasePhenotype == i,]
      
      sum.cells <- sum(cellNumbers_subset$freq)
      
      cell.types <- unique(cellNumbers_subset$CellGroupLabel)
      for(j in cell.types){
        cellNumbers[(cellNumbers$CellGroupLabel == j & cellNumbers$DiseasePhenotype == i), paste0("Proportion_", patient.ids[k])] <- (cellNumbers_subset[cellNumbers_subset$CellGroupLabel == j, "freq"] / sum.cells)
      }
    }
    
  }
  
  cell.types <- unique(cellNumbers$CellGroupLabel)
  cellNumbers <- cellNumbers[,!(names(cellNumbers) %in% "freq")]
  
  pdf(file = paste0(Sys.Date(), "_BarPlot_HighlightingCelltypeProportionInCLLvsLBCL.pdf"), width = 2, height =3)
  for(i in 1:length(cell.types)){
    
    # subset one cell type
    cellNumbers_subset <- cellNumbers[cellNumbers$CellGroupLabel == cell.types[i], ]
    rownames(cellNumbers_subset) <- cellNumbers_subset$DiseasePhenotype
    cellNumbers_subset <- cellNumbers_subset[, !(names(cellNumbers_subset) %in%"Proportion_All")]
    names(cellNumbers_subset) <- gsub("Proportion_", "", names(cellNumbers_subset))
    cellNumbers_subset <- cellNumbers_subset %>% t %>% as.data.frame
    cellNumbers_subset$Patient <- rownames(cellNumbers_subset)
    cellNumbers_subset <- cellNumbers_subset[!(rownames(cellNumbers_subset) %in% c("CellGroupLabel", "DiseasePhenotype")),]
    
    library(data.table)
    cellNumbers_long <- melt(setDT(cellNumbers_subset), id.vars = "Patient", variable.name = "DiseasePhenotype")
    cellNumbers_long$value <- as.numeric(cellNumbers_long$value)
    
    # plot
    print(
      ggplot(cellNumbers_long,
             aes(y = value, x = DiseasePhenotype, col=Patient)) +
        geom_point() +
        geom_line(aes(group = Patient), color = "gray", size = 0.2) +
        theme_classic() +
        scale_x_discrete(limits= c("CLL", "LBCL")) +
        xlab("") + ylab("Cell Proprotion") + ggtitle(cell.types[i]) +
        theme(legend.title=element_blank(), legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    )
  }
  dev.off()
  
  # calculate signature score based on published signatures and highlight score per cell group
  ## define cell groups of interest
  cell.groups <- Variables$signature.group
  
  ## get excel sheet number with genes of interest
  sheets <- Variables$signature.sheets

  ## create empty list to store genes
  gene.list <- list()
  
  ## extract genes from literature
  for(i in 1:length(sheets)){
    ### import file
    genes.interest_patient <- readxl::read_excel(Variables$signature.path, sheet = sheets[i]) %>% as.data.frame()
    
    ### get sample ID
    sampleID <- genes.interest_patient[2,1]
    sampleID <- gsub(".*clusters of ", "", sampleID)
    sampleID <- gsub(" ", "", sampleID)
    
    ### remove unnecessary rows
    names(genes.interest_patient) <- genes.interest_patient[4,]
    genes.interest_patient <- genes.interest_patient[5:nrow(genes.interest_patient),]
    
    ### keep significantly upregulated genes
    genes.interest_patient$p_val_adj <- as.numeric(genes.interest_patient$p_val_adj)
    genes.interest_patient <- genes.interest_patient[genes.interest_patient$p_val_adj < 0.05,]
    genes.interest_patient <- genes.interest_patient[genes.interest_patient$avg_log2FC >0,]
    
    for(j in 1:length(cell.groups)){
      cell.group.interest <- gsub(" ", "", cell.groups[j])
      
      ### extract list of marker genes for specific group
      genes <- genes.interest_patient[genes.interest_patient$annotation == cell.groups[j], "gene"]
      
      ### save list of marker genes for specific group
      gene.list[[paste0(cell.group.interest, "_", sampleID)]] <- genes
    }
    
    rm(genes.interest_patient); rm(genes); rm(cell.group.interest); rm(sampleID); gc()
  }
  
  ## calculate score per cell for each siganture
  for(i in 1:length(cell.groups)){
    cell.group.interest <- gsub(" ", "", cell.groups[i])
    
    ### extract gene list for each cell group
    gene.list_subset <- gene.list[grepl(cell.group.interest, names(gene.list))]
    names(gene.list_subset) <- gsub(".*_", "", names(gene.list_subset))
    
    ### remove empty lists
    gene.list_subset <- gene.list_subset[lapply(gene.list_subset,length)>0]
    
    ### Identify genes upregulated in cell group from multiple patients
    CompareGenes <- venn::venn(gene.list_subset, simplify = FALSE, intersections = TRUE, show.plot=TRUE, zcolor = "style")
    genes.intersect <- attr(CompareGenes, "intersections")
    names(genes.intersect) <- gsub(":", "_", names(genes.intersect))

    if(length(gene.list_subset)>1){
      
      ### Keep genes upregulated in same cell group in at least 2 samples
      min2samples.genes <- unlist(genes.intersect[grepl("_", names(genes.intersect))], recursive = FALSE)
      
      ### Calculate gene module score based on those genes
      integrated.data <- AddModuleScore(
        integrated.data,
        features = list(min2samples.genes),
        ctrl = 100,
        k = FALSE,
        assay = "RNA",
        name = paste0("Score_", cell.group.interest, "_2min"),
        seed = 123
      )
      
    }else{
      
      ### Keep genes upregulated in same cell group
      patient.genes <- genes.intersect[[1]]
      patient.ID <- names(genes.intersect)
      
      ### Calculate gene module score based on those genes
      integrated.data <- AddModuleScore(
        integrated.data,
        features = list(patient.genes),
        ctrl = 100,
        k = FALSE,
        assay = "RNA",
        name = paste0("Score_", cell.group.interest, "_", patient.ID),
        seed = 123
      )
    }
  }
  
  # Dot plot of signature score per cell groups
  module.score <- names(integrated.data@meta.data)[grepl("Score_", names(integrated.data@meta.data))]
  pdf(file = paste0(Sys.Date(), "_DotPlot_HighlightingSignatureScores.pdf"), width = 6, height = 6)
  print(DotPlot(integrated.data,
                group.by = "CellGroupLabel", 
                features = module.score,
                cols = c("yellow", "purple")) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          xlab("") +ylab("")
  )
  dev.off()
  
}

saveRDS(integrated.data, paste0(Sys.Date(), "_", Variables$AnalysisName, ".rds"))

rm(list = ls()); gc(); gc()