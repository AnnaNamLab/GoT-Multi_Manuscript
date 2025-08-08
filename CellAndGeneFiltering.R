#########################################
## Rscript for cell and gene filtering ##
#########################################

# This code aims at
# A : creating seurat objects from cellranger output files for each sequencing run
# B : filtering out low-complexity/quality cells as well as genes detected in too few cells
# C. Creating box plots of QC metrics after cell filtering.

#########################################################################################
# Step A. Open cellranger output files and create Seurat object for each sequencing run #
#########################################################################################

# Define variables
Variables <- list(
  # name of experiment or sample. Will be added to the seurat object as the project name and also to the rds output file name
  SampleName = "Primary_Experiment1", 
  
  # Type of data available
  ## can be GEX_HTO or GEX
  ## If GEX_HTO, the HTO data will be added as an assay to the seurat object
  DataType = "GEX",
  
  # path to the "filtered_feature_bc_matrix" directory of the cellranger output
  ## directory should contain the following files with the exact file name : barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz
  InputDir = "sample_filtered_feature_bc_matrix/",
  
  # path to directory for storing output files
  OutDir = "." 
)

# Load required packages
library(stringr)
library(Seurat)

# Set working directory
setwd(Variables$OutDir)

# Import data
DatasetX <- Read10X(paste0(Variables$InputDir,"/"), gene.column = 2, cell.column = 1, unique.features = TRUE, strip.suffix = FALSE)

# Create Seurat object
if(Variables$DataType == "GEX_HTO"){
  res.data <- CreateSeuratObject(counts = DatasetX$`Gene Expression`, project = Variables$SampleName)
  res.data[["HTO"]] <- CreateAssayObject(counts = DatasetX$`Antibody Capture`)
}else{
  if(Variables$DataType == "GEX"){
    res.data <- CreateSeuratObject(counts = DatasetX, project = Variables$SampleName)
  }
}

# Add prefix to cell barcode to distinguish cells from different experiments/samples with a meaningful ID
res.data <- RenameCells(res.data, add.cell.id=Variables$SampleName)

# Save Seurat object as RDS file
saveRDS(res.data, file = paste0(Sys.Date(), "_SeuratObject_GEX_CellsFrom" , Variables$SampleName, ".rds"))

# clear environment
rm(list = ls()); gc()

#####################################
## Step B. Cell and gene filtering ##
#####################################

# Define variables
Variables <- list(
  # path to directory with seurat object(s) from Step A
  InputDir = "./", 
  
  # pattern to look for in to get RDS files with seurat object(s) that distinguishes experiments/samples
  pattern = "FromPrimary_",
  
  # path to directory for storing output files
  OutDir = "./", 
  
  # suffix that will appear in output file name
  SampleNames = "Primary_Experiments1and2", 
  
  # name of analysis (will also appear in output file name)
  AnalysisName = "CellQualityAssessment", 
  
  # Type of data available
  ## can be GEX_HTO or GEX. If GEX_HTO, the HTO data will also be subsetted.
  DataType = "GEX_HTO", 
  
  # threshold for minimum number of genes to keep a cell in the dataset
  ## set final threshold for each experiment separately after assessing the QC metrics
  nGene_cutoff = 100,
  
  # threshold for maximum percentage of mitochondrial genes detected per cell to keep a cell in the dataset
  PercentMito = 20,
  
  # threshold for minimum number of cells in which a gene should be detected for it to be retained in the dataset
  GeneFilter = 3 
)

# Load required packages
library(Seurat)
library(ggplot2)
library(stringr)
library(gghighlight)

# Set working directory
setwd(Variables$OutDir)

# Get names of RDS file(s) containing Seurat object(s) with gene expression data
temp <- list.files(pattern = ".rds", path = Variables$InputDir)

# Load the Seurat object(s)
for(i in 1:length(temp)){
  assign(gsub(paste0(".*", Variables$pattern, "(.+).rds.*"), "\\1", temp[i]), readRDS(paste0(Variables$InputDir, temp[i])))
}

# Merge data if multiple experiments/samples are being analyzed
nameList <- gsub(paste0(".*", Variables$pattern, "(.+).rds.*"), "\\1", temp)
if(length(nameList) == 1){
  res.data <- get(nameList[1])
}else{
  res.data <- merge(x=get(nameList[1]), y = c(get(nameList[2])), project = "AllSamples")
  
  if(length(nameList) > 2){
    for(i in 3:length(nameList)){
      res.data <- merge(x=res.data, y = c(get(nameList[i])), project = "AllSamples")
    }
  }
}

# Assess percentage of mitochondrial genes detected  per cell
res.data[["percent.mt"]] <- PercentageFeatureSet(res.data, pattern = "^MT-")

# save cell metadata file with all cells and their relevant info before filtering
write.csv(res.data@meta.data, paste0(Sys.Date(), "_", Variables$SampleNames, "_CellMetadata_UnfilteredData.csv"))

# Visualize quality-control (QC) metrics (nCount_RNA, nFeature_RNA, percent.mt) per cell per experiment/sample
## nCount_RNA = total UMI counts detected per cell, nFeature_RNA = number of genes detected per cell
pdf(paste0(Sys.Date(), "_QCmetrics_", Variables$SampleNames, ".pdf"), width = 5, height = 5)
VlnPlot(res.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

ggplot(res.data@meta.data, aes(x=nCount_RNA, y = nFeature_RNA, color = orig.ident)) +
  geom_point(size=1.5) +
  gghighlight(use_group_by = FALSE, use_direct_label = FALSE) +  
  facet_wrap( ~ orig.ident, ncol=2) +
  theme(legend.position = "none")

VlnPlot(res.data, features = "nFeature_RNA", pt.size = 0) +
  theme(legend.position = "none")
res.data@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = Variables$nGene_cutoff)

VlnPlot(res.data, features = "nCount_RNA", pt.size = 0) +
  theme(legend.position = "none")
res.data@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density")

VlnPlot(res.data, features = "percent.mt", pt.size = 0) +
  theme(legend.position = "none")
res.data@meta.data %>% 
  ggplot(aes(color=orig.ident, x=percent.mt, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = Variables$PercentMito)

dev.off()

# Filter out low-complexity cells
## (1) Based on number of genes detected and percentage of mitochondrial genes detected
filtered_res.data <- subset(res.data, subset=(nFeature_RNA >= Variables$nGene_cutoff) & (percent.mt < Variables$PercentMito))

## (2) Based on number of transcripts (UMIs) detected
lower_limit <- mean(filtered_res.data$nCount_RNA) - 3*sd(filtered_res.data$nCount_RNA)
upper_limit <- mean(filtered_res.data$nCount_RNA) + 3*sd(filtered_res.data$nCount_RNA)
filtered_res.data <- subset(filtered_res.data, subset = nCount_RNA < upper_limit & nCount_RNA > lower_limit)

# Filter out genes detected in < X cells
counts <- GetAssayData(filtered_res.data, slot="counts", assay="RNA")   
genes.detection <- rowSums(counts!= 0 )
genes.filter <- names(genes.detection[genes.detection >= Variables$GeneFilter])
filtered_res.data <- subset(filtered_res.data, features = genes.filter)

# Extract HTO data and add to filtered seurat object (if applicable)
filtered_res.data1 <- subset(res.data, cells = rownames(filtered_res.data@meta.data))
if(Variables$DataType == "GEX_HTO"){
  filtered_res.data[["HTO"]] <- filtered_res.data1[["HTO"]]
}
rm(filtered_res.data1);gc(); gc()

# box plot of QC metrics
print(table(filtered_res.data@meta.data$orig.ident))
pdf(file = paste0(Sys.Date(), "_BoxPlot_HighlightingQCmetricsPerSample.pdf"), width = 3, height = 4)
print(ggplot(filtered_res.data@meta.data, aes(x=orig.ident, y=nCount_RNA)) + 
        geom_boxplot(outlier.size = 0.2, outlier.stroke = 0.2, size = 0.5, color = "black", outlier.color = "grey") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        xlab("") + ylab("Total UMI counts") +
        theme_classic() +
        coord_flip()
) 
print(ggplot(filtered_res.data@meta.data, aes(x=orig.ident, y=nFeature_RNA)) + 
        geom_boxplot(outlier.size = 0.2, outlier.stroke = 0.2, size = 0.5, color = "black", outlier.color = "grey") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        xlab("") + ylab("Number genes detected ") +
        theme_classic() +
        coord_flip()
)  
dev.off()

# Save Seurat object for downstream analyses
saveRDS(filtered_res.data, file = paste0(Sys.Date(), "_SeuratObject_GEX_CellsFrom" , Variables$SampleNames, "_FilteredCellsAndGenes.rds"))