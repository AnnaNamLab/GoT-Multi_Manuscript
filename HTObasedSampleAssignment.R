##############################################
## Assign sample identity based on HTO data ##
################################################

# This code aims at demultiplexing samples and phenotypes (CLL/LBCL) based on HTO information.
# It has been used for generating Figure S3B and S4B.

# Define variables
Variables <- list(
  # Path to seurat object after cell and gene filtering (see CellAndGeneFiltering.R code for cell and gene filtering, UMI counts as input)
  ## seurat object should contain data from one experiment
  InputPath = "SeuratObject_Primary.rds",

  # suffix that will appear in output file name
  SampleNames = "Primary_Experiment1", 
  
  # assign which HTOs are from which patient
  HTO.patient = list("CLL2" = c("HTO1_Other", "HTO2_LBCL"), "CLL4" = c("HTO3_Other", "HTO4_LBCL"))

)

# Load required packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)

# Load the Seurat object containing the GEX data (after filtering genes and cells)
res.data <- readRDS(Variables$InputPath)

# Extract HTO data from seurat object
HTO.data <- GetAssayData(res.data, slot="counts", assay="HTO")   
HTO.data <- as.data.frame(HTO.data)
HTO.data <- HTO.data[,rownames(res.data@meta.data)]
HTO.data <- as.data.frame(t(HTO.data))

# Add HTO data to metadata table
HTO.data <- HTO.data[rownames(res.data@meta.data),]
res.data@meta.data <- cbind(res.data@meta.data, HTO.data)

# log-transform HTO UMI values
HTO.col <- names(res.data@meta.data)[grepl("HTO", names(res.data@meta.data))]
HTO.col <- HTO.col[grepl("-", HTO.col)]
for(i in 1:length(HTO.col)){
  res.data@meta.data[, paste0(HTO.col[i], "_log10")] <- log10(res.data@meta.data[,HTO.col[i]] +1)
}

# Sum up UMI counts for HTOs from the same patient sample
patient.ids <- names(Variables$HTO.patient)
  
  for(i in 1:length(patient.ids)){
    
    # extract relevant HTO
    HTO.keep <- Variables$HTO.patient[[i]]
    
    res.data@meta.data[, paste0("TotalHTO_", patient.ids[i])] <- res.data@meta.data[, HTO.keep[1]] + res.data@meta.data[, HTO.keep[2]]
    rm(HTO.keep); gc()
  }
  

HTO.col <- names(res.data@meta.data)[grepl("TotalHTO_", names(res.data@meta.data))]
for(i in 1:length(HTO.col)){
  res.data@meta.data[, paste0(HTO.col[i], "_log10")] <- log10(res.data@meta.data[,HTO.col[i]] +1)
}

# Manually select sample1 versus sample2 (here, CLL2 versus CLL4 cells as an example)
  q1.table <- data.frame(x=res.data@meta.data[, paste0("TotalHTO_", patient.ids[1], "_log10")], y=res.data@meta.data[, paste0("TotalHTO_", patient.ids[2], "_log10")], row.names = rownames(res.data@meta.data))
  
  ## Mark cells from sample1
  X11()
  plot(q1.table, pch=16)
  Sample1.cells <- gatepoints::fhs(q1.table, mark = TRUE)
  Sample1.cells <- as.data.frame(Sample1.cells)
  names(Sample1.cells)[1] <- "CellID"
  Sample1.cells$HTO_SampleLabel <- patient.ids[1]
  
  ## Mark cells from sample2
  X11()
  plot(q1.table, pch=16)
  Sample2.cells <- gatepoints::fhs(q1.table, mark = TRUE)
  Sample2.cells <- as.data.frame(Sample2.cells)
  names(Sample2.cells)[1] <- "CellID"
  Sample2.cells$HTO_SampleLabel <- patient.ids[2]
  BothSamples.cells <- rbind(Sample1.cells, Sample2.cells)
  rm(q1.table)
  
  # extract cell barcodes from each sample
  metadata.Sample1 <- res.data@meta.data[rownames(res.data@meta.data) %in% Sample1.cells$CellID,]
  metadata.Sample2 <- res.data@meta.data[rownames(res.data@meta.data) %in% Sample2.cells$CellID,]
  rm(Sample1.cells); rm(Sample2.cells)

# Manually select CLL versus DLBCL cells from each patient sample
    ## extract relevant HTO for sample1
    HTO.keep <- Variables$HTO.patient[[1]]
    q1.table <- data.frame(x=res.data@meta.data[, paste0(HTO.keep[1], "_log10")], y=res.data@meta.data[, paste0(HTO.keep[2], "_log10")], row.names = rownames(res.data@meta.data))
    
    ## Mark CLL cells from sample1
    X11()
    plot(q1.table, pch=16)
    CLL.cells <- gatepoints::fhs(q1.table, mark = TRUE)
    CLL.cells <- as.data.frame(CLL.cells)
    names(CLL.cells)[1] <- "CellID"
    CLL.cells$HTO_Phenotype <- "CLL"
    
    ## Mark LBLC cells from sample1
    X11()
    plot(q1.table, pch=16)
    DLBCL_T.cells <- gatepoints::fhs(q1.table, mark = TRUE)
    DLBCL_T.cells <- as.data.frame(DLBCL_T.cells)
    names(DLBCL_T.cells)[1] <- "CellID"
    DLBCL_T.cells$HTO_Phenotype <- "DLBCL_T"
    CLL_DLBCL.cells <- rbind(CLL.cells, DLBCL_T.cells)
    rm(CLL.cells); rm(DLBCL_T.cells); rm(q1.table)
    
    ## extract relevant HTO for sample2
    HTO.keep <- Variables$HTO.patient[[2]]
    q1.table <- data.frame(x=res.data@meta.data[, paste0(HTO.keep[1], "_log10")], y=res.data@meta.data[, paste0(HTO.keep[2], "_log10")], row.names = rownames(res.data@meta.data))
    
    ## Mark CLL cells from sample2
    X11()
    plot(q1.table, pch=16)
    CLL.cells <- gatepoints::fhs(q1.table, mark = TRUE)
    CLL.cells <- as.data.frame(CLL.cells)
    names(CLL.cells)[1] <- "CellID"
    CLL.cells$HTO_Phenotype <- "CLL"
    
    ## Mark LBLC cells from sample2
    X11()
    plot(q1.table, pch=16)
    DLBCL_T.cells <- gatepoints::fhs(q1.table, mark = TRUE)
    DLBCL_T.cells <- as.data.frame(DLBCL_T.cells)
    names(DLBCL_T.cells)[1] <- "CellID"
    DLBCL_T.cells$HTO_Phenotype <- "DLBCL_T"
    CLL_DLBCL.cells <- rbind(CLL_DLBCL.cells, CLL.cells, DLBCL_T.cells)
    rm(CLL.cells); rm(DLBCL_T.cells); rm(q1.table)
    
    ## merge tables
    BothSamples.cells <- full_join(x=BothSamples.cells, y=CLL_DLBCL.cells, by="CellID")
    rownames(BothSamples.cells) <- BothSamples.cells$CellID
    rm(CLL_DLBCL.cells); rm(metadata.Sample1); rm(metadata.Sample2)
    
# Add HTO-based sample assignment to Seurat object
res.data <- AddMetaData(res.data, metadata = BothSamples.cells)
table(res.data@meta.data$HTO_SampleLabel, res.data@meta.data$HTO_Phenotype)

# Assign doublet status
res.data@meta.data$HTO_SampleLabel[is.na(res.data@meta.data$HTO_Phenotype)==TRUE] <- "doublet"
res.data@meta.data$HTO_Phenotype[is.na(res.data@meta.data$HTO_Phenotype)==TRUE] <- "doublet"
table(res.data@meta.data$HTO_SampleLabel, res.data@meta.data$HTO_Phenotype)

metadata.full <- res.data@meta.data

# plot HTO data after sample assignment
pdf(paste0(Sys.Date(), "_HTO UMI values after sample assignment.pdf"))

ggplot(res.data@meta.data, aes(x=res.data@meta.data[,paste0("TotalHTO_", patient.ids[1], "_log10")], y=res.data@meta.data[,paste0("TotalHTO_", patient.ids[2], "_log10")], col = HTO_SampleLabel)) + 
  geom_point() +
  xlab(HTO.keep[1]) + ylab(HTO.keep[2])

ggplot(res.data@meta.data, aes(x=res.data@meta.data[,paste0("TotalHTO_", patient.ids[1], "_log10")], y=res.data@meta.data[,paste0("TotalHTO_", patient.ids[2], "_log10")], col = HTO_Phenotype)) + 
  geom_point() +
  xlab(HTO.keep[1]) + ylab(HTO.keep[2])

HTO.keep <- Variables$HTO.patient[[1]]
ggplot(res.data@meta.data, aes(x=res.data@meta.data[,paste0(HTO.keep[1], "_log10")], y=res.data@meta.data[,paste0(HTO.keep[2], "_log10")], col = HTO_Phenotype)) + 
  geom_point() +
  xlab(HTO.keep[1]) + ylab(HTO.keep[2]) + ggtitle(names(Variables$HTO.patient)[1])

HTO.keep <- Variables$HTO.patient[[2]]
ggplot(res.data@meta.data, aes(x=res.data@meta.data[,paste0(HTO.keep[1], "_log10")], y=res.data@meta.data[,paste0(HTO.keep[2], "_log10")], col = HTO_Phenotype)) + 
  geom_point() +
  xlab(HTO.keep[1]) + ylab(HTO.keep[2]) + ggtitle(names(Variables$HTO.patient)[2])

dev.off()

write.csv(res.data@meta.data, paste0(Sys.Date(), "_CellMetadataWithHTObasedDoublets.csv"))

# Filter out HTO-based doublets
res.data <- subset(res.data, subset = HTO_SampleLabel != "doublet")
dim(res.data)
table(res.data@meta.data$HTO_SampleLabel, res.data@meta.data$HTO_Phenotype)

# keep only genes detected in at least 3 cells
selected_genes <- rownames(res.data)[Matrix::rowSums(res.data) > 3] 
res.data.filt <- subset(res.data, features = selected_genes)
dim(res.data.filt)

# Save Seurat object for downstream analyses
saveRDS(res.data.filt, file = paste0(Sys.Date(), "_SeuratObject_GEX_CellsFrom" , Variables$SampleNames, "_FilteredCellsAndGenesWithHTOsampleAssignment.rds"))