##############################################################################################################
## Differential expression (DE) analysis between mutant (MUT) and wildtype (WT) cells for each mutation     ##
##############################################################################################################

# This code aims at identifying genes deregulated between mutant and wild type cells for each mutation target using Logistic Regression (LR) method.
# Only cells from the sample in which the mutation is expected are retained in this analysis.
# The code has been used for generating
# (1) Figures 3B and 4H, and
# (2) Results presented in Tables S5 and S6.

# Define variables
Variables <- list(
  # path to file containing genotyping status per cell for each mutation target
  ## should have a column containing (1) genotypes, (2) sample in which mutation is expected (column name = "expected_sample")
  mut.path = "genotype.csv",

  # Path to seurat object with gene expression (GEX) data and cell group labels
  InputPath = "../",
  
  # path to genes of interest (e.g. Transcription factors from Lambert et al 2018)
  ## can be NULL or a path.
  ## if NULL, variable genes are used for the analysis. If a path is specified, only genes of interest are used.
  genes.interest = "TFs_Lambert2018.csv",
  
  # name of column containing cell group labels
  col.keep = "CellGroupLabel",
  
  # name of column containing genotypes 
  ## here, using after ML-refinement of IronThrone results
  col.geno = "best_pred_adj",
  
  # group(s) to exclude (here, non-malignant cells)
  ## can be NULL or a string of group labels
  group.exclude = c("Myeloid", "Tcells")
  
)

# Load packages
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(EnhancedVolcano)

# get name of files containing seurat object
temp <- list.files(pattern = ".rds", recursive = TRUE, path=Variables$InputPath)

# Import seurat object
res.data <- readRDS(paste0(Variables$InputPath, temp))
DefaultAssay(res.data) <- "RNA"

# Exclude cell groups not needed for DE analysis
if(length(Variables$group.exclude) != 0){
  metadata <- res.data@meta.data[!(res.data@meta.data[, Variables$col.keep] %in% Variables$group.exclude),]
  res.data <- subset(res.data, cells= rownames(metadata))
  rm(metadata); gc()
}

# import mutation status per cell
mut.data <- read.csv(Variables$mut.path, header = TRUE)
mut.data$GenoData <- mut.data[, Variables$col.geno]

# get names of all mutations
mut.id <- unique(mut.data$target)

# add genotype per cell per target to seurat object
for(j in 1:length(mut.id)){
  
  # subset for one target
  mut.data_target <- mut.data[mut.data$target == mut.id[j],]
  mut.data_target[, paste0(mut.id[j], "_MutStatus")] <- mut.data_target$GenoData 
  
  # keep cells in expected sample
  mut.data_target <- mut.data_target %>% filter(expected_sample==sample)
  rownames(mut.data_target) <- mut.data_target$cell_id
  
  # add to metadata
  res.data <- AddMetaData(res.data, metadata = mut.data_target)
  
  # assign unprofiled cells
  table(res.data@meta.data[, paste0(mut.id[j], "_MutStatus")])
  res.data@meta.data[, paste0(mut.id[j], "_MutStatus")][is.na(res.data@meta.data[, paste0(mut.id[j], "_MutStatus")]) ==TRUE] <- "Unprofiled"
  table(res.data@meta.data[, paste0(mut.id[j], "_MutStatus")])
  
}

# get all targets
mut.id <- names(res.data@meta.data)[grepl("_MutStatus", names(res.data@meta.data))]

for(mut.interest in 1:length(mut.id)){
  
  # get mutation name
  mut.name <- gsub("_MutStatus", "", mut.id[mut.interest])
  
  # subset to keep only cells profiled for that target in expected sample
  res.data@meta.data$Target <- res.data@meta.data[, mut.id[mut.interest]]
  
  # get number of genotyped cells
  nGeno <- res.data@meta.data[res.data@meta.data$Target != "Unprofiled", ]

  if(nrow(nGeno) > 20){ # minimum 20 cells in analysis
    
    # get number of MUT and WT cells
    nMUT <- nrow(nGeno[nGeno$Target == "MUT", ])
    nWT <- nrow(nGeno[nGeno$Target == "WT", ])
    
    if(nWT > 10 & nMUT > 10){
      
      # keep only genotyped cells
      res.data_subset <- subset(res.data, Target != "Unprofiled")

      # Filter out genes detected in < X cells
      counts <- GetAssayData(res.data_subset, slot="counts", assay="RNA")
      genes.detection <- rowSums(counts!= 0 )
      genes.filter <- names(genes.detection[genes.detection >= 3])  #select genes detected in at least X cells
      res.data_subset <- subset(res.data_subset, features = genes.filter)
      rm(counts); rm(genes.detection); rm(genes.filter); gc()
      
      # Define grouping column
      res.data_subset@meta.data$Target <- as.factor(res.data_subset@meta.data$Target)
      res.data_subset <- SetIdent(res.data_subset, value = "Target")
      
      if(length(Variables$genes.interest) != 0){
        # Import genes of interest
        genes.interest <- read.csv("/Users/mss/Desktop/DataUsedForAnalysis/GeneLists/2021-02-23_GeneMetadata_TFs_Lambert2018.csv", header = TRUE, row.names = 1)
        features.use <- rownames(res.data_subset)[rownames(res.data_subset) %in% genes.interest$HGNC.symbol]
      }else{
        # Identify variable genes
        res.data_subset <- FindVariableFeatures(res.data_subset, selection.method = "vst", nfeatures = 2000)
        features.use <- VariableFeatures(res.data_subset)
      }

      # Identify RPL/RPS genes
      rpl_genes <- features.use[startsWith(features.use, "RPL")]
      rps_genes <- features.use[startsWith(features.use, "RPS")]
      rps.rpl_genes <- c(rpl_genes, rps_genes)
      
      # Remove RPL/RPS genes
      features.use <- features.use[!features.use %in% rps.rpl_genes]
      
      # Differential expression analysis between MUT and WT cells across all cell states
      de_results <- FindMarkers(object=res.data_subset, ident.1 = "MUT", ident.2 = "WT", min.pct = 0.1,
                                          logfc.threshold = 0, test.use = "LR", assay = "RNA",
                                          group.by = "Target", features = features.use)
      
      # Adjust p-value
      de_results$BH_adjPval <- p.adjust(de_results$p_val, method = "BH")
      
      # Save DE results
      write.csv(de_results, paste0(Sys.Date(), "_GenesDistinguishing ", mut.name, " MUTfromWTcells_LRmethod.csv"))
      
      # get non-deregulated genes
      nondereg.genes <- rownames(de_results[de_results$BH_adjPval > 0.05, ])
      
      # get top 20 upregulated genes based on adj p-value
      upreg.genes <- de_results[de_results$BH_adjPval < 0.05 & de_results$avg_log2FC >0, ]
      upreg.genes <- upreg.genes[order(upreg.genes$BH_adjPval, decreasing = FALSE),]
      upreg.genes <- rownames(upreg.genes[1:20,])
      
      # get top 20 downregulated genes based on adj p-value
      downreg.genes <- de_results[de_results$BH_adjPval < 0.05 & de_results$avg_log2FC <0, ]
      downreg.genes <- downreg.genes[order(downreg.genes$BH_adjPval, decreasing = FALSE),]
      downreg.genes <- rownames(downreg.genes[1:20,])
      
      # label genes based on their being up- or down-regulated
      de_results$gene.plot[rownames(de_results) %in% nondereg.genes] <- "_"
      de_results$gene.plot[rownames(de_results) %in% upreg.genes] <- "UP"
      de_results$gene.plot[rownames(de_results) %in% downreg.genes] <- "DOWN"
      
      de_results$gene.label <- rownames(de_results)
      de_results$gene.label[!(de_results$gene.plot %in% c("UP", "DOWN"))] <- NA
      de_results$gene.plot[is.na(de_results$gene.plot)==TRUE] <- "Others"
      rm(upreg.genes); rm(downreg.genes); rm(nondereg.genes); gc()
      genes.plot <- de_results$gene.label
      genes.plot <- genes.plot[is.na(genes.plot)==FALSE]
      
      # define colors for gene groups
      col.palette <- c("UP" = "red2", "DOWN" = "blue2", "_"= "lightgrey", "Others"= "black")
      
      # Volcano Plot
      pdf(paste0(Sys.Date(), "_VolcanoPlot_GenesDistinguishing ", mut.name, " MUTfromWTcells_LRmethod.pdf"), width = 5, height = 5)
      print(ggplot(data = de_results, aes(x = avg_log2FC, y = -log10(BH_adjPval), col = gene.plot, label = gene.label)) +
              geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
              geom_point(size = 0.5) +
              scale_color_manual(values = col.palette) +
              scale_y_sqrt() +
              labs(x = "AveragelogFoldChange", y = "-log10(Adj Pval)", title = paste0(mut.name, " MUTfromWTcells")) +
              geom_text_repel(max.overlaps = Inf, position = "identity", direction = "both") +
              theme_classic() +
              theme(legend.position = "bottom")
      )
      dev.off()

    }
  }
}

rm(list = ls()); gc(); gc()
