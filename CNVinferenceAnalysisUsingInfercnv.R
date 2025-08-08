########################################################################
## Rscript for inferring Copy Number Variations (CNV) scRNA-seq data  ##
########################################################################

# This code aims at
# A. Preparing data (add cell group labels to seurat object containing gene expression (GEX) data and extract data for each patient)
# B. Inferring CNVs from cells' transcriptomes of a given patient

################################
##   Step A : Prepare data    ##
################################

# Define variables
Variables <- list(
  # Path to seurat object after cell and gene filtering containing data from multiple samples (UMI counts as input)
  GEX_path = "SeuratObject_Primary.rds",
  
  # Path to csv file containing cell group labels
  cell.type.annot = "CellMetadataWithCellGroupLabels.csv",
  
  # name of analysis (will also appear in output file name)
  AnalysisName = "GEXdata_Primary",
  
  # threshold for minimum number of cells in which a gene should be detected for it to be retained in the dataset
  GeneFilter = 3,
  
  # name of column containing cell group labels
  group.label.col = "CellGroupLabel",
  
  # group to exclude (e.g. cell group with too few cells per patient)
  ## can be NULL or a string of group names
  group.exclude = c("MyeloidCells")
  
)

# Load required packages
library(Seurat)
library(dplyr)

# Import GEX data
GEX.data <- readRDS(Variables$GEX_path)
GEX.data@meta.data$CellID <- rownames(GEX.data@meta.data)

# Import cell group labels (cell annotations)
metadata <- read.csv(Variables$cell.type.annot, header = TRUE, row.names = 1)
table(metadata[, Variables$group.label.col])

# Add only columns of interest to seurat object
metadata <- metadata[, c("CellID", Variables$group.label.col)]
GEX.data <- AddMetaData(GEX.data, metadata = metadata)
rm(metadata); gc()
GEX.data@meta.data[, Variables$group.label.col][is.na(GEX.data@meta.data[, Variables$group.label.col])==TRUE] <- "_"
table(GEX.data@meta.data[, Variables$group.label.col])
GEX.data$CellGroups <- GEX.data@meta.data[, Variables$group.label.col]

# remove non-annotated cells
GEX.data <- subset(GEX.data, subset = CellGroups != "_")
table(GEX.data$CellGroups)

# exclude specific cell groups
if(length(Variables$group.exclude) !=0){
  metadata <- GEX.data@meta.data[!(GEX.data@meta.data$CellGroups %in% Variables$group.exclude),]
  GEX.data <- subset(GEX.data, cells = rownames(metadata))
  rm(metadata); gc()
}

# identify patient names
samples <- unique(GEX.data@meta.data$HTO_SampleLabel)

for(i in 1:length(samples)){
  
  # subset object for one patient
  GEX.data_subset <- subset(GEX.data, subset = HTO_SampleLabel == samples[i])
  
  # Filter out genes detected in < X cells
  counts <- GetAssayData(GEX.data_subset, slot="counts", assay="RNA")   
  genes.detection <- rowSums(counts!= 0 )
  genes.filter <- names(genes.detection[genes.detection >= Variables$GeneFilter])  #select genes detected in at least X cells
  GEX.data_subset <- subset(GEX.data_subset, features = genes.filter)
  
  # set RNA assay as default assay
  DefaultAssay(GEX.data_subset) <- "RNA"
  
  # save seurat object
  saveRDS(GEX.data_subset, paste0(Sys.Date(), "_SeuratObject_GEXdata_", samples[i], ".rds"))
  
  rm(GEX.data_subset); rm(genes.filter); rm(genes.detection); rm(counts); gc(); gc()
}

rm(list = ls()); gc(); gc()

###############################################################
##   Step B : CNV inference based on transcriptomic data     ##
###############################################################

# Define variables
Variables <- list(
  # pattern indicating that the file contains the seurat object with GEX data (files generate in StepA)
  pattern = ".rds",

  # path to gene annotation file
  ## downloaded from https://data.broadinstitute.org/Trinity/CTAT/cnv/
  gene.ann_path = "gencode_v21_gen_pos.complete.txt",
  
  # length of sliding window for gene smoothing
  genewindowlength=25,
  
  # whether or not to save the matrix used for generating the heatmaps
  heatmap.table.save = TRUE, 
  
  # mode for running the analysis
  ## "subclusters" to partition cells into groups having consistent patterns of CNV or "samples" 
  analysis_mode = "samples", 
  
  # name of cell group(s) to use as reference groups
  ## reference groups = groups with no expected CNVs
  ref_group = c("Tcells"),
  
  # name of column containing cell group labels
  col.keep = "CellGroupLabel"
)

# load required packages
library(Seurat)
library(dplyr)
library(tidyr)
library(infercnv)
library(base)

# get name of files containing seurat objects with GEX data
temp <- list.files(pattern = Variables$pattern)

for(temp.patient in 1:length(temp)){
  
  # get name of patient under study
  patient.id <- gsub(Variables$pattern, "", temp[temp.patient])
  patient.id <- gsub(".*_", "", patient.id)
  
  # Load Seurat object containing GEX data
  seurat.object <- readRDS(temp[temp.patient])
  
  # Create directory for storing results and set it as the working directory
  dir.create(patient.id)
  setwd(patient.id)

  # Extract counts matrix and cell metadata for patient under study
  rawcount <- GetAssayData(object = seurat.object, slot = "counts")
  metadata <- seurat.object@meta.data[colnames(rawcount),]
  metadata$Full.cell_id <- rownames(metadata)
  
  # Get cell annotation as needed by infercnv
  cell.annotation <- metadata[, c("Full.cell_id", Variables$col.keep)]
  write.table(cell.annotation, paste(Sys.Date(), "_CellAnnotation_", patient.id, ".txt"), row.names = FALSE, col.names = FALSE)
  
  # Read the annotation file
  cell.annotation <- read.table(paste(Sys.Date(), "_CellAnnotation_", patient.id, ".txt"), row.names = 1)
  table(cell.annotation$V2)
  
  # Import the gene annotation file and extract relevant info
  genocode <- read.delim(Variables$gene.ann_path, sep="\t", header =FALSE)
  genocode <- genocode %>%
    separate(V1, c("gene", "ID"), sep = "([|])")
  row.names(genocode) <- make.names(genocode$gene, unique = TRUE)
  gene_list <- as.list(genocode$gene)
  
  genocode <- genocode %>%
    dplyr::select(-c(gene, ID))
  genocode$V2 <- as.factor(genocode$V2)
  colnames(genocode) <- NULL
  
  # remove unnecessary files from environment
  rm(list=setdiff(ls(), c("rawcount", "genocode", "cell.annotation", "Variables", "patient.id", "temp", "temp.patient")))
  
  # Create infercnv object
  infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = rawcount,
                                       gene_order_file = genocode,
                                       annotations_file = cell.annotation,
                                       ref_group_names = Variables$ref_group,
                                       chr_exclude = "chrM")
  
  rm(list=setdiff(ls(), c("infercnv_obj", "Variables", "patient.id", "temp", "temp.patient")))
  
  # Specify chromosome order 
  new_gene_order <- data.frame()
  for (chr_name in c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")) {
    new_gene_order <- rbind(new_gene_order, infercnv_obj@gene_order[which(infercnv_obj@gene_order[["chr"]] == chr_name) , , drop=FALSE])
  }
  names(new_gene_order) <- c("chr", "start", "stop")
  infercnv_obj@gene_order <- new_gene_order
  infercnv_obj@expr.data <- infercnv_obj@expr.data[rownames(new_gene_order), , drop=FALSE]
  gc()
  
  # Run infercnv analysis
  infercnv_obj <- infercnv::run(infercnv_obj,
                                cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                out_dir="infercnv_results", # dir for storing outputs
                                cluster_by_groups=TRUE,
                                window_length = Variables$genewindowlength, # length of sliding window for smoothing
                                num_threads = 16,
                                denoise=TRUE,
                                analysis_mode=Variables$analysis_mode,
                                HMM=TRUE, 
                                HMM_type = "i6",
                                write_expr_matrix = Variables$heatmap.table.save)
  
  # save infercnv results
  saveRDS(infercnv_obj, paste0(Sys.Date(), "_InfercnvResults_", patient.id, ".rds"))
  
  setwd("../")
  rm(list = setdiff(ls(), c("Variables", "temp", "temp.patient"))); gc(); gc()
  file.remove(temp[temp.patient])
}

rm(list = ls()); gc(); gc()