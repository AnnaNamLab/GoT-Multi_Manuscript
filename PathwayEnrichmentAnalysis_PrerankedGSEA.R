######################################################################
## Gene set enrichment analysis (GSEA) using a pre-ranked gene list ##
######################################################################

# This code aims at
# A. Evaluating whether certain pathways are more enriched among deregulated genes than among random genes, and
# B. Visualizing results as a heatmap
# The code has been used for generating Figure 3D, highlighting pathways in Figure 4H, and for generating results presented in Tables S5 and S6.

################################
##      Step A : Run GSEA     ##
################################

# Define variables
Variables <- list(

  ## p-value cutoff for significance level
  pval.threshold = 0.05,
  
  ## Type of pathways to analyze
  ## Options: Hallmark, GO.BP (gmt files for each was downloaded from MSigDB website (https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#H))
  pathways = c("GO.BP", "Hallmark"),
  
  ## path to folder containing the gmt files for each type of pathway to analyze
  gmt.path = "/Users/mss/Desktop/DataUsedForAnalysis/GeneLists/",
 
   ## path to file containing genes of interest (here, DE result table from DifferentialExpressionAnalysisUsingLR.R code)
  pattern = "GenesDistinguishing"
)

# load required packages
library(fgsea)
library(data.table)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)

# Get names of files containing DE results
temp <- list.files(path = "../", pattern = Variables$pattern)
temp.exclude <- temp[grepl("VolcanoPlot", temp)]
temp <- temp[!(temp %in% temp.exclude)]
rm(temp.exclude)

for(file in 1:length(temp)){
  # Import genes
  DE.res_full <- read.csv(paste0("../", temp[file]), header = TRUE, row.names = 1)
  DE.res_full$gene <- rownames(DE.res_full)
  
  # Convert gene symbols to Entrez IDs
  ## use mapIds method to obtain Entrez IDs
  DE.res_full$ENTREZID <- mapIds(org.Hs.eg.db, DE.res_full$gene, 'ENTREZID', 'SYMBOL')
  
  # Remove genes without Entrez ID
  DE.res <- DE.res_full[is.na(DE.res_full$ENTREZID)==FALSE,]

  # Calculate sign of FC * -log10Pval
  ## Get min positive p_val value
  min.p <- min(DE.res$p_val[DE.res$p_val > 0])
  
  ## Add min positive p_val value to each zero-value p_val value to avoid log10(0)
  ### Added only to the p_val=0 since adding it to all makes most of the (sign of FC * -log10Pval)=0
  for(j in 1:nrow(DE.res)){
    if(DE.res[j, "p_val"] == 0){
      DE.res[j, "RankingMeasureGSEA"] <- as.numeric(sign(DE.res[j, "avg_log2FC"]) * -log10(min.p))
    } else {
      DE.res[j, "RankingMeasureGSEA"] <- as.numeric(sign(DE.res[j, "avg_log2FC"]) * -log10(DE.res[j, "p_val"]))
    }
  }
  
  # Create ranks
  ranks <- DE.res$RankingMeasureGSEA
  names(ranks) <- DE.res$ENTREZID
  head(ranks)
  
  # Create prefix for output file names
  out.prefix <- gsub(".*GenesDistinguishing ", "", temp[file])
  out.prefix <- gsub("_LRmethod", "", out.prefix)
  out.prefix <- gsub(".csv", "", out.prefix)
  
  # # Plot the ranked genes
  # pdf(paste0(Sys.Date(), "_", out.prefix, "_RankedGenes.pdf"))
  # barplot(sort(ranks, decreasing = TRUE), axisnames = FALSE)
  # dev.off()
  
  for(pathway.analyze in Variables$pathways){
    # create directory to store results per pathway type
    dir.create(file.path("./", pathway.analyze), showWarnings = FALSE)
    setwd(file.path("./", pathway.analyze))
    
    # define pathways to analyze
    if(pathway.analyze == "Hallmark"){
      pathways.use <- gmtPathways(paste0(Variables$gmt.path, "AllHallmarkGenes_msigdb_v2023.2.Hs_Downloaded20240119fromMSigDB.entrez.gmt"))
      names(pathways.use) <- gsub("HALLMARK_", "", names(pathways.use))
    }
    if(pathway.analyze == "GO.BP"){
      pathways.use <- gmtPathways(paste0(Variables$gmt.path, "GeneOntology_BiologicalProcesses_c5.go.bp.v2023.2.Hs_Downloaded20240119fromMSigDB.entrez.gmt"))
      names(pathways.use) <- gsub("GOBP_", "", names(pathways.use))
    }

    # run GSEA
    set.seed(42)
    fgseaRes <- fgsea(pathways = pathways.use,
                      stats = ranks,
                      eps = 0,
                      minSize=5,
                      maxSize = length(ranks)-1)
    
    if(nrow(fgseaRes)>0){
      # make leading edge more human-readable
      fgseaRes[, leadingEdge := list(mapIdsList(
        x=org.Hs.eg.db,
        keys=leadingEdge,
        keytype="ENTREZID",
        column="SYMBOL"))]
      
      # Identify significantly enriched pathways
      ## padj = BH-adjusted p-value
      fgseaRes.sig <- fgseaRes[fgseaRes$padj < Variables$pval.threshold,] %>% as.data.frame()
      
      # save the results in a text format
      fwrite(fgseaRes, file=paste0(Sys.Date(), "_", out.prefix, "_", pathway.analyze, "_FullTable.txt"), sep="\t", sep2=c("", " ", ""))
    }
    
    setwd("../")
    rm(list = setdiff(ls(), c("Variables", "ranks", "DE.res", "DE.res_full", "pathway.analyze", "temp", "file", "out.prefix")))
  }
    
}

rm(list = ls()); gc(); gc()

####################################################
## Step B : Heatmap to visualize GSEA results     ##
####################################################

# Define variables
Variables <- list(
  # name of analysis (will also appear in output file name)
  analysis.type = "PrerankedGSEA",
  
  # path to directory containing results of pathway enrichment analysis (From StepA)
  input.dir = "./",
  
  # pattern of files containing results of pathway enrichment analysis
  pattern = ".txt", 
  
  # name of pathway database for which to do the heatmap
  pathway.keep = c("Hallmark"),
  
  # list specifying which mutation is expected in which patient sample
  expected.sample = list("CLL-RT3" = c("EMD", "B2M", "POU2F2_Mut1", "MCL1"),
                         "CLL-RT4" = c("PLCG2_Mut1", "PLCG2_Mut2", "PLCG2_Mut3"),
                         "CLL-RT5" = c("JAK1", "PRKDC", "SRRM2"),
                         "CLL-RT6" = c("JUNB", "IKZF3", "SF3B1", "IRF8_Mut1", "HIST1H1C"),
                         "CLL-RT7" =c("BTK")),
  
  # colors to assign to each patient
  sample.color = c("CLL-RT3" = "#6ba04d", "CLL-RT4" = "#8551a5", "CLL-RT5" = "#bda23c", "CLL-RT6" = "#6881d6", "CLL-RT7" = "#b86938")
  
)

# load required packages
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ComplexHeatmap)

# list files containing pathway enrichment analysis results
temp <- list.files(path = Variables$input.dir, pattern = Variables$pattern, recursive = TRUE)
names(temp) <- gsub("/.*", "", temp)
temp <- temp[(names(temp) %in% Variables$pathway.keep)]
temp.exclude <- temp[grepl("AnyMUTvsTripleNeg", temp)]
temp <- temp[!(temp %in% temp.exclude)]

# get mutation name
mut.names <- gsub(".*/", "", temp)
mut.names <- gsub(" MUT.*", "", mut.names)
mut.names <- gsub(".*_", "", mut.names)
names(temp) <- mut.names
mut.names <- unique(mut.names)

# Identify significantly deregulated pathways to keep
sig.pathways <- list()
for(i in 1:length(temp)){
  
  ## import data for one mutation
  pathways1 <- read.delim(paste0(Variables$input.dir, temp[i]), sep = "\t", row.names = 1)
  
  ## keep only significantly enriched pathways
  pathways1 <- pathways1[is.na(pathways1$padj)==FALSE,]
  pathways1 <- pathways1[pathways1$padj < 0.05,]
  
  ## store
  sig.pathways[[mut.names[i]]] <- rownames(pathways1)
}
sig.pathways.all <- unlist(sig.pathways) %>% unique

# Plot enrichment scores
## import data for one mutation
pathways.all <- read.delim(paste0(Variables$input.dir, temp[1]), sep = "\t", row.names = 1)

## keep all pathways of interest
pathways.all <- pathways.all[is.na(pathways.all$padj)==FALSE,]
pathways.all <- pathways.all[rownames(pathways.all) %in% sig.pathways.all, ]

## calculate signNES * -log10AdjPval
pathways.all$PlotValue <- sign(pathways.all$NES) * -(log10(pathways.all$padj))

## keep only columns of interest
pathways.all$pathway <- rownames(pathways.all)
pathways.all <- pathways.all[, c("pathway", "PlotValue")]
names(pathways.all)[2] <- mut.names[1]
rownames(pathways.all) <- NULL

# create duplicate table with only significant terms
mut.name <- names(pathways.all)[2]
sig.pathways.targets <- sig.pathways[[mut.name]]
pathways.sig <- pathways.all
for(j in 1:nrow(pathways.all)){
  if(!(pathways.sig[j, "pathway"] %in% sig.pathways.targets)){
    pathways.sig[j, mut.name] <- NA
  }
}

for(i in 2:length(mut.names)){
  
  ## import data for one mutation
  pathways1 <- read.delim(paste0(Variables$input.dir, temp[i]), sep = "\t", row.names = 1)
  
  ## keep all pathways of interest
  pathways1 <- pathways1[is.na(pathways1$padj)==FALSE,]
  pathways1 <- pathways1[rownames(pathways1) %in% sig.pathways.all, ]
  
  ## calculate signNES * -log10AdjPval
  pathways1$PlotValue <- sign(pathways1$NES) * -(log10(pathways1$padj))
  
  ## keep columns of interest
  pathways1$pathway <- rownames(pathways1)
  pathways1 <- pathways1[, c("pathway", "PlotValue")]
  names(pathways1)[2] <- mut.names[i]
  rownames(pathways1) <- NULL
  
  # merge results for all mutations
  pathways.all <- full_join(x=pathways.all, y=pathways1, by="pathway")
  
  # create duplicate table with only significant terms
  mut.name <- names(pathways1)[2]
  sig.pathways.targets <- sig.pathways[[mut.name]]
  pathways.sig2 <- pathways1
  for(j in 1:nrow(pathways.sig2)){
    if(!(pathways.sig2[j, "pathway"] %in% sig.pathways.targets)){
      pathways.sig2[j, mut.name] <- NA
    }
  }
  pathways.sig <- full_join(x=pathways.sig, y=pathways.sig2, by="pathway")
  
  rm(pathways.sig2); rm(pathways1); gc
}
rm(temp); rm(sig.pathways.targets); rm(mut.name); rm(j); rm(i)

rownames(pathways.all) <- pathways.all$pathway
rownames(pathways.sig) <- pathways.sig$pathway
pathways.all <- pathways.all[, !(names(pathways.all) %in% "pathway")]
pathways.sig <- pathways.sig[, !(names(pathways.sig) %in% "pathway")]

pathways.sig.stars <- pathways.sig
pathways.sig.stars[!(is.na(pathways.sig.stars))] <- "*"
pathways.sig.stars[(is.na(pathways.sig.stars))] <- ""

# create table for mutation annotation
mut.sample <- data.frame(mutation = names(pathways.all))
mut.sample$sample <- NA
for(i in 1:length(Variables$expected.sample)){
  mut.sample$sample[mut.sample$mutation %in% Variables$expected.sample[[i]]] <- names(Variables$expected.sample)[i]
}
rownames(mut.sample) <- mut.sample$mutation
mut.sample$Status <- "MUTvsWT"
mut.sample <- mut.sample[, c("sample", "Status")]
mut.sample <- mut.sample[names(pathways.all1),]

# define annotation colors
my_colour <- list(
  sample = Variables$sample.color,
  Status = c(MUTvsWT = "grey")
)

# change order of mutation
col.order <- unlist(Variables$expected.sample) %>% as.character()
col.order <- col.order[col.order %in% names(pathways.all1)]
pathways.all2 <- pathways.all1[,col.order]
mut.sample2 <- mut.sample[names(pathways.all2),]
pathways.sig.stars <- pathways.sig.stars[rownames(pathways.all2), col.order]

# prepare data for heatmap
pathways.all1 <- pathways.all
pathways.all1[is.na(pathways.all1)] <- 0
my_palette2 <- colorRampPalette(brewer.pal(9, "RdBu"))(12)
breaksList <- seq(-3, 3, by = 0.5)

# change order of mutation
col.order <- c("BTK", "PLCG2_Mut1", "PLCG2_Mut2", "SRRM2", "JUNB", "PRKDC",  "EMD", "POU2F2_Mut1", "IRF8_Mut1", "PLCG2_Mut3", "IKZF3")
pathways.all2 <- pathways.all1[,col.order]
mut.sample2 <- mut.sample[names(pathways.all2),]
pathways.sig.stars <- pathways.sig.stars[rownames(pathways.all2), col.order]

# heatmap
pdf(paste0(Sys.Date(), "_Heatmap_SignificantPathways_", Variables$analysis.type, ".pdf"), width = 5, height = 5)
print(pheatmap(as.matrix(pathways.all2),
               breaks = breaksList,
               color = rev(my_palette2),
               cluster_rows = TRUE,
               clustering_method = "ward.D", # ward.D, single, complete, average, mcquitty
               cluster_cols = TRUE,
               show_rownames = TRUE,
               treeheight_row = 0,
               treeheight_col = 0,
               border_color = "grey",
               na_col = "#FFFFFF",
               scale = "none",
               fontsize = 10,
               fontsize_row = 5,
               annotation_col = mut.sample2,
               annotation_colors = my_colour,
               display_numbers = as.matrix(pathways.sig.stars),
               fontsize_number = 20
    )
  )
dev.off()