# GoT-Multi_Manuscript
This repository contains code for analyses related to the manuscript by Pak M, Saurty-Seerunghen MS et al., "Co-mapping Clonal and Transcriptional Heterogeneity in Somatic Evolution via GoT-Multi."

## 1. CellAndGeneFiltering.R

This code aims at

A. creating seurat objects from cellranger output files for each sequencing run (Step A)

B. filtering out low-complexity/quality cells as well as genes detected in too few cells

C. Creating box plots of QC metrics after cell filtering.


## 2. GEXbasedDoubletPrediction.R

This code aims at identifying doublets based on aberrant gene expression profiles using two independent methods (DoubletFinder & scDblFinder) and retaining only the cells labelled as "doublets" by both methods as potential doublets.


## 3. HTObasedSampleAssignment.R

This code aims at demultiplexing samples and phenotypes (CLL/LBCL) based on HTO information. It has been used for generating Figures S3B and S4B.


## 4. CellClusteringAnalysisUsingSeurat.R

This code aims at

A. Running cell clustering analysis without any prior sample integration, and visualising the results on a UMAP

B. Generating dot plots of expression levels of markers of each cell group

C. Assigning cell group identities, regenerating dot plot of marker gene expression and assessing number of cells per cell group

D. Assessing cell proportion differences/similarities between samples

The code has been used for generating Figures 1C, 4B-C, S1C-D and S3C-E.


## 5. IntegratedCellClusteringAnalysisUsingSeurat.R

This code aims at

A. Running cell clustering analysis after integrating data from different samples, and visualising the results on a UMAP.

B. Generating dot plots of expression levels of markers of each cell group

C. Assigning cell group identities, regenerating dot plot of marker gene expression and assessing number of cells per cell group

D. Assessing cell proportion differences/similarities between samples/phenotypes

E. Evaluating signature scores per cell group (e.g. cell cycle, …)

The code has been used for generating Figures 2A-D, S4C-H.


## 6. CNVinferenceAnalysisUsingInfercnv.R

This code aims at

A. Preparing data (i.e. adding cell group labels to seurat object containing gene expression (GEX) data and extracting data for each patient)

B. Inferring CNVs from cells' transcriptomes of a given patient

Of note, the analysis is done in a loop for each patient separately.

The code has been used for generating Figures S6 (left panels) and S11.


## 7. DifferentialExpressionAnalysisUsingLR.R

This code aims at identifying genes deregulated between mutant and wild type cells for each mutation target using Logistic Regression (LR) method. Only cells from the sample in which the mutation is expected are retained in this analysis.

The code has been used for generating

(1) Figures 3B and 4H, Figure S8 and

(2) Results presented in Tables S5, S6 and S7.


## 8. PathwayEnrichmentAnalysis_PrerankedGSEA.R

This code aims at

A. Evaluating whether certain pathways are more enriched among deregulated genes than among random genes

B. Visualizing results as a heatmap

The code has been used for generating Figure 3D, highlighting pathways in Figure 4H, and for generating results presented in Tables S5 and S6.


## 9. Doublet_Detection.py

This code performs multi-method doublet detection (SOLO, DoubletDetection, Scrublet) on raw count AnnData objects, adds QC metrics (mitochondrial, ribosomal, optional hemoglobin), flags potential outliers using MAD thresholds, and provides a consensus doublet call (`doublet_intersected`).

## 10. GenotypingDataVisualization.py

This code generates (i) UMAP projections colored by genotype calls (MUT / WT / Unprofiled) and (ii) a scatter plot comparing per-cell genotyping (GEN) UMI counts versus transcriptomic (GEX) UMI counts.

## 11. Genotyping_Across_CellStates_Visualization.py

This code computes mutant vs wild-type cell proportions across cell states and cell types, normalizes per-target mutation distributions, and visualizes via heatmaps and per-target stacked MUT/WT barplots.

## 12. MutantCell_Frequency_Heatmap.py

This code constructs genotype frequency matrix for a selected sample, and generates a co-occurrence heatmap of mutant cell frequencies across targets.

## 13. ML_Model_Performance_Visualization.py

This code aggregates genotyping performance across multiple ML models (precision, recall, combined score) and visualizes model performance as heatmaps.


---

Note: For stacked barplots showing genotyping results (proportion distribution of MUT vs. WT) across samples for each target, refer to the "ml_genotyping.py" script of the GoT-Multi-ML repository: https://github.com/AnnaNamLab/GoT-Multi-ML

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17093294.svg)](https://doi.org/10.5281/zenodo.17093294)

