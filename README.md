# AraLightMeta

# Meta-analysis of high light (HL) response of A. thaliana using our integrative resource AraLightDEGs (https://www.sysbio.ru/aralightdegs/).
![image](https://github.com/user-attachments/assets/8ac4e44c-3627-4c7a-b565-0b0bb2ac304e)

This repository contains a single computation script named AraLightMeta.R and transcriptomic data for A. thaliana leaves and seedlings in high light conditions. 

The analysis starts with data clustering and ends with reconstruction of gene regulatory networks.



# Main objectives:
Comparative analysis of short-term and long-term responses to high light.

Identification of key transcription factors and signaling pathways involved in the HL response.

Detailed visualization of data and results: dendrograms, heatmaps, Sankey diagrams, and regulatory network graphs.

Functional annotation of genes via Gene Ontology (GO) enrichment analysis.

# Technologies and libraries:
The project is implemented in R and utilizes the following packages:
edgeR, dendextend, RColorBrewer, circlize, ggplot2, dplyr, clusterProfiler, org.At.tair.db, purrr, stringr, tidyverse, scales, tidyr, pheatmap, ggvenn, ggraph, tidygraph, igraph, reshape2, gridExtra, WGCNA, biomaRt, broom, GENIE3, DIANE

# Project structure:
base_dir - input data files:

Datasets:

aralightdegs_counts_metadata.csv

aralightdegs_all_degs.csv

ATH_GO_GOSLIM.txt

edges_test_(short/long)_(leaves/seedlings).rds

Main R script for running the full analysis pipeline:

AraLightMeta.R 

# Output folders and files: 
1_step_results folder:

Results of the first stage: initial processing and clustering

Aggregated and normalized expression data

Clustering of experimental conditions

GO enrichment for frequently DEGs


2_step_results_(leaves|seedlings) folder:

Results of the second stage: gene classification and tissue-specific analysis

Gene classification by response type (short-term vs. long-term)

Heatmaps of key functional groups (PhANGs, photomorphogenesis, redoxins, etc.)

Venn diagrams and detailed statistics for specific DEGs


3_step_results_(leaves|seedlings) folder:

Results of the third stage: co-expression and regulatory network reconstruction

WGCNA-based coexpression networks

GENIE3-derived regulatory interactions

Combined GRN tables with annotation

Network visualizations

# Installation and execution:
1. Install Dependencies
   
If required packages are not installed, run the following commands:

Install from CRAN

install.packages(c("edgeR", "dendextend", "RColorBrewer", "circlize",
                   "ggplot2", "pheatmap", "ggraph", "tidygraph", "igraph",
                   "WGCNA", "GENIE3", "tidyverse", "gridExtra", "reshape2",
                   "purrr", "stringr"))

Install Bioconductor packages

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "org.At.tair.db", "biomaRt"))

2. Open the main script (AraLightMeta.R) in RStudio or execute it in an R console.
   
Ensure that all paths are correctly set and that input files are located in their respective directories (only base_dir with input files should be set at line 78, and in this directory, output subdirectories will be created).

To set this analysis to run for our leaves/seedlings datasets, just set the variable analysis on lines 791 or 792. It will calculate steps 2-3 for the required dataset and will save all files in specific folders.
 
# Key results: 
 
Clustering of datasets and visualization of their characteristics.

GO enrichment analysis for core dataset of genes.

Heatmaps of gene expression for relevant genes.

Sankey diagrams for visualization of Gene Ontology categories and processes upregulated/downregulated in short/long high light.

Gene regulatory networks (GRNs) for short- and long-term responses. 

All output files are saved in the subdirectories: 1_step_results, 2_step_results_(leaves/seedlings), and 3_step_results_(leaves/seedlings).

# Key functions: 
 
plot_expression_heatmap()—generates normalized heatmaps of gene expression with a certain minimal percentage threshold of DEGs occurrence in groups. 

create_sankey_go()—constructs Sankey diagrams based on GO terms.

visualize_network()—visualizes regulatory and co-expression networks with metadata: regulation, clusters, and their GO-enriched terms.

run_go_analysis()—performs GO enrichment analysis for specified gene lists.

run_wgcna—modular analysis of gene expression with automatic selection of soft-thresholding parameters.

visualize_network—inference of gene regulatory networks using tree-based methods. 

# Output files examples: 

*.pdf—vector graphics for plots and diagrams. 

*.tsv—tab-separated tables with results from GO analysis, co-expression, and statistical summaries in tables. 

*.Rdata—R environment snapshots for reproducibility. 

For questions or suggestions regarding this project, please open an issue in the repository or contact the author (Aleksandr V. Bobrovskikh) directly via email (avb@bionet.nsc.ru).
