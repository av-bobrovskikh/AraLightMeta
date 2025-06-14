# AraLightMeta

🌱 Meta-analysis of high light (HL) response of A. thaliana using our integrative resource data AraLightDEGs (https://www.sysbio.ru/aralightdegs/).

This repository contains single computation script names AraLightMeta.R and transcriptomic data for A. thaliana leaves and seedlings in high light conditions. The analysis starts with data clustering and ends with reconstruction of gene regulatory networks.

🔍 Main Objectives
Comparative analysis of short-term and long-term responses to high light.
Identification of key transcription factors and signaling pathways involved in the HL response.
Detailed visualization of data and results: dendrograms, heatmaps, Sankey diagrams, and regulatory network graphs.
Functional annotation of genes via Gene Ontology (GO) enrichment analysis.

🧰 Technologies and Libraries
The project is implemented in R and utilizes the following packages:
edgeR, dendextend, RColorBrewer, circlize, ggplot2, dplyr, clusterProfiler, org.At.tair.db, purrr, stringr, tidyverse, scales, tidyr, pheatmap, ggvenn, ggraph, tidygraph, igraph, reshape2, gridExtra, WGCNA, biomaRt, broom, GENIE3, DIANE

📁 Project Structure
├── base_dir/                     
│   ├── Input data files            # Raw and precomputed input datasets:
│                                   # - aralightdegs_counts_metadata.csv
│                                   # - aralightdegs_all_degs.csv
│                                   # - ATH_GO_GOSLIM.txt
│                                   # - edges_test_(short/long)_(leaves/seedlings).rds
│
│   └── analysis.R / AraLightMeta.R # Main R script for running the full analysis pipeline
│
├── 1_step_results/                 # Results of the first stage: initial processing and clustering
│                                   # - Aggregated and normalized expression data
│                                   # - Clustering of experimental conditions
│                                   # - GO enrichment for frequently DEGs
│
├── 2_step_results_(leaves|seedlings)/  
│                                   # Results of the second stage: gene classification and tissue-specific analysis
│                                   # - Gene classification by response type (short-term vs. long-term)
│                                   # - Heatmaps of key functional groups (PhANGs, photomorphogenesis, redoxins, etc.)
│                                   # - Venn diagrams and detailed statistical summaries
│
├── 3_step_results_(leaves|seedlings)/  
│                                   # Results of the third stage: co-expression and regulatory network reconstruction
│                                   # - WGCNA-based co-expression networks
│                                   # - GENIE3-derived regulatory interactions
│                                   # - Combined GRN tables with annotation
│                                   # - Network visualizations in PDF format

