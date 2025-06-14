# Meta-analysis of Arabidopsis thaliana gene expression under high light conditions
#
# This script performs a comprehensive analysis of Arabidopsis gene expression data 
# under high light stress conditions in leaves and tissues. It includes:
# 1. Data loading and preprocessing
# 2. Expression data aggregation and normalization
# 3. Correlation analysis and clustering of experimental conditions
# 4. Gene Ontology enrichment analysis of frequently differentially expressed genes
# 5. Classification of long- and short-specific genes 
# 6. Gene regulatory network reconstruction combining both co-expression (WGCNA) and machine-learning approaches (GENIE3)
#
# Required input files in base directory (base_dir):
# - aralightdegs_counts_metadata.csv: Raw expression data with metadata
# - aralightdegs_all_degs.csv: Precomputed DEGs data from AraLightDEGs database
# - ATH_GO_GOSLIM.txt: Gene Ontology processes (https://www.arabidopsis.org/download/list?dir=Public_Data_Releases/TAIR_Data_20231231)
# - 4 files: edges_test_(long/short)_(leaves/seedlings).rds with precomputed edges of GENIE3 with p-value (DIANE)

# Computation results are saved in subdirectories:
# - 1_step_results/: Initial analysis results (clustering, GO enrichment)
# - 2_step_results_(seedlings/leaves)/: Gene classification and detailed analysis for leaves and seedlings datasets
# - 3_step_results_(seedlings/leaves)/: Co-expression and gene regulatory network analysis for leaves and seedlings datasets

# Check and install required packages
required_cran <- c("edgeR", "dendextend", "RColorBrewer", "circlize", "ggplot2", 
                  "dplyr", "purrr", "stringr", "tidyverse", "scales", "tidyr",
                  "pheatmap", "ggvenn", "ggraph", "tidygraph", "igraph", "reshape2",
                  "gridExtra", "WGCNA", "biomaRt", "broom")

required_bioc <- c("clusterProfiler", "org.At.tair.db", "GENIE3")

# Install CRAN packages
for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# Load libraries
library(edgeR)
library(dendextend)
library(RColorBrewer)
library(circlize)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.At.tair.db)
library(purrr)
library(stringr)
library(tidyverse)
library(scales)
library(tidyr)
library(pheatmap)
library(ggvenn)
library(ggraph)
library(tidygraph)
library(igraph)
library(reshape2)
library(gridExtra)
library(WGCNA)
library(biomaRt)
library(broom)
library(GENIE3)

# Enable multi-threading for WGCNA
WGCNA::enableWGCNAThreads()

# Define paths
base_dir <- "path/to/your/working/dir/with/requiring/files"
output_dir_1 <- file.path(base_dir, "1_step_results")

# Create output directories
dir.create(output_dir_1, showWarnings = FALSE)

# FIRST STEP OF ANALYSIS
# GENE CLASSIFICATION & DETAILED ANALYSIS OF GENE ONTOLOGY GROUPS

# DATA READING AND PREPARATION 

# Read and process metadata
data_file <- "aralightdegs_counts_metadata.csv"
data_path <- file.path(base_dir, data_file)

metadata_raw <- read.csv(data_path, header = FALSE, nrows = 6, stringsAsFactors = TRUE)
metadata <- as.data.frame(t(metadata_raw))
colnames(metadata) <- c("EXPERIMENT", "TISSUE", "CONDITION", "LIGHT_INTENSITY", "TIME", "AGE")
metadata <- metadata[-1, ]

# Process metadata
metadata$EXPERIMENT <- factor(metadata$EXPERIMENT)
metadata$TISSUE <- factor(metadata$TISSUE)
metadata$CONDITION <- factor(ifelse(grepl("CONTROL", metadata$CONDITION), "CONTROL", "HIGH_LIGHT"))
metadata$TIME <- as.numeric(metadata$TIME)
metadata$LIGHT_INTENSITY <- as.numeric(gsub("[^0-9]", "", metadata$LIGHT_INTENSITY))
metadata$AGE <- as.numeric(gsub("[^0-9]", "", metadata$AGE))

# Read expression data
gene_expression <- read.csv(data_path, skip = 6, header = TRUE)
rownames(gene_expression) <- gene_expression$GeneID
gene_counts <- gene_expression[, -1]
rownames(metadata) <- colnames(gene_counts)

# DATA AGGREGATION 

# Create FULL_DATA with underscores as separators
metadata$FULL_DATA <- paste(
  metadata$EXPERIMENT,
  metadata$TISSUE,
  ifelse(metadata$CONDITION == "CONTROL", "C", "HL"),
  paste0(metadata$LIGHT_INTENSITY, "PPFD"),
  paste0(metadata$TIME, "MIN"),
  paste0(metadata$AGE, "D"),
  sep = "_"
)

# Aggregate by identical conditions
unique_conditions <- unique(metadata$FULL_DATA)
unique_condition_expression <- data.frame(row.names = rownames(gene_counts))

for (cond in unique_conditions) {
  samples <- rownames(metadata)[metadata$FULL_DATA == cond]
  if (length(samples) > 1) {
    unique_condition_expression[[cond]] <- rowSums(gene_counts[, samples, drop = FALSE])
  } else {
    unique_condition_expression[[cond]] <- gene_counts[, samples]
  }
}

write.table(unique_condition_expression, file.path(output_dir_1, "1.Condition_expression_aggregated.tsv"), sep="\t", quote=FALSE)

# NORMALIZATION AND CORRELATION ANALYSIS

normalize_counts <- function(counts) {
  cpm <- cpm(counts)
  mean_cpm <- rowMeans(cpm)
  filtered_cpm <- cpm[mean_cpm >= 3, , drop = FALSE]
  return(filtered_cpm)
}

# CPM normalization
dge_unique <- DGEList(counts = unique_condition_expression)
dge_unique <- calcNormFactors(dge_unique)
cpm_unique_filtered <- normalize_counts(dge_unique)

write.table(round(as.data.frame(cpm_unique_filtered), 2), file.path(output_dir_1, "1.Condition_expression_aggregated_cpm_filtered.tsv"), sep="\t", quote=FALSE)

# Calculate correlation matrix
cor_matrix_unique <- cor(cpm_unique_filtered, method = "pearson")
write.table(round(as.data.frame(cor_matrix_unique), 2), file.path(output_dir_1, "1.Correlation_matrix_all_conditions.tsv"), sep="\t", quote=FALSE)

# Convert to dissimilarity
dissimilarity <- 1 - abs(cor_matrix_unique)
dist_dissimilarity <- as.dist(dissimilarity)

# CLUSTERING AND VISUALIZATION 

metadata$TISSUE_GROUPED <- ifelse(
  metadata$TISSUE %in% c("LOCAL LEAF", "SYSTEMIC LEAF"),
  "LEAF",
  as.character(metadata$TISSUE)
)

# Create dendrogram
hc <- hclust(dist_dissimilarity, method = "ward.D2")
dend <- as.dendrogram(hc)

# Get dendrogram leaf order
dend_order <- order.dendrogram(dend)
dend_labels <- labels(dend)

# Prepare metadata for annotations
condition_info_dendr <- metadata[match(dend_labels, metadata$FULL_DATA), 
                          c("TISSUE_GROUPED", "CONDITION", "TIME", "FULL_DATA")]

# Verify order
stopifnot(all(dend_labels == condition_info_dendr$FULL_DATA))

# Color schemes
tissue_colors <- c("SEEDLING" = "#1f78b4", "LEAF" = "#33a02c")

# Function for condition colors
sample_colors <- function(time, condition) {
  colors <- character(length(time))
  for (i in seq_along(time)) {
    if (condition[i] == "CONTROL") {
      colors[i] <- "black"
    } else if (is.na(time[i])) {
      colors[i] <- "white"
    } else if (time[i] <= 90) {
      colors[i] <- "#FFD700"     # orange - Short-term
    } else if (time[i] >= 120 & time[i] <= 960) {
      colors[i] <- "#FFA500"     # yellow - Mid-term
    } else if (time[i] >= 1440) {
      colors[i] <- "red"         # red - Long-term
    } else {
      colors[i] <- "white"
    }
  }
  return(colors)
}

# Create color data frame
sample_col <- sample_colors(condition_info_dendr$TIME, condition_info_dendr$CONDITION)
tissue_col <- tissue_colors[as.character(condition_info_dendr$TISSUE_GROUPED)]
color_df <- data.frame(
  tissue = tissue_col,
  condition = sample_col,
  row.names = condition_info_dendr$FULL_DATA
)

# Configure dendrogram
dend <- dend %>%
  set("labels_cex", 0.1) %>%
  set("branches_lwd", 1.2) %>%
  set("leaves_pch", 19) %>%
  set("leaves_cex", 0.25) %>%
  set("leaves_col", "gray40") %>%
  hang.dendrogram(hang_height = 0.35)

# Create PDF output
pdf(file.path(output_dir_1, "1.Clustering_dendrogram.pdf"), 
    width = 8, height = 6)

layout(matrix(c(1, 2), nrow = 2, byrow = TRUE), heights = c(1, 1))

# Top panel: dendrogram
par(mar = c(0, 4, 4, 1))
plot(dend,
     main = "Clustering of Experimental Conditions",
     xlab = "", ylab = "Ward distance",
     horiz = FALSE,
     axes = FALSE,
     labels = FALSE)

y_range <- par("yaxp")
y_labels <- seq(floor(y_range[1]), ceiling(y_range[2]), by = 1)
axis(side = 2, at = y_labels, labels = y_labels)

box()

# Legends
legend("right",
       legend = c("Seedlings", "Leaves"),
       fill = c("#1f78b4", "#33a02c"),
       title = "Tissue",
       cex = 0.685)

legend("topright",
       legend = c("Control", "Short-term HL", "Mid-term HL", "Long-term HL"),
       fill = c("black", "#FFD700", "#FFA500", "red"),
       title = "Conditions",
       cex = 0.685)

# Bottom panel: color annotations
par(mar = c(2, 4, 0, 1))
plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")

n <- length(dend_order)
x_coords <- seq(0, 1, length.out = n + 1)

# Tissue annotation
rect(xleft = x_coords[-(n+1)],
     xright = x_coords[-1],
     ybottom = 0.9, ytop = 1.0,
     col = color_df$tissue,
     border = NA)

# Condition annotation
rect(xleft = x_coords[-(n+1)],
     xright = x_coords[-1],
     ybottom = 0.8, ytop = 0.9,
     col = color_df$condition,
     border = NA)

# Sample labels
x_centers <- (x_coords[-1] + x_coords[-(n+1)]) / 2
text(x = x_centers,
     y = 0.8,
     labels = dend_labels,
     srt = 90,
     adj = 1,
     cex = 0.3,
     xpd = NA)

dev.off()

# METADATA VISUALIZATION
condition_info <- metadata[match(dend_labels, metadata$FULL_DATA), 
                          c("TISSUE_GROUPED", "CONDITION", "TIME", "FULL_DATA", "EXPERIMENT","LIGHT_INTENSITY", "AGE")]

# Calculate counts for all samples
metadata_counts_intensity_time <- condition_info %>%
  filter(!str_detect(EXPERIMENT, "PRJNA699408|GSE251796")) %>%
  filter(CONDITION == "HIGH_LIGHT") %>%
  group_by(LIGHT_INTENSITY, TIME) %>%
  summarise(count = n(), .groups = 'drop') %>%
  arrange(LIGHT_INTENSITY, TIME)

# Calculate counts only for LEAF samples
metadata_counts_intensity_time_leaf <- condition_info %>%
  filter(!str_detect(EXPERIMENT, "PRJNA699408|GSE251796")) %>%
  filter(CONDITION == "HIGH_LIGHT", TISSUE_GROUPED == "LEAF") %>%
  group_by(LIGHT_INTENSITY, TIME) %>%
  summarise(leaf_count = n(), .groups = 'drop') %>%
  arrange(LIGHT_INTENSITY, TIME)

# Merge the two dataframes
metadata_counts_combined <- metadata_counts_intensity_time %>%
  left_join(metadata_counts_intensity_time_leaf, by = c("LIGHT_INTENSITY", "TIME")) %>%
  mutate(leaf_count = ifelse(is.na(leaf_count), 0, leaf_count))

# Add time & intensity grouping
metadata_counts_combined <- metadata_counts_combined %>%
  mutate(
    intensity_group = case_when(
      LIGHT_INTENSITY < 1000 ~ "Moderate",
      LIGHT_INTENSITY >= 1000 & LIGHT_INTENSITY < 1500 ~ "Intensive",
      LIGHT_INTENSITY >= 1500 ~ "Severe"
    ),
    time_group = case_when(
      TIME <= 90 ~ "Short",
      TIME >= 120 & TIME <= 960 ~ "Medium",
      TIME >= 1440 ~ "Long-term"
    )
  )

# Data aggregation for heatmap
intensity_time_counts <- metadata_counts_combined %>%
  group_by(intensity_group, time_group) %>%
  summarise(total_count = sum(count), .groups = "drop") %>%
  pivot_wider(names_from = time_group, values_from = total_count, values_fill = list(total_count = 0))

# Factor arrangement
intensity_time_counts$intensity_group <-
  factor(intensity_time_counts$intensity_group,
         levels = c("Moderate", "Intensive", "Severe"))

# Long format for ggplot conversion
heatmap_data <- intensity_time_counts %>%
  pivot_longer(cols = -intensity_group, names_to = "time_group", values_to = "count") %>%
  mutate(time_group = factor(time_group,
                             levels = c("Short", "Medium", "Long-term")))

# Draw heatmap
heatmap_occur <- ggplot(heatmap_data, aes(x = time_group, y = intensity_group, fill = count)) +
  geom_tile(color = "black") +
  scale_fill_gradient(
    low = "white", 
    high = "#006d2c", 
    na.value = "gray",
    name = "Count"
  ) +
  geom_text(aes(label = count), size = 4) +
  labs(
    title = "",
    x = "Time group",
    y = "Light intensity group"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title.x = element_text(color = "black", size = 10),
    axis.title.y = element_text(color = "black", size = 10),
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),  
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(clip = "off") +
  guides(fill = guide_colourbar(title = "Count", barwidth = 0.5, barheight = 5))

ggsave(
  filename = file.path(output_dir_1, "1.Experimental_conditions_heatmap_Time_intensity.pdf"),
  plot = heatmap_occur,
  width = 4,
  height = 4,
  dpi = 300,
  limitsize = FALSE
)

# Create plot with both all samples (light gray) and leaf samples (black)
metadata_plot <- ggplot(metadata_counts_combined, aes(x = TIME, y = LIGHT_INTENSITY)) +
  # All samples - light gray points
  geom_point(aes(size = count), color = "lightgray", alpha = 0.7) +
  # Leaf samples - black points
  geom_point(data = subset(metadata_counts_combined, leaf_count > 0),
             aes(size = leaf_count), color = "black", alpha = 0.9) +
  scale_size_continuous(range = c(3, 10),
                       breaks = seq(min(metadata_counts_combined$count), 
                                   max(metadata_counts_combined$count), 
                                   by = 1)) +
  scale_x_log10(breaks = c(1, 10, 30, 60, 120, 240, 480, 720, 1440, 2880, 4320, 7200),
                labels = c("1", "10", "30", "60", "120", "240", "480", "720", "1440", "2880", "4320", "7200")) +
  scale_y_continuous(breaks = seq(0, max(metadata_counts_combined$LIGHT_INTENSITY), by = 500)) +
  labs(
    x = expression("TIME (minutes, log"[10]*" scale)"), 
    y = expression("LIGHT INTENSITY ("*mu*"mol*m"^-2*"*s"^-1*")"),
    size = "Number of\nexperimental groups") +
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Save plot
ggsave(
  filename = file.path(output_dir_1, "1.Experimental_conditions_counts.pdf"),
  plot = metadata_plot,
  width = 10,
  height = 8,
  dpi = 300,
  limitsize = FALSE
)

# Data filtering - excluding PRJNA699408 и GSE251796
filtered_unique_conditions <- unique_condition_expression[, !grepl("PRJNA699408|GSE251796", colnames(unique_condition_expression))]

# Creating classification function
classify_sample <- function(sample_name) {
  ppfd <- as.numeric(str_extract(sample_name, "(?<=_)\\d+(?=PPFD)"))
  
  time <- as.numeric(str_extract(sample_name, "(?<=_)\\d+\\.?\\d*(?=MIN)"))
  
  if (grepl("SEEDLING", sample_name)) {
    tissue <- "SEEDLING"
  } else if (grepl("LEAF", sample_name)) {
    tissue <- "LEAF"
  } else {
    return(NA)
  }
  
  if (ppfd < 200) {
    return(paste0("CONTROL_", tissue))
  } else if (ppfd >= 500) {
    if (time <= 90) {
      return(paste0("SHORT_HL_", tissue))
    } else if (time >= 120) {
      return(paste0("LONG_HL_", tissue))
    }
  }
  return(NA) 
}

sample_groups <- sapply(colnames(filtered_unique_conditions), classify_sample)
sample_groups_df <- enframe(sample_groups, name = "Sample", value = "Group")

# Samples grouping and merge counts
grouped_counts <- data.frame(GeneID = rownames(filtered_unique_conditions))

for (group in unique(sample_groups)) {
  group_samples <- names(sample_groups)[sample_groups == group]
  if (length(group_samples) > 0) {
    if (length(group_samples) > 1) {
      grouped_counts[[group]] <- rowSums(filtered_unique_conditions[, group_samples, drop = FALSE])
    } else {
      grouped_counts[[group]] <- filtered_unique_conditions[, group_samples]
    }
  }
}

# Remove GeneID and convert to matrix
count_matrix_groups <- as.matrix(grouped_counts[, -1])
rownames(count_matrix_groups) <- grouped_counts$GeneID

# CPM normalization
dge_counts_groups <- DGEList(counts = count_matrix_groups)
dge_counts_groups <- calcNormFactors(dge_counts_groups)
cpm_data_groups <- normalize_counts(dge_counts_groups)

# Correlation matrix
cor_matrix_groups <- cor(cpm_data_groups, method = "spearman")

# Ordering samples
desired_order <- c("CONTROL_LEAF", "CONTROL_SEEDLING", 
                   "SHORT_HL_LEAF", "SHORT_HL_SEEDLING", 
                   "LONG_HL_LEAF", "LONG_HL_SEEDLING")

desired_order <- desired_order[desired_order %in% colnames(cor_matrix_groups)]

# Short names
short_names <- c("C.LEAVES", "C.SEEDLINGS", "S.LEAVES", "S.SEEDLINGS", "L.LEAVES", "L.SEEDLINGS")
names(short_names) <- c("CONTROL_LEAF", "CONTROL_SEEDLING", 
                        "SHORT_HL_LEAF", "SHORT_HL_SEEDLING", 
                        "LONG_HL_LEAF", "LONG_HL_SEEDLING")

short_names <- short_names[desired_order]

# Rename corr. matrix
cor_matrix_groups <- cor_matrix_groups[desired_order, desired_order]
colnames(cor_matrix_groups) <- short_names[desired_order]
rownames(cor_matrix_groups) <- short_names[desired_order]

# Heatmap visualization 
cor_groups_htmp <-pheatmap(
  cor_matrix_groups,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "steelblue"))(100),
  main = "Correlation between groups",
  display_numbers = TRUE,
  number_format = "%.2f",
  number_color = "black",
  border_color = "black",
  fontsize = 14,
  fontsize_number = 25,
  angle_col = 0,
  legend = TRUE
)

ggsave(
  filename = file.path(output_dir_1, "1.Experimental_conditions_groups_correlations.pdf"),
  plot = cor_groups_htmp,
  width = 10,
  height = 8,
  dpi = 300,
  limitsize = FALSE
)

write.csv(file.path(output_dir_1, cor_matrix_groups), "1.Correlation_matrix_conditions.csv")
write.csv(file.path(output_dir_1, cpm_data_groups), "1.Filtered_cpm_data.csv")


metadata_counts_age_tissue <- condition_info %>%
  filter(!str_detect(EXPERIMENT, "PRJNA699408|GSE251796")) %>%
  filter(CONDITION == "HIGH_LIGHT") %>%
  mutate(
    TISSUE_GROUP = case_when(
      TISSUE_GROUPED == "SEEDLING" ~ "Seedling",
      TRUE ~ "Leaf"
    ),
    AGE_GROUP = case_when(
      AGE >= 7 & AGE <= 15 ~ "Juvenile",
      AGE >= 21 & AGE <= 30 ~ "Mature",
      AGE >= 32 & AGE <= 49 ~ "Old",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(AGE_GROUP, TISSUE_GROUP) %>%
  summarise(count = n(), .groups = 'drop') %>%
  filter(AGE_GROUP != "Other") %>%
  arrange(AGE_GROUP, TISSUE_GROUP)

# Long format for ggplot2
heatmap_data_age_tissue <- metadata_counts_age_tissue %>%
  mutate(
    AGE_GROUP = factor(AGE_GROUP, levels = c("Juvenile", "Mature", "Old")),
    TISSUE_GROUP = factor(TISSUE_GROUP, levels = c("Seedling", "Leaf")))
names(heatmap_data_age_tissue) <- c("age_group", "tissue_group", "count")

# heatmap age-tissue
heatmap_age_tissue <- ggplot(
  heatmap_data_age_tissue,
  aes(x = age_group, y = tissue_group, fill = count)
) +
  geom_tile(color = "black") +
  scale_fill_gradient(
    low = "white",
    high = "#006d2c",
    na.value = "gray",
    name = "Count"
  ) +
  geom_text(aes(label = count), size = 4) +
  labs(
    title = "",
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(clip = "off") +
  guides(fill = guide_colourbar(title = "Count", barwidth = 0.5, barheight = 5))

ggsave(
  filename = file.path(output_dir_1, "1.Experimental_conditions_heatmap_Age_Tissue.pdf"),
  plot = heatmap_age_tissue,
  width = 4,
  height = 2.5,
  dpi = 300,
  limitsize = FALSE
)

# 6. GENE ONTOLOGY ANALYSIS

# Load DEG data
aralight_data <- read.csv(file.path(base_dir, "aralightdegs_all_degs.csv"), stringsAsFactors = FALSE)

# Calculate median log2FC
genes_median_fc <- aralight_data %>%
  mutate(
    upregulated_fc = map(upregulated_fc, ~ {
      if (is.na(.x)) return(NA_real_)
      as.numeric(str_split(.x, ";")[[1]])
    }),
    downregulated_fc = map(downregulated_fc, ~ {
      if (is.na(.x) || .x == "") return(numeric(0))
      as.numeric(str_split(.x, ";")[[1]])
    }),
    median_fc = map2_dbl(upregulated_fc, downregulated_fc, ~ {
      all_values <- c(.x, .y)
      if (length(all_values) == 0) return(NA_real_)
      median(all_values, na.rm = TRUE)
    })
  ) %>%
  dplyr::select(GeneID = gene_id, MEDIAN_FC = median_fc, FREQUENCY = total_numbers)

# Filter frequent genes (50% of all conditions)
genes_num_threshold <- sum(metadata_counts_intensity_time$count)/2
frequent_genes <- genes_median_fc %>%
  filter(FREQUENCY >= genes_num_threshold)

# Filter by median fold change
log2FC_treshold = 0.5
upregulated_genes <- frequent_genes[frequent_genes$MEDIAN_FC >= log2FC_treshold, "GeneID"]
downregulated_genes <- frequent_genes[frequent_genes$MEDIAN_FC <= -log2FC_treshold, "GeneID"]
mixed_genes <- frequent_genes[frequent_genes$MEDIAN_FC > -log2FC_treshold & frequent_genes$MEDIAN_FC < log2FC_treshold, "GeneID"]

gene_groups <- list(
  Upregulated = upregulated_genes,
  Downregulated = downregulated_genes,
  Mixed = mixed_genes
)

# GO enrichment analysis
go_results_list <- list()

run_go_analysis <- function(genes, class_name) {
  cat("Running GO analysis for group:", class_name, "(", length(genes), "genes)\n")
  if (length(genes) < 3) {
    message("Not enough genes for ", class_name, " (min 3 required)")
    return(NULL)
  }

  ego <- enrichGO(
    gene = genes,
    OrgDb = org.At.tair.db,
    keyType = "TAIR",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )

  go_results_list[[class_name]] <<- ego

  if (length(ego) == 0) {
    message("No enriched GO terms found for ", class_name)
    return(NULL)
  }

  go_results <- as.data.frame(ego)
  go_results$class <- class_name

  gene_lists <- setNames(
    lapply(ego@geneSets[as.character(go_results$ID)], function(term_genes) {
      intersect(term_genes, genes)
    }),
    go_results$ID
  )

  go_results$geneIDs <- sapply(go_results$ID, function(id) {
    paste(gene_lists[[id]], collapse = ", ")
  })

  return(go_results)
}

# Run GO analysis for all groups
all_go_results <- lapply(names(gene_groups), function(class_name) {
  genes <- gene_groups[[class_name]]
  run_go_analysis(genes, class_name)
}) %>%
  bind_rows() %>%
  dplyr::select(class, ID, Description, GeneRatio, BgRatio, pvalue, qvalue, geneIDs)

# Filter and save results
filtered_go <- all_go_results %>%
  filter(qvalue < 1e-5) %>%
  mutate(
    ratio_split = str_split(GeneRatio, "/"),
    ratio_percent = as.numeric(sapply(ratio_split, `[`, 1)) / 
                    as.numeric(sapply(ratio_split, `[`, 2)) * 100
  ) %>%
  arrange(class, qvalue)

filtered_go$ratio_split <- NULL
write.table(
  filtered_go,
  file = file.path(output_dir_1, "1.GO_enrichment_top_genes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Prepare data for visualization
filtered_go <- filtered_go %>%
  separate(GeneRatio, into = c("n_genes", "total"), sep = "/") %>%
  mutate(
    n_genes = as.numeric(n_genes),
    total = as.numeric(total),
    ratio = n_genes / total,
    class_label = factor(
      class,
      levels = c("Upregulated", "Mixed", "Downregulated"),
      labels = c("UP", "MIXED", "DOWN")
    )
  ) %>%
  group_by(class_label) %>%
  arrange(ratio, .by_group = TRUE) %>%
  mutate(
    Description = factor(Description, levels = unique(Description))
  ) %>%
  ungroup()

# Create dot plot
dot_plot <- ggplot(filtered_go, aes(
  x = ratio,
  y = Description,
  size = n_genes,
  color = -log10(qvalue)
)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(q-value)") +
  scale_size_continuous(
    range = c(3, 8),
    breaks = pretty_breaks(n = 4),
    name = "Gene count"
  ) +
  scale_x_continuous(labels = percent) +
  facet_grid(
    class_label ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  labs(x = "Gene Ratio (%)", y = NULL) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10),
    strip.text.y.left = element_text(
      size = 12,
      face = "bold",
      angle = 0,
      hjust = 0
    ),
    strip.background = element_blank(),
    panel.spacing.y = unit(1, "lines"),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )

# Save plot
base_height <- 5
term_height <- 0.1
total_height <- base_height + term_height * nrow(filtered_go)

ggsave(
  filename = file.path(output_dir_1, "1.Gene_Ontology.pdf"),
  plot = dot_plot,
  width = 10,
  height = total_height,
  dpi = 300,
  limitsize = FALSE
)

# Save results
save.image(file.path(output_dir_1, "1.Initial_analysis_results.Rdata"))
message("First step of analysis completed successfully. Results saved to: ", output_dir_1)

# SECOND STEP OF ANALYSIS
# GENE CLASSIFICATION & DETAILED ANALYSIS OF GENE ONTOLOGY GROUPS

# BEFORE NEXT STEPS PLEASE SET UP THE DESIRED DATASET TO ANALYZE
# by setting analysis = "LEAVES" or analysis = "SEEDLINGS"

#analysis = "LEAVES"
analysis = "SEEDLINGS"

output_dir_2 <- file.path(base_dir, paste0("2_step_results_",tolower(analysis)))
output_dir_3 <- file.path(base_dir, paste0("3_step_results_",tolower(analysis)))

# Create output directories
dir.create(output_dir_2, showWarnings = FALSE)
dir.create(output_dir_3, showWarnings = FALSE)

make_count_tissue <- function(analysis_type = "LEAVES") {
  pattern <- if (analysis_type == "LEAVES") "LEAF" else "SEEDLING"
  
  function(x) {
    if (is.na(x) || x == "") return(0)
    parts <- unlist(str_split(x, ";"))
    sum(str_detect(parts, pattern))
  }
}

count_tissue <- make_count_tissue(analysis)

genes_tissue <- aralight_data %>%
  rowwise() %>%
  mutate(
    tissue_up   = count_tissue(upregulated_experiments),
    tissue_down = count_tissue(downregulated_experiments),
    tissue_total = tissue_up + tissue_down
  ) %>%
  filter(tissue_total >= 1) %>%  
  ungroup()

gene_data_tissue <- tibble(
  gene_id = genes_tissue$gene_id,
  freq = genes_tissue$tissue_total,
  is_tf = ifelse(genes_tissue$tf_family != "" & !is.na(genes_tissue$tf_family), 
                 "Transcription Factor", 
                 "Non-TF Gene")
)

freq_counts_tissue <- gene_data_tissue %>%
  count(freq, is_tf, name = "n") %>%
  spread(key = is_tf, value = n, fill = 0) %>%
  mutate(
    `Transcription Factor` = as.integer(`Transcription Factor`),
    `Non-TF Gene` = as.integer(`Non-TF Gene`),
    Total = `Transcription Factor` + `Non-TF Gene`
  ) %>%
  arrange(freq)

plot_data_freq_tissue <- freq_counts_tissue %>%
  pivot_longer(cols = c("Non-TF Gene", "Transcription Factor"),
               names_to = "GeneType",
               values_to = "Count")

freq_tissue_plot <- ggplot(plot_data_freq_tissue, aes(x = as.factor(freq), y = Count, fill = GeneType)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +  # Thin black outline
  labs(
    title = paste0("Number of ",tolower(analysis),"-related DEGs by their occurrence"),
    x = "Number of detections",
    y = "Number of genes",
    fill = "Gene type"
  ) +
  scale_fill_manual(values = c("Non-TF Gene" = "#69b3a2", "Transcription Factor" = "#404080")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(hjust = 1, vjust = 6, colour = "black", size = 14),
    axis.text.y = element_text(colour = "black", size = 14),
    legend.text = element_text(size = 14),      # Legend text size
    legend.title = element_text(size = 14),     # Legend title size
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.5, 0.7),
    axis.title.x = element_text(margin = margin(t = -10), size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5),
  ) +
  scale_x_discrete(breaks = c("1", "5", "10", "15", "20", "25"))  # Show only needed labels

ggsave(
  filename = file.path(output_dir_2, paste0("2.Frequency_DEGs_", tolower(analysis), ".pdf")),
  plot = freq_tissue_plot,
  width = 8,
  height = 6,
  dpi = 300,
  limitsize = FALSE
)

# Select genes with FREQUENCY >= threshold
threshold_freq <- c(LEAVES = 5, SEEDLINGS = 3)[[analysis]]

gene_data_tissue_selected <- gene_data_tissue[gene_data_tissue$freq >= threshold_freq,]

exclude_tissue <- if (analysis == "LEAVES") {
  "SEEDLING"
} else if (analysis == "SEEDLINGS") {
  "LEAF"
}


samples_to_keep_tissue <- rownames(metadata)[
    !metadata$EXPERIMENT %in% c("GSE251796", "PRJNA699408") & 
        metadata$TISSUE_GROUPED != exclude_tissue
]
gene_counts_tissue <- gene_counts[, samples_to_keep_tissue]
metadata_leaves <- metadata[samples_to_keep_tissue, ]

# Normalize and select DEGs relevant to leaves
norm_data_tissue <- normalize_counts(gene_counts_tissue)
selected_genes_leaves <- gene_data_tissue_selected$gene_id

# Select genes present in both tables
selected_genes_presented <- intersect(selected_genes_leaves, rownames(norm_data_tissue))
selected_data_tissue <- norm_data_tissue[selected_genes_presented,]

filtered_genes_tissue <- genes_tissue %>%
  filter(gene_id %in% selected_genes_presented)

# Filter metadata: only LEAF and HIGH_LIGHT
metadata_tissue_high_light <- metadata_leaves %>%
  filter(
    TISSUE_GROUPED != exclude_tissue,
    CONDITION == "HIGH_LIGHT"
  )

# Remove duplicates by FULL_DATA (keep one sample per condition)
metadata_unique_conditions_tissue <- metadata_tissue_high_light %>%
  distinct(FULL_DATA, .keep_all = TRUE)

# Add groups
metadata_condition_groups_tissue <- metadata_unique_conditions_tissue %>%
  mutate(
    intensity_group = case_when(
      LIGHT_INTENSITY < 1000 ~ "Moderate",
      LIGHT_INTENSITY >= 1000 & LIGHT_INTENSITY < 1500 ~ "Intensive",
      TRUE ~ "Severe"
    ),
    time_group = case_when(
      TIME <= 90 ~ "Short",
      TIME >= 120 & TIME <= 960 ~ "Medium",
      TIME >= 1440 ~ "Long-term"
    ),
    AGE_GROUP = case_when(
      AGE >= 7 & AGE <= 15 ~ "Juvenile",
      AGE >= 21 & AGE <= 30 ~ "Mature",
      AGE >= 31 & AGE <= 49 ~ "Old"
    ),
    time_group_s = case_when(
      TIME <= 90 ~ "Short",
      TIME >= 120 ~ "Long",
    ),
  )

# Aggregate comprehensive statistic
condition_stats_tissue <- metadata_condition_groups_tissue %>%
  group_by(intensity_group, time_group, AGE_GROUP) %>%
  summarise(
    condition_count = n(),
    conditions = paste(FULL_DATA, collapse = ", "),
    .groups = 'drop'
  ) %>%
  mutate(group_key = paste(intensity_group, time_group, AGE_GROUP, sep = "_"))

# Simplified time grouping for further metaanalysis
condition_stats_tissue_time <- metadata_condition_groups_tissue %>%
  group_by(time_group_s) %>%
  summarise(
    condition_count = n(),
    conditions = paste(FULL_DATA, collapse = ", "),
    .groups = 'drop'
  ) %>%
  mutate(group_key = paste(time_group_s))

# Function to classify experiments by time
classify_experiment_time <- function(exp, analysis = "LEAVES") {
  # Determine the type of tissue to analyze
  tissue_type <- if (analysis == "LEAVES") "LEAF" else if (analysis == "SEEDLINGS") "SEEDLING" else NA
  
  if (!grepl(tissue_type, exp)) return(NA_character_)
  
  # Change LOCAL_*/SYSTEMIC_ to simple LEAF
  exp <- gsub(paste0("(LOCAL|SYSTEMIC)_", tissue_type), tissue_type, exp)
  parts <- unlist(strsplit(exp, "_"))
  
  # Structure check of data
  if (length(parts) < 6) return(NA_character_)
  
  # Time data
  time_part <- paste(parts[5:length(parts)], collapse = "_")
  
  # Convert to minutes
  time_min <- case_when(
    grepl("seconds", time_part) ~ as.numeric(sub("_seconds", "", time_part)) / 60,
    grepl("minutes", time_part) ~ as.numeric(sub("_minutes", "", time_part)),
    grepl("hour", time_part) ~ {
      hours <- sub("_hour(s)?", "", time_part)
      hours <- gsub("_", ".", hours)
      60 * as.numeric(hours)
    },
    grepl("day", time_part) ~ 1440 * as.numeric(sub("_day(s)?", "", time_part)),
    grepl("minute", time_part) ~ as.numeric(sub("_minute(s)?", "", time_part)),
    TRUE ~ NA_real_
  )
  
  if (is.na(time_min)) return(NA_character_)
  
  # Time classification
  time_group <- ifelse(time_min <= 90, "Short", "Long")
  
  return(time_group)
}

# Function to process a column of experiments with progress tracking
classify_experiments_column <- function(exp_column, col_name) {
  total <- length(exp_column)
  result <- character(total)
  
  for (i in seq_along(exp_column)) {
    # Print progress every 100 genes
    if (i %% 100 == 0) {
      message(sprintf("Processed %d of %d genes for column %s", i, total, col_name))
    }
    
    x <- exp_column[i]
    if (is.na(x) || x == "") {
      result[i] <- NA_character_
      next
    }
    
    # Split multiple experiments
    exps <- unlist(strsplit(x, ";\\s*"))
    classifications <- sapply(exps, classify_experiment_time, analysis = analysis, USE.NAMES = FALSE)
    valid_classifications <- classifications[!is.na(classifications)]
    
    if (length(valid_classifications) == 0) {
      result[i] <- NA_character_
    } else {
      # Keep all classifications including duplicates
      result[i] <- paste(valid_classifications, collapse = "; ")
    }
  }
  
  message(sprintf("Completed processing column %s (%d genes)", col_name, total))
  return(result)
}

# Main processing with progress messages
message("Starting processing of upregulated_experiments...")
filtered_genes_tissue$classification_upregulated_exp <- classify_experiments_column(
  filtered_genes_tissue$upregulated_experiments, 
  "upregulated_experiments"
)

message("\nStarting processing of downregulated_experiments...")
filtered_genes_tissue$classification_downregulated_exp <- classify_experiments_column(
  filtered_genes_tissue$downregulated_experiments, 
  "downregulated_experiments"
)

message("Processing completed successfully!")

# Create all column names
up_columns <- paste0(condition_stats_tissue_time$group_key, "_Up")
down_columns <- paste0(condition_stats_tissue_time$group_key, "_Down")
all_new_columns <- c(up_columns, down_columns)

# Initialize new columns with 0
filtered_genes_tissue[all_new_columns] <- 0

# Modified processing function with progress tracking
process_condition_counts <- function(df) {
  total_genes <- nrow(df)
  message("Starting processing of ", total_genes, " genes...")
  
  for (i in 1:total_genes) {
    if (i %% 100 == 0) {
      message("Processed ", i, " of ", total_genes, " genes (", 
              round(i/total_genes*100, 1), "%)")
    }
    
    # Process upregulated
    up_class <- df$classification_upregulated_exp[i]
    if (!is.na(up_class) && nchar(up_class) > 0) {
      up_terms <- unlist(strsplit(up_class, ";\\s*"))
      up_counts <- table(up_terms)
      
      for (cond in names(up_counts)) {
        col_name <- paste0(cond, "_Up")
        if (col_name %in% all_new_columns) {
          df[i, col_name] <- up_counts[cond]
        }
      }
    }
    
    # Process downregulated
    down_class <- df$classification_downregulated_exp[i]
    if (!is.na(down_class) && nchar(down_class) > 0) {
      down_terms <- unlist(strsplit(down_class, ";\\s*"))
      down_counts <- table(down_terms)
      
      for (cond in names(down_counts)) {
        col_name <- paste0(cond, "_Down")
        if (col_name %in% all_new_columns) {
          df[i, col_name] <- down_counts[cond]
        }
      }
    }
  }
  
  message("Completed processing all ", total_genes, " genes")
  return(df)
}

# Apply the processing
filtered_genes_tissue <- process_condition_counts(filtered_genes_tissue)

# Verify the results
message("Successfully added new columns:")
print(head(filtered_genes_tissue[, all_new_columns]))

# Function to identify top 25% genes in each column
get_top_percent <- function(df, col_name, percent) {
  # Calculate 75% quantile to determine top 25% cutoff
  cutoff <- quantile(df[[col_name]], probs = 1 - percent, na.rm = TRUE)
  
  # Create new column name
  new_col_name <- paste0(col_name, "_top25")
  
  # Add column marking top 25%
  df[[new_col_name]] <- ifelse(df[[col_name]] >= cutoff, 
                              col_name, 
                              NA_character_)
  return(df)
}

# Apply function to all 4 columns
filtered_genes_tissue <- filtered_genes_tissue %>%
  get_top_percent("Long_Up", 0.25) %>%
  get_top_percent("Short_Up", 0.25) %>%
  get_top_percent("Long_Down", 0.25) %>%
  get_top_percent("Short_Down", 0.25)

# Create Venn diagram of gene intersections
top_genes_lists <- list(
    "Short Up" = filtered_genes_tissue$gene_id[!is.na(filtered_genes_tissue$Short_Up_top25)],
    "Long Up" = filtered_genes_tissue$gene_id[!is.na(filtered_genes_tissue$`Long_Up_top25`)],
    "Long Down" = filtered_genes_tissue$gene_id[!is.na(filtered_genes_tissue$`Long_Down_top25`)],
    "Short Down" = filtered_genes_tissue$gene_id[!is.na(filtered_genes_tissue$Short_Down_top25)]
)

gene_intersection_plot <- ggvenn(
    top_genes_lists,
    fill_color = c("#f0a7a7", "#a62b2b", "#1f4b85", "#8abdd2"),
    stroke_size = 0.5,
    set_name_size = 7,
    text_size = 5.5
)

print(gene_intersection_plot)

ggsave(
  filename = file.path(output_dir_2, paste0("2.Venn_diagram_specific_genes_", tolower(analysis), ".pdf")),
  plot = gene_intersection_plot,
  width = 12,
  height = 10,
  dpi = 300,
  limitsize = FALSE
)

# Calculate median fold change
filtered_genes_tissue <- filtered_genes_tissue %>%
  mutate(
    # Combine all FC values
    all_fc_values = map2(upregulated_fc, downregulated_fc, ~ {
      up_vals <- as.numeric(unlist(strsplit(.x, ";")))
      down_vals <- as.numeric(unlist(strsplit(.y, ";")))
      c(up_vals, down_vals)
    }),
    
    # Calculate median FC
    median_fc = map_dbl(all_fc_values, median, na.rm = TRUE),
    
    # Remove temporary column
    all_fc_values = NULL
  )

# Add direction column based on median_fc
filtered_genes_tissue <- filtered_genes_tissue %>%
  mutate(
    direction = case_when(
      median_fc >= 0.5 ~ "Up",
      median_fc <= -0.5 ~ "Down",
      TRUE ~ "Mixed"
    )
  )

# Get all possible intersections
get_all_intersections <- function(gene_lists) {
  all_genes <- unique(unlist(gene_lists))
  membership <- map_dfc(names(gene_lists), ~ tibble(!!.x := all_genes %in% gene_lists[[.x]]))
  
  # Create results list
  results <- tibble(GeneID = all_genes, classification = NA_character_)
  
  # Function to classify genes
  classify_gene <- function(gene_row) {
    groups <- names(gene_row)[gene_row]
    n_groups <- length(groups)
    
    if (n_groups == 0) {
      return(NA_character_)
    } else if (n_groups == 1) {
      return(groups[1])
    } else if (n_groups == 2) {
      return(paste(sort(groups), collapse = " & "))
    } else {
      return("Frequent")
    }
  }
  
  # Apply classification to each gene
  results$classification <- apply(membership, 1, classify_gene)
  
  # Remove NAs (genes not in any group)
  results <- results %>% filter(!is.na(classification))
  
  return(results)
}

# Get gene classification
gene_classification <- get_all_intersections(top_genes_lists)

# Recode classification names
gene_classification_fin <- gene_classification %>%
  mutate(classification = recode(
    classification,
    "Long Up & Short Up" = "Short–Long Up",
    "Long Down & Short Down" = "Short–Long Down",
    "Short Down & Short Up"  = "Short Up–Down",
    "Long Up & Short Down" = "Short Down–Long Up",
    "Long Down & Short Up" = "Short Up–Long Down",
    "Long Down & Long Up" = "Long Up–Down",
  ))

# Save classification
write.table(
  gene_classification_fin,
  file = file.path(output_dir_2, paste0("2.DEGs_", tolower(analysis),"_classification.tsv")),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Merge classification with main data
filtered_genes_tissue <- filtered_genes_tissue %>%
  left_join(
    gene_classification_fin,
    by = c("gene_id" = "GeneID")
  ) %>%
  # Replace NA with "not specific"
  mutate(classification = ifelse(is.na(classification), "Not specific", classification))

# Get gene symbols from Ensembl
mart <- biomaRt::useMart(
    host = "https://plants.ensembl.org",   
    biomart = "plants_mart",
    dataset = "athaliana_eg_gene"
)

# Get gene symbols for filtered_genes_tissue
gene_symbols <- tryCatch({
    biomaRt::getBM(
        attributes = c("ensembl_gene_id", "external_gene_name"),
        filters = "ensembl_gene_id",
        values = filtered_genes_tissue$gene_id,
        mart = mart
    )
}, error = function(e) {
    warning("Failed to retrieve gene symbols from BioMart")
    tibble(ensembl_gene_id = filtered_genes_tissue$gene_id, 
           external_gene_name = filtered_genes_tissue$gene_id)
})

# Rename columns for clarity
colnames(gene_symbols) <- c("gene_id", "symbol_id")

# Merge with main data
filtered_genes_tissue <- filtered_genes_tissue %>%
  left_join(gene_symbols, by = "gene_id") %>%
  # Move symbol_id after gene_id
  relocate(symbol_id, .after = gene_id)

# Fill missing symbols with gene_id
filtered_genes_tissue <- filtered_genes_tissue %>%
  mutate(symbol_id = ifelse(is.na(symbol_id) | symbol_id == "", gene_id, symbol_id))

# Save full data
write.table(
  filtered_genes_tissue,
  file = file.path(output_dir_2, paste0("2.DEGs_", tolower(analysis),"_full_data.tsv")),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Run GO analysis for all groups
all_go_results_tissue <- lapply(names(top_genes_lists), function(class_name) {
  genes <- top_genes_lists[[class_name]]
  run_go_analysis(genes, class_name)
}) %>%
  bind_rows() %>%
  dplyr::select(class, ID, Description, GeneRatio, BgRatio, pvalue, qvalue, geneIDs)

# Filter and save results
filtered_go_tissue <- all_go_results_tissue %>%
  filter(qvalue < 1e-5) %>%
  mutate(
    ratio_split = str_split(GeneRatio, "/"),
    ratio_percent = as.numeric(sapply(ratio_split, `[`, 1)) / 
                    as.numeric(sapply(ratio_split, `[`, 2)) * 100
  ) %>%
  filter(ratio_percent >= 3) %>% 
  arrange(class, qvalue)

# Sankey visualization of top GO processes
edges_tissue <- filtered_go_tissue %>%
  dplyr::select(Description, Group = class, geneIDs) %>%
  dplyr::mutate(
    gene_count = lengths(strsplit(geneIDs, ",")),
    Group_split = stringr::str_split(Group, " ∩ "),
    Group_split = stringr::str_remove(Group_split, "-specific"),
    Group_split = lapply(Group_split, function(x) {
      ifelse(length(x) == 0, x, trimws(x))
    })
  ) %>%
  tidyr::unnest(Group_split) %>%
  dplyr::mutate(
    source = Group_split,
    target = Description
  ) %>%
  dplyr::group_by(source, target) %>%
  dplyr::summarise(
    value = sum(gene_count),
    .groups = "drop"
  )

edges_tissue <- edges_tissue %>%
  mutate(
    source = as.character(source),
    source_group = gsub(" ", "", source)
  )

# Set up nodes
nodes_tissue <- data.frame(
  name = unique(c(edges_tissue$source, edges_tissue$target)),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    group = ifelse(name %in% unique(filtered_go_tissue$Description), "Process", "Category"),
    level = ifelse(group == "Process", 1, 0)
  )

# Add indexes for network
nodes_tissue$id <- 0:(nrow(nodes_tissue) - 1)

edges_tissue <- edges_tissue %>%
  dplyr::left_join(nodes_tissue %>% dplyr::select(name, id), by = c("source" = "name")) %>%
  dplyr::rename(source_idx = id) %>%
  dplyr::left_join(nodes_tissue %>% dplyr::select(name, id), by = c("target" = "name")) %>%
  dplyr::rename(target_idx = id)

# Color scale
colourScale_l <- networkD3::JS(
  paste0(
    'd3.scaleOrdinal()
     .domain(["ShortUp", "LongUp", "ShortDown", "LongDown"])
     .range(["#f0a7a7", "#a62b2b", "#8abdd2", "#1f4b85"])
     .unknown("#888")'
  )
)

# Create Sankey diagram
create_sankey_go <- function(edges, nodes, file_out, colourScale, LinkGroup = "source_group",
                             fontSize = 18, nodeWidth = 15, nodePadding = 12,
                             margin = list(top = 50, right = 100, bottom = 50, left = 400),
                             width = 1000, height = 1200) {

  sankey <- networkD3::sankeyNetwork(
    Links = edges,
    Nodes = nodes,
    Source = "source_idx",
    Target = "target_idx",
    Value = "value",
    NodeID = "name",
    NodeGroup = "group",
    LinkGroup = LinkGroup,
    colourScale = colourScale,
    fontSize = fontSize,
    nodeWidth = nodeWidth,
    nodePadding = nodePadding,
    margin = margin,
    sinksRight = TRUE,
    iterations = 64,
    width = width,
    height = height
  )

  sankey <- htmlwidgets::onRender(sankey, '
    function(el) {
      d3.select(el).selectAll(".node text")
        .attr("text-anchor", function(d) {
          return d.group === "Process" ? "start" : "end";
        })
        .attr("x", function(d) {
          return d.group === "Process" ? 25 : -10;
        })
        .style("pointer-events", "none");

      function updateLinks() {
        d3.selectAll(".link")
          .attr("d", function(d) {
            const sourceX = d.source.x + (d.source.group === "Category" ? 20 : 0);
            const targetX = d.target.x - (d.target.group === "Process" ? 10 : 0);
            const sourceY = d.source.y;
            const targetY = d.target.y;

            return `
              M${sourceX},${sourceY}
              C${(sourceX + targetX)/2},${sourceY}
               ${(sourceX + targetX)/2},${targetY}
               ${targetX},${targetY}`;
          })
          .lower();
      }

      function adjustCategories() {
        const categories = ["Age-specific", "Multi-factor", "Time-specific", "Intensity-specific"];
        const spacing = el.clientHeight / (categories.length + 1);

        categories.forEach((cat, i) => {
          const node = d3.select(`.node[id*="${cat}"]`).node();
          if (node) {
            const d = node.__data__;
            d.y = spacing * (i + 1);
            d3.select(node).attr("transform", `translate(${d.x},${d.y})`);
          }
        });
      }

      setTimeout(function() {
        adjustCategories();
        updateLinks();
      }, 100);
    }
  ')

  htmlwidgets::saveWidget(widget = sankey, file = file_out, selfcontained = TRUE)
  return(sankey)
}

sankey_tissue <- create_sankey_go(
  edges = edges_tissue,
  nodes = nodes_tissue,
  file_out = file.path(output_dir_2, paste0("2.Sankey_GO_", tolower(analysis),".html")),
  colourScale = colourScale_l,
  LinkGroup = "source_group",
  fontSize = 18,
  nodeWidth = 15,
  nodePadding = 12,
  margin = list(top = 50, right = 100, bottom = 50, left = 400),
  width = 1000,
  height = 1200
)

# Load GO data for parsing
go_data <- readr::read_tsv(
  file.path(base_dir, "ATH_GO_GOSLIM.txt"),
  skip = 4,
  col_names = c("GeneID", "locus_tag", "gene_symbol", "relationship", 
               "description", "GO_ID", "taxon_id", "aspect", "term", 
               "evidence_code", "with_from", "interpro_ids", "analysis_ref", 
               "assigned_by", "date")
)


# Function to plot expression heatmap
plot_expression_heatmap <- function(data_list, title, percentage_threshold = 50, output_file, width = 12, height = 8) {
  
  # Normalize data
  data_norm <- data_list %>%
    mutate(
      Long_Up_pct = (Long_Up / filter(condition_stats_tissue_time, time_group_s == "Long")$condition_count * 100),
      Short_Up_pct = (Short_Up / filter(condition_stats_tissue_time, time_group_s == "Short")$condition_count * 100),
      Long_Down_pct = - (Long_Down / filter(condition_stats_tissue_time, time_group_s == "Long")$condition_count * 100),
      Short_Down_pct = - (Short_Down / filter(condition_stats_tissue_time, time_group_s == "Short")$condition_count * 100)
    ) %>%
    dplyr::select(symbol_id, ends_with("_pct")) %>% 
    mutate(across(ends_with("_pct"), ~ replace_na(., 0))) %>%
    filter(if_any(ends_with("_pct"), ~ abs(.) >= percentage_threshold)) %>%
    arrange(desc(rowSums(abs(dplyr::select(., ends_with("_pct"))))))
  
  if (nrow(data_norm) == 0) {
    warning("No data to display after filtering.")
    return(NULL)
  }
  
  # Prepare matrix
  heatmap_matrix <- data_norm %>%
    as.data.frame() %>%
    column_to_rownames("symbol_id") %>%
    as.matrix()
  
  # Change column order
  desired_order <- c("Short_Up_pct", "Long_Up_pct", "Short_Down_pct", "Long_Down_pct")
  heatmap_matrix <- heatmap_matrix[, desired_order]
  
  # Column annotations
  column_annotation <- data.frame(
    Regulation = factor(rep(c("Up", "Down"), each = 2)),
    Time = rep(c("Short", "Long"), times = 2), 
    row.names = colnames(heatmap_matrix)
  )
  
  # Color schemes
  annotation_colors <- list(
    Regulation = c(Up = "#d73027", Down = "#4575b4"),
    Time = c(Short = "#c7e9c0", Long = "#238b45")
  )
  
  # Heatmap palette
  my_palette <- colorRampPalette(c("#4575b4", "white", "#d73027"))(100)
  breaks <- seq(-100, 100, length.out = 101)
  
  # Create heatmap
  pdf(output_file, width = width, height = height)
  
  pheatmap(
    heatmap_matrix,
    color = my_palette,
    breaks = breaks,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_col = column_annotation,
    annotation_row = NULL,
    annotation_colors = annotation_colors,
    border_color = "black",
    fontsize_row = 8,
    main = title,
    cellwidth = 30,
    cellheight = 12,
    gaps_row = NULL,
    annotation_legend = TRUE
  )
  
  dev.off()
  
  return(heatmap_matrix)
}

# Analyze antioxidant genes
antioxidant_genes <- filtered_genes_tissue %>%
  filter(
    str_detect(gene_description, regex("glutathione", ignore_case = TRUE)) |
    str_detect(gene_description, regex("ascorbate", ignore_case = TRUE)) |
    str_detect(gene_description, regex("catalase", ignore_case = TRUE)) |
    str_detect(gene_description, regex("superoxide", ignore_case = TRUE)) 
  )

aox_htmp <- plot_expression_heatmap(
  data_list = antioxidant_genes,
  title = paste("Antioxidant genes", if (analysis == "LEAVES") "(leaves)" else "(seedlings)"),
  percentage_threshold =  if (analysis == "SEEDLINGS") 50 else 50,
  output_file = file.path(output_dir_2, paste0("2.Antioxidant_heatmap_", tolower(analysis),".pdf")),
  width = 6,
  height = if (analysis == "SEEDLINGS") 8 else 5
) 

# Analyze transcription factors
tfs_tissue <- filtered_genes_tissue %>%
  filter(!is.na(tf_family) & tf_family != "")

tfs_htmp <- plot_expression_heatmap(
  data_list = tfs_tissue,
  title = paste("Transcription factors", if (analysis == "LEAVES") "(leaves)" else "(seedlings)"),
  percentage_threshold =  if (analysis == "SEEDLINGS") 100 else 80,
  output_file = file.path(output_dir_2, paste0("2.TFs_heatmap_", tolower(analysis),".pdf")),
  width = 6,
  height = if (analysis == "SEEDLINGS") 8 else 8
) 

# Analyze photoreceptor genes
photoreceptor_genes <- c(
  # Main photoreceptors
  "PHYA", "PHYB", "PHYC", "PHYD", "PHYE",
  "CRY1", "CRY2", "CRY3", "CRY-DASH",
  "PHOT1", "PHOT2", "UVR8",
  "ZTL", "FKF1", "LKP2", 
  
  # Signaling components
  "PIF1", "PIF3", "PIF4", "PIF5", "PIF7",
  "SPA1", "SPA2", "SPA3", "SPA4",
  "COP1", "COP9", "DET1",
  "HY5", "HYH",
  "FHY1", "FHY3", "FAR1",
  "GUN1", "GLK1", "GUN5", "GLK2", 
  
  # Phototropism-related genes
  "NPH3", "RPT2", "NPL1",
  "PKS1", "PKS2", "PKS4",
  "JAC1", "PMI1", "THR1",
  
  # Auxin signaling
  "ARF1", "ARF2", "ARF6", "ARF8", "ARF9", "ARF10", 
  "ARF11", "ARF16", "ARF17", "ARF18", "ARF19", 
  "PID", "WAG1", "WAG2",
  
  # Circadian clocks
  "ELF3", "ELF4", "LUX",
  "PRR5", "PRR7", "PRR9",
  "GI", "TOC1", "CCA1", "LHY"
)

photoreceptor_presented <- photoreceptor_genes[photoreceptor_genes %in% filtered_genes_tissue$symbol_id]

photoreceptor_genes_data <- filtered_genes_tissue %>%
  filter(symbol_id %in% photoreceptor_presented)

phot_htmp <- plot_expression_heatmap(
  data_list = photoreceptor_genes_data,
  title = paste("Photoreceptor-related", if (analysis == "LEAVES") "(leaves)" else "(seedlings)"),
  percentage_threshold = if (analysis == "SEEDLINGS") 50 else 50,
  output_file = file.path(output_dir_2, paste0("2.Photoreceptors_heatmap_", tolower(analysis),".pdf")),
  width = 6,
  height = if (analysis == "SEEDLINGS") 6 else 5
) 

# Analyze B-box genes
bbx_genes <- filtered_genes_tissue %>%
  filter(grepl("BBX", symbol_id, ignore.case = FALSE))

bbx_htmp <- plot_expression_heatmap(
  data_list = bbx_genes,
  title = paste("B-box genes", if (analysis == "LEAVES") "(leaves)" else "(seedlings)"),
  percentage_threshold = if (analysis == "SEEDLINGS") 50 else 50,
  output_file = file.path(output_dir_2, paste0("2.BBX_heatmap_", tolower(analysis),".pdf")),
  width = 6,
  height = if (analysis == "SEEDLINGS") 3 else 3
)

# Analyze photomorphogenesis genes
photomorpho_genes <- go_data %>%
  filter(
    grepl("photomorphogenesis", term, ignore.case = TRUE) |
      grepl("photomorphogenesis", description, ignore.case = TRUE)
  ) %>%
  distinct(GeneID, .keep_all = TRUE)  

photomorpho_genes_l <- photomorpho_genes %>%
  pull(GeneID) %>% 
  str_split(",\\s*") %>%            
  unlist() %>%    
  unique()

photomorpho_genes_data <- filtered_genes_tissue %>%
  filter(gene_id %in% photomorpho_genes_l)

photomorphogenesis_htmp <- plot_expression_heatmap(
  data_list = photomorpho_genes_data,
  title = paste("Photomorphogenesis", if (analysis == "LEAVES") "(leaves)" else "(seedlings)"),
  percentage_threshold = if (analysis == "SEEDLINGS") 50 else 50,
  output_file = file.path(output_dir_2, paste0("2.Photomorphogenesis_heatmap_", tolower(analysis),".pdf")),
  width = 6,
  height = if (analysis == "SEEDLINGS") 6 else 5
)

# Analyze photosynthesis-associated nuclear genes (PhANGs)
phang_genes <- c(
  "LHCA1", "LHCA2", "LHCA3", "LHCA4", "Lhca6", "LHCB2.1",
  "LHCB2.2", "LHCB2.3", "LHCB3", "LHCB4.1", "LHCB4.2",
  "LHCB4.3", "LHCB5", "LHCB6",
  "PSB28", "psbA", "psbB", "psbC", "psbD", "PSBO1", "PSBO2", 
  "PSBP-1", "PSBP-2", "PSBQ-2", "PSBQA", "PSBR", "PSBTN", 
  "PSBW", "PSBX", "PSBY",
  "RBCS1A","RBCS1B", "RBCS2B", "RBCS3B", "PGK", "GAPDH", "TPI", "PRK",
  "FBA1", "FBA4", "FBA5", "FBA6", "FBA7", "FBA8"
)

phang_presented <- phang_genes[phang_genes %in% filtered_genes_tissue$symbol_id]

phang_genes_data <- filtered_genes_tissue %>%
  filter(symbol_id %in% phang_presented)

phang_htmp <- plot_expression_heatmap(
  data_list = phang_genes_data,
  title = paste("PhANGs", if (analysis == "LEAVES") "(leaves)" else "(seedlings)"),
  percentage_threshold = if (analysis == "SEEDLINGS") 50 else 50,
  output_file = file.path(output_dir_2, paste0("2.PhANGs_heatmap_", tolower(analysis),".pdf")),
  width = 6,
  height = if (analysis == "SEEDLINGS") 7 else 4.5
) 

# Analyze chloroplast organization genes
chlor_org <- all_go_results_tissue %>%
  filter(Description == "chloroplast organization")

chlor_org_l <- chlor_org %>%
  pull(geneIDs) %>% 
  str_split(",\\s*") %>%            
  unlist() %>%    
  unique()

chlor_org_data <- filtered_genes_tissue %>%
  filter(gene_id %in% chlor_org_l)

chlor_org_htmp <- plot_expression_heatmap(
  data_list = chlor_org_data,
  title = paste("Chloroplast organization", if (analysis == "LEAVES") "(leaves)" else "(seedlings)"),
  percentage_threshold = if (analysis == "SEEDLINGS") 50 else 50,
  output_file = file.path(output_dir_2, paste0("2.Chloroplast_organization_heatmap_", tolower(analysis),".pdf")),
  width = 6,
  height = if (analysis == "SEEDLINGS") 9 else 9
)

# Analyze redoxins
redoxins <- filtered_genes_tissue %>%
  filter(str_detect(gene_description, regex("redoxin", ignore_case = TRUE)))

redoxins_htmp <- plot_expression_heatmap(
  data_list = redoxins,
  title = paste("Redoxins", if (analysis == "LEAVES") "(leaves)" else "(seedlings)"),
  percentage_threshold = if (analysis == "SEEDLINGS") 50 else 50,
  output_file = file.path(output_dir_2, paste0("2.Redoxins_heatmap_", tolower(analysis),".pdf")),
  width = 6,
  height = if (analysis == "SEEDLINGS") 7 else 6
)

save.image(file.path(output_dir_2, paste0("2.Genes_analysis_results_",tolower(analysis),".Rdata")))
# THIRD STEP OF ANALYSIS
# COEXPRESSION & GENE REGULATORY NETWORK RECONSTRUCTION

# Set global parameters
min_module_size <- 30   # For WGCNA

if (analysis == "LEAVES") {
  min_correlation <- 0.3
} else if (analysis == "SEEDLINGS") {
  min_correlation <- 0.55
}


exclude_tissue <- if (analysis == "LEAVES") {
  "SEEDLING"
} else if (analysis == "SEEDLINGS") {
  "LEAF"
}

# Select samples for short-term and long-term high light response
short_samples_tissue <- rownames(metadata)[
    !metadata$EXPERIMENT %in% c("GSE251796", "PRJNA699408") & 
        metadata$TISSUE_GROUPED != exclude_tissue &
        metadata$CONDITION == "HIGH_LIGHT" &
        metadata$TIME <= 90]

long_samples_tissue <- rownames(metadata)[
    !metadata$EXPERIMENT %in% c("GSE251796", "PRJNA699408") & 
        metadata$TISSUE_GROUPED != exclude_tissue &
        metadata$CONDITION == "HIGH_LIGHT" &
        metadata$TIME >= 120]

n_leaves_long <- length(long_samples_tissue)
n_leaves_short <- length(short_samples_tissue)

# Normalize counts
#norm_data_leaves <- normalize_counts(gene_counts)

gene_counts_short <- norm_data_tissue[, short_samples_tissue]
gene_counts_long <- norm_data_tissue[, long_samples_tissue]

# Filter genes for short-term response
selected_genes_short <- filtered_genes_tissue[
  !(filtered_genes_tissue$classification %in% c("Long Down", "Long Up", "Long Up–Down", "Not specific")), 
  ]

selected_genes_short_id <- selected_genes_short$gene_id

# Filter genes for long-term response
selected_genes_long <- filtered_genes_tissue[
  !(filtered_genes_tissue$classification %in% c("Short Down", "Short Up", "Short Up–Down", "Not specific")), 
  ]

selected_genes_long_id <- selected_genes_long$gene_id

# Get valid genes present in both tables
valid_genes_short <- intersect(selected_genes_short_id, rownames(gene_counts_short))
valid_genes_long <- intersect(selected_genes_long_id, rownames(gene_counts_long))

data_short_tissue <- gene_counts_short[valid_genes_short, ]
data_long_tissue <- gene_counts_long[valid_genes_long, ]

# Update direction classification for short-term genes
selected_genes_short <- selected_genes_short %>%
  mutate(direction = case_when(
    classification == "Short Up–Down" ~ "Mixed",
    classification == "Short Up–Long Down" ~ "Up",
    classification == "Short Down" ~ "Down",
    classification == "Short Down–Long Up" ~ "Down",
    classification == "Short–Long Up" ~ "Up",
    classification == "Short Up" ~ "Up",
    classification == "Short–Long Down" ~ "Down",
    TRUE ~ direction
  ))

# Update direction classification for long-term genes
selected_genes_long <- selected_genes_long %>%
  mutate(direction = case_when(
    classification == "Short Up–Long Down" ~ "Down",
    classification == "Short Down–Long Up" ~ "Up",
    classification == "Short–Long Up" ~ "Up",
    classification == "Short–Long Down" ~ "Down",
    classification == "Long Up–Down" ~ "Mixed",
    classification == "Long Up" ~ "Up",
    classification == "Long Down" ~ "Down",
    TRUE ~ direction
  ))

# Prepare data for WGCNA
datExpr_tissue_short <- t(data_short_tissue)
datExpr_tissue_long <- t(data_long_tissue) 

run_wgcna <- function(datExpr, output_dir_3, prefix, minModuleSize = 30) {
  
  # Check data quality
  gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 3)
  if (!gsg$allOK) {
    if (sum(!gsg$goodGenes) > 0) {
      message(paste("Removing genes:", 
                    paste(colnames(datExpr)[!gsg$goodGenes], collapse = ", ")))
    }
    if (sum(!gsg$goodSamples) > 0) {
      message(paste("Removing samples:", 
                    paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")))
    }
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  }
  
  # Stop if not enough genes remain
  num_genes <- ncol(datExpr)
  if (num_genes < minModuleSize) {
    stop("Too few genes remaining after filtering (", num_genes, ") to form modules (minModuleSize = ", minModuleSize, "). Reduce minModuleSize or check data quality.")
  }

  # Determine soft-thresholding power
  powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
  sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

  # Plot scale independence and mean connectivity
  grDevices::pdf(file.path(output_dir_3, paste0(prefix,"_",tolower(analysis),"_soft_threshold_plots.pdf")), width = 10, height = 6)
  par(mfrow = c(1, 2))
  cex1 <- 0.9

  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab = "Soft Threshold (power)", 
       ylab = "Scale Free Topology Model Fit, signed R^2",
       type = "n", main = "Scale independence")
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels = powers, cex = cex1, col = "red")
  abline(h = 0.90, col = "red")

  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab = "Soft Threshold (power)", 
       ylab = "Mean Connectivity", 
       type = "n", main = "Mean connectivity")
  text(sft$fitIndices[,1], sft$fitIndices[,5], 
       labels = powers, cex = cex1, col = "red")
  grDevices::dev.off()

  # Get power estimate
  power <- ifelse(is.na(sft$powerEstimate), 6, sft$powerEstimate)
  message("Using soft-thresholding power: ", power)

  # Adjust minModuleSize if needed
  adjusted_minModuleSize <- min(minModuleSize, floor(num_genes / 10))
  message("Adjusting minModuleSize from ", minModuleSize, " to ", adjusted_minModuleSize, " based on number of genes.")

  # Construct network
  net <- WGCNA::blockwiseModules(
    datExpr,
    power = power,
    TOMType = "unsigned",
    minModuleSize = adjusted_minModuleSize,
    reassignThreshold = 0,
    mergeCutHeight = 0.25,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs = TRUE,
    saveTOMFileBase = file.path(output_dir_3, paste0(prefix, "_",tolower(analysis), ".geneTOM")),
    verbose = 3
  )

  # Save module colors
  moduleColors <- WGCNA::labels2colors(net$colors)

  # Plot dendrogram with modules
  grDevices::pdf(file.path(output_dir_3, paste0(prefix, "_",tolower(analysis), "_module_dendrogram.pdf")), , width = 10, height = 6)
  WGCNA::plotDendroAndColors(
    net$dendrograms[[1]],
    moduleColors[net$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05
  )
  grDevices::dev.off()

  return(list(
    net = net,
    moduleColors = moduleColors,
    power = power
  ))
}

# Perform WGCNA analysis
wgcna_results_short <- run_wgcna(datExpr_tissue_short, output_dir_3, prefix="3.Short", minModuleSize = min_module_size)
wgcna_results_long <- run_wgcna(datExpr_tissue_long, output_dir_3, prefix="3.Long", minModuleSize = min_module_size)

extract_tom_network <- function(wgcna_results, datExpr, min_correlation = 0.3, output_path = NULL) {
  # Parsing WGCNA metadata
  moduleLabels <- wgcna_results$net$colors
  moduleColors <- wgcna_results$moduleColors
  MEs <- wgcna_results$net$MEs
  power <- wgcna_results$power
  
  # Create gene-module df
  gene_module_df <- data.frame(
    GeneID = colnames(datExpr),
    Module = moduleColors
  )
  
  # Build TOM
  TOM <- WGCNA::TOMsimilarityFromExpr(datExpr, power = power)
  rownames(TOM) <- colnames(TOM) <- colnames(datExpr)
  
  # Convert TOM to edgelist
  tom_df <- tom_to_edgelist(TOM, min_correlation = min_correlation)
  
  # Save to path
  if (!is.null(output_path)) {
    readr::write_tsv(tom_df, output_path)
  }
  
  # Return data
  list(
    tom_df = tom_df,
    gene_module_df = gene_module_df,
    power = power,
    MEs = MEs
  )
}

tom_to_edgelist <- function(TOM, min_correlation = min_correlation) {
  edgelist <- reshape2::melt(TOM) %>%
    dplyr::rename(Gene1 = Var1, Gene2 = Var2, Correlation = value) %>%
    dplyr::mutate(
      Gene1 = as.character(Gene1),
      Gene2 = as.character(Gene2)
    ) %>%
    dplyr::filter(Gene1 < Gene2 & Correlation >= min_correlation)
  return(edgelist)
}


tom_network_short <- extract_tom_network(wgcna_results_short, datExpr_tissue_short,
  min_correlation = min_correlation, output_path = file.path(output_dir_3, paste0("3.Short_",tolower(analysis),"_",min_correlation,"_TOM_df.tsv")))

tom_network_long <- extract_tom_network(wgcna_results_long, datExpr_tissue_long,
  min_correlation = min_correlation, output_path = file.path(output_dir_3, paste0("3.Long_",tolower(analysis),"_",min_correlation,"_TOM_df.tsv")))

# For calculation of GENIE3 edges you can run the following code below, 
# But DIANE step requiring at least 128 Gb of RAM and time to compute !!!
# You can just simply load into your working env TWO FILES edges_test_short.rds, edges_test_long.rds from base_dir to skip this step


genes_coexpressed_short <- unique(c(tom_network_short$tom_df$Gene1, tom_network_short$tom_df$Gene2))
genes_coexpressed_long <- unique(c(tom_network_long$tom_df$Gene1, tom_network_long$tom_df$Gene2))

#coexpressed_genes_data_short <- data_short_tissue[genes_coexpressed_short, , drop = FALSE]
#TF_genes_short <- filtered_genes_tissue %>%
#    filter(gene_id %in% genes_coexpressed_short) %>% 
#    filter(!is.na(tf_family), tf_family != "") %>% 
#    pull(gene_id)   

#coexpressed_genes_data_long <- data_long_tissue[genes_coexpressed_long, , drop = FALSE]
#TF_genes_long <- filtered_genes_tissue %>%
#    filter(gene_id %in% genes_coexpressed_long) %>% 
#    filter(!is.na(tf_family), tf_family != "") %>% 
#    pull(gene_id)

#weight_matrix_short <- GENIE3(coexpressed_genes_data_short, regulators = TF_genes_short, nCores = 30, verbose = TRUE)

#weight_matrix_long <- GENIE3(coexpressed_genes_data_long, regulators = TF_genes_long, nCores = 30, verbose = TRUE)

# For statistical evaluation we will use function test_edges from DIANE
# install.packages("https://cran.r-project.org/src/contrib/Archive/swfscMisc/swfscMisc_1.6.6.tar.gz",  repos = NULL, type = "source")
# install.packages("https://cran.r-project.org/src/contrib/Archive/rfPermute/rfPermute_2.5.4.tar.gz",  repos = NULL, type = "source")
# install.packages("HDInterval", "kknn", "modeest")
# remotes::install_github("OceaneCsn/DIANE")

# backup all before DIANE computation

#save.image(file.path(output_dir_3, "3.Before_permutations_test.Rdata"))

#library(DIANE)
#edges_test_short <- test_edges(
#   weight_matrix_short,
#   coexpressed_genes_data_short,
#   nrow(coexpressed_genes_data_short),
#   length(TF_genes_short),
#   density = 0.065,
#   nTrees = 1000,
#   nShuffle = 10000,
#   nCores = 30,
#   verbose = TRUE
#)

#saveRDS(edges_test_short, file.path(base_dir, "edges_test_short.rds"))

#edges_test_long <- test_edges(
#   weight_matrix_long,
#   coexpressed_genes_data_long,
#   nrow(coexpressed_genes_data_long),
#   length(TF_genes_long),
#   density = 0.065,
#   nTrees = 1000,
#   nShuffle = 10000,
#   nCores = 30,
#   verbose = TRUE
#)

#saveRDS(edges_test_long, file.path(base_dir, "edges_test_long.rds"))

edges_test_short <- readRDS(file.path(base_dir, paste0("edges_test_short_",tolower(analysis),".rds")))
edges_test_long <- readRDS(file.path(base_dir, paste0("edges_test_long_",tolower(analysis),".rds")))

filtered_edges_short <- edges_test_short$links[edges_test_short$links$fdr < 0.05, ]
filtered_edges_long <- edges_test_long$links[edges_test_long$links$fdr < 0.05, ]

get_common_edges <- function(
  filtered_edges,
  tom_df,
  correlation_threshold = 0.3
) {
  
  # Preparing GENIE3 edges
  genie_edges <- filtered_edges %>%
    dplyr::mutate(
      Gene1_lower = tolower(regulatoryGene),
      Gene2_lower = tolower(targetGene),
      pair1 = paste(pmin(Gene1_lower, Gene2_lower), pmax(Gene1_lower, Gene2_lower), sep = "_")
    ) %>%
    dplyr::select(regulatoryGene, targetGene, pval, fdr, pair1)
  
  # Preparing coexpression edges
  coexpr_pairs <- tom_df %>%
    dplyr::mutate(
      Gene1_lower = tolower(Gene1),
      Gene2_lower = tolower(Gene2),
      pair1 = paste(pmin(Gene1_lower, Gene2_lower), pmax(Gene1_lower, Gene2_lower), sep = "_")
    ) %>%
    dplyr::filter(abs(Correlation) >= correlation_threshold) %>%
    dplyr::select(pair1 = pair1, Correlation)
  
  # Intersection of both pairse and >= WGCNA correlation threshold
  common_edges <- genie_edges %>%
    dplyr::inner_join(coexpr_pairs, by = "pair1") %>%
    dplyr::mutate(source = "GENIE&WGCNA") %>%
    dplyr::rename(Gene1 = regulatoryGene, Gene2 = targetGene) %>%
    dplyr::select(Gene1, Gene2, Correlation, source)
  
  # Remove duplicate edges
  common_edges_unique <- common_edges %>%
    rowwise() %>%
    mutate(gene_a = min(Gene1, Gene2), gene_b = max(Gene1, Gene2)) %>%
    ungroup() %>%
    dplyr::distinct(gene_a, gene_b, .keep_all = TRUE) %>%
    dplyr::select(-gene_a, -gene_b) %>%
    dplyr::arrange(desc(Correlation))
  
  return(common_edges_unique)
}

common_edges_short <- get_common_edges(filtered_edges_short, tom_network_short$tom_df, correlation_threshold = 0.3)
common_edges_long <- get_common_edges(filtered_edges_long, tom_network_long$tom_df, correlation_threshold = 0.3)

write.table(common_edges_short, file.path(output_dir_3, paste0("3.WGCNA_GENIE3_",tolower(analysis),"_short.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
write.table(common_edges_long, file.path(output_dir_3, paste0("3.WGCNA_GENIE3_",tolower(analysis),"_long.tsv")), sep="\t", quote=FALSE, row.names=FALSE)


#extract genes associated with high light response according to Gene Ontology
known_genes <- go_data %>%
  dplyr::filter(term == "response to light stimulus") %>%
  dplyr::distinct(GeneID)

#Function to add metadata to network
annotate_interactions <- function(df, classification, count_columns, known_genes, photoreceptor_genes = NULL) {
  
  # Transform count_columns to symbols for across()
  count_vars <- rlang::syms(count_columns)
  
  # Counting the total expression occurencies in target groups 
  classification_ext <- classification %>%
    dplyr::mutate(
      total = rowSums(dplyr::select(., !!!count_vars), na.rm = TRUE)
    )
  
  # Main part of annotation
  result <- df %>%
    # Add data for Gene1
    dplyr::left_join(
      classification_ext %>% 
        dplyr::select(gene_id, classification, direction, total, tf_family),
      by = c("Gene1" = "gene_id")
    ) %>%
    dplyr::rename(
      Gene1_class = classification,
      Gene1_direction = direction,
      Gene1_total = total,
      Gene1_tf_family = tf_family
    ) %>%
    
    # Add data for Gene1
    dplyr::left_join(
      classification_ext %>% 
        dplyr::select(gene_id, classification, direction, total, tf_family),
      by = c("Gene2" = "gene_id")
    ) %>%
    dplyr::rename(
      Gene2_class = classification,
      Gene2_direction = direction,
      Gene2_total = total,
      Gene2_tf_family = tf_family
    ) %>%
    
    # Remove non-informative edges
    dplyr::filter(!is.na(Gene1_class), !is.na(Gene2_class)) %>%
    
    # Flags
    dplyr::mutate(
      Gene1_known = Gene1 %in% known_genes,
      Gene2_known = Gene2 %in% known_genes,
      Gene1_TF = !is.na(Gene1_tf_family),
      Gene2_TF = !is.na(Gene2_tf_family),
      Gene1_photoreceptor = if (!is.null(photoreceptor_genes)) Gene1 %in% photoreceptor_genes else FALSE,
      Gene2_photoreceptor = if (!is.null(photoreceptor_genes)) Gene2 %in% photoreceptor_genes else FALSE,
      
      # Labels
      Gene1_label = dplyr::case_when(
        Gene1_photoreceptor ~ paste("Photoreceptor,", Gene1_class),
        Gene1_TF ~ paste("TF,", Gene1_class),
        TRUE ~ Gene1_class
      ),
      Gene2_label = dplyr::case_when(
        Gene2_photoreceptor ~ paste("Photoreceptor,", Gene2_class),
        Gene2_TF ~ paste("TF,", Gene2_class),
        TRUE ~ Gene2_class
      ),
      
      # Defining interaction types
      interaction_type = dplyr::case_when(
        Gene1_photoreceptor & Gene2_photoreceptor ~ "photo-photo",
        Gene1_photoreceptor & Gene2_TF ~ "photo-TF",
        Gene1_TF & Gene2_photoreceptor ~ "TF-photo",
        Gene1_photoreceptor & !Gene2_TF ~ "photo-gene",
        !Gene1_TF & Gene2_photoreceptor ~ "gene-photo",
        Gene1_TF & !Gene2_TF ~ "TF-gene",
        !Gene1_TF & Gene2_TF ~ "gene-TF",
        Gene1_TF & Gene2_TF ~ "TF-TF",
        TRUE ~ "gene-gene"
      ),
      
      same_class = Gene1_class == Gene2_class,
      same_direction = Gene1_direction == Gene2_direction
    ) %>%
    
    # Remove gene-gene interactions (if any)
    dplyr::filter(interaction_type != "gene-gene") %>%
    
    # Final ordering
    dplyr::select(
      Gene1, Gene2, Correlation,
      Gene1_class, Gene2_class,
      Gene1_direction, Gene2_direction,
      Gene1_total, Gene2_total,
      Gene1_TF, Gene2_TF,
      Gene1_photoreceptor, Gene2_photoreceptor,
      Gene1_known, Gene2_known,
      interaction_type,
      same_class, same_direction,
      source, Gene1_label, Gene2_label,
      dplyr::everything()
    )
  
  return(result)
}

selected_genes_short$tf_family[selected_genes_short$tf_family == ""] <- NA
coexpr_short_annotated <- annotate_interactions(
  df = common_edges_short,
  classification = selected_genes_short,
  count_columns = c("Short_Up", "Short_Down"), #("Short_Up", "Short_Down"),
  known_genes = known_genes$GeneID,
  photoreceptor_genes = photoreceptor_genes_data$gene_id
)

selected_genes_long$tf_family[selected_genes_long$tf_family == ""] <- NA
coexpr_long_annotated <- annotate_interactions(
  df = common_edges_long,
  classification = selected_genes_long,
  count_columns = c("Long_Up", "Long_Down"), #("Short_Up", "Short_Down"),
  known_genes = known_genes$GeneID,
  photoreceptor_genes = photoreceptor_genes_data$gene_id
)

# Save annotated networks in table format
readr::write_tsv(coexpr_short_annotated, file.path(output_dir_3, paste0("3.WGCNA_GENIE3_short_",tolower(analysis),"_annotated.tsv")))
readr::write_tsv(coexpr_long_annotated, file.path(output_dir_3, paste0("3.WGCNA_GENIE3_long_",tolower(analysis),"_annotated.tsv")))

# Cluster analysis

analyze_and_annotate_clusters <- function(coexpr_annotated, filtered_genes_tissue, correlation_threshold = 0.3) {
    # Correlation threshold filtering
    coexpr_annotated <- coexpr_annotated %>%
        dplyr::filter(Correlation >= correlation_threshold)
    
    print("After correlation filtering")
    print(head(coexpr_annotated))
    
    # Selecting TF-gene edges and defining the clusters
    tf_edges <- coexpr_annotated %>%
        dplyr::filter(
            interaction_type == "TF-gene" | 
                (Gene1_TF & !Gene2_TF) | 
                (!Gene1_TF & Gene2_TF)
        )
    
    print("TF-gene edges:")
    print(head(tf_edges))
    
    if (nrow(tf_edges) == 0) {
        stop("No TF-target interactions found at the given correlation threshold")
    }
    
    # Creating TF-target pairs list
    tf_targets <- dplyr::bind_rows(
        tf_edges %>% 
            dplyr::filter(Gene1_TF) %>% 
            dplyr::select(TF = Gene1, target = Gene2),
        tf_edges %>% 
            dplyr::filter(Gene2_TF) %>% 
            dplyr::select(TF = Gene2, target = Gene1)
    ) %>%
        dplyr::distinct()
    
    print("List of TF-target pairs:")
    print(head(tf_targets))
    
    # Building the graph and clusters search
    tf_graph <- igraph::graph_from_data_frame(tf_targets, directed = FALSE)
    tf_clusters <- igraph::components(tf_graph)
    
    print("Clusters:")
    print(names(tf_clusters$membership))
    
    # Creating table of clusters-genes relationship
    cluster_membership <- data.frame(
        gene = names(tf_clusters$membership),
        cluster = tf_clusters$membership
    ) %>%
        dplyr::mutate(cluster = as.integer(cluster))
    
    print("Clusters membership table")
    print(head(cluster_membership))
    
    # Filtering out the small clusters
    cluster_sizes <- cluster_membership %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(size = n(), .groups = "drop")
    
    valid_clusters <- cluster_sizes %>%
        dplyr::filter(size >= 10) %>%
        dplyr::pull(cluster)
    
    if (length(valid_clusters) == 0) {
        stop("No clusters with at least 10 genes found")
    }
    
    cluster_membership <- cluster_membership %>%
        dplyr::filter(cluster %in% valid_clusters)
    
    print("Clusters after filtering")
    print(head(cluster_membership))
    
    # Adding clustering data to input data
    edges_with_clusters <- coexpr_annotated %>%
        dplyr::left_join(
            cluster_membership %>% dplyr::rename(Gene1 = gene, cluster1 = cluster),
            by = "Gene1"
        ) %>%
        dplyr::left_join(
            cluster_membership %>% dplyr::rename(Gene2 = gene, cluster2 = cluster),
            by = "Gene2"
        ) %>%
        dplyr::mutate(
            cluster = dplyr::coalesce(cluster1, cluster2)
        ) %>%
        dplyr::select(-cluster1, -cluster2)
    
    print("Information about clusters")
    print(head(edges_with_clusters))
    
    # Add genes without the clusters to separate ones (if their annotations is not TF-gene)
    all_genes <- unique(c(edges_with_clusters$Gene1, edges_with_clusters$Gene2))
    clustered_genes <- unique(cluster_membership$gene)
    unclustered_genes <- setdiff(all_genes, clustered_genes)
    
    if (length(unclustered_genes) > 0) {
        max_cluster <- max(cluster_membership$cluster, 0)
        new_clusters <- data.frame(
            gene = unclustered_genes,
            cluster = seq(max_cluster + 1, length.out = length(unclustered_genes)))
        
        all_cluster_membership <- dplyr::bind_rows(cluster_membership, new_clusters)
        
        # Add to main df
        edges_with_clusters <- edges_with_clusters %>%
            dplyr::left_join(
                all_cluster_membership %>% dplyr::rename(Gene1 = gene, cluster1_new = cluster),
                by = "Gene1"
            ) %>%
            dplyr::left_join(
                all_cluster_membership %>% dplyr::rename(Gene2 = gene, cluster2_new = cluster),
                by = "Gene2"
            ) %>%
            dplyr::mutate(
                cluster = dplyr::coalesce(cluster, cluster1_new, cluster2_new)
            ) %>%
            dplyr::select(-cluster1_new, -cluster2_new)
    }
    
    print("After adding step:")
    print(head(edges_with_clusters))
    
    # Again, filtering small clusters
    final_cluster_sizes <- edges_with_clusters %>%
        tidyr::pivot_longer(cols = c(Gene1, Gene2), names_to = "type", values_to = "gene") %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(size = n_distinct(gene), .groups = "drop")
    
    valid_final_clusters <- final_cluster_sizes %>%
        dplyr::filter(size >= 10) %>%
        dplyr::pull(cluster)
    
    edges_filtered <- edges_with_clusters %>%
        dplyr::filter(cluster %in% valid_final_clusters)
    
    print("Finally filtered data:")
    print(head(edges_filtered))
    
    # Rename N clusters in 1-N order
    unique_clusters <- sort(unique(edges_filtered$cluster))
    edges_filtered <- edges_filtered %>%
        dplyr::mutate(
            cluster = as.integer(factor(cluster, levels = unique_clusters))
        )
    
    print("Renamed clusters")
    print(head(edges_filtered))
    
    # Gene classes counting in clusters
    count_class_freq <- function(cluster_data) {
        all_genes <- unique(c(cluster_data$Gene1, cluster_data$Gene2))
        
        gene_classes <- filtered_genes_tissue %>%
            dplyr::filter(gene_id %in% all_genes) %>%
            dplyr::select(gene_id, classification) %>%
            dplyr::distinct()
        
        class_freq <- gene_classes %>%
            dplyr::count(classification, name = "count") %>%
            tidyr::complete(classification = unique(filtered_genes_tissue$classification), fill = list(count = 0))
        
        return(class_freq)
    }
    
    cluster_class_freq <- edges_filtered %>%
        dplyr::group_by(cluster) %>%
        dplyr::group_modify(~ count_class_freq(.x)) %>%
        dplyr::ungroup()
    
    print("Statistics on gene classes in clusters:")
    print(cluster_class_freq)
    
    # Background for Statistical Tests
    all_genes_in_clusters <- unique(c(edges_filtered$Gene1, edges_filtered$Gene2))
    
    background_class_counts <- filtered_genes_tissue %>%
        dplyr::filter(gene_id %in% all_genes_in_clusters) %>%
        dplyr::group_by(classification) %>%
        dplyr::summarise(total_in_background = n(), .groups = 'drop') %>%
        dplyr::filter(!is.na(classification))
    
    total_background <- length(all_genes_in_clusters)
    
    print("Background for statistical tests:")
    print(background_class_counts)
    
    # Fisher's Test for class enrichment
    perform_fisher_test <- function(class_count, total_in_cluster, class_background, total_background) {
        contingency <- matrix(
            c(
                class_count,
                total_in_cluster - class_count,
                class_background - class_count,
                total_background - total_in_cluster - (class_background - class_count)
            ),
            nrow = 2
        )
        
        test_result <- tryCatch({
            fisher.test(contingency)
        }, error = function(e) {
            return(list(estimate = NA, p.value = NA))
        })
        
        return(tibble::tibble(
            odds_ratio = as.numeric(test_result$estimate),
            p_value = test_result$p.value
        ))
    }
    
    # First, we create a temporary dataframe with the test results
    temp_stats <- cluster_class_freq %>%
        dplyr::left_join(background_class_counts, by = "classification") %>%
        dplyr::group_by(cluster, classification) %>%
        dplyr::mutate(
            total_in_cluster = sum(count)
        ) %>%
        dplyr::mutate(
            fisher_test = purrr::pmap(list(
                class_count = count,
                total_in_cluster = total_in_cluster,
                class_background = total_in_background,
                total_background = total_background
            ), perform_fisher_test)
        ) %>%
        tidyr::unnest_wider(fisher_test)
    
    # Intermediate output of temporary dataframe with results
    print("Intermediate dataframe with test results:")
    print(head(temp_stats))
    
    # Apply the Benjamini-Hochberg correction to all
    temp_stats <- temp_stats %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
            p_adjusted = p.adjust(p_value, method = "BH")
        )

    print("Dataframe with adjusted p-values:")
    print(head(temp_stats))
    
    # Forming the final result
    cluster_class_stats <- temp_stats %>%
        dplyr::select(cluster, classification, count, total_in_cluster, total_in_background, 
                     odds_ratio, p_value, p_adjusted) %>%
        dplyr::arrange(cluster, p_adjusted)
    
    print("Summary statistics for gene classes:")
    print(cluster_class_stats)
    
    # GO Analysis for Clusters
    if (!requireNamespace("org.At.tair.db", quietly = TRUE)) {
        warning("org.At.tair.db package not installed, skipping GO analysis")
        go_data <- NULL
        top_go_per_cluster <- NULL
    } else {
        cluster_gene_lists <- edges_filtered %>%
            dplyr::select(Gene1, Gene2, cluster) %>%
            tidyr::pivot_longer(cols = c(Gene1, Gene2), names_to = "gene_role", values_to = "gene_id") %>%
            dplyr::distinct(cluster, gene_id) %>%
            dplyr::group_by(cluster) %>%
            dplyr::summarise(genes = list(unique(gene_id)), .groups = "drop") %>%
            tibble::deframe()
        
        run_go_analysis <- function(genes, cluster_num) {
            ego <- tryCatch({
                clusterProfiler::enrichGO(
                    gene = genes,
                    OrgDb = org.At.tair.db::org.At.tair.db,
                    keyType = "TAIR",
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE
                )
            }, error = function(e) {
                return(NULL)
            })
            
            if (is.null(ego) || nrow(ego) == 0) return(NULL)
            
            go_results <- as.data.frame(ego)
            go_results$cluster <- cluster_num
            
            return(go_results)
        }
        
        all_go_results <- purrr::map_dfr(
            names(cluster_gene_lists),
            ~ run_go_analysis(cluster_gene_lists[[.x]], .x)
        )
        
        if (!is.null(all_go_results) && nrow(all_go_results) > 0) {
            go_data <- all_go_results %>%
                tidyr::separate(GeneRatio, into = c("gene_count", "total_in_cluster"), sep = "/") %>%
                dplyr::mutate(
                    across(c(gene_count, total_in_cluster), as.numeric),
                    percent_genes = (gene_count / total_in_cluster) * 100
                )
            
            top_go_per_cluster <- go_data %>%
                dplyr::group_by(cluster) %>%
                dplyr::arrange(p.adjust, .by_group = TRUE) %>%
                dplyr::slice_head(n = 3) %>%
                dplyr::ungroup()
        } else {
            go_data <- NULL
            top_go_per_cluster <- NULL
        }
    }
    
    print("Go analysis results")
    print(top_go_per_cluster)
    
    # Returning results
    return(list(
        annotated_edges = edges_filtered,
        cluster_stats = cluster_class_stats,
        go_all = go_data,
        go_top = top_go_per_cluster
    ))
}

if (analysis == "LEAVES") {
  correlation_threshold <- 0.5
} else if (analysis == "SEEDLINGS") {
  correlation_threshold <- 0.6
}

output_short_cl <- analyze_and_annotate_clusters(coexpr_short_annotated, selected_genes_short, correlation_threshold = correlation_threshold)
genie_coexpr_annot_with_cluster_short <- output_short_cl$annotated_edges
clusters_short <- output_short_cl$cluster_stats

output_long_cl <- analyze_and_annotate_clusters(coexpr_long_annotated, selected_genes_long)
genie_coexpr_annot_with_cluster_long <- output_long_cl$annotated_edges
clusters_long <- output_long_cl$cluster_stats

go_all_short <- output_short_cl$go_all
go_top_short <- output_short_cl$go_top

go_all_long <- output_long_cl$go_all
go_top_long <- output_long_cl$go_top

readr::write_tsv(go_all_short, file.path(output_dir_3, paste0("3.clusters_short_",tolower(analysis),"_all_enriched.tsv")))
readr::write_tsv(go_all_long, file.path(output_dir_3, paste0("3.clusters_long_",tolower(analysis),"_all_enriched.tsv")))

readr::write_tsv(genie_coexpr_annot_with_cluster_short, file = file.path(output_dir_3, paste0("3.WGCNA_GENIE3_short_",tolower(analysis),"_annotated_fin.tsv")))
readr::write_tsv(genie_coexpr_annot_with_cluster_long, file = file.path(output_dir_3, paste0("3.WGCNA_GENIE3_long_",tolower(analysis),"_annotated_fin.tsv")))

# Build and visualize network with metadata
visualize_network <- function(edge_list, output_file, layout_alg="nicely", nudge_text = 0.15, width = 40, height = 40, go_terms_df = NULL) {
    
    required_columns <- c("Gene1", "Gene2", "Correlation",
                          "Gene1_direction", "Gene2_direction",
                          "Gene1_total", "Gene2_total",
                          "Gene1_label", "Gene2_label",
                          "Gene1_tf_family", "Gene2_tf_family",
                          "Gene1_photoreceptor", "Gene2_photoreceptor",
                          "Gene1_known", "Gene2_known",
                          "source")
    
    missing <- required_columns[!required_columns %in% names(edge_list)]
    if (length(missing) > 0) {
        stop("Missing required columns in edge_list: ", paste(missing, collapse = ", "))
    }
    
    # Getting Gene Symbols from Ensembl
    mart <- biomaRt::useMart(
        host = "https://plants.ensembl.org", 
        biomart = "plants_mart",
        dataset = "athaliana_eg_gene"
    )
    
    gene_ids <- unique(c(edge_list$Gene1, edge_list$Gene2))
    gene_symbols <- tryCatch({
        biomaRt::getBM(
            attributes = c("ensembl_gene_id", "external_gene_name"),
            filters = "ensembl_gene_id",
            values = gene_ids,
            mart = mart
        )
    }, error = function(e) {
        warning("Failed to retrieve gene symbols from BioMart")
        tibble(ensembl_gene_id = gene_ids, external_gene_name = gene_ids)
    })
    
    edge_list <- edge_list %>%
        dplyr::mutate(
            Gene1_short = gene_symbols$external_gene_name[match(Gene1, gene_symbols$ensembl_gene_id)],
            Gene2_short = gene_symbols$external_gene_name[match(Gene2, gene_symbols$ensembl_gene_id)],
            Gene1_short = ifelse(is.na(Gene1_short), Gene1, Gene1_short),
            Gene2_short = ifelse(is.na(Gene2_short), Gene2, Gene2_short)
        )
    
    # Mapping gene names to unique IDs for tidygraph
    gene_names <- unique(c(edge_list$Gene1, edge_list$Gene2))
    node_mapping <- tibble::tibble(name = gene_names) %>%
        dplyr::mutate(id = row_number())
    
    # Collecting node attributes (without clusters)
    node_attrs <- dplyr::bind_rows(
        edge_list %>%
            dplyr::select(name = Gene1,
                          Class = Gene1_label,
                          Direction = Gene1_direction,
                          Frequency = Gene1_total,
                          Known = Gene1_known,
                          ShortName = Gene1_short,
                          TF_family = Gene1_tf_family,
                          photoreceptor = Gene1_photoreceptor),
        edge_list %>%
            dplyr::select(name = Gene2,
                          Class = Gene2_label,
                          Direction = Gene2_direction,
                          Frequency = Gene2_total,
                          Known = Gene2_known,
                          ShortName = Gene2_short,
                          TF_family = Gene2_tf_family,
                          photoreceptor = Gene2_photoreceptor)
    ) %>%
        dplyr::distinct(name, .keep_all = TRUE) %>%
        dplyr::mutate(
            Direction = as.character(Direction),
            TF_type = ifelse(!is.na(TF_family) & TF_family != "", "TF", "Gene")
        )
    
    # Add clusters if they are present in the data
    if ("cluster" %in% names(edge_list)) {
        node_clusters <- dplyr::bind_rows(
            edge_list %>% dplyr::select(name = Gene1, cluster),
            edge_list %>% dplyr::select(name = Gene2, cluster)
        ) %>%
            dplyr::group_by(name) %>%
            dplyr::summarise(cluster = min(cluster, na.rm = TRUE)) %>%
            dplyr::ungroup()
        
        node_attrs <- node_attrs %>%
            dplyr::left_join(node_clusters, by = "name")
    }
    
    # Preparing the edges
    edges_fixed <- edge_list %>%
        dplyr::mutate(
            edge_color = ifelse(Correlation >= 0.5, "#006400", "gray"),
            edge_width = 1
        ) %>%
        dplyr::left_join(node_mapping, by = c("Gene1" = "name")) %>%
        dplyr::rename(from = id) %>%
        dplyr::left_join(node_mapping, by = c("Gene2" = "name")) %>%
        dplyr::rename(to = id) %>%
        dplyr::select(from, to, edge_color, edge_width, Correlation, source)
    
    # Creating a graph
    graph_data <- tidygraph::as_tbl_graph(edges_fixed, directed = FALSE) %>%
        tidygraph::activate(nodes) %>%
        dplyr::mutate(
            id = 1:dplyr::n(),
            name = node_mapping$name[match(id, node_mapping$id)]
        ) %>%
        dplyr::left_join(node_attrs, by = "name")
    
    # Coordinates calculations
    layout <- ggraph::create_layout(graph_data, layout = layout_alg)
    
    # Preparing data for clusters
    if ("cluster" %in% names(node_attrs)) {
        node_coords <- layout %>%
            tibble::as_tibble() %>%
            dplyr::select(name, x, y, Direction, cluster) %>%
            dplyr::mutate(cluster = as.character(cluster))
        
        cluster_stats <- node_coords %>%
            dplyr::group_by(cluster) %>%
            dplyr::summarise(
                Up = sum(Direction == "Up", na.rm = TRUE),
                Down = sum(Direction == "Down", na.rm = TRUE),
                Mixed = sum(Direction == "Mixed", na.rm = TRUE)
            ) %>%
            tidyr::pivot_longer(cols = c(Up, Down, Mixed), names_to = "type", values_to = "count") %>%
            dplyr::filter(count > 0) %>%
            tidyr::unite("stat", c(type, count), sep = ": ") %>%
            dplyr::group_by(cluster) %>%
            dplyr::summarise(label_stat = paste(stat, collapse = ", "))

        if (!is.null(go_terms_df)) {
            go_terms <- go_terms_df %>%
                dplyr::group_by(cluster) %>%
                dplyr::slice_head(n = 3) %>%
                dplyr::summarise(go_labels = paste(Description, collapse = "\n"))
        } else {
            go_terms <- NULL
        }
        
        cluster_labels <- node_coords %>%
            dplyr::group_by(cluster) %>%
            dplyr::summarise(
                x = mean(x),
                y = min(y) - 2
            ) %>%
            dplyr::left_join(cluster_stats, by = "cluster") %>%
            dplyr::left_join(if (!is.null(go_terms)) go_terms else tibble(cluster = character(), go_labels = character()), 
                           by = "cluster") %>%
            dplyr::mutate(
                label = paste0("Cluster ", cluster),
                stats_label = ifelse(is.na(label_stat), label, paste0(label, "\n", label_stat)),
                full_label = ifelse(is.na(go_labels), stats_label, paste0(stats_label, "\n", go_labels))
            )
    } else {
        cluster_labels <- NULL
    }
    
    # Colors for expression visualization
    direction_colors <- c(
        "Down" = "#3366CC",
        "Mixed" = "#9966CC",
        "Up" = "#CC3333"
    )
    
    # Network visualization
    network_plot <- ggraph::ggraph(layout) + 

        ggraph::geom_edge_link(
            aes(width = edge_width, color = edge_color),
            alpha = 0.4,
            lineend = "round"
        ) +
        
        ggraph::geom_node_point(
            aes(fill = Direction, shape = TF_type, size = TF_type),
            color = "black",
            alpha = 0.4,
            stroke = 0.5
        ) +
        
        ggplot2::scale_size_manual(
            name = "Gene Type",
            values = c("TF" = 15, "Gene" = 8),
            guide = "none"
        ) +
        
        ggplot2::scale_shape_manual(
            name = "Gene Type",
            values = c("TF" = 24, "Gene" = 21),
            guide = guide_legend(override.aes = list(fill = "grey50"))
        ) +
        
        ggplot2::scale_fill_manual(
            name = "Expression Direction",
            values = direction_colors,
            na.value = "gray",
            breaks = c("Up", "Down", "Mixed")
        ) +
        
        ggraph::geom_node_text(
            aes(label = ShortName),
            size = 4.5,
            repel = TRUE,
            nudge_y = nudge_text
        ) +
        
        ggraph::scale_edge_color_identity(
            name = "Edge Type",
            guide = "legend",
            labels = c("Strong", "Weak"),
            breaks = c("#006400", "gray")
        ) +
        
        ggplot2::guides(
            shape = guide_legend(order = 1, override.aes = list(size = 6)),
            fill = guide_legend(order = 2, override.aes = list(shape = 21)),
            edge_color = guide_legend(order = 3),
            edge_width = "none"
        ) +
        
        ggplot2::labs(
            title = "Gene Coexpression Network",
            subtitle = "Node color = Expression Direction\nEdge color - power of correlation"
        ) +
        ggraph::theme_graph(base_family = "sans") +
        ggplot2::theme(
            legend.position = "bottom",
            legend.box = "vertical",
            legend.spacing.y = grid::unit(0.1, "cm"),
            plot.margin = grid::unit(c(0.2, 0.2, 2.5, 0.2), "cm")
        )
    
    # 11. Add cluster labels if any
    if (!is.null(cluster_labels)) {
        network_plot <- network_plot +
            ggrepel::geom_text_repel(
                data = cluster_labels,
                aes(x = x, y = y, label = full_label),
                size = 4,
                fontface = "bold",
                color = "black",
                inherit.aes = FALSE,
                lineheight = 0.8,
                hjust = 0,
                direction = "y",
                nudge_y = -0.5,
                segment.color = NA,
                box.padding = 0.5,
                xlim = c(-Inf, Inf),
                ylim = c(-Inf, Inf)
            )
    }
    
    # Save to PDF
    ggplot2::ggsave(
        output_file,
        network_plot,
        width = width,
        height = height,
        device = cairo_pdf,
        limitsize = FALSE
    )
}

# Visualize GRN for short- and long-time response data
visualize_network(
  edge_list = genie_coexpr_annot_with_cluster_short,
  layout_alg = "kk",
  nudge_text = 0.05,
  output_file = file.path(output_dir_3, paste0("3.GRN_genie_coexpr_short_",tolower(analysis),".pdf")),
  width = 25,
  height = 25,
  go_terms_df = go_top_short
)

visualize_network(
  edge_list = genie_coexpr_annot_with_cluster_long,
  layout_alg = "nicely",
  nudge_text = 0.05,
  output_file = file.path(output_dir_3, paste0("3.GRN_genie_coexpr_long_",tolower(analysis),".pdf")),
  width = 25,
  height = 25,
  go_terms_df = go_top_long
)

# Save all data!!!
# P.S. if you got very thick edges in visualized networks, please just load again your R environment
# with 3.Backup_full_analysis_results.Rdata, library(dplyr), library(ggplot2) and run just visualize_network code

save.image(file.path(output_dir_3, paste0("3.Backup_full_analysis_results_",tolower(analysis),".Rdata")))