#Hirarchical analysis for all the samples

library(readr)
library(dplyr)
library(stats)
library(gplots)  
library(ggplot2)
library(pheatmap)

# Set directory
setwd("/Users/user/Documents/Genestack")

# Read the significant_genes_cell_lines.csv file
significant_genes_cell_lines <- read.csv("significant_genes_cell_lines.csv")

# Specify all the columns for clustering
clustering_columns <- c(
  "ALT_lung_1", "ALT_lung_2", "ALT_lung_3", "ALT_lung_4", "ALT_lung_5", "ALT_lung_6",
  "ALT_liver_1", "ALT_liver_2", "Tel_ovarian_1", "Tel_ovarian_2", "Tel_bladder_1",
  "Tel_bladder_2", "Tel_cervical_1", "Tel_cervical_2", "Tel_fibosarcoma_1", "Tel_fibosarcoma_2",
  "hMSC_1", "hMSC_2", "hMSC_3", "hMSC_4", "hMSC_5", "hMSC_6", "hMSC_7", "hMSC_8",
  "normal_lung_1", "normal_lung_2", "normal_lung_3", "normal_lung_4",
  "hMSC_9", "hMSC_10", "hMSC_11", "hMSC_12", "hMSC_13", "hMSC_14", "hMSC_15", "hMSC_16",
  "Liposarcoma_ALT_1", "Liposarcoma_ALT_2", "Liposarcoma_ALT_3", "Liposarcoma_ALT_4", "Liposarcoma_ALT_5",
  "Liposarcoma_ALT_6", "Liposarcoma_ALT_7", "Liposarcoma_ALT_8", "Liposarcoma_ALT_9", "Liposarcoma_ALT_10",
  "Liposarcoma_Tel_1", "Liposarcoma_Tel_2", "Liposarcoma_Tel_3", "Liposarcoma_Tel_4", "Liposarcoma_Tel_5",
  "Liposarcoma_Tel_6", "Liposarcoma_Tel_7", "Liposarcoma_Tel_8"
)

# Extract the relevant columns
clustering_data <- significant_genes_cell_lines[, clustering_columns]

# Set row names
rownames(clustering_data) <- significant_genes_cell_lines$Genes

# Perform hierarchical clustering using Spearman correlation
distance_matrix <- as.dist(1 - cor(clustering_data, method = "spearman"))
cluster_result <- hclust(distance_matrix)

# Create a heatmap using pheatmap
heatmap <- pheatmap(clustering_data, cluster_rows = TRUE, cluster_cols = TRUE,
                    clustering_method = "complete", scale = "row",
                    main = "Hierarchical Clustering Heatmap", fontsize_row = 5)

# Save the heatmap as a PNG file
ggsave("heatmap_All_samples.png", heatmap, width = 10, height = 8, units = "in")