#Hierarchical analysis for Liposarcoma samples (clinical)

library(readr)
library(dplyr)
library(stats)
library(ggplot2)
library(pheatmap)

#Set directory
setwd("/Users/user/Documents/Genestack")

# Read the significant_genes_liposarcoma.csv file
significant_genes_liposarcoma <- read.csv("significant_genes_liposarcoma.csv")

# Select the columns of interest for clustering
liposarcoma_columns <- c(
  "Liposarcoma_ALT_1", "Liposarcoma_ALT_2", "Liposarcoma_ALT_3", "Liposarcoma_ALT_4", "Liposarcoma_ALT_5",
  "Liposarcoma_ALT_6", "Liposarcoma_ALT_7", "Liposarcoma_ALT_8", "Liposarcoma_ALT_9", "Liposarcoma_ALT_10",
  "Liposarcoma_Tel_1", "Liposarcoma_Tel_2", "Liposarcoma_Tel_3", "Liposarcoma_Tel_4", "Liposarcoma_Tel_5",
  "Liposarcoma_Tel_6", "Liposarcoma_Tel_7", "Liposarcoma_Tel_8"
)

# Extract the relevant columns
clustering_data_liposarcoma <- significant_genes_liposarcoma[, liposarcoma_columns]

# Set row names
rownames(clustering_data_liposarcoma) <- significant_genes_liposarcoma$Genes

# Perform hierarchical clustering using Spearman correlation
distance_matrix_liposarcoma <- as.dist(1 - cor(clustering_data_liposarcoma, method = "spearman"))
cluster_result_liposarcoma <- hclust(distance_matrix_liposarcoma)

# Create a heatmap using pheatmap
heatmap_liposarcoma <- pheatmap(clustering_data_liposarcoma, cluster_rows = TRUE, cluster_cols = TRUE,
                                clustering_method = "complete", scale = "row",
                                main = "Hierarchical Clustering Heatmap (Liposarcoma)", fontsize_row = 5)

# Save the heatmap as a PNG file
ggsave("heatmap_liposarcoma.png", heatmap_liposarcoma, width = 10, height = 8, units = "in")
