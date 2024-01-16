
#Hierarchical representation of the significantly different genes on the cell lines (ALT and Tel)

library(readr)
library(dplyr)
library(stats)
library(gplots)  
library(ggplot2)
library(pheatmap)

#Set directory
setwd("/Users/user/Documents/Genestack")

# Read the significant_genes_cell_lines.csv file
significant_genes_cell_lines <- read.csv("significant_genes_cell_lines.csv")

# Select the columns of interest for clustering
clustering_columns <- c(
  "ALT_lung_1", "ALT_lung_2", "ALT_lung_3", "ALT_lung_4", "ALT_lung_5", "ALT_lung_6",
  "ALT_liver_1", "ALT_liver_2", "Tel_ovarian_1", "Tel_ovarian_2", "Tel_bladder_1",
  "Tel_bladder_2", "Tel_cervical_1", "Tel_cervical_2", "Tel_fibosarcoma_1", "Tel_fibosarcoma_2"
)

# Extract the relevant columns
clustering_data <- significant_genes_cell_lines[, clustering_columns]

# Set row names
rownames(clustering_data) <- significant_genes_cell_lines$Genes

# Perform hierarchical clustering using Spearman correlation
distance_matrix <- as.dist(1 - cor(clustering_data, method = "spearman"))
cluster_result <- hclust(distance_matrix)

# Create a heatmap using pheatmap
pheatmap(clustering_data, cluster_rows = TRUE, cluster_cols = TRUE,
         clustering_method = "complete", scale = "row",
         main = "Hierarchical Clustering Heatmap", fontsize_row = 5)

heatmap <- pheatmap(clustering_data, cluster_rows = TRUE, cluster_cols = TRUE,
                    clustering_method = "complete", scale = "row",
                    main = "Hierarchical Clustering Heatmap", fontsize_row = 5)

# Save the heatmap as a PNG file
ggsave("heatmap_cell_lines.png", heatmap, width = 10, height = 8, units = "in")
