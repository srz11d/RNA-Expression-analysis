library(readr)
library(dplyr)
library(stats)
library(gplots)  
library(ggplot2)
library(pheatmap)

# Set directory
setwd("/Users/user/Documents/Genestack")

# Read the combined_genes.csv file
main_lines <- read.csv("combined_genes.csv")

# Define the order of columns
column_order <- c(
  "Tel_ovarian_1_Cell_Lines", "Tel_ovarian_2_Cell_Lines", 
  "Tel_bladder_1_Cell_Lines", "Tel_bladder_2_Cell_Lines", 
  "Tel_cervical_1_Cell_Lines", "Tel_cervical_2_Cell_Lines", 
  "Tel_fibosarcoma_1_Cell_Lines", "Tel_fibosarcoma_2_Cell_Lines", 
  "ALT_lung_1_Cell_Lines", "ALT_lung_2_Cell_Lines", "ALT_lung_3_Cell_Lines", 
  "ALT_lung_4_Cell_Lines", "ALT_lung_5_Cell_Lines", "ALT_lung_6_Cell_Lines",
  "ALT_liver_1_Cell_Lines", "ALT_liver_2_Cell_Lines",
  "hMSC_1_Cell_Lines", "hMSC_2_Cell_Lines", "hMSC_3_Cell_Lines", 
  "hMSC_4_Cell_Lines", "hMSC_5_Cell_Lines", "hMSC_6_Cell_Lines", 
  "hMSC_7_Cell_Lines", "hMSC_8_Cell_Lines", 
  "hMSC_9_Cell_Lines", "hMSC_10_Cell_Lines", 
  "hMSC_11_Cell_Lines", "hMSC_12_Cell_Lines", "hMSC_13_Cell_Lines", 
  "hMSC_14_Cell_Lines", "hMSC_15_Cell_Lines", "hMSC_16_Cell_Lines"
)

# Extract the relevant columns in the specified order
clustering_data <- main_lines[, column_order]

# Set row names
rownames(clustering_data) <- main_lines$Genes

# Perform hierarchical clustering using Spearman correlation
distance_matrix <- as.dist(1 - cor(clustering_data, method = "spearman"))
cluster_result <- hclust(distance_matrix)

# Create a heatmap using pheatmap, disabling column clustering
heatmap <- pheatmap(clustering_data, cluster_rows = TRUE, cluster_cols = FALSE,
                    clustering_method = "complete", scale = "row",
                    main = "Hierarchical Clustering Heatmap ALT - Tel - hMSC", fontsize_row = 4)

# Save the heatmap as a PNG file
ggsave("heatmap_overlap_all.png", heatmap, width = 10, height = 8, units = "in")
