#Spearman correlation to create a network graph

# Load required libraries
library(readr)
library(Hmisc)

#Set directory
setwd("/Users/user/Documents/Genestack")

#The data comes from two files. In this step, the specific genes are not relevant, we 
#want to find the correlations based on the Spearman coefficient

# Read the gene expression data for REST and ALT samples
rest_data <- read.csv("Rest_genes_sign.csv", header = TRUE, row.names = 1)
alt_data <- read.csv("ALT_genes_sign.csv", header = TRUE, row.names = 1)

# Normalize the data (convert to proportions)
rest_data_normalized <- rest_data / rowSums(rest_data)
alt_data_normalized <- alt_data / rowSums(alt_data)

# Convert normalized data frames to numeric matrices
rest_data_numeric <- as.matrix(rest_data_normalized)
alt_data_numeric <- as.matrix(alt_data_normalized)

# Get the number of columns in the data matrices
num_cols_rest <- ncol(rest_data_numeric)
num_cols_alt <- ncol(alt_data_numeric)

# Initialize matrices to store the correlation coefficients and p-values
cor_matrix <- matrix(NA, nrow = num_cols_alt, ncol = num_cols_rest)
p_values <- matrix(NA, nrow = num_cols_alt, ncol = num_cols_rest)

# Calculate Spearman correlations and p-values for ALT samples vs. REST genes
for (i in 1:num_cols_alt) {
  for (j in 1:num_cols_rest) {
    cor_test_result <- cor.test(alt_data_numeric[, i], rest_data_numeric[, j], method = "spearman")
    cor_matrix[i, j] <- cor_test_result$estimate
    p_values[i, j] <- cor_test_result$p.value
  }
}

# Apply a significance threshold (e.g., p < 0.05)
significant_cor_matrix <- cor_matrix * (p_values < 0.05)

# Export the significant correlation matrix as a CSV file
write.csv(significant_cor_matrix, "alt_rest_significant_correlations.csv")



##NOTE: This code creates the matrix with the significant correlations (with standard names, let's work on this)
#It is important to change the format if the plan is to create a new graph (network)