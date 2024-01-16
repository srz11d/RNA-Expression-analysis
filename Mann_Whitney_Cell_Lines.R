
#Mann-Whitney U Test for multiple rows with Median and Z value

#First we need to find the significant differences betweeen ALT and Tel in cell lines

library(readr)
library(dplyr)
library(stats)

#Set directory
setwd("/Users/user/Documents/Genestack")

# Load the cleaned data from Python
telomere_data_cleaned <- read_csv("telomere_data_cleaned.csv")

# Select the columns for Indoor and Outdoor samples
alt_columns <- c("ALT_lung_1", "ALT_lung_2", "ALT_lung_3", "ALT_lung_4", "ALT_lung_5", "ALT_lung_6", "ALT_liver_1", "ALT_liver_2")
tel_columns <- c("Tel_ovarian_1", "Tel_ovarian_2", "Tel_bladder_1", "Tel_bladder_2", "Tel_cervical_1", "Tel_cervical_2", "Tel_fibosarcoma_1", "Tel_fibosarcoma_2")

# Initialize vectors to store p-values, medians, and z-values
p_values <- numeric()
medians_alt <- numeric()
medians_tel <- numeric()
z_values <- numeric()

# Loop through each row and perform Mann-Whitney U test
for (i in 1:nrow(telomere_data_cleaned)) {
  Alt <- unlist(telomere_data_cleaned[i, alt_columns])
  Tel <- unlist(telomere_data_cleaned[i, tel_columns])
  
  # Perform Mann-Whitney U test
  mw_test <- wilcox.test(Alt, Tel, alternative = "two.sided")
  
  # Append the p-value to the vector
  p_values <- c(p_values, mw_test$p.value)
  
  # Calculate medians for indoor and outdoor samples
  median_alt <- median(Alt)
  median_tel <- median(Tel)
  
  # Append the medians to the vectors
  medians_alt <- c(medians_alt, median_alt)
  medians_tel <- c(medians_tel, median_tel)
  
  # Calculate the "z" value
  z_value <- mw_test$statistic
  z_values <- c(z_values, z_value)
}

# Add p-values, medians, and z-values to the Full_COG_MW.csv dataframe
telomere_data_cleaned$p_value <- p_values
telomere_data_cleaned$median_alt <- medians_alt
telomere_data_cleaned$median_tel <- medians_tel
telomere_data_cleaned$z_value <- z_values

# Save the updated dataframe to a new CSV file
write.csv(telomere_data_cleaned, "Cell_lines_MW.csv")

#Extract the significant genes

# Read the updated CSV file
telomere_data_mw <- read.csv("Cell_lines_MW.csv")

# Filter rows with p-values less than 0.05
significant_rows <- telomere_data_mw %>%
  filter(p_value < 0.05)

# Save the filtered dataframe to a new CSV file
write.csv(significant_rows, "significant_genes_cell_lines.csv", row.names = FALSE)

