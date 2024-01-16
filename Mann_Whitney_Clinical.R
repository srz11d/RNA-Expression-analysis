library(readr)
library(dplyr)
library(stats)

#Finding significant differences in Liposarcoma cells (clinical)

# Set directory
setwd("/Users/user/Documents/Genestack")

# Load the cleaned data from Python
telomere_data_cleaned <- read_csv("telomere_data_cleaned.csv")

# Select the columns for Liposarcoma_ALT and Liposarcoma_Tel samples
liposarcoma_alt_columns <- c("Liposarcoma_ALT_1", "Liposarcoma_ALT_2", "Liposarcoma_ALT_3", "Liposarcoma_ALT_4", "Liposarcoma_ALT_5", "Liposarcoma_ALT_6", "Liposarcoma_ALT_7", "Liposarcoma_ALT_8", "Liposarcoma_ALT_9", "Liposarcoma_ALT_10")
liposarcoma_tel_columns <- c("Liposarcoma_Tel_1", "Liposarcoma_Tel_2", "Liposarcoma_Tel_3", "Liposarcoma_Tel_4", "Liposarcoma_Tel_5", "Liposarcoma_Tel_6", "Liposarcoma_Tel_7", "Liposarcoma_Tel_8")

# Initialize vectors to store p-values, medians, and z-values
p_values <- numeric()
medians_liposarcoma_alt <- numeric()
medians_liposarcoma_tel <- numeric()
z_values <- numeric()

# Loop through each row and perform Mann-Whitney U test
for (i in 1:nrow(telomere_data_cleaned)) {
  Liposarcoma_ALT <- unlist(telomere_data_cleaned[i, liposarcoma_alt_columns])
  Liposarcoma_Tel <- unlist(telomere_data_cleaned[i, liposarcoma_tel_columns])
  
  # Perform Mann-Whitney U test
  mw_test <- wilcox.test(Liposarcoma_ALT, Liposarcoma_Tel, alternative = "two.sided")
  
  # Append the p-value to the vector
  p_values <- c(p_values, mw_test$p.value)
  
  # Calculate medians for Liposarcoma_ALT and Liposarcoma_Tel samples
  median_liposarcoma_alt <- median(Liposarcoma_ALT)
  median_liposarcoma_tel <- median(Liposarcoma_Tel)
  
  # Append the medians to the vectors
  medians_liposarcoma_alt <- c(medians_liposarcoma_alt, median_liposarcoma_alt)
  medians_liposarcoma_tel <- c(medians_liposarcoma_tel, median_liposarcoma_tel)
  
  # Calculate the "z" value
  z_value <- mw_test$statistic
  z_values <- c(z_values, z_value)
}

# Add p-values, medians, and z-values to the Liposarcoma_MW.csv dataframe
telomere_data_cleaned$p_value_liposarcoma <- p_values
telomere_data_cleaned$median_liposarcoma_alt <- medians_liposarcoma_alt
telomere_data_cleaned$median_liposarcoma_tel <- medians_liposarcoma_tel
telomere_data_cleaned$z_value_liposarcoma <- z_values

# Save the updated dataframe to a new CSV file
write.csv(telomere_data_cleaned, "Liposarcoma_MW.csv")

# Read the updated CSV file
telomere_data_liposarcoma_mw <- read.csv("Liposarcoma_MW.csv")

# Filter rows with p-values less than 0.05
significant_liposarcoma_rows <- telomere_data_liposarcoma_mw %>%
  filter(p_value_liposarcoma < 0.05)

# Save the filtered dataframe to a new CSV file
write.csv(significant_liposarcoma_rows, "significant_genes_liposarcoma.csv", row.names = FALSE)
