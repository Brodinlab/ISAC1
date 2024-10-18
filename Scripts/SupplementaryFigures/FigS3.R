# Script to generate supplementary figure 3:

# Load packages
library(tidyverse)
library(corrplot)

# Read blood immune cell data
cell <- read.csv("Data/cytof_freq_lineage_baseline.csv", check.names = TRUE)

# Read plasma protein data
protein <- read.csv("Data/olink_npx_baseline.csv", check.names = TRUE)

# Read tumor cell frequency data
tumor_nb <- read.csv("Data/immune_infiltration_cells_for_NB.csv")
tumor_wilms <- read.csv("Data/immune_infiltration_cells_for_wilms.csv")

# Function to calculate correlations
calulateCor <- function(tumor_data, cell_data, protein_data) {
  
  # Merge blood immune cell and protein data
  data_blood <- data.frame(cytof_barcode_id = tumor_data$cytof_barcode_id, olink_id = tumor_data$olink_id) %>%
    left_join(cell_data) %>%
    left_join(protein_data)
  
  # Calculate correlations
  cor_matrix <- matrix(NA, nrow = ncol(tumor_data) - 3, ncol = ncol(data_blood) - 2)
  cor_p_matrix <- matrix(NA, nrow = ncol(tumor_data) - 3, ncol = ncol(data_blood) - 2)
  
  rownames(cor_matrix) <- colnames(tumor_data[,-c(1:3)])
  colnames(cor_matrix) <- colnames(data_blood[,-c(1:2)])
  rownames(cor_p_matrix) <- colnames(tumor_data[,-c(1:3)])
  colnames(cor_p_matrix) <- colnames(data_blood[,-c(1:2)])
  
  for (tumor_feature in rownames(cor_matrix)) {
    for (blood_feature in colnames(cor_matrix)) {
      cor_test <- cor.test(tumor_data[[tumor_feature]], data_blood[[blood_feature]])
      cor_matrix[tumor_feature, blood_feature] <- cor_test$estimate
      cor_p_matrix[tumor_feature, blood_feature] <- cor_test$p.value
    }
  }
  
  list(cor_matrix = cor_matrix, cor_p_matrix = cor_p_matrix)
}

# Supplementary figure 3a -----
nb_results <- calulateCor(tumor_nb, cell, protein)
corrplot(nb_results$cor_matrix, tl.col = "black", addgrid.col = NA, col = paletteer::paletteer_c("grDevices::ArmyRose", n = 100, direction = -1))

# Supplementary figure 3b -----
wilms_results <- calulateCor(tumor_wilms, cell, protein)
corrplot(wilms_results$cor_matrix, tl.col = "black", addgrid.col = NA, col = paletteer::paletteer_c("grDevices::ArmyRose", n = 100, direction = -1))
