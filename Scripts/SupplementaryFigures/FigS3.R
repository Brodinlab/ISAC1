# Script to generate supplementary figure 3:

# Load packages
library(tidyverse)
library(corrplot)

# Read metadata
meta_all <- read.csv2("Data/metadata.csv", stringsAsFactors = FALSE) %>% 
  filter(tumor_occassion != "Relapse" | is.na(tumor_occassion))

# Read blood cell frequency data
cell <- read.csv("Data/cytof_freq_lineage_baseline.csv", check.names = TRUE)

# Read plasma protein NPX data
protein <- read.csv("Data/olink_npx_baseline.csv", check.names = TRUE)

# Read tumor cell frequency data
tumor_nb <- read.csv("Data/immune_infiltration_cells_for_NB.csv")
tumor_wilms <- read.csv("Data/immune_infiltration_cells_for_wilms.csv")

# Helper function to process tumor data and calculate correlations
process_tumor_data <- function(tumor_data, meta, cell_data, protein_data) {
  tumor_data <- tumor_data %>%
    mutate(study_id = as.integer(substring(study_id, 5, 7))) %>%
    left_join(meta %>% select(study_id, cytof_barcode_id, olink_id)) %>%
    na.omit() %>%
    select(study_id, cytof_barcode_id, olink_id, B, CAF, CD4, CD8, Endothelial, Macrophage, NK)
  
  # Rename tumor columns
  colnames(tumor_data)[4:10] <- paste("tumor", colnames(tumor_data)[4:10], sep = "_")
  
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
nb_results <- process_tumor_data(tumor_nb, meta_all, cell, protein)
corrplot(nb_results$cor_matrix, tl.col = "black", addgrid.col = NA, col = paletteer::paletteer_c("grDevices::ArmyRose", n = 100, direction = -1))

# Supplementary figure 3b -----
wilms_results <- process_tumor_data(tumor_wilms, meta_all, cell, protein)
corrplot(wilms_results$cor_matrix, tl.col = "black", addgrid.col = NA, col = paletteer::paletteer_c("grDevices::ArmyRose", n = 100, direction = -1))
