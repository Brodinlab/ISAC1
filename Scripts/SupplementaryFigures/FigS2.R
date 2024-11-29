# Script to generate supplementary figure 2:

# Load packages
library(tidyverse)

# Read functions
source("Scripts/functions/RandomForest.R")

# Read metadata
meta_all <- read.csv2("Data/metadata.csv", stringsAsFactors = FALSE) %>% filter(tumor_occassion != "Relapse" | is.na(tumor_occassion))

# Read Olink protein NPX data
protein <- read.csv("Data/olink_npx_baseline.csv", check.names = TRUE)

# Supplementary figure 2a -----

# Data preparation
meta <- meta_all %>% filter(tumor_level1 == "Neuroblastoma" & !is.na(tmb_group) & !is.na(olink_id)) %>% 
  mutate(tmb_group = factor(tmb_group, levels = c("low", "high")))
data <- protein %>% filter(olink_id %in% meta$olink_id)
identical(meta$olink_id, data$olink_id)
data$olink_id <- NULL
group <- meta$tmb_group

# Lollipop plot
data %>%
  summarise(across(everything(), ~ mean(.[group == "high"]) - mean(.[group == "low"]))) %>%
  pivot_longer(cols = everything(), names_to = "feature", values_to = "log2fc") %>%
  mutate(feature = fct_reorder(feature, log2fc)) %>% 
  arrange(desc(feature)) %>%
  ggplot(aes(x = feature, y = log2fc)) +
  geom_segment(aes(x = feature, xend = feature, y = 0, yend = log2fc), alpha = 0.5, show.legend = FALSE) +
  geom_point(size = 4) +
  coord_flip() +
  labs(x = "Feature", y = "Log2 fold change") +
  theme_bw() +
  theme(text = element_text(color = "black"), 
        panel.background = element_blank(),
        panel.grid.major = element_blank())

# Supplementary figure 2b -----
data %>% 
  select(IL6, MUC.16, MMP12) %>%
  cbind(meta %>% select(study_id, tmb_group)) %>%
  gather(feature, value, -study_id, -tmb_group) %>%
  ggplot(aes(x = tmb_group, y = value, group = tmb_group, fill = tmb_group)) +
  geom_boxplot(color = "black", alpha = 0.2, outlier.shape = NA) +
  geom_jitter(alpha = 0.8, shape = 21) +
  facet_wrap(~feature, scales = "free") +
  scale_color_manual(values = c("grey", "#BC80BDFF")) +
  scale_fill_manual(values = c("grey", "#BC80BDFF")) +
  theme(text = element_text(color = "black"), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Supplementary figure 2c -----
## please find the codes at Scripts/Figure4/Fig4.R section Fig 4d & 4e and related supplementary figures


# Supplementary figure 2d and 2e

# Load packages
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

# Supplementary figure 2d -----
nb_results <- calulateCor(tumor_nb, cell, protein)
corrplot(nb_results$cor_matrix, tl.col = "black", addgrid.col = NA, col = paletteer::paletteer_c("grDevices::ArmyRose", n = 100, direction = -1))

# Supplementary figure 2e -----
wilms_results <- calulateCor(tumor_wilms, cell, protein)
corrplot(wilms_results$cor_matrix, tl.col = "black", addgrid.col = NA, col = paletteer::paletteer_c("grDevices::ArmyRose", n = 100, direction = -1))



