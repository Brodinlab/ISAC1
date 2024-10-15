# Script to generate figure 2:

# Load packages
library(tidyverse)
library(robCompositions)
library(paletteer)
library(relaimpo)
library(ggrepel)

# Read metadata
meta_all <- read.csv2("Data/metadata.csv", stringsAsFactors = FALSE) %>% filter(tumor_occassion != "Relapse" | is.na(tumor_occassion))

# Read CyTOF cell frequency data
cell_cluster <- read.csv("Data/cytof_freq_cluster_baseline.csv", check.names = TRUE)
cell_lineage <- read.csv("Data/cytof_freq_lineage_baseline.csv", check.names = TRUE)

# Read Olink protein NPX data
protein <- read.csv("Data/olink_npx_baseline.csv", check.names = TRUE)

# Read healthy reference cell frequency data
cell_healthy <- read.csv("Data/cytof_freq_lineage_healthy.csv", check.names = TRUE)

# Figure 2b -----

# Data preparation
meta <- meta_all %>% filter(!is.na(cytof_barcode_id))
data <- cell_cluster

# Check sample order consistency in data and metadata
identical(meta$cytof_barcode_id, data$cytof_barcode_id)
data$cytof_barcode_id <- NULL

# Compute Aitchison distances between samples
aitchison_dist <- aDist(data)

# Perform classical multidimensional scaling
mds <- cmdscale(aitchison_dist)

# Add metadata
mds_df <- cbind(mds, meta)

# Generate MDS plot
ggplot(mds_df, aes(x = `1`, y = `2`, color = age_month)) +
  geom_point(size = 3) +
  labs(x = "MDS1", y = "MDS2") +
  coord_fixed() +
  scale_color_paletteer_c("grDevices::Purple-Orange", direction = -1) +
  theme(text = element_text(color = "black"), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Figure 2c -----

# Data preparation
meta <- meta_all %>% filter(!is.na(olink_id))
data <- protein

# Check sample order consistency in data and metadata
identical(meta$olink_id, data$olink_id)
data$olink_id <- NULL

# Perform principal component analysis
pca <- prcomp(data, scale. = TRUE) 

# Generate PCA plot
autoplot(pca, data = meta, color = "age_month", size = 3) +
  coord_fixed() +
  scale_color_paletteer_c("grDevices::Purple-Orange", direction = -1) +
  theme(text = element_text(color = "black"), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Figure 2d -----

# Data preparation
meta <- meta_all %>% filter(!is.na(cytof_barcode_id)) %>% mutate(age_scale = scale(age_month))
data <- cell_lineage

# Set factor levels
meta$cytof_batch <- as.factor(meta$cytof_batch)
meta$tumor_level1 <- as.factor(meta$tumor_level1)
meta$sex <- as.factor(meta$sex)

# Check sample order consistency in data and metadata
identical(meta$cytof_barcode_id, data$cytof_barcode_id)
data$cytof_barcode_id <- NULL
data_df <- cbind(meta, data)

# Fit a linear model for each major cell population
relimp_cell <- list()
for (cell in c("B.cells", "CD4.T", "CD8.T", "gdT", "Monocytes", "Neutrophils", "NK.cells")) {
  formula <- as.formula(paste(cell, "~ cytof_batch + tumor_level1 + age_scale + sex + age_scale*sex"))
  lm_cell <- lm(formula, data = data_df)
  # Calculate relative importance metrics of each variable in the model
  relimp_cell[[cell]] <- calc.relimp(lm_cell, type = "lmg", rela = TRUE)
}

# Get R-square table
relimp_cell_r2 <- do.call(rbind, lapply(relimp_cell, function(x) x@lmg))

# Regress out the contribution of batch effects
relimp_cell_r2_rm_batch <- relimp_cell_r2[,-1]

# Recalculate R2 to sum up to 1
relimp_cell_r2_rm_batch <- relimp_cell_r2_rm_batch/rowSums(relimp_cell_r2_rm_batch)

# Plot relative contribution of each factor in each cell type
relimp_cell_r2_rm_batch %>%
  as.data.frame() %>%
  rownames_to_column("feature") %>%
  gather(factor, R2, -feature) %>%
  mutate(factor = factor(factor, levels = c("age_scale:sex", "sex", "age_scale", "tumor_level1"))) %>%
  ggplot(aes(feature, R2, fill = factor)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  labs(x = "Cell type", y = "Fraction of variance", 
       title = "Variance explained by each factor in each cell type") +
  scale_fill_manual(values = c("#E69E9CFF", "#CB74ADFF", "#BF9BDDFF", "#D8AEDDFF")) +
  theme(text = element_text(color = "black"), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Figure 2e -----

protein_var <- data.frame(Proteins = names(protein[,-1]),
                          AverageExpression = colMeans(protein[,-1]),
                          StandardizedVariance = apply(protein[,-1], 2, function(x) stats::var(x)/mean(x)))
ggplot(protein_var, aes(x = StandardizedVariance, y = AverageExpression, 
                        fill = StandardizedVariance, color = StandardizedVariance)) +
  geom_point(aes(size = StandardizedVariance), alpha = 0.5) +
  xlim(0,1) +
  ylim(0,10) +
  geom_text_repel(aes(label = ifelse(StandardizedVariance > mean(protein_var$StandardizedVariance), Proteins, "")), size = 3, color = "black") +
  geom_vline(xintercept = mean(protein_var$StandardizedVariance), linetype = "dashed", color = "black") +
  labs(y = "Average Expression", x = "Standardized Variance") +
  paletteer::scale_color_paletteer_c("viridis::plasma", direction = -1) + 
  paletteer::scale_fill_paletteer_c("viridis::plasma", direction = -1) +
  scale_size(range = c(.1, 10)) +
  guides(color = guide_legend(), fill = guide_legend(), size = guide_legend()) + 
  theme(text = element_text(color = "black"), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Figure 2f -----

# Data preparation
meta <- meta_all %>% filter(!is.na(olink_id)) %>% mutate(age_scale = scale(age_month))
data <- protein[,c("olink_id", protein_var$Proteins[which(protein_var$StandardizedVariance >= mean(protein_var$StandardizedVariance))])] 

# Set factor levels
meta$olink_batch <- as.factor(meta$olink_batch)
meta$tumor_level1 <- as.factor(meta$tumor_level1)
meta$sex <- as.factor(meta$sex)

# Check sample order consistency in data and metadata
identical(meta$olink_id, data$olink_id)
data$olink_id <- NULL
data_df <- cbind(meta, data)

# Fit a linear model for each protein
relimp_protein <- list()
for (protein in colnames(data)) {
  formula <- as.formula(paste(protein, "~ olink_batch + tumor_level1 + age_scale + sex + age_scale*sex"))
  lm_protein <- lm(formula, data = data_df)
  # Calculate relative importance metrics of each variable in the model
  relimp_protein[[protein]] <- calc.relimp(lm_protein, type = "lmg", rela = TRUE)
}

# Get R-square table
relimp_protein_r2 <- do.call(rbind, lapply(relimp_protein, function(x) x@lmg))

# Regress out the contribution of batch effects
relimp_protein_r2_rm_batch <- relimp_protein_r2[,-1]

# Recalculate R2 to sum up to 1
relimp_protein_r2_rm_batch <- relimp_protein_r2_rm_batch/rowSums(relimp_protein_r2_rm_batch)

# Plot relative contribution of each factor in each cell type
relimp_protein_r2_rm_batch %>%
  as.data.frame() %>%
  rownames_to_column("feature") %>%
  gather(factor, R2, -feature) %>%
  mutate(factor = factor(factor, levels = c("age_scale:sex", "sex", "age_scale", "tumor_level1"))) %>%
  ggplot(aes(feature, R2, fill = factor)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  labs(x = "Protein", y = "Fraction of variance", 
       title = "Variance explained by each factor in each protein") +
  scale_fill_manual(values = c("#E69E9CFF", "#CB74ADFF", "#BF9BDDFF", "#D8AEDDFF")) +
  theme(text = element_text(color = "black"), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Figure 2g-2l -----

cell_healthy %>%
  mutate(age_group = factor(age_group, levels = c("Cord blood", "<1yr", "1-5yr", "5-10yr", "10-15yr", "15-18yr"))) %>%
  gather(feature, value, -group, -subject_id, -sex, -age_days, -age_months, -age_years, -sample_id, -age_group) %>%
  group_by(feature) %>%
  arrange(desc(value)) %>%
  slice(-1:-3) %>%
  ungroup() %>%
  ggplot(aes(x = age_group, y = value, fill = group)) +
  geom_boxplot(position = position_dodge(0.8), outlier.shape = NA) +
  geom_jitter(aes(color = group), 
              position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2, seed = 1109), 
              size = 0.3, alpha = 0.9) +
  facet_wrap(~feature, scales = "free") +
  scale_color_manual(values = c("ISAC" = "#654321", "Healthy reference" = "black")) +
  scale_fill_manual(values = c("Healthy reference" = "grey", "ISAC" = "orange")) +
  labs(x = "Age group", y = "Fraction of cells") + 
  theme(text = element_text(color = "black"), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
