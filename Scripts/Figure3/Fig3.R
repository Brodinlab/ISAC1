# Script to generate figure 3:

# Load packages
library(tidyverse)
library(pheatmap)
library(caret)
library(pROC)

# Read functions
source("Scripts/functions/RandomForest.R")

# Read metadata
meta_all <- read.csv2("Data/metadata.csv", stringsAsFactors = FALSE) %>% filter(tumor_occassion != "Relapse" | is.na(tumor_occassion))

# Read CyTOF cell frequency data
cell <- read.csv("Data/cytof_freq_cluster_baseline.csv", check.names = TRUE)

# Read Olink protein NPX data
protein <- read.csv("Data/olink_npx_baseline.csv", check.names = TRUE)

# Figure 3a & 3b -----

# Data preparation
meta <- meta_all %>%
  filter(!is.na(cytof_barcode_id) & !is.na(olink_id) & tumor_level0 != "Others") %>%
  mutate(group = factor(tumor_level0, levels = c("Extra Cranial Tumor", "Intra Cranial Tumor")))

# Filter and combine cell and protein data
data.cell <- cell %>%
  filter(cytof_barcode_id %in% meta$cytof_barcode_id) %>%
  select(-cytof_barcode_id)
data.protein <- protein %>%
  filter(olink_id %in% meta$olink_id) %>%
  select(-olink_id)

# Combine and scale the data
data.ml <- cbind(data.cell, data.protein, age = meta$age_month) %>%
  scale() %>%
  replace(is.na(.), 0)  # Replace NAs with 0

# Convert group to binary
group <- as.numeric(meta$group) - 1  # "Extra Cranial Tumor" becomes 0, "Intra Cranial Tumor" becomes 1

# Random forest model fit and variable selection
res.rf <- rf_fit(data.ml, group)

# Plot AUC curve
plot(res.rf$roc_curve, print.auc = TRUE)

# Evaluate performance by nested CV
res.rf <- rf_evaluate(data.ml, group, res.rf)

# Establish significance of original performance (nested CV)
res.rf <- rf_significance(data.ml, group, res.rf)

# Plot selected features
data.ml %>%
  as.data.frame() %>%
  summarise(across(everything(), ~ mean(.[group == 1]) - mean(.[group == 0]))) %>%
  pivot_longer(cols = everything(), names_to = "feature", values_to = "mean_diff") %>%
  mutate(
    direction = ifelse(mean_diff > 0, "up", "down"),
    type = case_when(
      feature == "age" ~ "age",
      feature %in% colnames(data.cell) ~ "cell",
      feature %in% colnames(data.protein) ~ "protein"
    )
  ) %>%
  left_join(data.frame(feature = res.rf$names.selected.features, imp = res.rf$importances[res.rf$index.selected.features]), by = "feature") %>%
  na.omit() %>%
  mutate(feature = fct_reorder(feature, imp)) %>%
  arrange(desc(feature)) %>%
  ggplot(aes(x = feature, y = imp)) +
  geom_segment(aes(xend = feature, y = 0, yend = imp, color = direction), alpha = 0.5, show.legend = FALSE) +
  geom_point(aes(shape = type, color = direction, fill = direction), size = 4, alpha = 0.7) +
  scale_color_manual(values = c("down" = "grey30", "up" = "#BC80BDFF")) +
  scale_fill_manual(values = c("down" = "grey30", "up" = "#BC80BDFF")) +
  scale_shape_manual(values = c("cell" = 21, "protein" = 23, "age" = 22)) +
  coord_flip() +
  labs(x = "Feature", y = "Permutation importance") +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_line(linetype = "dashed"),
    panel.grid.major = element_line(linetype = "dashed"),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 10)
  )

# Figure 3c -----

protein %>%
  gather(feature, value, -olink_id) %>%
  left_join(meta_all %>% select(olink_id, tumor_level1)) %>%
  filter(!tumor_level1 %in% c("Ependymal tumors", "Rare tumors", "Germ cell tumors", "Choroid plexus tumor")) %>%
  filter(feature %in% c("CAIX", "CXCL13", "IFN.gamma", "IL12", "TNF")) %>%
  mutate(tumor_level1 = factor(tumor_level1, 
                               levels = c("Embryonal tumors", "Meningioma", "Liver tumor", "Kidney tumor", 
                                          "Neuroblastoma", "Lymphoma", "Retinoblastoma", "Bone tumor", 
                                          "Diffuse astrocytic and oligodenroglial tumors",
                                          "Neuronal and mixed neuronal-glial tumors", "Soft tissue sarcoma"))) %>%
  ggplot(aes(x = tumor_level1, y = value, color = tumor_level1, fill = tumor_level1)) +
  facet_wrap(~feature) +
  geom_jitter(alpha = 0.8, shape = 16) +
  geom_boxplot(color = "black", alpha = 0.5, outlier.shape = NA) +
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, by = 2.5)) +  
  scale_color_manual(values = c("#8DD3C7FF", "#FFFFB3FF", "#BEBADAFF", "#FB8072FF", 
                                "#80B1D3FF", "#FDB462FF", "#B3DE69FF", "#FCCDE5FF", 
                                "#D9D9D9FF", "#BC80BDFF", "#CCEBC5FF")) +
  scale_fill_manual(values = c("#8DD3C7FF", "#FFFFB3FF", "#BEBADAFF", "#FB8072FF", 
                               "#80B1D3FF", "#FDB462FF", "#B3DE69FF", "#FCCDE5FF", 
                               "#D9D9D9FF", "#BC80BDFF", "#CCEBC5FF")) +
  labs(y = "NPX") +
  theme(text = element_text(color = "black"), 
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Figure 3d -----

# Compute distances between tumor types
dist <- protein %>%
  left_join(meta_all %>% select(olink_id, tumor_level1)) %>%
  group_by(tumor_level1) %>%
  filter(n() > 2) %>%  # Keep tumor types with more than 2 samples
  gather(feature, value, -olink_id, -tumor_level1) %>%
  group_by(tumor_level1, feature) %>%
  summarise(median = median(value), .groups = 'drop') %>%
  spread(feature, median) %>%
  column_to_rownames("tumor_level1") %>%
  dist(method = "euclidean") %>%
  as.matrix()

# Reshape distance table into a dataframe
dist_df <- as.data.frame(as.table(dist))
colnames(dist_df) <- c("From", "To", "Distance")

# Plot heatmap and hierarchical tree
heatmap <- pheatmap(dist, cluster_rows = TRUE, cluster_cols = TRUE,
                    angle_col = 315, fontsize = 15, display_numbers = TRUE)

# Reorder tumor types according to hierarchical tree
dist_df <- dist_df %>%
  mutate(From = factor(From, levels = heatmap$tree_row$labels[heatmap$tree_row$order]),
         To = factor(To, levels = rev(heatmap$tree_row$labels[heatmap$tree_row$order])))

# Plot dotplot
dotplot <- ggplot(dist_df, aes(x = From, y = To, color = -Distance, size = -Distance)) +
  geom_point() +
  scale_color_gradient(low = "white", high = "gray30") +
  scale_size_continuous(range = c(1, 15)) +
  scale_y_discrete(position = "right") +
  labs(x = "", y = "", title = "Pairwise distances between tumor types") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0), text = element_text(color = "black"),
        plot.title = element_text(size = 10), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Add legend to dotplot
combined_legend <- guides(
  color = guide_legend(title = "-Distance", override.aes = list(shape = 16)),
  size = guide_legend(title = "-Distance", override.aes = list(shape = 16))
)

dotplot + combined_legend

# Figure 3e -----

# Data preparation
meta <- meta_all %>%
  filter(!is.na(cytof_barcode_id) & !is.na(olink_id) & tumor_level1 %in% 
           c("Kidney tumors", "Neuroblastoma", "Lymphoma", "Retinoblastoma", "Bone tumors", 
             "Diffuse astrocytic and oligodenroglial tumors", 
             "Neuronal and mixed neuronal-glial tumors", "Soft tissue sarcoma")) %>%
  mutate(group = ifelse(tumor_level1 %in% c("Kidney tumors", "Neuroblastoma", "Lymphoma", "Retinoblastoma"), "Cluster2", "Cluster3")) %>%
  mutate(group = factor(group, levels = c("Cluster3", "Cluster2")))

# Filter and combine cell and protein data
data.cell <- cell %>%
  filter(cytof_barcode_id %in% meta$cytof_barcode_id) %>%
  select(-cytof_barcode_id)
data.protein <- protein %>%
  filter(olink_id %in% meta$olink_id) %>%
  select(-olink_id)

# Combine and scale the data
data.ml <- cbind(data.cell, data.protein, age = meta$age_month) %>%
  scale() %>%
  replace(is.na(.), 0)  # Replace NAs with 0

# Convert group to binary
group <- as.numeric(meta$group) - 1  # "Extra Cranial Tumor" becomes 0, "Intra Cranial Tumor" becomes 1

# Random forest model fit and variable selection
res.rf <- rf_fit(data.ml, group)

# Plot AUC curve
plot(res.rf$roc_curve, print.auc = TRUE)

# Evaluate performance by nested CV
res.rf <- rf_evaluate(data.ml, group, res.rf)

# Establish significance of original performance (nested CV)
res.rf <- rf_significance(data.ml, group, res.rf)

# Plot selected features
data.ml %>%
  as.data.frame() %>%
  summarise(across(everything(), ~ mean(.[group == 1]) - mean(.[group == 0]))) %>%
  pivot_longer(cols = everything(), names_to = "feature", values_to = "mean_diff") %>%
  mutate(
    direction = ifelse(mean_diff > 0, "up", "down"),
    type = case_when(
      feature == "age" ~ "age",
      feature %in% colnames(data.cell) ~ "cell",
      feature %in% colnames(data.protein) ~ "protein"
    )
  ) %>%
  left_join(data.frame(feature = res.rf$names.selected.features, imp = res.rf$importances[res.rf$index.selected.features]), by = "feature") %>%
  na.omit() %>%
  mutate(feature = fct_reorder(feature, imp)) %>%
  arrange(desc(feature)) %>%
  ggplot(aes(x = feature, y = imp)) +
  geom_segment(aes(xend = feature, y = 0, yend = imp, color = direction), alpha = 0.5, show.legend = FALSE) +
  geom_point(aes(shape = type, color = direction, fill = direction), size = 4, alpha = 0.7) +
  scale_color_manual(values = c("down" = "grey30", "up" = "#BC80BDFF")) +
  scale_fill_manual(values = c("down" = "grey30", "up" = "#BC80BDFF")) +
  scale_shape_manual(values = c("cell" = 21, "protein" = 23, "age" = 22)) +
  coord_flip() +
  labs(x = "Feature", y = "Permutation importance") +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_line(linetype = "dashed"),
    panel.grid.major = element_line(linetype = "dashed"),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 10)
  )
