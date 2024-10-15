# Script to generate figure 3:

# Load packages
library(tidyverse)

# Read functions
source("../functions/RandomForest.R")

# Read metadata
meta_all <- read.csv2("../../Data/metadata.csv", stringsAsFactors = FALSE)

# Read CyTOF cell frequency data
cell <- read.csv("../../Data/cytof_freq_cluster_baseline.csv", check.names = TRUE)

# Read Olink protein NPX data
protein <- read.csv("../../Data/olink_npx_baseline.csv", check.names = TRUE)

# Figure 3a & 3b -----

# Data preparation
meta <- meta_all %>% filter(!is.na(cytof_barcode_id) & !is.na(olink_id) & tumor_level0 != "Others") %>% mutate(group = tumor_level0) 
data.cell <- cell %>% filter(cytof_barcode_id %in% meta$cytof_barcode_id) %>% select(-cytof_barcode_id)
data.protein <- protein %>% filter(olink_id %in% meta$olink_id) %>% select(-olink_id)
data.ml <- cbind(data.cell, data.protein, age = meta$age_month) %>% scale()
data.ml[is.na(data.ml)] <- 0
meta$group <- factor(meta$group, levels = c("Extra Cranial Tumor", "Intra Cranial Tumor"))
group <- ifelse(as.numeric(meta$group) == 1, 0, 1)

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
  mutate(direction = ifelse(mean_diff > 0, "up", "down")) %>%
  mutate(type = case_when(feature == "age" ~ "age",
                          feature %in% colnames(data.cell) ~ "cell",
                          feature %in% colnames(data.protein) ~ "protein")) %>%
  left_join(data.frame(feature = res.rf$names.selected.proteins, imp = res.rf$importances[res.rf$index.selected.proteins])) %>%
  mutate(feature = fct_reorder(feature, imp)) %>% 
  arrange(desc(feature)) %>%
  ggplot(aes(x = feature, y = imp)) +
  geom_segment(aes(x = feature, xend = feature, y = 0, yend = imp, color = direction), alpha = 0.5, show.legend = FALSE) +
  geom_point(aes(shape = type, color = direction, fill = direction), size=4, alpha=0.7) +
  scale_color_manual(values = c("down" = "grey30", "up" = "#BC80BDFF")) +
  scale_fill_manual(values = c("down" = "grey30", "up" = "#BC80BDFF")) +
  scale_shape_manual(values = c("cell" = 21, "protein" = 23, "age" = 22)) +
  coord_flip() +
  labs(x = "Feature", y = "Permutation importance", color = "Direction", shape = "Type", title = title) +
  guides(fill = "none") +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_line(linetype = "dashed"), 
    panel.grid.major = element_line(linetype = "dashed"), 
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(color = rev(c("down" = "black", "up" = "#BC80BDFF")[imp.df$direction])),
    plot.title = element_text(size = 10)
  )

# Figure 3c -----



