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

