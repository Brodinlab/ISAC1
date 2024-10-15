# Script to generate figure 1:

# Load packages
library(tidyverse)
library(webr)
library(ggridges)

# Read functions
source("../functions/PieDonutCustom.R")

# Read meta data
meta_all <- read.csv2("../../Data/metadata.csv", stringsAsFactors = FALSE)

# Set factor levels
meta_all$sex <- factor(meta_all$sex, levels = c("Female", "Male"))
meta_all$tumor_type_grouped <- factor(meta_all$tumor_type_grouped, 
                                      levels = c("Lymphoma", "Brain tumor", "Neuroblastoma", 
                                                 "Retinoblastoma", "Kidney tumor", "Liver tumor", 
                                                 "Bone tumor", "Soft tissue sarcoma", 
                                                 "Germ cell tumors", "Others"))

# Figure 1a -----

meta_all %>%
  count(sex, age_year) %>%
  ggplot(aes(age_year, n, fill = sex)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  scale_x_continuous(breaks = seq(0, 18, by = 1)) +
  scale_y_continuous(breaks = seq(0, 20, by = 1)) +
  labs(y = "Number of patients", x = "Age (Years)") +
  scale_fill_manual(values = c("Male"="#D9D9D9FF", "Female"="#BEBADAFF")) +
  theme(text = element_text(color = "black"), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Figure 1b -----

meta_all %>% 
  group_by(tumor_type_grouped) %>% 
  mutate(tumor_type_grouped_n = paste0(tumor_type_grouped, " (n= ", n(), ")")) %>%
  group_by(tumor_type_grouped, tumor_level2_fig1b) %>%
  mutate(tumor_level2_n = paste0(tumor_level2_fig1b, " (n =", n(), ")")) %>%
  ungroup() %>%
  group_by(tumor_type_grouped_n, tumor_level2_n) %>%
  tally() %>%
  PieDonutCustom(aes(tumor_type_grouped_n, tumor_level2_n, count = n),
                 ratioByGroup = FALSE,
                 showRatioThreshold = 0,
                 labelpositionThreshold = 0.3,
                 showDonutName = FALSE,
                 showPieName = FALSE,
                 showRatioDonut = FALSE,
                 showRatioPie = FALSE,
                 palette = c("#BEBADAFF", "#BC80BDFF", "#80B1D3FF", "#8DD3C7FF", 
                             "#B3DE69FF", "#CCEBC5FF", "#FFFFB3FF", "#FFED6FFF", 
                             "#FDB462FF", "#FB8072FF"))

# Figure 1c -----

# Reorder tumor types by median age
p_df <- meta_all %>%
  group_by(tumor_level1_fig1c) %>%
  mutate(median_age = median(age_month)) %>%
  ungroup() %>%
  arrange(desc(median_age)) %>%
  mutate(tumor_level1_fig1c = factor(tumor_level1_fig1c, levels = unique(tumor_level1_fig1c)))
ggplot(p_df, aes(x = age_month, y = tumor_level1_fig1c, fill = tumor_level1_fig1c)) + 
  geom_density_ridges(alpha = 0.8, quantile_lines = TRUE, quantiles = 2, scale = 1,
                      jittered_points = TRUE, point_shape = 21, point_alpha = 1) +
  geom_segment(data = subset(p_df, tumor_level1_fig1c == "Choroid plexus tumor"), 
               aes(x = median_age, y = "Choroid plexus tumor", xend = median_age, yend = "Retinoblastoma")) +
  geom_segment(data = subset(p_df, tumor_level1_fig1c == "Rare tumors"), 
               aes(x = median_age, y = "Rare tumors", xend = median_age, yend = "Meningioma")) +
  geom_point(data = subset(p_df, tumor_level1_fig1c %in% c("Choroid plexus tumor", "Rare tumors")), 
             position = position_jitter(width = 0, height = 0.1, seed = 1), size = 1.5, shape = 21) +
  labs(x = "Age in months", y = "Tumor type") +
  scale_fill_manual(values = c("Bone tumor"="#BEBADAFF", 
                               "Brain tumor"="#BC80BDFF", 
                               "Germ cell tumors"="#80B1D3FF",       
                               "Kidney tumor"="#8DD3C7FF",     
                               "Liver tumor"="#B3DE69FF",        
                               "Lymphoma"="#CCEBC5FF",         
                               "Neuroblastoma"="#FFFFB3FF",          
                               "Others"="#FFED6FFF",
                               "Retinoblastoma"="#FDB462FF", 
                               "Soft tissue sarcoma"="#FB8072FF",
                               "Rare tumors"="#E31A1CFF",   
                               "Diffuse astrocytic and oligodenroglial tumors"="#FCCDE5FF",
                               "Embryonal tumors"="#D9D9D9FF",
                               "Meningioma"="#A6CEE3FF",
                               "Neuronal and mixed neuronal-glial tumors"="#FB9A99FF",
                               "Choroid plexus tumor"="#6A3D9AFF", 
                               "Ependymal tumors"="#B2DF8AFF",
                               "Pilocytic astrocytoma" = "#CAB2D6FF")) +
  theme(text = element_text(color = "black"), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none")