rm(list=ls())

library(maftools)
library(tidyverse)
library(ggsignif)
library(data.table)
library(patchwork)

setwd("~/isac1/")

#Fig 4a----------------- 

meta_data <- read.csv2("Data/metadata.csv") %>% dplyr::select(study_id,tumor_level2, tmb, Tumor_Sample_Barcode, tumor_occassion)

meta_data$isac_id <- paste("ISAC", formatC(meta_data$study_id, width = 3, format = "d", flag = "0"), sep = "")




## read maf files


maf_merge <- read.maf("Data/isac1_filtered_all_135_samples_maf_maftools.maf", clinicalData = meta_data, verbose = F) ## in total, 135 patients.

source("Scirpts/Fig4/tcgacompare.R")

isac.mutload = tcgaCompare_test(maf = maf_merge, cohortName = 'ISAC', logscale = TRUE,primarySite = F, capture_size = 35.8, 
                                rm_zero = FALSE, medianCol = "#7F3F98",bg_col = c("#EDF8B1", "white"),cohortFontSize = 1, axisFontSize = 1)



#Fig 4b----------------- 

tmb_exom <- tmb(maf_merge, captureSize = 35.8)
meta_table_subset <- maf_merge@clinical.data
tmb_merge <- merge(tmb_exom,meta_table_subset, by.x = "Tumor_Sample_Barcode", by.y = "Tumor_Sample_Barcode")

table_median <- tmb_merge %>% group_by(tumor_level2) %>% summarise(median = median(total_perMB)) %>% arrange(median)

tmb_merge$tumor_level2 <- factor(tmb_merge$tumor_level2, levels = table_median$tumor_level2)
tmb_merge[tmb_merge$Tumor_Sample_Barcode %in% c("810T", "P13714_104", "P26101_102", "P13714_126", "P25501_136"),]$tumor_occassion <- "Primary_0"
tmp <- dplyr::select(tmb_merge, tumor_level2, total_perMB, tumor_occassion, Tumor_Sample_Barcode, study_id)
bbb <- c("Neuroblastoma", 0.028, "Primary", "1149T", "ISAC046") # TMB of this sample is 0, a small value makes the plot looks better.
bbb <- as.data.frame(t(bbb))
colnames(bbb) <- colnames(tmp)
tmp <- rbind(tmp, bbb)
tmp$total_perMB <- as.numeric(tmp$total_perMB)
ggplot(tmp, mapping = aes(x=tumor_level2, y=total_perMB)) + geom_boxplot(fill = "#91928E", lwd = 0.05, outlier.color=NA) +
  geom_point(position = position_jitter(width = 0.2, height = 0), aes(color=tumor_occassion),size =2) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        panel.border = element_blank(), 
        panel.grid.major  = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15, colour = "black")) +   
  scale_color_manual(values=c("#394A89", "black", "#D6A449", "#4D2E6A"),

  ) +
  xlab("Tumor type") + ylab("TMB (per MB)") +
  labs(color = NULL) +
  scale_y_log10()


## Fig 4c-----------------------------

meta_data <- read.csv2("Data/metadata.csv") %>% dplyr::select(study_id,tumor_level2, tmb, Tumor_Sample_Barcode, tumor_occassion)
meta_data$isac_id <- paste("ISAC", formatC(meta_data$study_id, width = 3, format = "d", flag = "0"), sep = "")
meta_data <- meta_data %>% 
  mutate(tumor_occassion = ifelse(is.na(tumor_occassion), "Primary", tumor_occassion))
meta_data$study_id_new <- if_else(meta_data$tumor_occassion == "Relapse", paste(meta_data$isac_id,"_" , meta_data$tumor_occassion, sep = ""), meta_data$isac_id)
meta_data[meta_data$study_id_new == "ISAC060_Relapse",]$study_id_new <- "ISAC060_Metastasis"

tpm_table_new <- read.table("Data/count_to_tpm_remove_batch_deseq2_normalized_new_meta20240903.csv", sep = ",", header = T)

epic <- immunedeconv::deconvolute(tpm_table_new, "epic")
epic <- t(epic) %>% data.frame()
colnames(epic) <- epic[1,]
epic <- epic[-1,]
epic$id <- if_else(endsWith(row.names(epic), "_Primary"), word(row.names(epic),sep = fixed("_"), start = 1), row.names(epic))


merge_table <- merge(epic, meta_data, by.x = "id", by.y = "study_id_new")
colnames(merge_table)[1] <- "study_id_new"

values_color <- c(paletteer_d("ggsci::category20b_d3"), paletteer_d("ggsci::default_nejm"), paletteer_d("ggsci::default_jama"))
values_color_bar <- c(paletteer_d("ggsci::default_jama"), paletteer_d("ggsci::default_nejm"))

merged_table_nb <- merge_table
merged_table_nb <- merged_table_nb[,-9] 
merged_table_nb[, 2:8] <- lapply(merged_table_nb[, 2:8], as.numeric)
merged_table_nb[, 2:8] <- t(apply(merged_table_nb[, 2:8], 1, function(x) x / sum(x)))

meta_tmb <- tmb_merge %>% dplyr::select(Tumor_Sample_Barcode, total_perMB)
merged_table_nb <- merge(merged_table_nb, meta_tmb, by.x = "Tumor_Sample_Barcode", by.y = "Tumor_Sample_Barcode", all.x = TRUE) 
merged_table_nb$total_perMB <- as.numeric(merged_table_nb$total_perMB)
merged_table_nb <- dplyr::arrange(merged_table_nb, tumor_level2, total_perMB)

sample_order <- merged_table_nb$study_id_new

long_data <- merged_table_nb %>%
  pivot_longer(cols = 3:9,          
               names_to = "cell_type",    
               values_to = "value")       
long_data$value <- as.numeric(long_data$value)

values_color_bar <- c("#374E55FF","#DF8F44FF","#00A1D5FF","#FFDC91FF","#79AF97FF","#6A6599FF","#80796BFF","#BC3C29FF")

bar_plot <- ggplot(long_data) +
  geom_bar(aes(fill = cell_type, y = value, x = factor(study_id_new, levels = unique(sample_order))), 
           position = "stack", 
           stat = "identity") +
  geom_point(aes(x = factor(study_id_new, levels = unique(sample_order)), 
                 y = (log10(total_perMB) + 2) * 0.2),  
             color="black", group = 1) + 
  theme(
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.text = element_text(size = 8),  
    legend.key.size = unit(0.1, "cm"),  
    plot.margin = margin(15, 0, 0, 0, "pt"),
    axis.ticks.length.x = unit(0, "pt"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank()  
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    sec.axis = sec_axis(~ 10^((. / 0.2) - 2),  
                        name = "TMB (per MB)", 
                        breaks = c(0.1, 1, 10, 100, 500, 1000),  
                        labels = scales::label_number())  
  ) +
  scale_fill_manual(values = values_color_bar) + ylab("Cell Fractions")

merge_merge_table_select_tile <- dplyr::select(merged_table_nb, study_id_new,tumor_level2)

merge_merge_table_select_tile_long <- merge_merge_table_select_tile  %>%
  pivot_longer(cols = c(2),           
               names_to = "variable",     
               values_to = "value")      

List <- split(merge_merge_table_select_tile_long, merge_merge_table_select_tile_long$variable)

# Function for plots
myfun <- function(x, last_plot = FALSE) {
  G <- ggplot(x, aes(x = factor(study_id_new, levels = sample_order), y = variable, fill = value)) +
    geom_tile() +
    coord_fixed(ratio = 1) +   # Ensure square tiles
    theme(
      axis.text.x = if (last_plot) element_text(angle = 90, vjust = 0.5, hjust = 1) else element_blank(),
      #axis.ticks.x = if (last_plot) element_line() else element_blank(),
      axis.line.y = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = margin(0, 0, 0, 0, "pt"),
      axis.ticks.length.x = unit(0, "pt"),
      # legend.position = "none",
      panel.grid = element_blank(),
      legend.key.size = unit(0.1, "cm"),  # Adjust the symbol size
      
      legend.text = element_text(size = 8),  # Adjust the text size
      panel.background = element_blank(),
      plot.background = element_blank()
    ) + 
    scale_fill_manual(values = values_color)
  return(G)
}

List2 <- lapply(seq_along(List), function(i) myfun(List[[i]], last_plot = i == length(List)))

tile_plot <- patchwork::wrap_plots(List2, `+`, ncol = 1)

tile_plot

bar_plot + tile_plot + plot_layout(guides ="collect",ncol = 1,heights = c(2, 0.1))


## Fig 4d & 4e-----------------------------

merged_table_nb_dot <- merged_table_nb
colnames(merged_table_nb_dot)[3:9] <- c("B", "CAF", "CD4", "CD8", "Endothelial", "Macrophage", "NK")
merged_table_nb_dot <- merged_table_nb_dot %>%
  filter(!is.na(total_perMB)) 
merged_table_nb_dot$logtmb <- log10(merged_table_nb_dot$total_perMB)
values_color <- c(paletteer_d("ggsci::category20b_d3"), paletteer_d("ggsci::default_nejm"), paletteer_d("ggsci::default_jama"))

tumor_type = "Wilms tumor" # 4e, for 4d, tumor_type =  Neuroblastoma
merged_table_nb_dot <- dplyr::filter(merged_table_nb_dot,tumor_level2 == tumor_type)

coef_cor <- sapply(colnames(merged_table_nb_dot)[3:9], function(cell) {
  cor_test_result <- cor(merged_table_nb_dot[[cell]], merged_table_nb_dot$logtmb)
}) %>% data.frame()
colnames(coef_cor) <- c("coef")
coef_cor$coefab <- abs(coef_cor$coef)
long_data_dot <- merged_table_nb_dot %>%
  pivot_longer(cols = 3:9,          
               names_to = "cell_type",     
               values_to = "value")   
long_data_dot$value <- as.numeric(long_data_dot$value)
long_data_dot <- merge(long_data_dot, coef_cor, by.x = "cell_type", by.y ="row.names")

values_color_bar <- c("#374E55FF","#DF8F44FF","#00A1D5FF","#FFDC91FF","#79AF97FF","#6A6599FF","#80796BFF","#BC3C29FF")
long_data_dot$shape <- ifelse(long_data_dot$cell_type != "CD8", 16, 1) 
long_data_dot$size <- ifelse(long_data_dot$cell_type != "CD8", 1, 2)
cor_test <- cor.test(merged_table_nb_dot$CD8, merged_table_nb_dot$logtmb)
rvalue <- cor_test$estimate
pvalue <- cor_test$p.value

ggplot(data = long_data_dot) +
  geom_point(mapping = aes(x = logtmb, y = value, color = cell_type, size = factor(size), shape = factor(shape)), 
             stroke = ifelse(long_data_dot$shape == 16, 1, 0),
             alpha = 0.7) +  
  geom_smooth(
    data = subset(long_data_dot, cell_type == "CD8"),
    mapping = aes(x = logtmb, y = value),
    method = "lm",   
    se = TRUE, span = 2, fill = "#BC80BDFF", color = "black"
  ) +
  geom_text(
    data = subset(long_data_dot, cell_type == "CD8"),
    mapping = aes(x = logtmb, y = value, label = study_id_new),
    vjust = -0.5, hjust = 0.5, size = 3, color = "black"
  ) +
  scale_shape_manual(values = c(16, 1)) +  
  theme_minimal() +
  #ylim(0, 1) +
  scale_color_manual(values = values_color) +  
  #scale_x_log10() +  
  labs(y = "Cell Composition", x = "log10(Total per MB)", title = paste("Comparison between TMB and Cell Composition in ", tumor_type, sep = ""), color = "Cell type"
  ) +  
  guides(
    color = guide_legend(override.aes = list(size = 5)),  
    shape = guide_legend(override.aes = list(size = 5))  
  ) +
  theme(legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 14)  
  ) + 
  annotate("text", x = -0.3, y = 0.9, label = paste("r = ", as.numeric(rvalue), "\np= ", as.numeric(pvalue), sep = ""), 
           size = 5, color = "black", hjust = 0)
