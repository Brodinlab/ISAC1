rm(list=ls())

library(maftools)
library(tidyverse)
library(ggsignif)
library(data.table)
library(patchwork)

#Fig 5a--------------------

setwd("~/isac1/")

source("Scirpts/Fig5/trust4_metric_functions.R")

isac_vector <- c("ISAC031", "ISAC062", "ISAC077", "ISAC100", "ISAC112", "ISAC125", "ISAC134", "ISAC141")
rna_seq <- read.delim("Data/tumor_tissue_tcr_from_trust_merge.tsv", sep = "\t") 


files <- unique(rna_seq$sample_name)
tcr_type <- "A"  ## changed to A B D G
results <- data.frame()
notcr_sample <- c()
cdr3_table <- data.frame()

for (isac_name in isac_vector) {
  file_name <- files[grepl(paste("^", isac_name, ".*report\\.tsv$",sep = ""), files)]
  cdr3 <- dplyr::filter(rna_seq, sample_name == file_name)
  cdr3 <- cdr3[-1,-1]
  cdr3 <- subset(cdr3, count > 0) %>% 
    mutate(V = as.character(V), J = as.character(J), C = as.character(C), CDR3aa = as.character(CDR3aa)) 
  cdr3$is_complete <- sapply(cdr3$CDR3aa, function(x) ifelse(x != "partial" && x != "out_of_frame" && !grepl("^_",x) && !grepl("^\\?", x),"Y","N"))
  cdr3 <- dplyr::filter(cdr3, is_complete == "Y")
  cdr3.tcr <- subset(cdr3, grepl(paste("^TR", tcr_type, sep = ""),V) | grepl(paste("^TR", tcr_type, sep = ""),J) |grepl(paste("^TR", tcr_type, sep = ""),C))
  if (length(unique(cdr3.tcr$CDR3aa)) != length(cdr3.tcr$CDR3aa)) {
    print(paste(isac_name, "-----------")) 
  }
  cdr3_table <- rbind(cdr3_table, cdr3.tcr)
  if(dim(cdr3.tcr)[1] != 0){
    cdr3.tcr$count <- as.numeric(cdr3.tcr$count)
    cdr3.tcr$frequency <- as.numeric(cdr3.tcr$frequency)
    cdr3.tcr <- cdr3.tcr %>% mutate(lib.size = sum(count)) 
    cdr3.tcr$sample <- isac_name 
    single_sample_tcr_clonality <- getClonalityTCR(isac_name,cdr3.tcr) 
    cd3_df <- as.data.frame(single_sample_tcr_clonality) %>% t() %>% data.frame()
    results <- rbind(results, cd3_df)
  }else{
    notcr_sample <- c(notcr_sample, isac_name)
  }
}
write.table(results, file = paste("clonality_of_select_samples_V", tcr_type, ".tsv", sep = ""), sep = "\t", row.names = F)

#Fig 5b-----------------------------

tcr_seq <- readr::read_csv("Data/scTCR_data_merge_paired_chain_from_tan_20240415.csv.gz")

tcr_seq$isac_id <- word(tcr_seq$Sample_Name, sep = fixed("_"), end = 1)
tcr_seq$isac_id <- ifelse(tcr_seq$isac_id == "ISAC", "ISAC134", tcr_seq$isac_id)
tcr_seq$isac_id <- sprintf("ISAC%03d", as.numeric(sub("ISAC", "", tcr_seq$isac_id)))

rna_seq <- read.delim("Data/tumor_tissue_tcr_from_trust_merge.tsv", sep = "\t") 

rna_seq_t_cells <- rna_seq[!grepl("^IG", rna_seq$V), ]
rna_seq_t_cells <- rna_seq_t_cells[grepl("^TR", rna_seq_t_cells$V), ]
rna_seq_t_cells <- dplyr::filter(rna_seq_t_cells, CDR3aa != "out_of_frame")

rna_seq_t_cells$count <- as.numeric(rna_seq_t_cells$count)

rna_seq_t_cells$CDR3aa_trunct <- str_sub(rna_seq_t_cells$CDR3aa, 2, -2)

rna_seq_t_cells$isac_id <- word(rna_seq_t_cells$sample_name, sep = fixed("_"), end = 1)

intersect_id <- intersect(rna_seq_t_cells$isac_id, tcr_seq$isac_id)

tcr_seq <- dplyr::filter(tcr_seq, isac_id %in% intersect_id)

tcr_seq <- tcr_seq %>%
  distinct(isac_id, TCR_Alpha_Gamma_CDR3_Translation_Dominant, .keep_all = TRUE)

rna_seq_t_cells <- dplyr::filter(rna_seq_t_cells, isac_id %in% intersect_id)

rna_seq_t_cells <- rna_seq_t_cells %>%
  distinct(isac_id, CDR3aa_trunct, .keep_all = TRUE)

filtered_rna_seq <- rna_seq_t_cells
fraction_bar_table <- rna_seq_t_cells

third_letters <- substr(filtered_rna_seq$V, 3, 3) 
filtered_rna_seq_A <- filtered_rna_seq[third_letters == "A", ]
filtered_rna_seq_G <- filtered_rna_seq[third_letters == "G", ]
filtered_rna_seq_B <- filtered_rna_seq[third_letters == "B", ]
filtered_rna_seq_D <- filtered_rna_seq[third_letters == "D", ]

tcr_seq_AB <- tcr_seq[substr(tcr_seq$TCR_Alpha_Gamma_V_gene_Dominant, 3, 3) == "A", ]
tcr_seq_GD <- tcr_seq[substr(tcr_seq$TCR_Alpha_Gamma_V_gene_Dominant, 3, 3) == "G", ]

result_A <- inner_join(tcr_seq_AB, filtered_rna_seq_A,  
                       by = "isac_id") %>%
  filter(TCR_Alpha_Gamma_CDR3_Translation_Dominant == CDR3aa_trunct) %>% 
  dplyr::select(isac_id, TCR_Alpha_Gamma_CDR3_Translation_Dominant,CDR3aa_trunct, CDR3aa, Sample_Name, CDR3aa_concat,sample_name, TCR_Alpha_Gamma_Molecule_Count, count, cell_type)

result_B <- inner_join(tcr_seq_AB, filtered_rna_seq_B, 
                       by = "isac_id") %>%
  filter(TCR_Beta_Delta_CDR3_Translation_Dominant == CDR3aa_trunct) %>% 
  dplyr::select(isac_id, TCR_Beta_Delta_CDR3_Translation_Dominant,CDR3aa_trunct, CDR3aa, Sample_Name, CDR3aa_concat,sample_name, TCR_Beta_Delta_Molecule_Count, count, cell_type)

result_D <- inner_join(tcr_seq_GD, filtered_rna_seq_D, 
                       by = "isac_id") %>%
  filter(TCR_Beta_Delta_CDR3_Translation_Dominant == CDR3aa_trunct) %>% 
  dplyr::select(isac_id, TCR_Beta_Delta_CDR3_Translation_Dominant,CDR3aa_trunct, CDR3aa, Sample_Name, CDR3aa_concat,sample_name, TCR_Beta_Delta_Molecule_Count, count, cell_type)

result_G <- inner_join(tcr_seq_GD, filtered_rna_seq_G, 
                       by = "isac_id") %>%
  filter(TCR_Alpha_Gamma_CDR3_Translation_Dominant == CDR3aa_trunct) %>% 
  dplyr::select(isac_id, TCR_Alpha_Gamma_CDR3_Translation_Dominant,CDR3aa_trunct, CDR3aa, Sample_Name, CDR3aa_concat,sample_name, TCR_Alpha_Gamma_Molecule_Count, count, cell_type)

colnames(result_A)[c(2,8)] <- c("CDR33a_tcr_table", "CDR33a_tcr_table_count")
colnames(result_B)[c(2,8)] <- c("CDR33a_tcr_table", "CDR33a_tcr_table_count")
colnames(result_D)[c(2,8)] <- c("CDR33a_tcr_table", "CDR33a_tcr_table_count")
colnames(result_G)[c(2,8)] <- c("CDR33a_tcr_table", "CDR33a_tcr_table_count")

result_merge <- rbind(result_A, result_B, result_D, result_G)


table_plot <- data.frame()
for (i in intersect_id) {
  patient_id <- i
  sub_table <- dplyr::filter(result_merge, isac_id == patient_id)
  sampleid <- rep(patient_id, 5) 
  cell_type <- c("CD4","CD8","gdT","Others", "No-overlap")
  values <- c(sum(sub_table[sub_table$cell_type == "CD4T",]$count), sum(sub_table[sub_table$cell_type == "CD8T",]$count),
              sum(sub_table[sub_table$cell_type == "gdT",]$count), sum(sub_table[sub_table$cell_type == "others",]$count),sum(fraction_bar_table[fraction_bar_table$isac_id == patient_id,]$count) - sum(sub_table$count))
  tmp_table <- data.frame(sampleid, cell_type, values)
  tmp_table$cell_type <- factor(tmp_table$cell_type, levels = rev(cell_type))
  if (sum(tmp_table$values) == sum(fraction_bar_table[fraction_bar_table$isac_id == patient_id,]$count)) {
    print(paste(patient_id, "YES"))
  }else{
    print(paste(patient_id, "NO"))
  }
  table_plot <- rbind(table_plot, tmp_table)
}


ggplot(table_plot, aes(fill=cell_type, y=values, x=factor(sampleid))) + 
  geom_bar(position="fill", stat="identity", width=0.7) +  # Adjust bar width
  scale_fill_manual(
    values = rev(c("#D39234", "#61A2DE", "#428D64", "#BA5222", "#878787")), 
    labels = rev(c("OL with blood CD4", "OL with blood CD8", "OL with blood gdT", "Others", "Not-Overlapped"))
  ) +  # Custom legend labels
  labs(y = "Freq. in tumor (%)", x = "Patient ID") + 
  theme_minimal() +  # Minimal theme
  theme(panel.grid = element_blank(), legend.title = element_blank()) +  # Removes grid lines and legend title
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Optional: rotate x-axis labels if needed



#Fig 5c----------------- 


rna_seq_blood <- read.delim("Data/blood_tissue_tcr_from_trust_merge.tsv", sep = "\t") 
colnames(rna_seq_blood)[2] <- "count"
rna_seq_blood_t_cells <- rna_seq_blood[!grepl("^IG", rna_seq_blood$V), ]
rna_seq_blood_t_cells <- rna_seq_blood_t_cells[grepl("^TR", rna_seq_blood_t_cells$V), ]
rna_seq_blood_t_cells <- dplyr::filter(rna_seq_blood_t_cells, CDR3aa != "out_of_frame")

rna_seq_blood_t_cells$count <- as.numeric(rna_seq_blood_t_cells$count)
rna_seq_blood_t_cells1 <- rna_seq_blood_t_cells %>%
  group_by(Filename, CDR3aa) %>%
  summarize(
    count = sum(count, na.rm = TRUE),  
    across(everything(), ~ first(na.omit(.)), .names = "first_{.col}"),
    .groups = 'drop'
  )


rna_seq_tissue <- read.delim("Data/tumor_tissue_tcr_from_trust_merge.tsv", sep = "\t") 
colnames(rna_seq_tissue)[2] <- "count"
rna_seq_tissue_t_cells <- rna_seq_tissue[!grepl("^IG", rna_seq_tissue$V), ]
rna_seq_tissue_t_cells <- rna_seq_tissue_t_cells[grepl("^TR", rna_seq_tissue_t_cells$V), ]
rna_seq_tissue_t_cells <- dplyr::filter(rna_seq_tissue_t_cells, CDR3aa != "out_of_frame")
rna_seq_tissue_t_cells$isac_id <- word(rna_seq_tissue_t_cells$sample_name, end = 1, sep = fixed("_"))

rna_seq_tissue_t_cells$count <- as.numeric(rna_seq_tissue_t_cells$count)
rna_seq_tissue_t_cells1 <- rna_seq_tissue_t_cells %>%
  group_by(sample_name, CDR3aa) %>%
  summarize(
    count = sum(count, na.rm = TRUE),  
    across(everything(), ~ first(na.omit(.)), .names = "first_{.col}"),
    .groups = 'drop'
  )

names_intersection_seq <- word(unique(rna_seq_blood_t_cells1$Filename), end = 2, sep = fixed("_"))
blood_ngi_data <- read.csv("Data/ISAC_baseline_NGI_id_20241001.csv", sep = ",")
blood_ngi_data$isac_id <- sprintf("ISAC%03d", as.numeric(sub("ISAC", "", blood_ngi_data$study_id)))
blood_ngi_data_select <- dplyr::filter(blood_ngi_data, NGI_id %in% names_intersection_seq)
names_intersection_isac <- sprintf("ISAC%03d", as.numeric(sub("ISAC", "", blood_ngi_data_select$study_id)))
rna_seq_tissue_t_cells1_inuse <- rna_seq_tissue_t_cells1 %>%
  filter(sapply(sample_name, function(x) any(grepl(paste(names_intersection_isac, collapse = "|"), x))))

rna_seq_tissue_t_cells1_inuse$index <- paste(rna_seq_tissue_t_cells1_inuse$first_isac_id, rna_seq_tissue_t_cells1_inuse$CDR3aa,sep = "_")


blood_ngi_data <- blood_ngi_data %>% dplyr::select(isac_id, NGI_id)
rna_seq_blood_t_cells1$seq_id <- word(rna_seq_blood_t_cells1$Filename, end = 2, sep = fixed("_"))
rna_seq_blood_t_cells1_inuse <- merge(rna_seq_blood_t_cells1, blood_ngi_data, by.x = "seq_id", by.y = "NGI_id")
rna_seq_blood_t_cells1_inuse$index <- paste(rna_seq_blood_t_cells1_inuse$isac_id, rna_seq_blood_t_cells1_inuse$CDR3aa, sep = "_")


third_letters_blood <- substr(rna_seq_blood_t_cells1_inuse$first_V, 3, 3)
third_letters_tissue <- substr(rna_seq_tissue_t_cells1_inuse$first_V, 3, 3)



tcr_chain <- "A" ## change to A B D G

filtered_rna_seq_tissue <- rna_seq_tissue_t_cells1_inuse[third_letters_tissue == tcr_chain, ]
filtered_rna_seq_blood <- rna_seq_blood_t_cells1_inuse[third_letters_blood == tcr_chain, ]

unique_cdr3aa <- union(filtered_rna_seq_tissue$index, filtered_rna_seq_blood$index)

filtered_rna_seq_tissue1 <-  filtered_rna_seq_tissue %>% dplyr::select(first_isac_id, CDR3aa, count, index)
colnames(filtered_rna_seq_tissue1) <- c("isac_id_tissue", "CDR3aa_tissue", "count_tissue", "index")
filtered_rna_seq_blood1 <- filtered_rna_seq_blood %>% dplyr::select(isac_id, CDR3aa, count, index)
colnames(filtered_rna_seq_blood1) <- c("isac_id_blood", "CDR3aa_blood", "count_blood", "index")

tissue_blood_merge <- merge(filtered_rna_seq_tissue1, filtered_rna_seq_blood1, by = "index", all = T) %>% 
  dplyr::select(index, count_tissue, count_blood) %>%
  mutate_all(~replace(., is.na(.), 0))

tissue_blood_merge$isac_id <- word(tissue_blood_merge$index, end = 1, sep = fixed("_"))
tissue_blood_merge$tcr <- word(tissue_blood_merge$index, start = 2, end = -1,sep = fixed("_"))
tissue_blood_merge <- tissue_blood_merge[,-1] 

meta_table <- read.csv2("~/OneDrive - Karolinska Institutet/personal_files/PhD/Rawdata/wgs/metadata_from_qi_20240903.csv", header = T)
meta_table$isac_id <- paste("ISAC", formatC(meta_table$study_id, width = 3, format = "d", flag = "0"), sep = "")
meta_table <- meta_table %>% dplyr::select(isac_id, tumor_level2)

tissue_blood_merge <- merge(tissue_blood_merge, meta_table, by = "isac_id", all.x = T)

library(dplyr)
library(ggplot2)

ggplot(data = tissue_blood_merge, aes(x = log(count_blood), y = log2(count_tissue))) + 
  geom_jitter(size = 3) +  
  theme_bw() + 
  xlab("log2(TCR reads in Blood)") + 
  ylab("log2(TCR reads in Tissues)") + 
  labs(title = paste(tcr_chain, " chain", sep = ""), fill = "Tumor types of samples with TCRs existing in blood and tissues")  



