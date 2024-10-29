library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
source("Scripts/functions/clone_expansion_plots.R")

clone_id_map <- read_csv("Data/clone_id_map.csv.gz")
data_scTCR <- read_csv("Data/scTCR_data_merge_paired_chain_from_tan_20240415.csv.gz") %>%
    mutate(seurat_clusters = as.character(seurat_clusters))

# construct sample list
sample_names <- data_scTCR %>%
    select(Sample_Name) %>%
    distinct() %>%
    filter(grepl("1$", Sample_Name)) %>%
    pull()

# donut plot
fig_dir <- file.path("Figures", "Figure5")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
for (sub_name in c("CD4T", "CD8T", "gdT")) {
    clone_exp_sub <- get_clone_exp_sub(data_scTCR, sub_name) %>%
        inner_join(clone_id_map, by = "CDR3_concat")
    lapply(sample_names, function(x) {
        clone_expansion_donut(x, clone_exp_sub)
    }) %>%
        ggarrange(plotlist = ., ncol = 3, nrow = 5) %>%
        ggexport(
            filename = file.path(fig_dir, paste0("clone_expansion_", sub_name, ".pdf")),
            width = 10, height = 20
        )
}
