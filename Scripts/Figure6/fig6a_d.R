library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(EnhancedVolcano)

source("Scripts/functions/clone_expansion_plots.R")

clone_id_map <- read_csv("Data/clone_id_map.csv.gz")
data_scTCR <- read_csv("Data/scTCR_data_merge_paired_chain_from_tan_20240415.csv.gz") %>%
    mutate(seurat_clusters = as.character(seurat_clusters))
virus_specific_id <- read_csv(file = "Data/virus_specific_clone_id.csv") %>% pull()
fig_dir <- file.path("Figures", "Figure6")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
sub_name <- "CD8T"
patient_name <- "ISAC02"

# Fig 6a
sample_names <- data_scTCR %>%
    dplyr::select(Sample_Name) %>%
    distinct() %>%
    filter(grepl(patient_name, Sample_Name)) %>%
    pull() %>%
    gtools::mixedsort()
clone_exp_sub <- get_clone_exp_sub(data_scTCR, sub_name) %>%
    inner_join(clone_id_map, by = "CDR3_concat")
lapply(sample_names, function(x) {
    clone_expansion_donut(x, clone_exp_sub)
}) %>%
    ggarrange(plotlist = ., ncol = 3, nrow = 5) %>%
    ggexport(
        filename = file.path(fig_dir, paste0(patient_name, "_clone_expansion_", sub_name, ".pdf")),
        width = 10, height = 20
    )

# Fig 6b
clone_exp_sub <- get_clone_exp_sub(data_scTCR, sub_name) %>%
    inner_join(clone_id_map, by = "CDR3_concat")
clone_expansion_alluvium(patient_name, clone_exp_sub, clone_label = FALSE) +
    labs(title = paste0(patient_name, "_", sub_name)) + guides(fill = "none", color = "none") +
    theme_bw()
ggsave(file.path(fig_dir, paste0(patient_name, "_clone_alluvium_", sub_name, ".pdf")))

# Fig 6c
trend_determination_plot(patient_name, clone_exp_sub, top_mod = "union")
increase_clone <- c(48795, 37805, 2890, 45479, 53017, 33890)
decrease_clone <- c(
    7918, 30260, 51808, 32056, 4293, 41385, 18961, 50670, 41552,
    47305, 20424, 18818, 47227, 32052, 43590, 32214, 55081
)
df <- clone_exp_alluvium_preparation(patient_name, clone_exp_sub, top_mod = "union")
df <- df %>% mutate(
    trend = case_when(
        clone_id %in% virus_specific_id ~ "virus_specific",
        clone_id %in% increase_clone ~ "increase",
        clone_id %in% decrease_clone ~ "decrease",
        TRUE ~ "other"
    )
)

ggplot(df, aes(
    x = time_point, stratum = clone_id, alluvium = clone_id,
    # y= clone_ratio,
    y = relative_clone_ratio,
    fill = trend, color = clone_id
)) +
    geom_alluvium() +
    geom_stratum(size = 0.1) +
    guides(color = "none") + # hide the legend of 'color'
    # scale_fill_manual(values = wes_palette("Rushmore1", length(unique(df$trend)), type = "continuous"))+
    # theme(legend.position = "none")+
    labs(title = patient_name) +
    theme_bw()
ggsave(file.path(fig_dir, paste0(patient_name, "_clone_alluvium_", sub_name, "_color_by_trend.pdf")))

# Fig 6d

res <- read_csv("Data/Volcano compare trends in ISAC02.csv")
rownames(res) <- res$gene
EnhancedVolcano(res,
    x = "avg_log2FC", y = "p_val_adj", lab = rownames(res),
    pCutoff = 1e-04, FCcutoff = 1, drawConnectors = TRUE,
    selectLab = res %>% filter(abs(avg_log2FC) > 2 | p_val < 1e-04) %>% rownames(),
    title = paste0(sub_name, " of ", patient_name),
    subtitle = "increase vs decrease", xlim = c(-1.5, 1.5)
)
ggsave(file.path(fig_dir, paste0(patient_name, "_", sub_name, "_expanding vs contracting volcano.pdf")), width = 10, height = 10)
