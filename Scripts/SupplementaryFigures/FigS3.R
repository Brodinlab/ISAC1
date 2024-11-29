library(dplyr)
library(readr)
library(igraph)
source("Scripts/functions/GLIPH_related_plots.R")


# Supplementary figure 3a -----

## please find the codes at Scripts/Figure5/fig5a_c.R section Fig 5c


# Supplementary figure 3b -----

data_TRA <- read_csv("Data/GLIPH2_TRA_annotated.csv")
data_TRB <- read_csv("Data/GLIPH2_TRB_annotated.csv")
data_scTCR <- read_csv("Data/scTCR_data_merge_paired_chain_from_tan_20240415.csv.gz")

nt2aa <- data_scTCR %>%
    select(CDR3_concat, CDR3aa_concat) %>%
    distinct()
clone_id_map <- data_scTCR %>%
    select(CDR3_concat, clone_id) %>%
    unique()

clone_exp <- data_scTCR %>%
    group_by(CDR3_concat, Sample_Name) %>%
    summarise(clone_count = n()) %>%
    ungroup() %>%
    inner_join(clone_id_map, by = "CDR3_concat") %>%
    mutate(
        clone_id = as.character(clone_id),
        patient_id = sub("_.*$", "", Sample_Name)
    )

fig_dir <- file.path("Figures", "SupplementaryFigures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# network for baseline samples
pdf(file.path(fig_dir, "FigS4b.pdf"))
plot.gliph(data_TRA %>% filter(condition %in% c("1", "v1")),
    data_TRB %>% filter(condition %in% c("1", "v1")),
    clone_exp %>% filter(grepl("(_v1$)|(_1$)", Sample_Name)),
    nt2aa,
    clone_thre = 1
)
title(main = "baseline samples")
dev.off()
