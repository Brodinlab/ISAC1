library(dplyr)
library(readr)
library(Seurat)
library(SeuratDisk)

obj <- LoadH5Seurat("single cell data/seurat_results_update2.h5Seurat")
virus_specific_id <- read_csv(file = "Data/virus_specific_clone_id.csv") %>% pull()
increase_clone <- c(48795, 37805, 2890, 45479, 53017, 33890)
decrease_clone <- c(
    7918, 30260, 51808, 32056, 4293, 41385, 18961, 50670, 41552,
    47305, 20424, 18818, 47227, 32052, 43590, 32214, 55081
)
sub_name <- "CD8T"
patient_name <- "ISAC02"
obj_sub <- subset(obj, subset = patient_id == patient_name &
    clone_id %in% union(increase_clone, decrease_clone) &
    !clone_id %in% virus_specific_id &
    cell_type == sub_name)
tmp <- obj_sub[[]] %>%
    mutate(
        trend = case_when(
            clone_id %in% virus_specific_id ~ "virus_specific",
            clone_id %in% increase_clone ~ "increase",
            clone_id %in% decrease_clone ~ "decrease",
            TRUE ~ "other"
        )
    ) %>%
    select(trend)
rownames(tmp) <- obj_sub[[]]$unique_index
obj_sub$trend <- tmp
Idents(obj_sub) <- "trend"
res <- FindMarkers(obj_sub,
    ident.1 = "increase", ident.2 = "decrease",
    logfc.threshold = 0
)
res$gene <- rownames(res)
write_csv(res, "Data/Volcano compare trends in ISAC02.csv")
