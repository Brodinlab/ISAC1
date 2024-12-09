library(dplyr)
library(readr)
library(Seurat)
library(SeuratDisk)
library(rstatix)
library(parallel)

# resampling DE analysis for Fig 7
## ranked list in Wu -----
# cannot import to Seurat because the data is not raw counts
df_wu <- read_csv("single cell data/adult data from public/blood mCD8T from dataset_Wu.csv.gz")
panel_wu <- colnames(df_wu)[1:3734]

calculate_DE_tmp <- function(df, panel) {
    tmp <- df %>%
        group_by(expansion) %>%
        summarise(across(all_of(panel), mean))
    od <- tmp$expansion
    tmp <- t(tmp[, panel])
    colnames(tmp) <- od
    res_fc <- tmp %>%
        as.data.frame() %>%
        mutate(logFC = expanded - non_expanded)
    res_fc$gene <- rownames(res_fc)
    tmp <- df %>% tidyr::pivot_longer(values_to = "level", names_to = "gene", cols = all_of(panel))
    res_p <- tmp %>%
        group_by(gene) %>%
        t_test(level ~ expansion) %>%
        adjust_pvalue(method = "fdr") %>%
        add_significance("p.adj")
    res_rank <- left_join(res_fc, res_p, by = "gene") %>%
        mutate(rank_metric = logFC * -log10(p.adj)) %>%
        arrange(desc(rank_metric)) %>%
        rename(p_val_adj = p.adj)
    return(res_rank)
}
# resampling (bootstrapping)
rank_wu_list <- mclapply(paste("resample", 1:100), function(x) {
    print(x)
    df_wu_re <- df_wu[sample(rownames(df_wu), size = 0.9 * length(rownames(df_wu)), replace = FALSE), ]
    rank_wu_re <- calculate_DE_tmp(df_wu_re, panel_wu)
    rank_wu_re$resample <- x
    return(rank_wu_re)
}, mc.cores = 4)
common_genes <- lapply(rank_wu_list, function(x) x$gene) %>% Reduce(f = intersect)
rank_wu_mean <- lapply(rank_wu_list, function(x) {
    x %>%
        filter(gene %in% common_genes) %>%
        select(gene, resample, rank_metric, logFC)
}) %>%
    do.call(what = rbind) %>%
    group_by(gene) %>%
    summarise(
        avg_rank_metric = mean(rank_metric),
        avg_log2FC = mean(logFC)
    ) %>%
    ungroup() %>%
    arrange(desc(avg_rank_metric))
write_csv(rank_wu_mean, "Data/DE_resample_results/adults_Wu.csv")

## ranked list in Guo -----
obj_guo_memt <- LoadH5Seurat("single cell data/adult data from public/dataset_Guo_memt.h5Seurat")

# resampling (bootstrapping)
rank_guo_list <- mclapply(paste("resample", 1:100), function(x) {
    obj_guo_memt_re <- obj_guo_memt[, sample(colnames(obj_guo_memt), size = 0.9 * length(colnames(obj_guo_memt)), replace = FALSE)]
    Idents(obj_guo_memt_re) <- "Clone.Status"
    rank_guo_re <- FindMarkers(obj_guo_memt_re,
        ident.1 = "Clonal", ident.2 = "NoClonal",
        logfc.threshold = 0, min.pct = 0.5
    ) %>%
        mutate(rank_metric = avg_log2FC * -log10(p_val_adj)) %>%
        arrange(desc(rank_metric))
    rank_guo_re$gene <- rownames(rank_guo_re)
    rank_guo_re$resample <- x
    return(rank_guo_re)
}, mc.cores = 4)
common_genes <- lapply(rank_guo_list, function(x) x$gene) %>% Reduce(f = intersect)
rank_guo_mean <- lapply(rank_guo_list, function(x) {
    x %>%
        filter(gene %in% common_genes) %>%
        select(gene, resample, rank_metric, avg_log2FC)
}) %>%
    do.call(what = rbind) %>%
    group_by(gene) %>%
    summarise(
        avg_rank_metric = mean(rank_metric),
        avg_log2FC = mean(avg_log2FC)
    ) %>%
    ungroup() %>%
    arrange(desc(avg_rank_metric))
write_csv(rank_guo_mean, "Data/DE_resample_results/adults_Guo.csv")

## ranked list in Zhang --------
obj_zhang_memt <- LoadH5Seurat("single cell data/adult data from public/dataset_Zhang_memt.h5Seurat")

# resampling (bootstrapping)
rank_zhang_list <- mclapply(paste("resample", 1:100), function(x) {
    obj_zhang_memt_re <- obj_zhang_memt[, sample(colnames(obj_zhang_memt), size = 0.9 * length(colnames(obj_zhang_memt)), replace = FALSE)]
    Idents(obj_zhang_memt_re) <- "Clonal.status"
    rank_zhang_re <- FindMarkers(obj_zhang_memt_re,
        ident.1 = "Clonal", ident.2 = "NoClonal",
        logfc.threshold = 0, min.pct = 0.5
    ) %>%
        mutate(rank_metric = avg_log2FC * -log10(p_val_adj)) %>%
        arrange(desc(rank_metric))
    rank_zhang_re$gene <- rownames(rank_zhang_re)
    rank_zhang_re$resample <- x
    return(rank_zhang_re)
}, mc.cores = 4)
common_genes <- lapply(rank_zhang_list, function(x) x$gene) %>% Reduce(f = intersect)
rank_zhang_mean <- lapply(rank_zhang_list, function(x) {
    x %>%
        filter(gene %in% common_genes) %>%
        select(gene, resample, rank_metric, avg_log2FC)
}) %>%
    do.call(what = rbind) %>%
    group_by(gene) %>%
    summarise(
        avg_rank_metric = mean(rank_metric),
        avg_log2FC = mean(avg_log2FC)
    ) %>%
    ungroup() %>%
    arrange(desc(avg_rank_metric))
write_csv(rank_zhang_mean, "Data/DE_resample_results/adults_Zhang.csv")

## ranked list in ISAC -----------
# exclude virus specific
obj <- LoadH5Seurat("single cell data/seurat_results_update2.h5Seurat") %>%
    subset(subset = cell_type == "CD8T" &
        virus_specific == "non_virus_specific" &
        TCR_Paired_Chains == 1)
for (i in unique(obj$tumor_type)) {
    obj_sub <- subset(obj,
        subset = tumor_type == i
    )
    # resampling (bootstrapping)
    rank_sub_list <- mclapply(paste("resample", 1:100), function(x) {
        obj_sub_re <- obj_sub[, sample(colnames(obj_sub), size = 0.9 * length(colnames(obj_sub)), replace = FALSE)]
        Idents(obj_sub_re) <- "expansion"
        rank_sub_re <- FindMarkers(obj_sub_re,
            ident.1 = "expanded", ident.2 = "non-expanded",
            logfc.threshold = 0, min.pct = 0
        ) %>%
            mutate(rank_metric = avg_log2FC * -log10(p_val_adj)) %>%
            arrange(desc(rank_metric))
        rank_sub_re$gene <- rownames(rank_sub_re)
        rank_sub_re$resample <- x
        return(rank_sub_re)
    }, mc.cores = 4)
    common_genes <- lapply(rank_sub_list, function(x) x$gene) %>% Reduce(f = intersect)
    rank_sub_mean <- lapply(rank_sub_list, function(x) {
        x %>%
            filter(gene %in% common_genes) %>%
            select(gene, resample, rank_metric, avg_log2FC)
    }) %>%
        do.call(what = rbind) %>%
        group_by(gene) %>%
        summarise(
            avg_rank_metric = mean(rank_metric),
            avg_log2FC = mean(avg_log2FC)
        ) %>%
        ungroup() %>%
        arrange(desc(avg_rank_metric))
    write_csv(rank_sub_mean, paste0("Data/DE_resample_results/ISAC_", i, ".csv"))
}

# calculate expansion ratio
df <- rbind(
    df_wu %>%
        group_by(patient, expansion, tumor_type) %>%
        tally() %>%
        group_by(patient) %>%
        mutate(sum = sum(n), fraction = n / sum),
    obj_guo_memt[[]] %>%
        group_by(Patient, Clone.Status) %>%
        tally() %>%
        group_by(Patient) %>%
        mutate(sum = sum(n), fraction = n / sum) %>%
        rename(patient = Patient, expansion = Clone.Status) %>%
        mutate(tumor_type = "NSCLC"),
    obj_zhang_memt[[]] %>%
        group_by(Patient_ID, Clonal.status) %>%
        tally() %>%
        group_by(Patient_ID) %>%
        mutate(sum = sum(n), fraction = n / sum) %>%
        rename(patient = Patient_ID, expansion = Clonal.status) %>%
        mutate(tumor_type = "Colorectal")
) %>%
    filter(expansion %in% c("expanded", "Clonal")) %>%
    mutate(group = "adults")

df <- rbind(
    df,
    obj[[]] %>%
        group_by(patient_id, expansion, tumor_type) %>%
        tally() %>% group_by(patient_id) %>%
        mutate(sum = sum(n), fraction = n / sum, group = "isac") %>%
        rename(patient = patient_id) %>%
        filter(expansion == "expanded")
)
write_csv(df, "Data/adult_isac_expansion_ratio.csv")
