# Script to generate supplementary figure 1:

# Load packages
library(pheatmap)

# Read cell cluster median marker expression data
marker <- read.csv("Data/cytof_MFI_zscore_cluster.csv")

# Supplementary figure 1 -----
pheatmap(marker[,4:35],
         color = paletteer::paletteer_d("rcartocolor::PurpOr", n = 100, type = 'continuous'), 
         scale = "column",
         labels_row = paste(marker$lineage, marker$cluster, sep = "_"),
         display_numbers = FALSE,
         angle_col = 45,
         breaks = seq(0, 2, 2/90),
         main = "Median marker expression per cluster (z-score transformed)")
