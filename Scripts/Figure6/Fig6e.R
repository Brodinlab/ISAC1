# Script to generate figure 6e:

# Load packages
library(tidyverse)
library(ggrepel)

# Read protein data
protein <- read.csv("Data/olink_npx_isac002.csv", row.names = 1, check.names = F) 

# Calculate log2 fold change
df <- protein %>%
  rownames_to_column("protein") %>%
  mutate(log2fc = `9.month` - `5.month`)

# Figure 6e -----

ggplot() +
  geom_vline(linetype = "dotted", xintercept = 0, color = "grey") +
  geom_point(data = df, aes(x = log2fc, y = reorder(protein, log2fc)), color = "grey", size = 4, shape = 1) +
  geom_point(data = subset(df, log2fc > 1), aes(x = log2fc, y = reorder(protein, log2fc)), color = "#16AB9B", size = 4) +
  geom_point(data = subset(df, log2fc < -1), aes(x = log2fc, y = reorder(protein, log2fc)), color = "#F3C142", size = 4) +
  geom_text_repel(data = subset(df, abs(log2fc) > 1), aes(x = log2fc, y = reorder(protein, log2fc), label = protein), box.padding = 0.5, max.overlaps = 100) +
  labs(x = "Log2fc(9 month/5 month)", y = "Proteins") +
  theme(text = element_text(color = "black"), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
