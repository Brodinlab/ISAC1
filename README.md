# ISAC1
Title：Systems-level immunomonitoring in children with solid tumors

# Description
This repo is used for reproducing the figures in the paper (~~add link here~~). 

# Contents
- ```Data/``` contains all the processed tables needed for reproducing the figures in the paper.
- ```Scripts/``` contains the codes for reproducing the figures in the paper.

# Data availablility
### Single cell RNA sequencing
Data available in ~~add link here~~ with .h5Seurat format. To regenerate the processed tables in ```Data/```, put the data uder ```single cell data/``` and use the ```preparation.R``` scripts within the folders. 
### Single cell TCR sequencing
Processed single cell table is ```Data/scTCR_data_merge_paired_chain_from_tan_20240415.csv.gz```. scTCR and scRNA tables can be linked via the column *unique_index*.
The raw sequencing files are considered sensitive and are available upon reasonable request. 
### Adult data from public datasets 
Original data are from ~~xxxxxx~~. The data of peripheral memory CD8+ T cells used in this paper are available in ~~add link here~~. To regenerate the processed tables in ```Data/```, put the data uder ```single cell data/``` and use the ```preparation.R``` scripts within the folders. 
### CyTOF
Raw FCS files are avaiable in http://flowrepository.org/id/FR-FCM-Z8C8.
### Plasma protein (Olink)
Processed Normalised Protein eXpression (NPX) value of plasma proteins is ```Data/olink_npx_baseline.csv```.
### Meta data of patients
~~Which file?~~
### Tissue sequencing ~~or~~
~~not avaiable? reasons?~~
### ~~Any more?~~ 

# Dependencies
```
R version 4.2.3
dplyr 1.1.1, readr 2.1.4, ggplot2 3.4.1, ggpubr 0.6.0, ggalluvial 0.12.5, ggrepel 0.9.3, wesanderson 0.3.6,
ComplexHeatmap 2.14.0, EnhancedVolcano 1.16.0, Seurat 4.3.0, SeuratDisk 0.0.0.9020, rstatix 0.7.2,
parallel 4.2.3, igraph 1.4.1
```
~~add more~~

# Contact
For any questions regarding data analyses, please contact *petter.brodin@ki.se*.
