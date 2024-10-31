# ISAC1
**Title**: Systems-level Immunomonitoring in Children with Solid Tumors

## Description
This repository contains scripts and data needed to reproduce the figures in the paper [*Systems-level Immunomonitoring in Children with Solid Tumors*](~~add link here~~).

## Repository Structure
- **`Data/`**: Contains all processed tables required for reproducing the figures in the paper.
- **`Scripts/`**: Contains code to regenerate the figures presented in the paper.

## Data Availability

### Single-Cell RNA Sequencing
Data is available in [~~add link here~~] in `.h5Seurat` format. To regenerate the processed tables in the `Data/` folder, place the data in `single cell data/` and run the `preparation.R` script located in the folder.

### Single-Cell TCR Sequencing
The processed single-cell TCR sequencing table is available as `Data/scTCR_data_merge_paired_chain_from_tan_20240415.csv.gz`. scTCR and scRNA tables can be linked via the column *unique_index*. The raw sequencing files are sensitive and available upon reasonable request.

### Adult Data from Public Datasets
Original data from adult datasets is from [~~xxxxxx~~]. The data for peripheral memory CD8+ T cells used in this paper can be accessed at [~~add link here~~]. To regenerate the processed tables in the `Data/` folder, place the data in `single cell data/` and run the `preparation.R` scripts.

### CyTOF
Raw FCS files are available at [FlowRepository (ID: FR-FCM-Z8C8)](http://flowrepository.org/id/FR-FCM-Z8C8).

### Plasma Protein (Olink)
Processed Normalized Protein eXpression (NPX) data is available as `Data/olink_npx_baseline.csv`.

### Metadata of Patients
Patient metadata can be found in the file `Data/metadata.csv`.

### Tissue Sequencing
~~Data availability status: not available? (Reasons if needed)~~

## Dependencies
To reproduce the figures, the following R version and packages are required:

- **R version**: 4.2.3
- **Packages**: dplyr 1.1.2, readr 2.1.4, ggplot2 3.4.4, ggpubr 0.6.0, ggalluvial 0.12.5, ggrepel 0.9.3, wesanderson 0.3.6, ComplexHeatmap 2.14.0, pheatmap 1.0.12, EnhancedVolcano 1.16.0, Seurat 4.3.0, SeuratDisk 0.0.0.9020, rstatix 0.7.2, parallel 4.2.3, igraph 1.4.1, pROC 1.18.5, corrplot 0.92, caret 6.0-94, lattice 0.22-5, relaimpo 2.2-7, paletteer 1.6.0, robCompositions 2.4.1, data.table 1.14.8, pls 2.8-3, ggridges 0.5.6, webr 0.1.5, lubridate 1.9.3, forcats 1.0.0, stringr 1.5.1, purrr 1.0.2, tidyr 1.3.0, tibble 3.2.1, tidyverse 2.0.0, mitools 2.4, survey 4.2-1, survival 3.5-7, Matrix 1.5-1, boot 1.3-28.1, MASS 7.3-60.

To install the required packages, run the following command in R:

```r
install.packages(c("dplyr", "readr", "ggplot2", "ggpubr", "ggalluvial", "ggrepel", 
                   "wesanderson", "ComplexHeatmap", "pheatmap", "EnhancedVolcano", 
                   "Seurat", "SeuratDisk", "rstatix", "parallel", "igraph", "pROC",
                   "corrplot", "caret", "lattice", "relimpo", "paletteer", "robCompositions",
                   "data.table", "pls", "ggridges", "webr", "lubridate", "forcats", 
                   "stringr", "purrr", "tidyr", "tibble", "tidyverse", "mitools", 
                   "survey", "survival", "Matrix", "boot", "MASS"))


# Contact
For any issues or questions regarding data analyses, please contact *petter.brodin@ki.se*.
