### Load libs required
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
### Load Data
All_merge <- readRDS("~/all_merge.rds") # rds data created in DataPrep script
All_Luminal_obj <- subset(All_merge, idents = c("0", "8", "2", "12", "1", "11", "6", "14", "5"))

# Generate Dim Plot
All_Luminal_Dim_plot <- DimPlot(All_Luminal_obj, reduction = "umap")

# Save Dim Plot
ggsave("~/ImageExports/All_Luminal_Dim_plot.eps", plot = All_Luminal_Dim_plot)
