### Load libs required
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
### Load Data
All_merge <- readRDS("~/all_merge.rds") # rds data created in DataPrep script
# Set colours
plot_cols <- c('Basal'='#A6CDE3','HR+ Luminal'='#FA9999','HR- Luminal'='#B1DE8A','Basal cycling'='#1F78B4','HR- Luminal cycling'='#33A02B')
# Generate Dim Plot
All_merge_Dim_plot <- DimPlot(All_merge, reduction = "umap", label = FALSE, pt.size = 0.5, cols = plot_cols)
# Save Dim Plot
ggsave("~/ImageExports/All_merge_Dim_plot.eps", plot = All_merge_Dim_plot)
