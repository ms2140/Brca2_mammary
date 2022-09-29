### Load libs required
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
### Load Data
All_merge <- readRDS("~/all_merge.rds") # rds data created in DataPrep script
All_Luminal_obj <- subset(All_merge, idents = c("0", "8", "2", "12", "1", "11", "6", "14", "5"))
All_Luminal_P0_obj <- subset(All_Luminal_obj, subset = (Condition == "WT_HU_P0" | Condition == "B2_HU_P0" | Condition == "WT_P0" | Condition == "B2_P0"))
All_Luminal_P4_obj <- subset(All_Luminal_obj, subset = (Condition == "WT_HU_P4" | Condition == "B2_HU_P4" | Condition == "WT_P4" | Condition == "B2_P4"))

# Find markers
All_Luminal_P4_markers_all_obj <- FindAllMarkers(object = All_Luminal_P4_obj)
All_Luminal_P4_markers_top20 <- All_Luminal_P4_markers_all_obj %>% group_by(cluster) %>% top_n(20, avg_log2FC)

# Generate HeatMap
All_Luminal_P4_Heatmap <- DoHeatmap(object = All_Luminal_P4_obj, features = All_Luminal_P4_markers_top20$gene)  + scale_fill_gradientn(colors = c("blue", "white", "red"))

# Save Dim Plot
ggsave("~/ImageExports/All_Luminal_P4_Heatmap.eps", plot = All_Luminal_P4_Heatmap)
ggsave("~/ImageExports/All_Luminal_P4_Heatmap.jpeg", plot = All_Luminal_P4_Heatmap)
