### Load libs required
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
### Load Data
All_merge <- readRDS("~/all_merge.rds") # rds data created in DataPrep script
# Set colours
plot_cols <- c('Basal'='#A6CDE3','HR+ Luminal'='#FA9999','HR- Luminal'='#B1DE8A','Basal cycling'='#1F78B4','HR- Luminal cycling'='#33A02B')
# Generate Dot Plot
All_merge_Dot_Plot <- DotPlot(object = All_merge, features = c("Lalba", "Mfge8", "Cited2", "Aqp5",  "Cd14", "Krt18", "Lmo4", "Elf5", "Kit", "Cebpb",  "Stat3", "Birc5", "Ccnb2", "Cdc25c", "Cdkn2d", "Ccnd1", "Krt14", "Acta2", "Id4", "Runx1", "Pdpn", "Trp63", "Vim", "Ctnnb1", "Areg", "Prlr", "Wnt4", "Foxa1", "Pgr", "Ar", "Erbb3", "Epcam")) + theme(axis.text.x = element_text(angle = 90))
# Save Dim Plot
ggsave("~/ImageExports/All_merge_Dot_Plot", plot = All_merge_Dot_Plot)
