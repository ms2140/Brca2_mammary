### Load libs required
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)

### Load Data
All_merge <- readRDS("~/all_merge.rds") # rds data created in DataPrep script
data_L_Alv <- read.csv("~/Saeki_L_Alv.csv")
data_LHor <- read.csv("~/Saeki_LHor.csv")
data_Bas <- read.csv("~/Saeki_Basal.csv")

# Set colours
plot_cols <- c('Basal'='#A6CDE3','HR+ Luminal'='#FA9999','HR- Luminal'='#B1DE8A','Basal cycling'='#1F78B4','HR- Luminal cycling'='#33A02B')
# Generate Full VlnPlots

# Luminal Alveolar
L_Alv <- AddModuleScore(All_merge, features = list(data_L_Alv$Gene),  ctrl = 20, name= "L_Alv")
L_Alv_All_plot <- VlnPlot(L_Alv, features = c("L_Alv1"),  pt.size = 0.0, cols = plot_cols) + geom_boxplot(width=0.15, fill = "white") +  geom_hline(yintercept=0, linetype="dashed", color="black")
L_Alv_All_plot$data$ident <- factor(x = L_Alv_All_plot$data$ident, levels = c("Basal", "Basal cycling", "HR- Luminal", "HR- Luminal cycling", "HR+ Luminal"))

# Luminal Hormone
L_Hor <- AddModuleScore(All_merge, features = list(data_LHor$Gene),  ctrl = 20, name= "L_Hor")
L_Hor_All_plot <- VlnPlot(L_Hor, features = c("L_Hor1"), pt.size = 0.0, cols = plot_cols) + geom_boxplot(width=0.15, fill = "white") +  geom_hline(yintercept=0, linetype="dashed", color="black")
L_Hor_All_plot$data$ident <- factor(x = L_Hor_All_plot$data$ident, levels = c("Basal", "Basal cycling", "HR- Luminal", "HR- Luminal cycling", "HR+ Luminal"))

# Basal
Bas <- AddModuleScore(All_merge, features = list(data_Bas$Gene),  ctrl = 20, name= "Bas")
Bas_All_plot <- VlnPlot(Bas, features = c("Bas1"), pt.size = 0.0, cols = plot_cols) + geom_boxplot(width=0.15, fill = "white") +  geom_hline(yintercept=0, linetype="dashed", color="black")
Bas_All_plot$data$ident <- factor(x = Bas_All_plot$data$ident, levels = c("Basal", "Basal cycling", "HR- Luminal", "HR- Luminal cycling", "HR+ Luminal"))

##### Save vlnPlots
#L_Alv_All_plot
ggsave("~/ImageExports/L_Alv_All_plot.eps", plot = L_Alv_All_plot)
ggsave("~/ImageExports/L_Alv_All_plot.jpeg", plot = L_Alv_All_plot)
#L_Hor_All_plot
ggsave("~/ImageExports/L_Hor_All_plot.eps", plot = L_Hor_All_plot)
ggsave("~/ImageExports/L_Hor_All_plot.jpeg", plot = L_Hor_All_plot)
#Bas_All_plot
ggsave("~/ImageExports/Bas_All_plot.eps", plot = Bas_All_plot)
ggsave("~/ImageExports/Bas_All_plot.jpeg", plot = Bas_All_plot)
