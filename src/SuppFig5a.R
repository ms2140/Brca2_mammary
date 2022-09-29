### Load libs required
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)

### Load Data
All_merge <- readRDS("~/all_merge.rds") # rds data created in DataPrep script
All_Luminal_Cluster_obj <- subset(All_merge, subset = (RNA_clusters == c("0") | RNA_clusters == c("1") | RNA_clusters == c("2") | RNA_clusters == c("5") | RNA_clusters == c("6") | RNA_clusters == c("8") | RNA_clusters == c("11") | RNA_clusters == c("12") | RNA_clusters == c("14") ))
BarChart_B2_Hu_P4_obj <- subset(All_Luminal_Cluster_obj, subset = (Condition == "B2_HU_P4"))
BarChart_B2_P4_obj <- subset(All_Luminal_Cluster_obj, subset = (Condition == "B2_P4"))
BarChart_WT_HU_P4_obj <- subset(All_Luminal_Cluster_obj, subset = (Condition == "WT_HU_P4"))
BarChart_WT_P4_obj <- subset(All_Luminal_Cluster_obj, subset = (Condition == "WT_P4"))

# Prep Data
b2_hu_p4 <- table(Idents(BarChart_B2_Hu_P4_obj), BarChart_B2_Hu_P4_obj$orig.ident)
b2_hu_p4 <- as.data.frame(b2_hu_p4)
b2_hu_p4$Var1 <- as.character(b2_hu_p4$Var1)
b2_hu_p4$Var2 <- "b2_hu_p4"

b2_p4 <- table(Idents(BarChart_B2_P4_obj), BarChart_B2_P4_obj$orig.ident)
b2_p4 <- as.data.frame(b2_p4)
b2_p4$Var1 <- as.character(b2_p4$Var1)
b2_p4$Var2 <- "b2_p4"

wt_hu_p4 <- table(Idents(BarChart_WT_HU_P4_obj), BarChart_WT_HU_P4_obj$orig.ident)
wt_hu_p4 <- as.data.frame(wt_hu_p4)
wt_hu_p4$Var1 <- as.character(wt_hu_p4$Var1)
wt_hu_p4$Var2 <- "wt_hu_p4"

wtp4 <- table(Idents(BarChart_WT_P4_obj), BarChart_WT_P4_obj$orig.ident)
wtp4 <- as.data.frame(wtp4)
wtp4$Var1 <- as.character(wtp4$Var1)
wtp4$Var2 <- "WTP4"

wtp4$id <- "WTP4"
wt_hu_p4$id <- "WTHUP4"
b2_p4$id <- "B2P4"
b2_hu_p4$id <- "B2HUP4"
grouped_data <- rbind(wtp4, wt_hu_p4, b2_p4, b2_hu_p4)

# Generate Plot
luminal_p4_cluster_bar_plot <- ggplot(grouped_data, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) + 
  scale_fill_manual(
    name = "Var1",
    labels = c("0", "1", "11", "12", "14", "2", "5", "6","8"),
    values = c("#F7766D", "#D49201", "#609DFF", "#DB72FC", "#FF61C3", "#94A900", "#0BBA39", "#0CC19F", "#00B9E3")
  )

# Save Plot
ggsave("~/ImageExports/All_Luminal_P4_Bargraph.eps", plot = luminal_p4_cluster_bar_plot)
ggsave("~/ImageExports/All_Luminal_P4_Bargraph.jpeg", plot = luminal_p4_cluster_bar_plot)
