### Load libs required
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)

### Load Data
All_merge <- readRDS("~/all_merge.rds") # rds data created in DataPrep script

### Subset data
WTP0_obj <- subset(All_merge, subset = Condition == "WT_P0")
WTP4_obj <- subset(All_merge, subset = Condition == "WT_P4")
WT_HU_P0_obj <- subset(All_merge, subset = Condition == "WT_HU_P0")
WT_HU_P4_obj <- subset(All_merge, subset = Condition == "WT_HU_P4")
B2P0_obj <- subset(All_merge, subset = Condition == "B2_P0")
B2P4_obj <- subset(All_merge, subset = Condition == "B2_P4")
B2_HU_P0_obj <- subset(All_merge, subset = Condition == "B2_HU_P0")
B2_HU_P4_obj <- subset(All_merge, subset = Condition == "B2_HU_P4")

# Set colours
plot_cols <- c('Basal'='#A6CDE3','HR+ Luminal'='#FA9999','HR- Luminal'='#B1DE8A','Basal cycling'='#1F78B4','HR- Luminal cycling'='#33A02B')

# Generate Plots
#Wild Type P0
wtp0 <- table(Idents(WTP0_obj), WTP0_obj$orig.ident)
wtp0 <- as.data.frame(wtp0)
wtp0$Var1 <- as.character(wtp0$Var1)
wtp0$Var2 <- "WTP0"
WT_P0_PLOT <- ggplot(wtp0, aes(x = "WT_P0", y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.position="none")

#Wild Type P4
wtp4 <- table(Idents(WTP4_obj), WTP4_obj$orig.ident)
wtp4 <- as.data.frame(wtp4)
wtp4$Var1 <- as.character(wtp4$Var1)
wtp4$Var2 <- "WTP4"
WT_P4_PLOT <- ggplot(wtp4, aes(x = "WT_P4", y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.position="none")

#Wild Type Drug Treated P0
wt_hu_p0 <- table(Idents(WT_HU_P0_obj), WT_HU_P0_obj$orig.ident)
wt_hu_p0 <- as.data.frame(wt_hu_p0)
wt_hu_p0$Var1 <- as.character(wt_hu_p0$Var1)
wt_hu_p0$Var2 <- "wt_hu_p0"
WT_HU_P0_PLOT <- ggplot(wt_hu_p0, aes(x = "WT_HU_P0", y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.position="none")

#Wild Type Drug Treated P4
wt_hu_p4 <- table(Idents(WT_HU_P4_obj), WT_HU_P4_obj$orig.ident)
wt_hu_p4 <- as.data.frame(wt_hu_p4)
wt_hu_p4$Var1 <- as.character(wt_hu_p4$Var1)
wt_hu_p4$Var2 <- "wt_hu_p4"
WT_HU_P4_PLOT <- ggplot(wt_hu_p4, aes(x = "WT_HU_P4", y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.position="none")

#Brca 2 P0
b2_p0 <- table(Idents(B2P0_obj), B2P0_obj$orig.ident)
b2_p0 <- as.data.frame(b2_p0)
b2_p0$Var1 <- as.character(b2_p0$Var1)
b2_p0$Var2 <- "b2_p0"
B2_P0_PLOT <- ggplot(b2_p0, aes(x = "B2_P0", y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.position="none")

#Brca 2 P4
b2_p4 <- table(Idents(B2P4_obj), B2P4_obj$orig.ident)
b2_p4 <- as.data.frame(b2_p4)
b2_p4$Var1 <- as.character(b2_p4$Var1)
b2_p4$Var2 <- "b2_p4"
B2_P4_PLOT <- ggplot(b2_p4, aes(x = "B2_P4", y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.position="none")

#Brca 2 Drug Treated P0
b2_hu_p0 <- table(Idents(B2_HU_P0_obj), B2_HU_P0_obj$orig.ident)
b2_hu_p0 <- as.data.frame(b2_hu_p0)
b2_hu_p0$Var1 <- as.character(b2_hu_p0$Var1)
b2_hu_p0$Var2 <- "b2_hu_p0"
B2_HU_P0_PLOT <- ggplot(b2_hu_p0, aes(x = "B2_HU_P0", y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.position="none")

#Brca 2 Drug Treated P4
b2_hu_p4 <- table(Idents(B2_HU_P4_obj), B2_HU_P4_obj$orig.ident)
b2_hu_p4 <- as.data.frame(b2_hu_p4)
b2_hu_p4$Var1 <- as.character(b2_hu_p4$Var1)
b2_hu_p4$Var2 <- "b2_hu_p4"
B2_HU_P4_PLOT <- ggplot(b2_hu_p4, aes(x = "B2_HU_P4", y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  scale_fill_manual(values = brewer.pal(12, "Paired"))

### Output graphs
bar_chart_grpah <- plot_grid(WT_P0_PLOT, WT_P4_PLOT, WT_HU_P0_PLOT, WT_HU_P4_PLOT,B2_P0_PLOT, B2_P4_PLOT, B2_HU_P0_PLOT, B2_HU_P4_PLOT, nrow = 1)
