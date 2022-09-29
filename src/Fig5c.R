### Load libs required
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)

### Load Data
All_merge <- readRDS("~/all_merge.rds") # rds data created in DataPrep script
All_Luminal_obj <- subset(All_merge, idents = c("0", "8", "2", "12", "1", "11", "6", "14", "5"))
All_Luminal_P4_obj <- subset(All_Luminal_obj, subset = (Condition == "WT_HU_P4" | Condition == "B2_HU_P4" | Condition == "WT_P4" | Condition == "B2_P4"))

# GO term values

balanced = TRUE
logfc.threshold = 0.25
assay = "RNA"
max.genes = 250
test.use = 'wilcox'
p.val.cutoff = 0.001
cols = NULL
enrich.database = "GO_Biological_Process_2021"
num.pathway = 5
return.gene.list = FALSE

ident.1 = "B2_HU_P4"
ident.2 = c("WT_HU_P4", "WT_P4", "B2_P4")

#Enrich find markers - Cluster 6

All_Luminal_P4_markers_C6_obj <- FindMarkers(
  object = All_Luminal_P4_obj,
  slot = "data",
  subset.ident = c("6"),
  group.by = "Condition",
  ident.1 = ident.1,
  ident.2 = ident.2,
  only.pos = FALSE,
  logfc.threshold = logfc.threshold,
  test.use = test.use,
  assay = assay
)
# Enrichr for a single cluster

enrich.all.markers <- All_Luminal_P4_markers_C6_obj
enrich.pos.markers <- enrich.all.markers[enrich.all.markers[, 2] > logfc.threshold & enrich.all.markers[, 1] < p.val.cutoff, , drop = FALSE]
enrich.pos.markers.list <- rownames(x = enrich.pos.markers)[1:min(max.genes, nrow(x = enrich.pos.markers))]
enrich.pos.er <- enrichR::enrichr(genes = enrich.pos.markers.list, databases = enrich.database)
enrich.pos.er <- do.call(what = cbind, args = enrich.pos.er)
enrich.pos.er$log10pval <- -log10(x = enrich.pos.er[, paste(enrich.database, sep = ".", "P.value")])
enrich.pos.er$term <- enrich.pos.er[, paste(enrich.database, sep = ".", "Term")]
enrich.pos.er <- enrich.pos.er[1:num.pathway, ]
enrich.pos.er$group <- 6
enrich.pos.er$label <- "6"
enrich.pos.er$term <- factor(x = enrich.pos.er$term, levels = enrich.pos.er$term[order(enrich.pos.er$log10pval)])
#enrich.gene.list <- list(pos = enrich.pos.er)
enrich.pos.c6 <- enrich.pos.er

# Generate Plot
enrich.pos.c6_plot <- ggplot(enrich.pos.c6, aes_string(x = "term", y = "log10pval")) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), aes(fill = log10pval, group = group))  +
  ylim(0,15) +
  coord_flip() + xlab("Pathway") +
  scale_fill_gradient(low = "#bed1f7", high = "#2a70fa") +
  ylab("-log10(pval)") +
  ggtitle(paste(enrich.database, ident.1, sep = "_", "C6 positive markers")) +
  theme_classic() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

enrich.pos.c6_plot_list <- list(pos = enrich.pos.c6)

# Save Plot
write.csv(enrich.pos.c6_plot_list, file = "~/Mona3/GO_Lum_P4_c6.csv")
ggsave("~/ImageExports/GO_Lum_P4_c6.eps", plot = enrich.pos.c6_plot)
ggsave("~/ImageExports/GO_Lum_P4_c6.jpeg", plot = enrich.pos.c6_plot)
