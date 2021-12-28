setwd("/data/gaos2/tmp/NHGRI/batch2/BMSeurat")
####R 3.4
#Note : We did this because we do not use correct parameters. Too many cells were filtered out 2500 is too small
library(Seurat)
library(Matrix)
rm(list=ls())
load("BM.combined_CCA_step1_withbatch2.RData")
# t-SNE and Clustering
S.combined_cca <- RunTSNE(S.combined_cca, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = T)
S.combined_cca <- FindClusters(S.combined_cca, reduction.type = "cca.aligned", resolution = 0.6, dims.use = 1:20)
# Visualization
png("tSNE_caa.png", width=2000, height=1000, res=200)
p1 <- TSNEPlot(S.combined_cca, do.return = T, pt.size = 0.5, group.by = "group")
p2 <- TSNEPlot(S.combined_cca, do.label = T, do.return = T, pt.size = 0.5)
print(plot_grid(p1, p2))
dev.off()
save(S.combined, S.combined_cca, vargene, file="BM.combined_CCA_step2_withbatch2.RData")

#nk.markers <- FindConservedMarkers(immune.combined, ident.1 = 7, grouping.var = "stim", print.bar = FALSE)
#head(nk.markers)
