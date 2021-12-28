setwd("/data/gaos2/tmp/NHGRI/batch2/BMSeurat")
####R 3.4
#Note : We did this because we do not use correct parameters. Too many cells were filtered out 2500 is too small
library(Seurat)
library(Matrix)
rm(list=ls())
load("BM.combined_CCA_step2_withbatch2.RData")
cells<-rownames(S.combined_cca@meta.data)
set.seed(666)
cellsSelected<-cells[sample(1:length(cells),20000)]
S.combined_cca_sub<-SubsetData(S.combined_cca, cells=cellsSelected)
S.combined_cca_sub@raw.data<-S.combined_cca@raw.data[, cellsSelected]
S.combined_cca_sub<-SetAllIdent(S.combined_cca_sub, id = "res.0.6")
for(cluster in 0:length(unique(S.combined_cca@meta.data$res.0.6))){
  cat(cluster)
  markers <- FindMarkers(S.combined_cca_sub, ident.1 = cluster,print.bar = TRUE, logfc.threshold = 0.1)
  write.csv(markers, file=paste("BM.aligned_cluster_", cluster, ".csv",sep=""))
}

