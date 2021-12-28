setwd("/data/gaos2/tmp/NHGRI/batch2/BMSeurat")
####R 3.4
library(Seurat)
rm(list=ls())
load("BM.combined_CCA_step2_withbatch2.RData")
####gET THE META DATA
metaData<-S.combined_cca@meta.data
tSNE<-S.combined_cca@dr$tsne@cell.embeddings
metaDatatSNE<-cbind(metaData, tSNE[rownames(metaData),])
write.csv(metaDatatSNE, file="metaDatatSNE.csv")
library(ggplot2)
png("tSNE_caa_plot.png", width=4000, height=4000, res=100)
print(ggplot(data=metaDatatSNE, aes(x=tSNE_1, y=tSNE_2, label=group, color=group))+geom_text())
dev.off()

png("tSNE_caa_plot_cluster.png", width=4000, height=4000, res=100)
print(ggplot(data=metaDatatSNE, aes(x=tSNE_1, y=tSNE_2, label=res.0.6, color=celltype))+geom_text())
dev.off()

