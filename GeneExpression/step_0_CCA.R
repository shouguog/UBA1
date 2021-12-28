####R 4.02
#library(Seurat)
#rm(list=ls())
###Let us read the new data
#PT6<-Read10X("../210517_A00941_0505_AHW5VFDMXX/cellranger/UBA1_GRCh38_VEXAS_PT6_KOWJ_BM/outs/filtered_feature_bc_matrix/")
#PT7<-Read10X("../210517_A00941_0505_AHW5VFDMXX/cellranger/UBA1_GRCh38_VEXAS_PT7_CIRE_GEX_BM/outs/filtered_feature_bc_matrix/")
#PT8<-Read10X("../210517_A00941_0505_AHW5VFDMXX/cellranger/UBA1_GRCh38_VEXAS_PT8_SALL_GEX_BM/outs/filtered_feature_bc_matrix/")
#PT9<-Read10X("../210517_A00941_0505_AHW5VFDMXX/cellranger/UBA1_GRCh38_VEXAS_PT9_DAND_GEX_BM/outs/filtered_feature_bc_matrix/")
#save(PT6, PT7, PT8, PT9, file="PT6_9.RData", version="2.6")

####R 3.4
setwd("/data/gaos2/tmp/NHGRI/batch2/BMSeurat")
library(Seurat)
rm(list=ls())
load("../../BMSeurat8runNoNorm/countsMeta.RData")
load("PT6_9.RData")
colnames(PT6)<-gsub("-1", "-10", colnames(PT6))
colnames(PT7)<-gsub("-1", "-11", colnames(PT7))
colnames(PT8)<-gsub("-1", "-12", colnames(PT8))
colnames(PT9)<-gsub("-1", "-13", colnames(PT9))
S.combined.data<-cbind(S.combined.data, PT6)
S.combined.data<-cbind(S.combined.data, PT7)
S.combined.data<-cbind(S.combined.data, PT8)
S.combined.data<-cbind(S.combined.data, PT9)
rm(PT6, PT7, PT8, PT9)
S.combined <- CreateSeuratObject(raw.data = S.combined.data)
cells<-rownames(metaData)
S.combined.1<-SubsetData(S.combined, cells = cells[grepl("-1", cells)])
S.combined.2<-SubsetData(S.combined, cells = cells[grepl("-2", cells)])
S.combined.3<-SubsetData(S.combined, cells = cells[grepl("-3", cells)])
S.combined.4<-SubsetData(S.combined, cells = cells[grepl("-4", cells)])
S.combined.5<-SubsetData(S.combined, cells = cells[grepl("-5", cells)])
S.combined.6<-SubsetData(S.combined, cells = cells[grepl("-6", cells)])
S.combined.7<-SubsetData(S.combined, cells = cells[grepl("-7", cells)])
S.combined.8<-SubsetData(S.combined, cells = cells[grepl("-8", cells)])
S.combined.9<-SubsetData(S.combined, cells = cells[grepl("-9", cells)])
S.combined.1@raw.data<-S.combined@raw.data[,cells[grepl("-1", cells)]]
S.combined.2@raw.data<-S.combined@raw.data[,cells[grepl("-2", cells)]]
S.combined.3@raw.data<-S.combined@raw.data[,cells[grepl("-3", cells)]]
S.combined.4@raw.data<-S.combined@raw.data[,cells[grepl("-4", cells)]]
S.combined.5@raw.data<-S.combined@raw.data[,cells[grepl("-5", cells)]]
S.combined.6@raw.data<-S.combined@raw.data[,cells[grepl("-6", cells)]]
S.combined.7@raw.data<-S.combined@raw.data[,cells[grepl("-7", cells)]]
S.combined.8@raw.data<-S.combined@raw.data[,cells[grepl("-8", cells)]]
S.combined.9@raw.data<-S.combined@raw.data[,cells[grepl("-9", cells)]]

cells<-colnames(S.combined.data)
rm(S.combined.data)
S.combined.10<-SubsetData(S.combined, cells = cells[grepl("-10", cells)])
S.combined.11<-SubsetData(S.combined, cells = cells[grepl("-11", cells)])
S.combined.12<-SubsetData(S.combined, cells = cells[grepl("-12", cells)])
S.combined.13<-SubsetData(S.combined, cells = cells[grepl("-13", cells)])

S.combined.10@raw.data<-S.combined@raw.data[,cells[grepl("-10", cells)]]
S.combined.11@raw.data<-S.combined@raw.data[,cells[grepl("-11", cells)]]
S.combined.12@raw.data<-S.combined@raw.data[,cells[grepl("-12", cells)]]
S.combined.13@raw.data<-S.combined@raw.data[,cells[grepl("-13", cells)]]

#####Let us normalize the data
S.combined <- NormalizeData(S.combined)
S.combined <- FindVariableGenes(S.combined, do.plot = T, y.cutoff = 0.6)
####union
vargene<-S.combined@var.genes
#rm(S.combined)

S.combined.1@meta.data$group<-"PT"
S.combined.2@meta.data$group<-"PT"
S.combined.3@meta.data$group<-"PT"
S.combined.4@meta.data$group<-"PT"
S.combined.5@meta.data$group<-"PT"
S.combined.6@meta.data$group<-"HD"
S.combined.7@meta.data$group<-"HD"
S.combined.8@meta.data$group<-"HD"
S.combined.9@meta.data$group<-"HD"
S.combined.10@meta.data$group<-"PT"
S.combined.11@meta.data$group<-"PT"
S.combined.12@meta.data$group<-"PT"
S.combined.13@meta.data$group<-"PT"

S.combined.1@meta.data$subject<-"PT1"
S.combined.2@meta.data$subject<-"PT2"
S.combined.3@meta.data$subject<-"PT3"
S.combined.4@meta.data$subject<-"PT4"
S.combined.5@meta.data$subject<-"PT5"
S.combined.6@meta.data$subject<-"HD1"
S.combined.7@meta.data$subject<-"HD2"
S.combined.8@meta.data$subject<-"HD3"
S.combined.9@meta.data$subject<-"HD4"
S.combined.10@meta.data$subject<-"PT6"
S.combined.11@meta.data$subject<-"PT7"
S.combined.12@meta.data$subject<-"PT8"
S.combined.13@meta.data$subject<-"PT9"

###scale data
S.combined.1<-ScaleData(S.combined.1, genes.use=vargene)
S.combined.2<-ScaleData(S.combined.2, genes.use=vargene)
S.combined.3<-ScaleData(S.combined.3, genes.use=vargene)
S.combined.4<-ScaleData(S.combined.4, genes.use=vargene)
S.combined.5<-ScaleData(S.combined.5, genes.use=vargene)
S.combined.6<-ScaleData(S.combined.6, genes.use=vargene)
S.combined.7<-ScaleData(S.combined.7, genes.use=vargene)
S.combined.8<-ScaleData(S.combined.8, genes.use=vargene)
S.combined.9<-ScaleData(S.combined.9, genes.use=vargene)
S.combined.10<-ScaleData(S.combined.10, genes.use=vargene)
S.combined.11<-ScaleData(S.combined.11, genes.use=vargene)
S.combined.12<-ScaleData(S.combined.12, genes.use=vargene)
S.combined.13<-ScaleData(S.combined.13, genes.use=vargene)


S.combined_cca <- RunMultiCCA(list(S.combined.1,S.combined.2,S.combined.3,S.combined.4,S.combined.5
                                   ,S.combined.6,S.combined.7,S.combined.8,S.combined.9
                                   ,S.combined.10,S.combined.11,S.combined.12,S.combined.13), num.ccs = 30,genes.use=vargene)
S.combined_cca <- AlignSubspace(S.combined_cca,reduction.type = "cca", grouping.var = "subject", dims.align = 1:20)
save(S.combined, S.combined_cca, vargene, file="BM.combined_CCA_step1_withbatch2.RData")
