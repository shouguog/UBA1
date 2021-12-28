###R4.0.2
rm(list=ls())
cellsSeurat<-rownames(read.csv("../../../../../mywork/NHGRI/batch2/BMSeurat/metaDatatSNE.csv", row.names = 1))
for(ID in 1:5){
  con <- file(paste0("../../PT", ID, "_BM.sam"), "r")
  lines<-readLines(con) # empty
  close(con)
  cells<-c()
  for(ii in 1:length(lines)){
    cells<-c(cells,strsplit(strsplit(lines[ii], "CB:Z:")[[1]][2], "\t")[[1]][1])
  }
  cells<-gsub("-1", paste0("-", ID), cells)
  zz<-file(paste0("PT", ID, "_BM_barcode.tsv"), "w")
  for(cell in unique(intersect(cells, cellsSeurat))){
    cat(gsub(paste0("-", ID),"-1",cell), file=zz);cat("\n", file = zz)
  }
  close(zz)
}

for(ID in 6:9){
  con <- file(paste0("../../PT", ID, "_BM.sam"), "r")
  lines<-readLines(con) # empty
  close(con)
  cells<-c()
  for(ii in 1:length(lines)){
    cells<-c(cells,strsplit(strsplit(lines[ii], "CB:Z:")[[1]][2], "\t")[[1]][1])
  }
  cells<-gsub("-1", paste0("-", ID+4), cells)
  zz<-file(paste0("PT", ID, "_BM_barcode.tsv"), "w")
  for(cell in unique(intersect(cells, cellsSeurat))){
    cat(gsub(paste0("-", ID+4),"-1",cell), file=zz);cat("\n", file = zz)
  }
  close(zz)
}

###R4.0.2
rm(list=ls())
cellsSeurat<-rownames(read.csv("../../../../../mywork/NHGRI/batch2/CD34Seurat/metaDatatSNE.csv", row.names = 1))
for(ID in 1:5){
  con <- file(paste0("../../PT", ID, "_CD34.sam"), "r")
  lines<-readLines(con) # empty
  close(con)
  cells<-c()
  for(ii in 1:length(lines)){
    cells<-c(cells,strsplit(strsplit(lines[ii], "CB:Z:")[[1]][2], "\t")[[1]][1])
  }
  cells<-gsub("-1", paste0("-", ID), cells)
  zz<-file(paste0("PT", ID, "_CD34_barcode.tsv"), "w")
  for(cell in unique(intersect(cells, cellsSeurat))){
    cat(gsub(paste0("-", ID),"-1",cell), file=zz);cat("\n", file = zz)
  }
  close(zz)
}

for(ID in 6:9){
  con <- file(paste0("../../PT", ID, "_CD34.sam"), "r")
  lines<-readLines(con) # empty
  close(con)
  cells<-c()
  for(ii in 1:length(lines)){
    cells<-c(cells,strsplit(strsplit(lines[ii], "CB:Z:")[[1]][2], "\t")[[1]][1])
  }
  cells<-gsub("-1", paste0("-", ID+4), cells)
  zz<-file(paste0("PT", ID, "_CD34_barcode.tsv"), "w")
  for(cell in unique(intersect(cells, cellsSeurat))){
    cat(gsub(paste0("-", ID+4),"-1",cell), file=zz);cat("\n", file = zz)
  }
  close(zz)
}
