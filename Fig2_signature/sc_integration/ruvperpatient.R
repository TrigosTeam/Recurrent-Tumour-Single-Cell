library(ruvIIInb)
library(SingleCellExperiment)
library(scater)
library(scran)
library(scuttle)
library(edgeR)
library(SingleR)
library(celldex)
library(hrbrthemes)
library(tidyverse)
library(ggplot2)
library(uwot)
library(scMerge)
library(Seurat)
library(randomcoloR)
library(dittoSeq)
library(pheatmap)
library(gridExtra)
library(igraph)
library(DelayedArray)
library(simspec)
library(simspec)
library(dplyr)
library(patchwork)
library(RColorBrewer)

data(Hs.schk)

# tumor normal --------------------------------------------------------------
paths <- system("realpath ~/Fig2_signature/sc_integration/perpatient/tumor_normal/*simspec_srt.Rds", intern = T)
patients <- substr(unlist(lapply(strsplit(paths, split = "/"), function(x) x[8])), 1, 6)
paths <- setNames(paths, patients)

for (i in patients){
test <- readRDS(paths[i])

samples <- test$sample
temp <- structure(seq(length(unique(samples))), names = unique(samples))
batch <- temp[samples]
ctl <- as.character(Hs.schk)
ind<-rownames(test) %in% ctl

M <- matrix(0,ncol(test),length(unique(test$seurat_clusters)))
cl<- sort(unique(as.numeric(unique(test$seurat_clusters))))
for(CL in cl){
  M[which(as.numeric(test$seurat_clusters)==CL),CL] <- 1
}

ruv3nb_out<-try(ruvIIInb::fastruvIII.nb(Y=DelayedArray(as.array(test@assays$RNA@counts)), # count matrix with genes as rows and cells as columns
                                        M=M, #Replicate matrix constructed as above
                                        ctl=ind, #A vector denoting control genes
                                        k=2, # dimension of unwanted variation factors
                                        ncores = 8,
                                        batch = batch
))

#saveRDS(ruv3nb_out, paste0("~/Fig2_signature/sc_integration/perpatient/tumor_normal/", i, "_ruvout.Rds"))
print("tumor only ruv finished")

logPAC <- as(log(get.res(ruv3nb_out,
                         type = "quantile", batch = ruv3nb_out$batch) + 1), "sparseMatrix")
rownames(logPAC) <- rownames(ruv3nb_out$counts)
colnames(logPAC) <- colnames(ruv3nb_out$counts)

saveRDS(logPAC, paste0("~/Fig2_signature/sc_integration/perpatient/tumor_normal/", i, "_logPAC.Rds"))

rm( list = "ruv3nb_out")
gc()

test[["ruv3"]]<- CreateAssayObject(data = logPAC)
DefaultAssay(test) <-"ruv3"
test[["ruv3"]]@scale.data <- as.matrix(test[["ruv3"]]@data)
test <- FindVariableFeatures(test)
test <- test %>% 
  RunPCA(features = VariableFeatures(test)) %>% 
  FindNeighbors(dims = 1:20) %>% 
  FindClusters(resolution = 1) %>% 
  RunUMAP(dims = 1:20)

sites<- test$sample
sites[grep("brain", sites)] <- "brain"
sites[grep("dura", sites)] <- "dura"
sites[grep("prostate", sites)] <- "prostate"
sites[grep("liver", sites)] <- "liver"
sites[grep("fat", sites)] <- "fat"
sites[grep("lymph|hilar", sites)] <- "LN"
sites[grep("rib|vertebra", sites)] <- "bone"
sites[grep("lung", sites)] <- "lung"
sites[grep("abdomen", sites)] <- "abdomen"
sites[grep("bladder", sites)] <- "bladder"
test$site <- sites
saveRDS(test, paste0("~/Fig2_signature/sc_integration/perpatient/tumor_normal/", i, "_ruvsrt.Rds"))

pdf(paste0("~/Fig2_signature/sc_integration/perpatient/tumor_normal/", i, "_ruvUMAP.pdf"), width = 12, height = 10)
p <- DimPlot(test, label = T, raster = F)
print(p)
p <- DimPlot(test, group.by = "patient", raster = F, label = T, repel = T) +NoLegend()
print(p)
p <- DimPlot(test, group.by = "sample", raster = F, label = T, repel = T)+NoLegend()
print(p)
p <- DimPlot(test, group.by = "Cell.Identity", raster = F, label = T)+NoLegend()
print(p)
p <- DimPlot(test, group.by = "site", raster = F, label = T)+NoLegend()
print(p)
for (g in  c("AR", "ASCL1", "FOLH1")[c("AR", "ASCL1", "FOLH1") %in% rownames(test)]){
  p <- FeaturePlot(test, g, cols = rev(brewer.pal(10, "RdYlBu")), raster = F)
  print(p)
}
dev.off()
}
rm(list = ls())
gc()


# tumor_only --------------------------------------------------------------

paths <- system("realpath ~/Fig2_signature/sc_integration/perpatient/tumor_only/*simspec_srt.Rds", intern = T)
patients <- substr(unlist(lapply(strsplit(paths, split = "/"), function(x) x[8])), 1, 6)
paths <- setNames(paths, patients)

for (i in patients){
  test <- readRDS(paths[i])
  
  samples <- test$sample
  temp <- structure(seq(length(unique(samples))), names = unique(samples))
  batch <- temp[samples]
  ctl <- as.character(Hs.schk)
  ind<-rownames(test) %in% ctl
  
  M <- matrix(0,ncol(test),length(unique(test$seurat_clusters)))
  cl<- sort(unique(as.numeric(unique(test$seurat_clusters))))
  for(CL in cl){
    M[which(as.numeric(test$seurat_clusters)==CL),CL] <- 1
  }
  
  ruv3nb_out<-try(ruvIIInb::fastruvIII.nb(Y=DelayedArray(as.array(test@assays$RNA@counts)), # count matrix with genes as rows and cells as columns
                                          M=M, #Replicate matrix constructed as above
                                          ctl=ind, #A vector denoting control genes
                                          k=2, # dimension of unwanted variation factors
                                          ncores = 8,
                                          batch = batch
  ))
  
  #saveRDS(ruv3nb_out, paste0("~/Fig2_signature/sc_integration/perpatient/tumor_only/", i, "_ruvout.Rds"))
  print("tumor only ruv finished")
  
  logPAC <- as(log(get.res(ruv3nb_out,
                           type = "quantile", batch = ruv3nb_out$batch) + 1), "sparseMatrix")
  rownames(logPAC) <- rownames(ruv3nb_out$counts)
  colnames(logPAC) <- colnames(ruv3nb_out$counts)
  
  saveRDS(logPAC, paste0("~/Fig2_signature/sc_integration/perpatient/tumor_only/", i, "_logPAC.Rds"))
  
  rm( list = "ruv3nb_out")
  gc()
  
  test[["ruv3"]]<- CreateAssayObject(data = logPAC)
  DefaultAssay(test) <-"ruv3"
  test[["ruv3"]]@scale.data <- as.matrix(test[["ruv3"]]@data)
  test <- FindVariableFeatures(test)
  test <- test %>% 
    RunPCA(features = VariableFeatures(test)) %>% 
    FindNeighbors(dims = 1:20) %>% 
    FindClusters(resolution = 1) %>% 
    RunUMAP(dims = 1:20)
  
  sites<- test$sample
  sites[grep("brain", sites)] <- "brain"
  sites[grep("dura", sites)] <- "dura"
  sites[grep("prostate", sites)] <- "prostate"
  sites[grep("liver", sites)] <- "liver"
  sites[grep("fat", sites)] <- "fat"
  sites[grep("lymph|hilar", sites)] <- "LN"
  sites[grep("rib|vertebra", sites)] <- "bone"
  sites[grep("lung", sites)] <- "lung"
  sites[grep("abdomen", sites)] <- "abdomen"
  sites[grep("bladder", sites)] <- "bladder"
  test$site <- sites
  saveRDS(test, paste0("~/Fig2_signature/sc_integration/perpatient/tumor_only/", i, "_ruvsrt.Rds"))
  
  pdf(paste0("~/Fig2_signature/sc_integration/perpatient/tumor_only/", i, "_ruvUMAP.pdf"), width = 12, height = 10)
  p <- DimPlot(test, label = T, raster = F)
  print(p)
  p <- DimPlot(test, group.by = "patient", raster = F, label = T, repel = T) +NoLegend()
  print(p)
  p <- DimPlot(test, group.by = "sample", raster = F, label = T, repel = T)+NoLegend()
  print(p)
  p <- DimPlot(test, group.by = "Cell.Identity", raster = F, label = T)+NoLegend()
  print(p)
  p <- DimPlot(test, group.by = "site", raster = F, label = T)+NoLegend()
  print(p)
  for (g in  c("AR", "ASCL1", "FOLH1")[c("AR", "ASCL1", "FOLH1") %in% rownames(test)]){
    p <- FeaturePlot(test, g, cols = rev(brewer.pal(10, "RdYlBu")), raster = F)
    print(p)
  }
  dev.off()
}