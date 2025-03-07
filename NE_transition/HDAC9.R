library(Seurat)
library(scran)
library(data.table)
library(dplyr)
library(patchwork)
library(tidyverse)
library(MAST)
library(pals)
#include new sample

srt <- readRDS("~/integration/2024_06/tumor_only_ruvsrt.Rds")

pdf("~/CASCADEpaper/paper/NE_transition/HDAC9_integrated_plot.pdf", width = 8, height = 4)
DimPlot(srt, label = T, raster = F, group.by = "patient")
FeaturePlot(srt, c("AR", "ASCL1"),  raster = F)
FeaturePlot(srt, c("HDAC9"),  raster = F) +  FeaturePlot(srt, "HDAC9", cols =  rev(brewer.rdylbu(10)),raster = F)+plot_layout(ncol = 2)
dev.off()


paths <- system("realpath /trigos_team/CASCADE/Analysis/230113_seurat_intron/*/*/*non_tumour_removed.Rds", intern = T)
paths_split <- unlist(lapply(strsplit(paths, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
paths2 <- system("realpath /trigos_team/CASCADE/Analysis/240603_seurat_intron/CA0027/dura*/*non_tumour_removed_v4.Rds", intern = T)
paths_split2 <- unlist(lapply(strsplit(paths2, split = "/"), function(x) paste(x[6], x[7], sep = "_")))

paths <- c(paths, paths2)
paths_split <- c(paths_split, paths_split2)

paths <- setNames(paths, paths_split)

pdf("~/CASCADEpaper/paper/NE_transition/HDAC9_persample_plot.pdf", width = 8, height = 4)
for(i in paths_split){
  srt <- readRDS(paths[i])
  
  DefaultAssay(srt) <- "RNA"
  lib.median <- median(srt$nCount_RNA)
  target_pseudocount = 1
  set.seed(100)
  poisson_fit <- modelGeneVarByPoisson(srt@assays$RNA@counts, dispersion = 0.2, size.factors = target_pseudocount * srt$nCount_RNA / lib.median)
  residuals  <- poisson_fit$total - poisson_fit@metadata$trend(poisson_fit$mean)
  names(residuals ) <- rownames(srt )
  top_genes  <- rownames(srt)[order(residuals , decreasing = TRUE)]
  
  srt <- srt %>% 
    RunPCA(features = top_genes[1:2000]) %>% 
    FindNeighbors(dims = 1:20) %>% 
    FindClusters(resolution = 1) %>% 
    RunUMAP(dims = 1:20)
  
  p <- try(FeaturePlot(srt, "HDAC9")+ VlnPlot(srt, features = "HDAC9")  + plot_annotation(title = i))+plot_layout(ncol = 2)
  print(p)
  print(i)
}
dev.off()