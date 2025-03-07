library(Seurat)
library(scran)
library(data.table)
library(dplyr)
library(patchwork)

paths <- system("realpath /trigos_team/CASCADE/Analysis/230113_seurat_intron/*/*/*non_tumour_removed.Rds", intern = T)
paths_split <- unlist(lapply(strsplit(paths, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
paths2 <- system("realpath /trigos_team/CASCADE/Analysis/240603_seurat_intron/CA0027/dura*/*non_tumour_removed_v4.Rds", intern = T)
paths_split2 <- unlist(lapply(strsplit(paths2, split = "/"), function(x) paste(x[6], x[7], sep = "_")))

paths <- c(paths, paths2)
paths_split <- c(paths_split, paths_split2)

paths <- setNames(paths, paths_split)

ind <- c(8, 19,20,27,32, 21:23)
paths <- paths[ind]

pdf("~/CASCADEpaper/paper/NE_transition/KO_targets.pdf", width = 8, height = 12)
for(i in paths_split[ind]){
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
  
  p <- FeaturePlot(srt, c("HDAC9", "HLA-A", "HLA-B", "HLA-E", "ASCL1", "IRF7", "IRF9", "STAT3", "TGFB1"), ncol = 3)
  f <- VlnPlot(srt, c("HDAC9", "HLA-A", "HLA-B", "HLA-E", "ASCL1", "IRF7", "IRF9", "STAT3", "TGFB1"), ncol = 3) + plot_layout(axes = "collect")
  print((p/f)+plot_layout(nrow  = 2) + plot_annotation(title = i))
  p <- FeaturePlot(srt, c("HDAC9","NFKB1A", "FOLH1", "FOXA1", "FOXA2", "SYP", "CHGA", "POU3F2","SOX2"), ncol = 3)
  f <- VlnPlot(srt, c("HDAC9","NFKB1A", "FOLH1", "FOXA1", "FOXA2", "SYP", "CHGA", "POU3F2","SOX2"), ncol = 3) + plot_layout(axes = "collect")
  print((p/f)+plot_layout(nrow  = 2) + plot_annotation(title = i))
  p <- FeaturePlot(srt, c("HDAC9","AR", "KLK3", "FOLH1", "DLL3", "CDH1"), ncol = 3)
  f <- VlnPlot(srt, c("HDAC9","AR", "KLK3", "FOLH1", "DLL3", "CDH1"), ncol = 3) + plot_layout(axes = "collect")
  print((p/f)+plot_layout(nrow  = 2) + plot_annotation(title = i))
  
  print(i)
  
}
dev.off()



pdf("~/CASCADEpaper/paper/NE_transition/targets_integrated.pdf", width = 8, height = 7)
for (i in c("STEAP1", "STEAP2", "DLL3", "IRF2", "SOX2", "HNF4G","NFKB1A")){
p <- try(FeaturePlot(test, "STEAP1", raster = F))
print(p)
p <- try(FeaturePlot(test, "STEAP1", cols = rev(brewer.pal(10, "RdYlBu")), raster = F))
print(p)
}
dev.off()


