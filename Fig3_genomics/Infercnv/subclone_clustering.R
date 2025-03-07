library(infercnv)
library(Seurat)
library(dplyr)
library(patchwork)
library(pals)
paths <- system("realpath /trigos_team/CASCADE/Analysis/230113_seurat_intron/*/*/*non_tumour_removed.Rds", intern = T)
paths_split <- unlist(lapply(strsplit(paths, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
paths <- setNames(paths, paths_split)
paths2 <- system("realpath /trigos_team/CASCADE/Analysis/240603_seurat_intron/CA0027/dura*/*non_tumour_removed.Rds", intern = T)
paths_split2 <- unlist(lapply(strsplit(paths2, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
paths2 <- setNames(paths2, paths_split2)
paths <- c(paths, paths2)

inf_path <- list.files("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/RNA2/unclustered")
subclone_anno <- readRDS("~/CASCADEpaper/paper/Fig3/Infercnv/subclone_anno.Rds")
subclone_clist <- list()

pdf("~/CASCADEpaper/paper/Fig3/Infercnv/subclone_cluster_test.pdf", width = 12, height = 8)
for (i in inf_path[1:3]){
  anno <- readRDS(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/anno/",i,".Rds"))
  inf_obj <- readRDS(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/subclone/", i, "/preliminary.infercnv_obj" ))
  mat <- inf_obj@expr.data
  mat <- read.delim()
  subclones <- subclone_anno[[i]]
  srt <- readRDS(paths[i])
  srt$subclones <- subclones[colnames(srt)]
  inf_srt <- CreateSeuratObject(mat, project = "SeuratProject", assay = "Infercnv",
                                min.cells = 0, min.features = 0, names.field = 1,
                                names.delim = "_", meta.data = NULL)
  
  inf_srt$rna_cluster <- anno[colnames(inf_srt), 1]
  inf_srt <- ScaleData(inf_srt , do.scale = F, do.center = F)
  inf_srt <- RunPCA(inf_srt, features = rownames(inf_srt))
  inf_srt<- FindNeighbors(inf_srt, dims = 1:20) %>% RunUMAP(dims = 1:20) %>% FindClusters(resolution = c(0.5, 1, 1.5))

  subclone_cluster <- inf_srt$Infercnv_snn_res.1.5
  srt$subclone_cluster <- subclone_cluster[colnames(srt)]
  
  g1 <- DimPlot(inf_srt, label = T, group.by = "Infercnv_snn_res.0.5") +DimPlot(inf_srt, label = T, group.by = "Infercnv_snn_res.1")+DimPlot(inf_srt, label = T, group.by = "Infercnv_snn_res.1.5")+NoLegend() 
  g2 <- DimPlot(inf_srt, label = T, group.by = "rna_cluster") + NoLegend() 
  g3 <-DimPlot(srt, group.by = "subclone_cluster", label = T, cols = as.character(polychrome(36)))
  g4 <- DimPlot(srt, group.by = "subclones", label = T)
  print(g1 + g2 +g3 +g4+ plot_layout(nrow = 2, ncol = 3) + plot_annotation(title = i))
  
  subclone_clist[[i]] <- inf_srt@meta.data
}
dev.off()
saveRDS(subclone_clist, "~/CASCADEpaper/paper/Fig3/Infercnv/subclone_cluster_meta.Rds")


