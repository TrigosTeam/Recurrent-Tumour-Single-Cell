library(Seurat)
library(scran)
library(simspec)
library(dplyr)
library(patchwork)
library(scCustomize)
library(ggplot2)
library(SingleR)
refs <- readRDS("~/integration/tumor_normal_seperation/atlasref.Rds")


normal <- readRDS("~/CASCADEpaper/paper/Fig3/Infercnv/normal_cells.Rds")




paths <- system("realpath /trigos_team/CASCADE/Analysis/230113_seurat_intron/*/*/*non_tumour_removed.Rds", intern = T)
paths_split <- unlist(lapply(strsplit(paths, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
paths <- setNames(paths, paths_split)

for(i in paths_split){
  srt <- readRDS(paths[i])
  tumor_normal <- merge(srt, normal)
  counts <- tumor_normal[["RNA2"]]@counts
  tumor_normal$cell_id[is.na(tumor_normal$cell_id)] <- tumor_normal$seurat_clusters[is.na(tumor_normal$cell_id)]
  anno <- as.data.frame(tumor_normal$cell_id)
  colnames(anno) <- NULL
  saveRDS(counts, paste0("~/CASCADEpaper/paper/Fig3/Infercnv/count_matrix/RNA2/", i, ".Rds"))
  saveRDS(anno, paste0("~/CASCADEpaper/paper/Fig3/Infercnv/anno/RNA2/", i, ".Rds"))
}

paths2 <- system("realpath /trigos_team/CASCADE/Analysis/240603_seurat_intron/CA0027/dura*/*non_tumour_removed.Rds", intern = T)
paths_split2 <- unlist(lapply(strsplit(paths2, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
paths2 <- setNames(paths2, paths_split2)

for(i in paths_split2){
  srt <- readRDS(paths2[i])
  # pred <- SingleR(srt[["RNA"]]@data, refs, labels = list(refs[[1]]$label.main,refs[[2]]$label.main))
  # srt[["Cell.Identity"]] <- pred[colnames(srt), "labels"]
  tumor_normal <- merge(srt, normal)
  counts <- tumor_normal[["RNA2"]]@counts
  tumor_normal$cell_id[is.na(tumor_normal$cell_id)] <- tumor_normal$seurat_clusters[is.na(tumor_normal$cell_id)]
  anno <- as.data.frame(tumor_normal$cell_id)
  colnames(anno) <- NULL
  saveRDS(counts, paste0("~/CASCADEpaper/paper/Fig3/Infercnv/count_matrix/RNA2/", i, ".Rds"))
  saveRDS(anno, paste0("~/CASCADEpaper/paper/Fig3/Infercnv/anno/RNA2/", i, ".Rds"))
}