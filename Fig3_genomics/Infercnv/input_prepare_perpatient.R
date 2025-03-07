library(Seurat)
library(scran)
library(dplyr)
library(patchwork)
library(scCustomize)
library(ggplot2)

normal <- readRDS("~/CASCADEpaper/paper/Fig3/Infercnv/normal_cells.Rds")



paths <- system("realpath /trigos_team/CASCADE/Analysis/230113_seurat_intron/*/*/*non_tumour_removed.Rds", intern = T)
paths_split <- unlist(lapply(strsplit(paths, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
paths <- setNames(paths, paths_split)
patients <- unique(substr(paths_split, 1, 6))

subclone_anno <- readRDS("~/CASCADEpaper/paper/Fig3/Infercnv/subclone_anno.Rds")
for (i in names(subclone_anno)){
  name <- paste(i, names(subclone_anno[[i]]), sep = "_") 
  subclone_anno[[i]] <- paste(i, subclone_anno[[i]], sep = "-") 
  names(subclone_anno[[i]]) <- name
}
names(subclone_anno) <- NULL
subclone_anno  <- unlist(subclone_anno)

for(i in patients[8:9]){
  samples <- grep(i, paths_split, value = T)
  srt <- lapply(paths[samples], readRDS)
  srt <- lapply(srt, function(x) {
    DefaultAssay(x) <- "RNA"
    return(x)
    })
  names(srt) <- samples
  srt <- Merge_Seurat_List(srt, project = i, merge.data = F, add.cell.ids = samples)
  srt$subclone <- "no subclone"
  srt$subclone[colnames(srt) %in% names(subclone_anno)] <- subclone_anno[intersect(colnames(srt), names(subclone_anno))]
  tumor_normal <- merge(srt, normal)
  counts <- tumor_normal[["RNA"]]@counts

  tumor_normal$cell_id[is.na(tumor_normal$cell_id)] <- tumor_normal$subclone[is.na(tumor_normal$cell_id)]
  anno <- as.data.frame(tumor_normal$cell_id)
  colnames(anno) <- NULL
  saveRDS(counts, paste0("~/CASCADEpaper/paper/Fig3/Infercnv/count_matrix/", i, ".Rds"))
  saveRDS(anno, paste0("~/CASCADEpaper/paper/Fig3/Infercnv/anno/", i, ".Rds"))
  print(i)
}