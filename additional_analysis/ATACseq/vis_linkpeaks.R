library(chromVAR)
library(motifmatchr)
library(Seurat)
library(scran)
library(Signac)
library(Matrix)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(patchwork)
## motif analysis 
library(JASPAR2020)
library(TFBSTools)
library(future)
library(dorothea)
library(ComplexHeatmap)




clean_module <- readRDS("~/CASCADEpaper/paper/Fig5_archetype/clean_module.Rds")
names(clean_module) <- c("AR", "Inflammation", "NE1", "NE2", "Cycling", "Glycolysis")


paths_split <- unlist(lapply(strsplit(paths, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
paths_split3 <- unlist(lapply(strsplit(paths3, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
paths_split <- c(paths_split, paths_split3)


srt <- readRDS(paste0( i, "/", i, "_multi_srt.Rds"))
DefaultAssay(srt) <- "RNA"
plotlist <- mapply(plot_fun, srt = list(srt), 
                   reduction = c("umap", "umap.atac", "weighted.umap", "umap.peak","weighted.peak.umap"), 
                   title = c("RNA UMAP", "ATAC UMAP", "Joint UMAP", "Peak UMAP", "Joint Peak UMAP"),
                   SIMPLIFY = F)
g <- wrap_plots(plotlist, ncol = 1)
ggsave(paste0(i ,"/Dim_featureplot.pdf"), width = 18, height = 12)

DefaultAssay(srt) <- "peaks"
srt <- LinkPeaks(
  object = srt,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = clean_module[["AR"]]
)
