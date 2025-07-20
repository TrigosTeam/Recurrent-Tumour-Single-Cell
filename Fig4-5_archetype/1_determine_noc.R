
reticulate::use_condaenv("reticulate_PCHA", conda = "auto",required = TRUE)
# devtools::install_version("matrixStats", version= "1.1.0", lib = "/trigos_team/R_Library/ParetoTI")

library(dbplyr) ## need to be version 2.3.4
library(matrixStats, lib.loc = "/trigos_team/R_Library/ParetoTI") # version 1.1.0
library(ParetoTI)
library(Seurat)
library(dplyr)
library(scater)
library(scales)
library(ggpubr)
library(tidyverse)
library(patchwork)
library(scran)
library(msigdbr)

paths <- system("realpath ~/230113_seurat_intron/*/*/*non_tumour_removed.Rds", intern = T)
paths_split <- unlist(lapply(strsplit(paths, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
paths <- setNames(paths, paths_split)
paths2 <- system("realpath ~/240603_seurat_intron/CA0027/dura*/*non_tumour_removed.Rds", intern = T)
paths_split2 <- unlist(lapply(strsplit(paths2, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
paths2 <- setNames(paths2, paths_split2)
paths <- c(paths, paths2)

Hall <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
KEGG <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")%>%
  dplyr::select(gs_name, gene_symbol)
Reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME")%>%
  dplyr::select(gs_name, gene_symbol)
C6 <- msigdbr(species = "Homo sapiens", category = "C6") %>%
  dplyr::select(gs_name, gene_symbol)
  
all <- rbind(Hall, KEGG, Reactome,C6)

markerlist <- split(all, all$gs_name) %>% lapply(., function(x) x[["gene_symbol"]])


for (i in paths_split){
  srt <- readRDS(paths[i])
  DefaultAssay(srt) <- "RNA"
  lib.median <- median(srt$nCount_RNA)
  target_pseudocount = 1
  set.seed(100)
  poisson_fit <- modelGeneVarByPoisson(srt[["RNA"]]@counts, dispersion = 0.2, size.factors = target_pseudocount * srt$nCount_RNA / lib.median)
  residuals  <- poisson_fit$total - poisson_fit@metadata$trend(poisson_fit$mean)
  names(residuals ) <- rownames(srt )
  top_genes  <- rownames(srt)[order(residuals , decreasing = TRUE)]

  srt <- srt %>%
    RunPCA(features = top_genes[1:2000]) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 1) %>%
    RunUMAP(dims = 1:20)
  # plotReducedDim(sce, ncomponents = 3, dimred = "PCA",colour_by = "seurat_clusters")
  PCs4arch= t(srt[["pca"]]@cell.embeddings)
  arc_ks = k_fit_pch(PCs4arch, ks = 2:8, check_installed = T,
                     bootstrap = T, bootstrap_N = 200, maxiter = 100,
                     bootstrap_type = "m", seed = 2543, 
                     volume_ratio = "t_ratio", # set to "none" if too slow
                     delta=0, conv_crit = 1e-03, order_type = "align",
                     sample_prop = 0.65)
  pdf(paste0("~/Fig4-5_archetype/1_noc_plots/", i, ".pdf"), width = 8, height = 6)
  f <- DimPlot(srt, label = T)
  print(f)
  f <- plot_arc_var(arc_ks, type = "varexpl", point_size = 2, line_size = 1.5) + theme_bw()
  f2 <- plot_arc_var(arc_ks, type = "res_varexpl", point_size = 2, line_size = 1.5) + theme_bw()
  p <- plot_arc_var(arc_ks, type = "total_var", point_size = 2, line_size = 1.5) +
    theme_bw() +
    ylab("Mean variance in position of vertices")
  
  # Show t-ratio
  g <- plot_arc_var(arc_ks, type = "t_ratio", point_size = 2, line_size = 1.5) + theme_bw()+
    ylab("t ratio")
  
  print(f+f2+p+g+plot_annotation(title = i) + plot_layout(ncol = 2, nrow = 2))
  dev.off()
  saveRDS(arc_ks,paste0("~/Fig4-5_archetype/1_arc_ks/", i, ".Rds"))
  
  srt <- AddModuleScore(srt, features = markerlist)
  modules <- srt@meta.data[, grep("Cluster", colnames(srt@meta.data))]
  colnames(modules) <- names(all)
  saveRDS(modules, paste0("~/Fig4-5_archetype/1_modulescore_msigdb/", i, ".Rds"))
  
  print(i)
}
