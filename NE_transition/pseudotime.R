library(Seurat)
library(scran)
library(data.table)
library(dplyr)
library(patchwork)
library(tidyverse)
library(MAST)
library(pals)

paths <- system("realpath /trigos_team/CASCADE/Analysis/230113_seurat_intron/*/*/*non_tumour_removed.Rds", intern = T)
paths_split <- unlist(lapply(strsplit(paths, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
paths2 <- system("realpath /trigos_team/CASCADE/Analysis/240603_seurat_intron/CA0027/dura*/*non_tumour_removed_v4.Rds", intern = T)
paths_split2 <- unlist(lapply(strsplit(paths2, split = "/"), function(x) paste(x[6], x[7], sep = "_")))

paths <- c(paths, paths2)
paths_split <- c(paths_split, paths_split2)
paths <- setNames(paths,paths_split)

srt<- readRDS(paths[1])
DefaultAssay(srt) <- "RNA"
lib.median <- median(srt$nCount_RNA)
target_pseudocount <- 1
srt <- srt %>% 
  NormalizeData(scale.factor = lib.median / target_pseudocount) %>% 
  ScaleData(do.scale = F)

set.seed(100)
poisson_fit2 <- modelGeneVarByPoisson(srt@assays$RNA@counts, dispersion = 0.2, size.factors = target_pseudocount * srt$nCount_RNA / lib.median)
residuals2 <- poisson_fit2$total - poisson_fit2@metadata$trend(poisson_fit2$mean)
names(residuals2) <- rownames(srt)
top_genes2 <- rownames(srt)[order(residuals2, decreasing = TRUE)]

srt <- srt %>% 
  RunPCA(features = top_genes2[1:2000]) %>% 
  FindNeighbors(dims = 1:20) %>% 
  FindClusters(resolution = 1) %>% 
  RunUMAP(dims = 1:20)

library(scuttle)
meta <- perCellQCMetrics(srt[["RNA"]]@counts)
r <- perCellQCFilters(meta)
colSums(as.matrix(r))
clean_module <- readRDS("~/regulon/new_regulon/archetype/2024_06/clean_module.Rds")

library(SeuratWrappers)
srt <- AddModuleScore(srt, features  = clean_module)
names(clean_module)
colnames(srt@meta.data)[grep("Cluster", colnames(srt@meta.data))] <- names(clean_module)

DimPlot(srt)
FeaturePlot(srt, names(clean_module))
library(monocle3)
cds <- as.cell_data_set(srt)
cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition =F, close_loop = F,learn_graph_control = list(orthogonal_proj_tip= TRUE ))
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = WhichCells(srt, ident = "2"))

plot_cells(
  cds =cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)

pseudo <- cds@principal_graph_aux@listData$UMAP$pseudotime
srt$monocle <- pseudo[colnames(srt)]
FeaturePlot(srt, "monocle")

plot_cells(cds, color_cells_by = "seurat_clusters", label_leaves=FALSE,
           label_branch_points=FALSE, show_trajectory_graph = F,label_groups_by_cluster=FALSE)
cds <- cluster_cells(cds, resolution=1e-5)
plot_cells(cds, label_leaves=FALSE,color_cells_by="partition", group_cells_by="partition",
           label_branch_points=FALSE, show_trajectory_graph = F)
cds <- learn_graph(cds)

plot_cells(cds,
           label_groups_by_cluster=T,
           label_leaves=T,
           label_branch_points=T)

