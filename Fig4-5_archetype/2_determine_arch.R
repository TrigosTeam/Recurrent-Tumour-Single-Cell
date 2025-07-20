reticulate::use_condaenv("reticulate_PCHA", conda = "auto",required = TRUE)
# devtools::install_version("matrixStats", version= "1.1.0")

library(dbplyr) ## need to be version 2.3.4
library(matrixStats) # version 1.1.0
library(ParetoTI)
library(Seurat)
library(MAST)
library(dplyr)
library(scater)
library(scales)
library(ggpubr)
library(tidyverse)
library(msigdbr)
library(scran)


paths <- system("realpath ~/230113_seurat_intron/*/*/*non_tumour_removed.Rds", intern = T)
paths_split <- unlist(lapply(strsplit(paths, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
paths <- setNames(paths, paths_split)
paths2 <- system("realpath ~/240603_seurat_intron/CA0027/dura*/*non_tumour_removed.Rds", intern = T)
paths_split2 <- unlist(lapply(strsplit(paths2, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
paths2 <- setNames(paths2, paths_split2)
paths <- c(paths, paths2)

# which noc is closest to 0.05 t-ratio
noc  = setNames(rep(4, 34), paths_split)
noc[c(1,3,6,11, 18, 19, 21, 24,25, 26, 31, 32)] <- 5

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
options(Seurat.object.assay.version = "v3")

for (i in paths_split){
  srt <- readRDS(paths[i])
  DefaultAssay(srt) <- "RNA"

  lib.median <- median(srt$nCount_RNA)
  target_pseudocount = 1
  set.seed(100)
  poisson_fit <- modelGeneVarByPoisson(srt[["RNA"]]$counts, dispersion = 0.2, size.factors = target_pseudocount * srt$nCount_RNA / lib.median)
  residuals  <- poisson_fit$total - poisson_fit@metadata$trend(poisson_fit$mean)
  names(residuals ) <- rownames(srt )
  top_genes  <- rownames(srt)[order(residuals , decreasing = TRUE)]

  srt <- srt %>%
    RunPCA(features = top_genes[1:2000]) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 1) %>%
    RunUMAP(dims = 1:20)

  PCs4arch= t(srt[["pca"]]@cell.embeddings)
  arc = fit_pch_bootstrap(PCs4arch, n = 200, sample_prop = 0.7, seed = 235,
                          noc = noc[i], delta = 0, conv_crit = 1e-03, type = "m")
  arc_1 = average_pch_fits(arc)
  saveRDS(arc_1, paste0("~/Fig4-5_archetype/2_arc/", i, ".Rds"))
  
  
  activ = readRDS(paste0("~/Fig4-5_archetype/1_modulescore_msigdb/", i, ".Rds"))
  colnames(activ) <- names(markerlist)
  activ$cells <- rownames(activ)
  
  # Merge distances, gene expression and gene set activity into one matrix
  data_attr = merge_arch_dist(arc_data = arc_1, data = PCs4arch, 
                              feature_data = as.matrix(srt[["RNA"]]$data),
                              colData = activ, 
                              dist_metric = c("euclidean", "arch_weights")[1],
                              colData_id = "cells", rank = F)
  saveRDS(data_attr, paste0("~/Fig4-5_archetype/2_data_attr/", i, ".Rds"))
  
  # Use Wilcox test to find genes maximally expressed in 10% closest to each vertex
  enriched_genes = find_decreasing_wilcox(data_attr$data, data_attr$arc_col,
                                          features = data_attr$features_col,
                                          bin_prop = 0.1, method = "BioQC")
  
  enriched_sets = find_decreasing_wilcox(data_attr$data, data_attr$arc_col,
                                         features = data_attr$colData_col,
                                         bin_prop = 0.1, method = "BioQC")
  
  # Take a look at top genes and functions for each archetype
  labs = get_top_decreasing(summary_genes = enriched_genes, summary_sets = enriched_sets, 
                            cutoff_genes = 0.01, cutoff_sets = 0.05, 
                            cutoff_metric = "wilcoxon_p_val", 
                            p.adjust.method = "fdr",
                            order_by = "mean_diff", order_decreasing = T,
                            min_max_diff_cutoff_g = 0.4, min_max_diff_cutoff_f = 0.03)
  
  saveRDS(list(gene = enriched_genes, set = enriched_sets, lab = labs), paste0("~/Fig4-5_archetype/2_enrichment/", i, ".Rds"))
  
  arc_cells <- bin_cells_by_arch(data_attr$data, data_attr$arc_col, return_names = T, bin_prop = 0.1)
  arc_ind <- unlist(lapply(names(arc_cells), function(x) structure(rep(x, length(arc_cells[[x]])), names = arc_cells[[x]])))
  
  
  srt$arc <- "no id"
  srt$arc[intersect(colnames(srt), names(arc_ind))] <- arc_ind[match(intersect(colnames(srt), names(arc_ind)), names(arc_ind))]
  Idents(srt) <- srt$arc
  pdf(paste0("~/Fig4-5_archetype/2_srt_plots/", i, ".pdf"), height = 4, width = 4)
  p <- DimPlot(srt, group.by = "arc", reduction = "umap", cols = c(hue_pal()(length(arc_cells)), "lightgrey"))
  print(p)
  for(n in names(arc_cells)){
    f <- DimPlot(srt, cells.highlight = arc_cells[[n]], reduction = "umap", sizes.highlight = 0.7)
    print(f)
  }
  dev.off()
  
  markers <- FindAllMarkers(srt, test.use = "MAST")
  saveRDS(markers, paste0("~/Fig4-5_archetype/2_marker/", i, ".Rds"))
  saveRDS(srt@meta.data,paste0("~/Fig4-5_archetype/2_meta/", i, ".Rds") )
  
  print(i)
}

# 
# dist = arch_dist(data = PCs4arch, arc_1$XC, dist_metric = "euclidean")
# temp <- as.data.frame(apply(dist, 1, function(x) names(which.min(x)))) %>% `colnames<-`("arc") %>% rownames_to_column("cell")
# temp <- split(temp, temp$arc)
# 
# temp <- lapply(temp , function(x){ 
#   x[,"dist"] <- dist[x[,"cell"], unique(x[,"arc"])]
#   return(x)
# })
# arc_cells <- lapply(temp, function(x) x$cell[which(ntile(x$dist, 5)%in% 1:4)])
# arc_ind <- unlist(lapply(names(arc_cells), function(x) structure(rep(x, length(arc_cells[[x]])), names = arc_cells[[x]])))
# srt$arc <- "no id"
# srt$arc[intersect(colnames(srt), names(arc_ind))] <- arc_ind[match(intersect(colnames(srt), names(arc_ind)), names(arc_ind))]
# 
# Idents(srt) <- srt$arc
# markers <- FindAllMarkers(srt)
# DimPlot(srt, cells.highlight = arc_cells$archetype_1, reduction = "umap")+
#   DimPlot(srt, cells.highlight = arc_cells$archetype_2, reduction = "umap")+
#   DimPlot(srt, cells.highlight = arc_cells$archetype_3, reduction = "umap")+
#   DimPlot(srt, cells.highlight = arc_cells$archetype_4, reduction = "umap")+
#   DimPlot(srt, cells.highlight = arc_cells$archetype_5, reduction = "umap")+
#   DimPlot(srt, cells.highlight = arc_cells$archetype_6, reduction = "umap")+
#   DimPlot(srt, cells.highlight = arc_cells$archetype_7, reduction = "umap")+
#   DimPlot(srt, cells.highlight = arc_cells$archetype_8, reduction = "umap")
# 
# 
# srt$group <- mta$cell_group_2
# DimPlot(srt, group.by = "arc", reduction = "umap", cols = c(hue_pal()(ncol(dist)), "lightgrey"))+
#   DimPlot(srt, group.by = "group", reduction = "umap")
# 
# library(plotly)
# srt <- RunUMAP(srt, dims = 1:20, n.components = 3)
# umap <- Embeddings(srt, "umap")
# plot_ly(as.data.frame(umap), x = ~umap_1, y = ~umap_2, z = ~umap_3, color = ~srt$seurat_clusters, size = 0.8 )
# 
# tb <- as.data.frame(table(markers$cluster, markers$gene))
# tb2 <- as.data.frame(table(marker$cluster, marker$gene))
# 
# markers  %>% filter(gene %in% intersect(marker$gene, markers$gene)) %>% filter(cluster != "no id") -> test
# 
# test <- test %>% filter(gene %in% names(which(table(test$gene) ==1)))
# 
# test2 <- marker %>% filter(gene %in% test$gene)%>% filter(gene %in% names(which(table(marker$gene) ==1)))
# 
# test <- test %>% filter(gene %in% test2$gene)
# testf <- merge(test, test2, by = "gene", all.x = T, all.y = F)
# 
# ggplot(testf,aes(x = avg_log2FC.x, y = avg_log2FC.y)) + geom_point() + geom_hline(yintercept = 0) + 
#   geom_vline(xintercept = 0)+
#   facet_grid(cluster.x ~ cluster.y) + xlim(c(-2.5, 2.5)) + ylim(c(-2.5, 2.5))