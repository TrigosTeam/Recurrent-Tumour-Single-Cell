library(Seurat)
library(dplyr)
library(harmony)
library(data.table)
library(SingleR)
library(patchwork)
library(ggplot2)
library(scran)
refs <- readRDS("~/Fig2_signature/sc_integration/tumor_normal_seperation/atlasref.Rds")
srt <- readRDS("~/Fig2_signature/sc_integration/tumor_normal_ruvsrt.Rds")
tumor_only_meta <-  readRDS("~/Fig4-5_archetype/5_integrated_final_module_meta.Rds")

DimPlot(srt, label = T, raster = F)
DimPlot(srt, cells = rownames(tumor_only_meta))
DimPlot(srt, cells = setdiff(colnames(srt), rownames(tumor_only_meta)), label = T, group.by = "Cell.Identity")


normal1 <- subset(srt, idents = as.character(c(20, 23, 22, 15, 29, 17, 32, 29,30)))
normal1 <- subset(normal1, cells = setdiff(colnames(normal1), rownames(tumor_only_meta)))

DimPlot(normal1)
rm(list = "srt")
gc()
saveRDS(normal1, "normal_cluster_srt.Rds")

f <- DimPlot(normal1, label = T)
clean_normal <- CellSelector(f)
umap <- as.data.frame(normal_cluster_srt@reductions$umap@cell.embeddings)

normal1$old_cluster <- normal1$seurat_clusters
normal1 <- FindVariableFeatures(normal1)
normal1 <- normal1 %>% 
  RunPCA(features = VariableFeatures(normal1)) %>% 
  FindNeighbors(dims = 1:20) %>% 
  FindClusters(resolution = 1) %>% 
  RunUMAP(dims = 1:20)

DimPlot(normal1, label = T)
DimPlot(normal1, label = T, group.by = "Cell.Identity") + NoLegend()
FeaturePlot(normal1, "AR")



# SingleR -----------------------------------------------------------------
pred <- SingleR(normal1[["ruv3"]]@data, refs, labels = list(refs[[1]]$label.main,refs[[2]]$label.main))
ids <- pred[colnames(normal1), "labels"]
ids <- gsub("_", " ", fixed = T, ids)
ids <- gsub("-", " ", fixed = T, ids)
ids <- gsub("s$", "", ids)
normal1[["Cell.Identity.new"]] <- ids
DimPlot(normal1, label = T, group.by = "Cell.Identity.new") + NoLegend()




# recluster with RNA ------------------------------------------------------
scATOMIC_result <- readRDS("~/normal_cell_annotation/scATOMIC_result.Rds")
head(scATOMIC_result)
normal1$scATOMIC_anno <- scATOMIC_result$scATOMIC_pred
DimPlot(normal1, label = F, group.by = "scATOMIC_anno")
DimPlot(normal1, group.by = "sample", label = T)
DimPlot(normal1, cells.highlight = colnames(normal1)[normal1$scATOMIC_anno == "PRAD"])
pdf("test.pdf", width = 50, height = 3)
DimPlot(normal1, split.by = "scATOMIC_anno")
dev.off()

normal_cluster_srt$anno <- scATOMIC_result$scATOMIC_pred
DimPlot(normal_cluster_srt,label = F, group.by = "anno")
DimPlot(normal_cluster_srt, cells.highlight = colnames(normal_cluster_srt)[normal_cluster_srt$anno == "Any Cell"])
DimPlot(normal_cluster_srt, cells= setdiff(colnames(normal_cluster_srt), colnames(normal_cluster_srt)[normal_cluster_srt$anno == "Any Cell"]))
DimPlot(normal_cluster_srt, cells.highlight = colnames(normal_cluster_srt)[normal_cluster_srt$anno == "PRAD"])

DimPlot(normal1, cells= setdiff(colnames(normal_cluster_srt), colnames(normal_cluster_srt)[normal_cluster_srt$anno%in% c("Any Cell", "PRAD")]))


normal2 <- subset(normal1, cells= setdiff(colnames(normal_cluster_srt), colnames(normal_cluster_srt)[normal_cluster_srt$anno%in% c("Any Cell", "PRAD")]))
DefaultAssay(normal2) <- "RNA"
lib.median <- median(normal2$nCount_RNA)
target_pseudocount <- 1
normal2 <- normal2 %>% 
  NormalizeData(scale.factor = lib.median / target_pseudocount) %>% 
  ScaleData(do.scale = F)

set.seed(100)
poisson_fit2 <- modelGeneVarByPoisson(normal2@assays$RNA@counts, dispersion = 0.2, size.factors = target_pseudocount * normal2$nCount_RNA / lib.median)
residuals2 <- poisson_fit2$total - poisson_fit2@metadata$trend(poisson_fit2$mean)
names(residuals2) <- rownames(normal2)
top_genes2 <- rownames(normal2)[order(residuals2, decreasing = TRUE)]

normal2 <- normal2 %>% 
  RunPCA(features = top_genes2[1:2000]) %>% 
  FindNeighbors(dims = 1:20) %>% 
  FindClusters(resolution = 1) %>% 
  RunUMAP(dims = 1:20)

DimPlot(normal2, label = T)
DimPlot(normal2, label =F, group.by = "Cell.Identity.new") 
DimPlot(normal2, label =F, group.by = "scATOMIC_anno") 

pred <- SingleR(normal2[["RNA"]]@data, refs, labels = list(refs[[1]]$label.main,refs[[2]]$label.main))
ids <- pred[colnames(normal2), "labels"]
ids <- gsub("_", " ", fixed = T, ids)
ids <- gsub("-", " ", fixed = T, ids)
ids <- gsub("s$", "", ids)
normal2[["Cell.Identity.new"]] <- ids

df <- normal2@meta.data
df <- df %>% group_by(seurat_clusters, Cell.Identity.new) %>% summarise(n = n()) %>% mutate(prop = n/sum(n))
df <- df%>% group_by(seurat_clusters) %>% filter(prop ==max(prop)) %>% mutate(anno = paste(Cell.Identity.new, collapse = "/"))
anno <- setNames(df$anno, df$seurat_clusters)
anno <- c(anno, c(`2` = "Fibroblast", `12` = "Fibroblast"))

normal2$cluster_anno <- anno[normal2$seurat_clusters]
DimPlot(normal2, group.by = "cluster_anno", cols = pals::brewer.set2(12))
DimPlot(normal2, group.by = "patient")
saveRDS(normal2, "normal_clean_srt.Rds")
# integration -------------------------------------------------------------
library(ruvIIInb)
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
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(harmony)

data(Hs.schk)

samples <- normal2$sample
temp <- structure(seq(length(unique(samples))), names = unique(samples))
batch <- temp[samples]
ctl <- as.character(Hs.schk)
ind<-rownames(normal2) %in% ctl

M <- matrix(0,ncol(normal2),length(unique(normal2$cluster_anno)))
colnames(M) <- unique(normal2$cluster_anno)
cl<- unique(normal2$cluster_anno)
for(CL in cl){
  M[which(normal2$cluster_anno==CL),CL] <- 1
}

ruv3nb_out<-try(ruvIIInb::fastruvIII.nb(Y=DelayedArray(as.array(normal2@assays$RNA@counts)), # count matrix with genes as rows and cells as columns
                                        M=M, #Replicate matrix constructed as above
                                        ctl=ind, #A vector denoting control genes
                                        k=2, # dimension of unwanted variation factors
                                        ncores = 1,
                                        batch = batch
))
#saveRDS(ruv3nb_out, paste0("~/integration/2024_06/tumor_normal_ruvout.Rds"))
logPAC <- as(log(get.res(ruv3nb_out,
                         type = "quantile", batch = ruv3nb_out$batch) + 1), "sparseMatrix")
rownames(logPAC) <- rownames(ruv3nb_out$counts)
colnames(logPAC) <- colnames(ruv3nb_out$counts)


normal2[["ruv3"]]<- CreateAssayObject(data = logPAC)
DefaultAssay(normal2) <-"ruv3"
normal2[["ruv3"]]@scale.data <- as.matrix(normal2[["ruv3"]]@data)
normal2 <- FindVariableFeatures(normal2)
normal2 <- normal2 %>%
  RunPCA(features = VariableFeatures(normal2)) %>%
  FindNeighbors(dims = 1:20, graph.name = c("ruv_nn", "ruv_snn")) %>%
  FindClusters(resolution = 1, graph.name = "ruv_snn") %>%
  RunUMAP(dims = 1:20, reduction.name = "ruv_umap")
saveRDS(normal2, "normal_ruvsrt.Rds")
DimPlot(normal2,label = T)
DimPlot(normal2, group.by = "cluster_anno", cols = pals::brewer.set2(12))
DimPlot(normal2, group.by = "patient")


# simspec integration ----------------------------------------------------------------
DefaultAssay(normal2) <- "RNA"
lib.median <- median(normal2$nCount_RNA)
target_pseudocount <- 1
normal2 <- normal2 %>% 
  NormalizeData(scale.factor = lib.median / target_pseudocount) %>% 
  ScaleData(do.scale = F)

set.seed(100)
poisson_fit2 <- modelGeneVarByPoisson(normal2@assays$RNA@counts, dispersion = 0.2, size.factors = target_pseudocount * normal2$nCount_RNA / lib.median)
residuals2 <- poisson_fit2$total - poisson_fit2@metadata$trend(poisson_fit2$mean)
names(residuals2) <- rownames(normal2)
top_genes2 <- rownames(normal2)[order(residuals2, decreasing = TRUE)]

normal2 <- normal2 %>% 
  RunPCA(features = top_genes2[1:2000]) %>% 
  FindNeighbors(dims = 1:20) %>% 
  FindClusters(resolution = 1) %>% 
  RunUMAP(dims = 1:20)

normal2 <- cluster_sim_spectrum(object = normal2, label_tag = "sample")
normal2 <- FindNeighbors(normal2, reduction = "css", dims = 1:ncol(Embeddings(normal2, "css")), graph.name = c("css_nn", "css_snn" ))
normal2 <- FindClusters(normal2, resolution = 1, graph.name = "css_snn")
normal2 <- RunUMAP(normal2, reduction = "css", dims = 1:ncol(Embeddings(normal2, "css")), reduction.name = "css_umap")
DimPlot(normal2, reduction = "css_umap", label = T)
DimPlot(normal2, group.by = "cluster_anno", cols = pals::brewer.set2(12), reduction = "css_umap", label = T)
DimPlot(normal2, group.by = "patient", reduction = "css_umap", label = T)

normal2 <- RunHarmony(normal2, "sample") %>% 
  RunUMAP(reduction = "harmony", dims = 1:30,reduction.name = "harmony_umap") %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30, graph.name = c("harmony_nn","harmony_snn")) %>% 
  FindClusters(resolution = 1, graph.name = "harmony_snn")

DimPlot(normal2, reduction = "harmony_umap", label = T)
DimPlot(normal2, group.by = "cluster_anno", cols = pals::brewer.set2(12), reduction = "harmony_umap", label = T)
DimPlot(normal2, group.by = "patient", reduction = "harmony_umap", label = T)
saveRDS(normal2, "normal_integratedsrt.Rds")


DotPlot(normal2, features = all_marker)+ theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
DotPlot(normal2, features = all_marker, group.by = "cell_id2")+ theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Azimuth subclassification -----------------------------------------------
# if (!requireNamespace('remotes', quietly = TRUE)) {
#   install.packages('remotes')
# }
remotes::install_github('satijalab/azimuth', ref = 'v0.3.2') ## should be compile with Seurat v4
.libPaths("/trigos_team/R_Library/4.2/Seuratv5/")
remotes::install_github('satijalab/azimuth', ref = 'master')
normal_integratedsrt <- readRDS("~/normal_cell_annotation/normal_integratedsrt.Rds")
DotPlot(normal_integratedsrt, features = all_marker)+ theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
DotPlot(normal_integratedsrt, features = all_marker, group.by = "cluster_anno")+ theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

library(Azimuth)
library(patchwork)
library(pals)

preprocess_srt <- function(normal2){
  DefaultAssay(normal2) <- "RNA"
  lib.median <- median(normal2$nCount_RNA)
  target_pseudocount <- 1
  normal2 <- normal2 %>% 
    NormalizeData(scale.factor = lib.median / target_pseudocount) %>% 
    ScaleData(do.scale = F)
  
  set.seed(100)
  poisson_fit2 <- modelGeneVarByPoisson(normal2@assays$RNA@counts, dispersion = 0.2, size.factors = target_pseudocount * normal2$nCount_RNA / lib.median)
  residuals2 <- poisson_fit2$total - poisson_fit2@metadata$trend(poisson_fit2$mean)
  names(residuals2) <- rownames(normal2)
  top_genes2 <- rownames(normal2)[order(residuals2, decreasing = TRUE)]
  normal2 <- normal2 %>% 
    NormalizeData(scale.factor = lib.median / target_pseudocount) %>% 
    ScaleData(do.scale = F)
  
  set.seed(100)
  poisson_fit2 <- modelGeneVarByPoisson(normal2@assays$RNA@counts, dispersion = 0.2, size.factors = target_pseudocount * normal2$nCount_RNA / lib.median)
  residuals2 <- poisson_fit2$total - poisson_fit2@metadata$trend(poisson_fit2$mean)
  names(residuals2) <- rownames(normal2)
  top_genes2 <- rownames(normal2)[order(residuals2, decreasing = TRUE)]
  
  normal2 <- normal2 %>% 
    RunPCA(features = top_genes2[1:2000]) %>% 
    FindNeighbors(dims = 1:20) %>% 
    FindClusters(resolution = 1) %>% 
    RunUMAP(dims = 1:20)
  normal2$RNA_cluster <- normal2$seurat_clusters
  
  normal2 <- cluster_sim_spectrum(object = normal2, label_tag = "sample")
  normal2 <- FindNeighbors(normal2, reduction = "css", dims = 1:ncol(Embeddings(normal2, "css")), graph.name = c("css_nn", "css_snn" ))
  normal2 <- FindClusters(normal2, resolution = 1, graph.name = "css_snn")
  normal2 <- RunUMAP(normal2, reduction = "css", dims = 1:ncol(Embeddings(normal2, "css")), reduction.name = "css_umap")
  normal2$css_cluster <- normal2$seurat_clusters
  
  normal2 <- RunHarmony(normal2, "sample") %>% 
    RunUMAP(reduction = "harmony", dims = 1:30,reduction.name = "harmony_umap") %>% 
    FindNeighbors(reduction = "harmony", dims = 1:30, graph.name = c("harmony_nn","harmony_snn")) %>% 
    FindClusters(resolution = 1, graph.name = "harmony_snn")
  normal2$harmony_cluster <- normal2$seurat_clusters
  Idents(normal2) <- normal2$css_cluster
  
  return(normal2)
}


f <- DimPlot(normal_integratedsrt, reduction = "css_umap", label = T)
DimPlot(normal_integratedsrt, reduction = "css_umap", label = T)+
DimPlot(normal_integratedsrt, group.by = "cluster_anno", cols = pals::brewer.set2(12), reduction = "css_umap", label = T)
DimPlot(normal_integratedsrt, group.by = "patient", reduction = "css_umap", label = T)

immunes_cell <- CellSelector(f)
immunes <- subset(normal_integratedsrt, cells = immunes_cell)
DimPlot(immunes, reduction = "harmony_umap")

p <- DimPlot(immunes, reduction = "ruv_umap")
junk <- CellSelector(p)
DimPlot(immunes, reduction = "css_umap")

immunes <- preprocess_srt(immunes)

DimPlot(immunes, reduction = "css_umap", cells = junk,group.by = "cluster_anno")+DimPlot(immunes, reduction = "css_umap", group.by = "cluster_anno")
DimPlot(immunes, reduction = "css_umap", label = T)+DimPlot(immunes, reduction = "css_umap", group.by = "cluster_anno")



predictions <- read.delim('azimuth_pred.tsv', row.names = 1)
immunes <- AddMetaData(
  object = immunes,
  metadata = predictions)
DimPlot(immunes, reduction = "css_umap", group.by = "predicted.celltype.l2")

DimPlot(immunes, reduction = "css_umap", cells = rownames(predictions)[predictions$predicted.celltype.l2 %in% names(which(table(predictions$predicted.celltype.l2)>10))])
df <- immunes@meta.data
df <- df %>% group_by(seurat_clusters, predicted.celltype.l2) %>% summarise(n = n()) %>% mutate(prop = n/sum(n))
df2 <- df%>% group_by(seurat_clusters) %>% filter(prop >0.4) %>% mutate(anno = paste(predicted.celltype.l2, collapse = "/"))
anno <- setNames(df2$anno, df2$seurat_clusters)
immunes$azimuth_cluster <- as.character(anno[as.character(immunes$seurat_clusters)])
DimPlot(immunes, reduction = "css_umap",label = T) +DimPlot(immunes, reduction = "css_umap", group.by = "azimuth_cluster")

Bcells <- subset(immunes, idents = c("1", "8", "9"))
Bcells <- preprocess_srt(Bcells)
Tcells <- subset(immunes, idents = as.character(c(0, 2, 6, 10)))
Tcells <- preprocess_srt(Tcells)
Macro <- subset(immunes, idents = as.character(c(3, 4, 5, 7, 11)))
Macro <- preprocess_srt(Macro)

saveRDS(Bcells, "~/normal_cell_annotation/subtype/Bcells.Rds")
saveRDS(Tcells, "~/normal_cell_annotation/subtype/Tcells.Rds")
saveRDS(Macro, "~/normal_cell_annotation/subtype/Macrophage.Rds")

stromal <- subset(normal_integratedsrt, cells = setdiff(colnames(normal_integratedsrt), immunes_cell))
saveRDS(stromal@assays$RNA@counts, "stromal_countmat.Rds")
stromal$old_cluster <- stromal$css_cluster
stromal <- preprocess_srt(stromal)
DimPlot(stromal, reduction = "css_umap", label = T, group.by = "css_cluster")+DimPlot(stromal, reduction = "css_umap", group.by = "cluster_anno", cols = brewer.set2(10))


predictions <- read.delim('azimuth_pred_stromal.tsv', row.names = 1)
stromal <- AddMetaData(
  object = stromal,
  metadata = predictions)
DimPlot(stromal, reduction = "css_umap", group.by = "predicted.celltype.l2") # liver ref
DimPlot(stromal, reduction = "css_umap", group.by = "predicted.celltype.l1", cols = as.character(brewer.set3(15))) #  ref

df <- stromal@meta.data
df <- df %>% group_by(seurat_clusters, predicted.celltype.l1) %>% summarise(n = n()) %>% mutate(prop = n/sum(n))
df2 <- df%>% group_by(seurat_clusters) %>% filter(prop >0.4) %>% mutate(anno = paste(predicted.celltype.l1, collapse = "/"))
anno <- setNames(df2$anno, df2$seurat_clusters)
stromal$azimuth_cluster <- as.character(anno[as.character(stromal$seurat_clusters)])
DimPlot(stromal, reduction = "css_umap",label = T) +DimPlot(stromal, reduction = "css_umap", group.by = "azimuth_cluster")
Idents(stromal) <- stromal$css_cluster
Epi <- preprocess_srt(subset(stromal, idents = as.character(c(2, 6))))
Endo <- preprocess_srt(subset(stromal, idents = as.character(c(14, 7, 1, 16))))
Fibro <- preprocess_srt(subset(stromal, idents = as.character(c(8, 9, 5, 4, 13, 0, 3, 12, 10))))
Neuro <- subset(stromal, idents = as.character(c(15)))
Hepo <- preprocess_srt(subset(stromal, idents = as.character(c(11))))


saveRDS(Epi, "~/normal_cell_annotation/subtype/Epi.Rds")
saveRDS(Endo, "~/normal_cell_annotation/subtype/Endo1.Rds")
saveRDS(Fibro, "~/normal_cell_annotation/subtype/Fibro.Rds")
saveRDS(Neuro, "~/normal_cell_annotation/subtype/Neuro.Rds")
saveRDS(Hepo, "~/normal_cell_annotation/subtype/Hepo.Rds")






