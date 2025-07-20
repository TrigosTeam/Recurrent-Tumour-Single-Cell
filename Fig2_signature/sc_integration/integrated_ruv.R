library(ruvIIInb)
library(SingleCellExperiment)
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
library(simspec)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(scales)
library(pals)

data(Hs.schk)
# # cycling cell ------------------------------------------------------------
# cycling_cells <- readRDS("~/regulon/cycling_srt/cycling_cells.Rds")
# cycling_cells <- lapply(names(cycling_cells), function(x){
#   temp <- paste(x, cycling_cells[[x]], sep = "_")
#   temp
# })
# cycling_cells <- unlist(cycling_cells)


# tumor only --------------------------------------------------------------

test <- readRDS("~/Fig2_signature/sc_integration/tumor_only_RNA_simspec_srt.Rds")
samples <- test$sample
temp <- structure(seq(length(unique(samples))), names = unique(samples))
batch <- temp[samples]
ctl <- as.character(Hs.schk)
ind<-rownames(test) %in% ctl

M <- matrix(0,ncol(test),length(unique(test$seurat_clusters)))
cl<- sort(unique(as.numeric(unique(test$seurat_clusters))))
for(CL in cl){
  M[which(as.numeric(test$seurat_clusters)==CL),CL] <- 1
}

ruv3nb_out<-try(ruvIIInb::fastruvIII.nb(Y=DelayedArray(as.array(test@assays$RNA@counts)), # count matrix with genes as rows and cells as columns
                                        M=M, #Replicate matrix constructed as above
                                        ctl=ind, #A vector denoting control genes
                                        k=2, # dimension of unwanted variation factors
                                        ncores = 1,
                                        batch = batch
))

#saveRDS(ruv3nb_out, "~/Fig2_signature/sc_integration/tumor_only_ruvout.Rds")
print("tumor only ruv finished")

logPAC <- as(log(get.res(ruv3nb_out,
                         type = "quantile", batch = ruv3nb_out$batch) + 1), "sparseMatrix")
rownames(logPAC) <- rownames(ruv3nb_out$counts)
colnames(logPAC) <- colnames(ruv3nb_out$counts)

saveRDS(logPAC, "~/Fig2_signature/sc_integration/tumor_only_logPAC.Rds")

rm( list = "ruv3nb_out")
gc()


test[["ruv3"]]<- CreateAssayObject(data = logPAC)
DefaultAssay(test) <-"ruv3"
test[["ruv3"]]@scale.data <- as.matrix(test[["ruv3"]]@data)
test <- FindVariableFeatures(test)
test <- test %>%
  RunPCA(features = VariableFeatures(test)) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 1) %>%
  RunUMAP(dims = 1:20)

sites<- test$sample
sites[grep("brain", sites)] <- "brain"
sites[grep("dura", sites)] <- "dura"
sites[grep("prostate", sites)] <- "prostate"
sites[grep("liver", sites)] <- "liver"
sites[grep("fat", sites)] <- "fat"
sites[grep("lymph|hilar", sites)] <- "LN"
sites[grep("rib|vertebra", sites)] <- "bone"
sites[grep("lung", sites)] <- "lung"
sites[grep("abdomen", sites)] <- "abdomen"
sites[grep("bladder", sites)] <- "bladder"
test$site <- sites
saveRDS(test, "~/Fig2_signature/sc_integration/tumor_only_ruvsrt.Rds")

pdf("~/Fig2_signature/sc_integration/tumor_only_ruv_srtplots.pdf", width = 12, height = 10)
DimPlot(test, label = T, raster = F)
DimPlot(test, group.by = "patient", raster = F, label = T, repel = T) +NoLegend()
DimPlot(test, group.by = "sample", raster = F, label = T, repel = T)+NoLegend()
DimPlot(test, group.by = "Cell.Identity", raster = F, label = T)+NoLegend()
DimPlot(test, group.by = "site", raster = F, label = T)+NoLegend()
FeaturePlot(test, "AR", cols = rev(brewer.pal(10, "RdYlBu")), raster = F)
FeaturePlot(test, "ASCL1", cols = rev(brewer.pal(10, "RdYlBu")), raster = F)
FeaturePlot(test, "FOLH1", cols = rev(brewer.pal(10, "RdYlBu")), raster = F)
dev.off()

rm(list = ls())
gc()


# tumor normal ------------------------------------------------------------
#15/07/2024
test <- readRDS("~/Fig2_signature/sc_integration/tumor_normal_RNA_simspec_srt.Rds")
final_normalsrt <- readRDS("~/normal_cell_annotation/subtype/final_normalsrt.Rds")
tumor_only_meta <-  readRDS("~/Fig4-5_archetype/5_integrated_final_module_meta.Rds")
test <- subset(test, cells = c(colnames(final_normalsrt), rownames(tumor_only_meta)))
DefaultAssay(test) <- "RNA"
lib.gmean3 <- exp(mean(log(colSums(test[["RNA"]]@counts))))
target_pseudocount = 1
test <- test %>%
  NormalizeData(scale.factor = lib.gmean3 / target_pseudocount) %>%
  ScaleData(do.scale = FALSE)
set.seed(100)
poisson_fit3 <- modelGeneVarByPoisson(test[["RNA"]]@counts, dispersion = 0.1, size.factors = colSums(test[["RNA"]]@counts)/lib.gmean3)
residuals3 <- poisson_fit3$total - poisson_fit3@metadata$trend(poisson_fit3$mean)
names(residuals3) <- rownames(test)
top_genes3 <- rownames(test)[order(residuals3, decreasing = TRUE)]
test <- test %>%
  RunPCA(features = top_genes3[1:2000]) %>%
  FindVariableFeatures(nfeatures = 2000)
test <- cluster_sim_spectrum(object = test, label_tag = "sample")
test <- RunUMAP(test, reduction = "css", dims = 1:ncol(Embeddings(test, "css")))
test <- FindNeighbors(test, reduction = "css", dims = 1:ncol(Embeddings(test, "css")))
test <- FindClusters(test, resolution = 0.5)
test$patient <- substr(test$sample, 1, 6)
test$cell_anno <- test$patient
test$cell_anno[colnames(final_normalsrt)] <- as.character(final_normalsrt$cell_anno)
test$cell_anno <- factor(test$cell_anno, levels = unique(test$cell_anno))
saveRDS(test, "~/Fig2_signature/sc_integration/tumor_normal_RNA_simspec_newsrt.Rds")
# saveRDS(test, "/trigos_team/CASCADE/Analysis/integration/objects/tumor_merged_27_34_35_43_46_58_76_83_90_srt.Rds")
# test <- readRDS("/trigos_team/CASCADE/Analysis/integration/objects/tumor_merged_27_34_35_43_46_58_76_83_90_srt.Rds")
pdf( "~/Fig2_signature/sc_integration/tumor_normal_RNA_simspec_newUMAP.pdf", width = 20, height = 10)
f <- DimPlot(test, group.by = "sample", label = T)
print(f)
f <- DimPlot(test, label = T)
print(f)
f <- DimPlot(test, group.by = "patient", label = T)
print(f)
f <- DimPlot(test, group.by = "patient")
print(f)

f <- DimPlot(test, group.by = "cell_anno", cols = c(brewer.set1(14), hue_pal()(9)))
print(f)


f <- DimPlot(test, group.by = "cell_anno", label = T, repel = T, cols = c(brewer.set1(14), hue_pal()(9)))
print(f)

dev.off()




samples <- test$sample
temp <- structure(seq(length(unique(samples))), names = unique(samples))
batch <- temp[samples]
ctl <- as.character(Hs.schk)
ind<-rownames(test) %in% ctl

M <- matrix(0,ncol(test),length(unique(test$seurat_clusters)))
cl<- sort(unique(as.numeric(unique(test$seurat_clusters))))
for(CL in cl){
  M[which(as.numeric(test$seurat_clusters)==CL),CL] <- 1
}

ruv3nb_out<-try(ruvIIInb::fastruvIII.nb(Y=DelayedArray(as.array(test@assays$RNA@counts)), # count matrix with genes as rows and cells as columns
                                        M=M, #Replicate matrix constructed as above
                                        ctl=ind, #A vector denoting control genes
                                        k=2, # dimension of unwanted variation factors
                                        ncores = 1,
                                        batch = batch
))
#saveRDS(ruv3nb_out, paste0("~/Fig2_signature/sc_integration/tumor_normal_ruvout.Rds"))
print("tumor normal ruv finished")

logPAC <- as(log(get.res(ruv3nb_out,
                         type = "quantile", batch = ruv3nb_out$batch) + 1), "sparseMatrix")
rownames(logPAC) <- rownames(ruv3nb_out$counts)
colnames(logPAC) <- colnames(ruv3nb_out$counts)

saveRDS(logPAC, "~/Fig2_signature/sc_integration/tumor_normal_logPAC.Rds")

rm( list = "ruv3nb_out")
gc()

logPAC <- readRDS( "~/Fig2_signature/sc_integration/tumor_normal_logPAC.Rds")
test[["ruv3"]]<- CreateAssayObject(data = logPAC)
DefaultAssay(test) <-"ruv3"
test[["ruv3"]]@scale.data <- as.matrix(test[["ruv3"]]@data)
test <- FindVariableFeatures(test)
test <- test %>%
  RunPCA(features = VariableFeatures(test)) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 1) %>%
  RunUMAP(dims = 1:20)
sites<- test$sample
sites[grep("brain", sites)] <- "brain"
sites[grep("dura", sites)] <- "dura"
sites[grep("prostate", sites)] <- "prostate"
sites[grep("liver", sites)] <- "liver"
sites[grep("fat", sites)] <- "fat"
sites[grep("lymph|hilar", sites)] <- "LN"
sites[grep("rib|vertebra", sites)] <- "bone"
sites[grep("lung", sites)] <- "lung"
sites[grep("abdomen", sites)] <- "abdomen"
sites[grep("bladder", sites)] <- "bladder"
test$site <- sites
saveRDS(test, "~/Fig2_signature/sc_integration/tumor_normal_ruvsrt.Rds")

pdf("~/Fig2_signature/sc_integration/tumor_normal_ruv_srtplots.pdf", width = 12, height = 10)
DimPlot(test, label = T, raster = F)
DimPlot(test, group.by = "patient", raster = F, label = T, repel = T) +NoLegend()
DimPlot(test, group.by = "patient", raster = F) 
DimPlot(test, group.by = "sample", raster = F, label = T, repel = T)+NoLegend()
DimPlot(test, group.by = "Cell.Identity", raster = F, label = T)+NoLegend()
DimPlot(test, group.by = "site", raster = F, label = T)+NoLegend()
DimPlot(test, group.by = "cell_anno", raster = F, label = T, cols = c(brewer.set1(14), hue_pal()(9)))+NoLegend()
DimPlot(test, group.by = "cell_anno", raster = F, cols = c(brewer.set1(14), hue_pal()(9))) 
FeaturePlot(test, "AR", raster = F)
FeaturePlot(test, "FOLH1", raster = F)
FeaturePlot(test, "ASCL1", raster = F)
FeaturePlot(test, "AR", cols = rev(brewer.pal(10, "RdYlBu")), raster = F)
FeaturePlot(test, "ASCL1", cols = rev(brewer.pal(10, "RdYlBu")), raster = F)
FeaturePlot(test, "FOLH1", cols = rev(brewer.pal(10, "RdYlBu")), raster = F)
dev.off()
