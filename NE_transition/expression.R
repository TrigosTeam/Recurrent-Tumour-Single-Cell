library(Seurat)
library(ggsci)
library(patchwork)
library(ggplot2)
library(pals)
source("~/CASCADEpaper/paper/functions.R", echo=TRUE)
#### per-patient expression ---------
CA0046_ruvsrt <- readRDS("~/integration/2024_06/perpatient/tumor_only/CA0046_ruvsrt.Rds")
CA0090_ruvsrt <- readRDS("~/integration/2024_06/perpatient/tumor_only/CA0090_ruvsrt.Rds")
CA0058_ruvsrt <- readRDS("~/integration/2024_06/perpatient/tumor_only/CA0058_ruvsrt.Rds")
CA0027_ruvsrt <- readRDS("~/integration/2024_06/perpatient/tumor_only/CA0027_ruvsrt.Rds")

clean_module <- readRDS("~/CASCADEpaper/paper/Fig5_archetype/clean_module.Rds")
names(clean_module)

CA0046_ruvsrt <- AddModuleScore(CA0046_ruvsrt, clean_module)
CA0090_ruvsrt <- AddModuleScore(CA0090_ruvsrt, clean_module)
CA0058_ruvsrt <- AddModuleScore(CA0058_ruvsrt, clean_module)
CA0027_ruvsrt <- AddModuleScore(CA0027_ruvsrt, clean_module)



plotall <- function(srt, patient){
  srt$labs <- gsub(paste0(patient, "_"), "", srt$sample)
  f1 <- rmaxes(DimPlot(srt, group.by = "labs", cols = ggsci::pal_bmj()(6)) +
    ggtitle(gsub("00", "", patient)) + NoLegend())
  f2 <- rmaxes(FeaturePlot(srt, "Cluster3", cols = rev(brewer.rdylbu(10)))+NoLegend()+ggtitle("NE1"))
  f3 <- coneraxes(FeaturePlot(srt, "Cluster4", cols = rev(brewer.rdylbu(10)))+NoLegend()+ggtitle("NE2"))
  p <- f1+f2+f3+plot_layout(nrow = 3)
  return(p)
}

p1 <- plotall(CA0046_ruvsrt, "CA0046")
plotall <- function(srt, patient){
  srt$labs <- gsub(paste0(patient, "_"), "", srt$sample)
  f1 <- DimPlot(srt, group.by = "labs", cols = ggsci::pal_bmj()(6)) +
    ggtitle(gsub("00", "", patient)) + NoLegend()+ 
    theme(axis.title = element_blank(), text = element_text(size = 9),plot.margin = margin(t = 0, l = 0.5, r = 0, b = 0.2, "cm"))+
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL)
  f2 <- FeaturePlot(srt, "Cluster3", cols = rev(brewer.rdylbu(10)))+NoLegend()+ggtitle("NE1")+
    theme(axis.title = element_blank(), text = element_text(size = 9),plot.margin = margin(t = 0, l = 0.5, r = 0, b = 0.2, "cm"))+
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL)
  f3 <- FeaturePlot(srt, "Cluster4", cols = rev(brewer.rdylbu(10)))+NoLegend()+ggtitle("NE2")+
    theme(axis.title = element_blank(), text = element_text(size = 9),plot.margin = margin(t = 0, l = 0.5, r = 0, b = 0.2, "cm"))+
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL)
  p <- f1+f2+f3+plot_layout(nrow = 3)
  return(p)
}
p2 <- plotall(CA0090_ruvsrt, "CA0090")
p3 <- plotall(CA0058_ruvsrt, "CA0058")
p4 <- plotall(CA0027_ruvsrt, "CA0027")

p <- p1|p2|p3|p4
ggsave(p, filename =  "moduleexpression.pdf", width = 14, height = 9,   bg = "transparent" )

"#2E2A2BFF" "#CF4E9CFF" "#8C57A2FF" "#358DB9FF" "#82581FFF" "#2F509EFF"
"#2A6EBBFF" "#F0AB00FF" "#C50084FF" "#7D5CC6FF" "#E37222FF" "#69BE28FF"

######   primiary tumour CA90 expression ------------
unloadNamespace("Seurat")
unloadNamespace("SeuratObject")
library("SeuratObject", lib.loc = "/trigos_team/R_Library/4.2/Seuratv5")
library("Seurat", lib.loc = "/trigos_team/R_Library/4.2/Seuratv5")
library(dplyr)
library(scran)
library(ggplot2)
CA90_13 <- readRDS("/trigos_team/CASCADE/Analysis/240731_seurat_flex/CA0090/p13/CA0090_p13_srt_new_filtering_doublets_non_tumour_removed.Rds")
meta <- CA90_13@meta.data
data <- CA90_13[["RNA"]]$data
scale.data <- CA90_13[["RNA"]]$scale.data
umap <- CA90_13@reductions$umap
options(Seurat.object.assay.version = 'v3')
srt3 <- CreateSeuratObject(counts = CA90_13[["RNA"]]$counts, meta.data = meta)

lib.median <- median(srt3$nCount_RNA)
target_pseudocount <- 1
srt3 <- srt3 %>% 
  NormalizeData(scale.factor = lib.median / target_pseudocount) %>% 
  ScaleData(do.scale = F, scale.max = Inf)
set.seed(100)
srt3 <- srt3 %>% 
  RunPCA(features = rownames(srt3), weight.by.var = T) %>% 
  FindNeighbors(dims = 1:30) %>% 
  FindClusters() %>% 
  RunUMAP(dims = 1:30)

DimPlot(srt3)
srt3 <- AddModuleScore(srt3, clean_module)
f1 <- FeaturePlot(srt3, "Cluster3")+ggtitle("NE1")+NoLegend()+scale_color_gradientn(colors = brewer.oranges(10))
f2 <- FeaturePlot(srt3, "Cluster4")+ggtitle("NE2")+NoLegend()+scale_color_gradientn(colors = brewer.oranges(10))
f3 <- FeaturePlot(srt3, "AR")+NoLegend()
f4 <- FeaturePlot(srt3, "ASCL1")+NoLegend()
f5 <- DimPlot(srt3, label = T)+NoLegend()
colnames(srt3@meta.data)[colnames(srt3@meta.data) %in%  c("Cluster3", "Cluster4")] <- c("NE1", "NE2")
v1 <- VlnPlot(srt3, features = c("NE1", "NE2"))
coneraxes(f3)+rmaxes(f4)+rmaxes(f1)+rmaxes(f2)+rmaxes(f5)+v1 + 
  plot_annotation(title = "CA90 primary p13")+plot_layout(nrow=1, ncol= 6, widths = c(1, 1, 1, 1, 1,3))
ggsave(filename = "primaryCA90p13.pdf", width = 18, height = 4, bg = "transparent")

FeaturePlot(srt3, c("Cluster3", "Cluster4"), blend = T)

saveRDS(srt3, "/trigos_team/CASCADE/Analysis/240731_seurat_flex/CA0090/p13/CA0090_p13_srt_new_filtering_doublets_non_tumour_removed_v4.Rds")
srt3 <- readRDS("/trigos_team/CASCADE/Analysis/240731_seurat_flex/CA0090/p13/CA0090_p13_srt_new_filtering_doublets_non_tumour_removed_v4.Rds")
normal <- readRDS("~/CASCADEpaper/paper/normal_cells_202406/subtype/final_normalsrt.Rds")
normal[["RNA"]] <- normal[["RNA2"]]
tumor_normal <- merge(srt3, normal)
counts <- tumor_normal[["RNA"]]@counts
tumor_normal$cell_anno[is.na(tumor_normal$cell_anno)] <- tumor_normal$seurat_clusters[is.na(tumor_normal$cell_anno)]
anno <- as.data.frame(tumor_normal$cell_anno)
colnames(anno) <- NULL
saveRDS(counts, paste0("infercnv/p13/CA90_p13_counts.Rds"))
saveRDS(anno, paste0("infercnv/p13/CA90_p13_anno.Rds"))

CA90_18 <- readRDS("/trigos_team/CASCADE/Analysis/240731_seurat_flex/CA0090/p18/CA0090_p18_srt_new_filtering_doublets_non_tumour_removed.Rds")
meta <- CA90_18@meta.data
data <- CA90_18[["RNA"]]$data
scale.data <- CA90_18[["RNA"]]$scale.data
umap <- CA90_18@reductions$umap
options(Seurat.object.assay.version = 'v3')
srt3 <- CreateSeuratObject(counts = CA90_18[["RNA"]]$counts, meta.data = meta)

lib.median <- median(srt3$nCount_RNA)
target_pseudocount <- 1
srt3 <- srt3 %>% 
  NormalizeData(scale.factor = lib.median / target_pseudocount) %>% 
  ScaleData(do.scale = F, scale.max = Inf)
set.seed(100)
srt3 <- srt3 %>% 
  RunPCA(features = rownames(srt3), weight.by.var = T) %>% 
  FindNeighbors(dims = 1:30) %>% 
  FindClusters() %>% 
  RunUMAP(dims = 1:30)

DimPlot(srt3)
srt3 <- AddModuleScore(srt3, clean_module)
f1 <- FeaturePlot(srt3, "Cluster3")+ggtitle("NE1")+NoLegend()+scale_color_gradientn(colors = brewer.oranges(10))
f2 <- FeaturePlot(srt3, "Cluster4")+ggtitle("NE2")+NoLegend()+scale_color_gradientn(colors = brewer.oranges(10))
f3 <- FeaturePlot(srt3, "AR")+NoLegend()
f4 <- FeaturePlot(srt3, "ASCL1")+NoLegend()
f5 <- DimPlot(srt3, label = T)+NoLegend()
colnames(srt3@meta.data)[colnames(srt3@meta.data) %in%  c("Cluster3", "Cluster4")] <- c("NE1", "NE2")
v1 <- VlnPlot(srt3, features = c("NE1", "NE2"))
coneraxes(f3)+rmaxes(f4)+rmaxes(f1)+rmaxes(f2)+rmaxes(f5)+v1 + 
  plot_annotation(title = "CA90 primary p18")+plot_layout(nrow=1, ncol= 6, widths = c(1, 1, 1, 1, 1,3))

ggsave(filename = "primaryCA90p18.pdf", width = 18, height = 4, bg = "transparent")

FeaturePlot(srt3, c("Cluster3", "Cluster4"), blend = T)

saveRDS(srt3, "/trigos_team/CASCADE/Analysis/240731_seurat_flex/CA0090/p18/CA0090_p18_srt_new_filtering_doublets_non_tumour_removed_v4.Rds")
normal <- readRDS("~/CASCADEpaper/paper/normal_cells_202406/subtype/final_normalsrt.Rds")
normal[["RNA"]] <- normal[["RNA2"]]
tumor_normal <- merge(srt3, normal)
counts <- tumor_normal[["RNA"]]@counts
tumor_normal$cell_anno[is.na(tumor_normal$cell_anno)] <- tumor_normal$seurat_clusters[is.na(tumor_normal$cell_anno)]
anno <- as.data.frame(tumor_normal$cell_anno)
colnames(anno) <- NULL
saveRDS(counts, paste0("infercnv/p18/CA90_p18_counts.Rds"))
saveRDS(anno, paste0("infercnv/p18/CA90_p18_anno.Rds"))

#### other NE marker expression ----------------
library(Seurat)
srt <- readRDS("~/integration/2024_06/tumor_only_ruvsrt.Rds")
DimPlot(srt, group.by = "patient", label = T)
FeaturePlot(srt, c("MYCN", "NEUROD1", "NFACT2", "ONECUT2", "OLIG3", "SOX2", "SOX9", 'SNAIL3', "ASCL2", "POU3F2"))+
VlnPlot(srt, "ONECUT2", group.by = "patient")
ggsave("expression of ONECUT2.pdf", height = 7, width = 5)
