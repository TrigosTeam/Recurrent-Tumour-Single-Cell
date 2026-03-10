library(Seurat)
library(Signac)
library(qs)
library(dplyr)
library(ggplot2)
library(pals)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(data.table)
library(JASPAR2024)
library(TFBSTools)
source("~/ATAC/functions.R")
setwd("~/ATAC")
## get individual sample ids and cluster based on TF only res = 0.5
clean_module <- readRDS("~/CASCADEpaper/paper/Fig5_archetype/clean_module.Rds")
names(clean_module) <- c("AR", "Inflammation", "NE1", "NE2", "Cycling", "Glycolysis")
# plan("multicore", workers = 4)
# options(future.globals.maxSize = 8000 * 1024^2)

module_col <- setNames( c("dodgerblue3", brewer.dark2(6)[2:6]),names(clean_module))
module_col <- c(module_col, setNames(c("grey90", "grey30", "red"), c("Background","Transition", "NE1&NE2")))
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
JASPAR2024 <- JASPAR2024()
pfm <- getMatrixSet(
  x = JASPAR2024@db,
  opts = list(collection = "CORE", species ='Homo sapiens')
)
# 
# combined_obj <- qread("integration/ruv_srt.qs")
# RNAmat <- readRDS("~/integration/2024_06/tumor_only_logPAC.Rds")
# combined_obj[["RNA"]] <- CreateAssayObject(data = RNAmat[, intersect(colnames(combined_obj), colnames(RNAmat))])
# DefaultAssay(combined_obj) <- "RNA"
# combined_obj[["RNA"]]@scale.data <- as.matrix(combined_obj[["RNA"]]@data)
# combined_obj <- FindVariableFeatures(combined_obj)
# combined_obj <- RunPCA(combined_obj,npc = 30)  
# 
# 
# 
# combined_obj <- FindMultiModalNeighbors(
#   object = combined_obj,
#   reduction.list = list("pca", "lsi.ruv3"), 
#   dims.list = list(1:20, 2:30),
#   modality.weight.name = "RNA.weight",
#   verbose = TRUE
# )
# 
# # build a joint UMAP visualization
# combined_obj <- RunUMAP(
#   object = combined_obj,
#   nn.name = "weighted.nn",
#   assay = "RNA",
#   reduction.name = "weighted.umap",
#   verbose = TRUE
# )
# # find clusters based on combined weighted graph of RNA and ATAC
# combined_obj <- FindClusters(combined_obj, graph.name = "wsnn")
# 
# f2 <- DimPlot(combined_obj, group.by = "patient", reduction = "weighted.umap", label = T, raster = F)
# f2 <-  coneraxes(f2, print = T, labelx = "Weighted UMAP1", labely = "Weighted UMAP2")
# ggsave("integration/pdf/RNAweighted.umap.bypatient.pdf", plot = f2, width = 6, height = 5.5 )
# 
# cat("==========Joint assay built =============\n")
# 
# 
# combined_obj <- AddModuleScore(combined_obj, features = clean_module, assay = "RNA")
# colnames(combined_obj@meta.data)[grep("Cluster",colnames(combined_obj@meta.data))] <- paste(names(clean_module), "Module")
# 
# module_anno <- give_module_anno(combined_obj@meta.data[, paste(names(clean_module), "Module")])
# combined_obj$module_group <- module_anno
# combined_obj$simple_group <- ifelse(grepl("&", combined_obj$module_group), "Transition", combined_obj$module_group )
# combined_obj$simple_group[combined_obj$module_group == "NE1&NE2"] <- "NE1&NE2"
# Idents(combined_obj) <- combined_obj$simple_group
# combined_obj <- SortIdents(combined_obj)
# DefaultAssay(combined_obj) <- "peaks"
# combined_obj <- RegionStats(combined_obj, genome = BSgenome.Hsapiens.UCSC.hg38)
# 
# #load("/trigos_team/CASCADE/Analysis/signatures_analysis/objects/all_gene_sets.rda")
# #genelist <- lapply(all_gene_sets, function(x) unique(unlist(x)))
# 
# 
# Annotation(combined_obj[["peaks"]]) <- annotation
# 
# 
# p <- DimPlot(combined_obj, group.by = "simple_group",cols = module_col, pt.size = 0.3, reduction = "weighted.umap", raster = F)+ ggtitle("Simple Module Group")
# p <-  coneraxes(p, print = T, labelx = "Weighted UMAP1", labely = "Weighted UMAP2")
# ggsave(paste0("integration/pdf/integrated_module_group.pdf"), width = 5, height = 5, plot = p)
# 
# 
# # link peaks to genes
# combined_obj <- LinkPeaks(
#   object = combined_obj,
#   peak.assay = "peaks",
#   expression.assay = "RNA",
#   genes.use = unlist(clean_module, use.name = F)
# )
# 
# p1 <- CoveragePlot(
#   object = combined_obj,
#   
#   region = "FOLH1",
#   features = "FOLH1",
#   expression.assay = "RNA",
#   group.by = "simple_group",
#   extend.upstream = 500,
#   extend.downstream = 500
# )
# ggsave("integration/pdf/integrated_module_group_FOLH1_peak.pdf", width = 5, height = 5.5, plot = p1)
# 
# # links <- Links(combined_obj)
# # combined_obj <- AddModuleScore(combined_obj, assay = "ruv3", features = list(links$peak[2:7]))
# # combined_obj$FOLH1_peak <- combined_obj$Cluster1
# # f1 <-FeaturePlot(combined_obj, features = "FOLH1_peak", reduction = "weighted.umap", raster = F)
# # f2 <- FeaturePlot(combined_obj, features = "FOLH1", reduction = "weighted.umap", raster = F)
# cat("==========Differential peaks between module group=============\n")
# DefaultAssay(combined_obj) <- "ruv3"
# Annotation(combined_obj[["ruv3"]]) <- annotation
# DE_peak <- FindAllMarkers(combined_obj, test.use = "LR", latent.vars = "nCount_peaks", group.by = "simple_group",
#                           only.pos = TRUE,min.pct = 0.05, logfc.threshold = 0.25)
# 
# 
# DE_peak$closest_gene <- ClosestFeature(combined_obj, DE_peak$gene)$gene_name
# DE_peak$distance <- ClosestFeature(combined_obj, DE_peak$gene)$distance
# 
# saveRDS(DE_peak, paste0("/integration/module_DE_peaks.Rds")) #save 
# 
# cat("==========Try find enriched motifs=============\n")
# 
# combined_obj <- AddMotifs(combined_obj, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
# qsave(combined_obj, "integration/integrated_srt.qs")
# meta.feature <- GetAssayData(combined_obj, assay = "peaks", layer = "meta.features")  
# # find peaks open in Pvalb or Sst cells
# open.peaks <- AccessiblePeaks(combined_obj, assay = "peaks")
# top.da.peak <- unique(DE_peak$gene[DE_peak$p_val_adj < 0.05&DE_peak$avg_log2FC>0.25])
# peaks.matched <- MatchRegionStats(
#   meta.feature = meta.feature[open.peaks, ],
#   query.feature = meta.feature[top.da.peak, ],
#   features.match = "GC.percent",
#   n = 50000, 
#   verbose = T
# )
# 
# sig.motifs <- lapply(unique(DE_peak$cluster), function(clus){
#   cat(clus, "\n")
#   top.da.peak <- DE_peak$gene[DE_peak$p_val_adj < 0.05&DE_peak$avg_log2FC>0.25&DE_peak$cluster == clus]
#   # compare the da.peak and all the open peaks
#   enriched.motifs  <- FindMotifs(
#     object = combined_obj,
#     features = top.da.peak, 
#     background = peaks.matched
#   )
#   sig.motif <- enriched.motifs[enriched.motifs$p.adjust < 0.05,]
#   return(sig.motif)
# })
# names(sig.motifs) <- unique(DE_peak$cluster)
# saveRDS(sig.motifs, paste0("integration/sig.motifs.Rds")) #save 
# 


combined_obj <- qread("integration/integrated_srt.qs")

combined_obj <- RunChromVAR(
  object = combined_obj,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

#change to motif assay
DefaultAssay(combined_obj) <- "chromvar"
motif_name <- sapply(row.names(combined_obj), function(x) {name(getMatrixByID(JASPAR2024@db, ID = x))})[rownames(combined_obj)]
combined_obj[["chromvar"]]@meta.features <- as.data.frame(motif_name)


de.motif.activity <- FindAllMarkers(
  object = combined_obj,
  fc.name = "avg_diff", 
  group.by = "simple_group"
)


motifs_names <- sapply(row.names(de.motif.activity), function(x) {name(getMatrixByID(JASPAR2024@db, ID = x))})
de.motif.activity$motif.name <- motifs_names
saveRDS(de.motif.activity, paste0("integration/DE_motifs_module.Rds")) #save

qsave(combined_obj, "integration/integrated_srt.qs")

cat("All finished!")