args <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(args[1])
library(qs)
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
library(JASPAR2024)
library(TFBSTools)
library(future)
library(dorothea)
library(ComplexHeatmap)
source("functions.R")
library(data.table)
library(pals)
setwd("~/ATAC")
JASPAR2024 <- JASPAR2024()
pfm <- getMatrixSet(
  x = JASPAR2024@db,
  opts = list(collection = "CORE", species ='Homo sapiens')
)

net <- dorothea_hs %>% dplyr::filter(confidence != "D") %>% dplyr::select(tf, target)

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

paths2 <- system("realpath ~/CASCADE/SingleCell/CA*/*/multiome/*/outs/raw_feature_bc_matrix.h5", intern = T)
#paths2.2 <- system("realpath ~/CASCADE/SingleCell/CA*/*/*/outs/raw_feature_bc_matrix.h5", intern = T)
#paths2 <- c(paths2, paths2.2)
paths_split <- unlist(lapply(strsplit(paths2, split = "/"), function(x) paste(x[5], x[6], sep = "_")))
paths2 <- setNames(paths2, paths_split)

fragpath <- paste0(unlist(lapply(strsplit(paths2, split = "/"), function(x) paste(x[1:(length(x)-1)], collapse  = "/"))), "/atac_fragments.tsv.gz")
fragpath <- setNames(fragpath, paths_split)

paths <- system("realpath ~/CASCADE/Analysis/230113_seurat_intron/*/*/*non_tumour_removed.Rds", intern = T)
paths_split <- unlist(lapply(strsplit(paths, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
paths3 <- system("realpath ~/CASCADE/Analysis/240603_seurat_intron/CA0027/dura*/*non_tumour_removed_v4.Rds", intern = T)
paths_split3 <- unlist(lapply(strsplit(paths3, split = "/"), function(x) paste(x[6], x[7], sep = "_")))

paths <- c(paths, paths3)
paths_split <- c(paths_split, paths_split3)
paths <- setNames(paths, paths_split)



## get individual sample ids and cluster based on TF only res = 0.5
clean_module <- readRDS("~/CASCADEpaper/paper/Fig5_archetype/clean_module.Rds")
names(clean_module) <- c("AR", "Inflammation", "NE1", "NE2", "Cycling", "Glycolysis")
# plan("multicore", workers = 4)
# options(future.globals.maxSize = 8000 * 1024^2)

module_col <- setNames( c("dodgerblue3", brewer.dark2(6)[2:6]),names(clean_module))
module_col <- c(module_col, setNames(c("grey90", "grey30", "red"), c("Background","Transition", "NE1&NE2")))

i = names(paths2)[arg]

  if (!dir.exists(i)) {
    dir.create(i)
  }
cat("==============Start", i, "processing============\n")
  
  # read in RNA srt as baseline
  srt <- readRDS(paths[i])
  test_mat <- Read10X_h5(paths2[i])
  if(!all(colnames(srt) %in% colnames(test_mat$Peaks))){
    srt <- subset(srt, cells = intersect(colnames(srt), colnames(test_mat$Peaks)))
    cat("cell id unmatched in", i, "\n")
  }
  
  srt <- preprocess_srt(srt)
  ## calculate expression of module signatures
  srt <- AddModuleScore(srt, features = clean_module, assay = "RNA")
  colnames(srt@meta.data)[grep("Cluster",colnames(srt@meta.data))] <- paste(names(clean_module), "Module")
  
  # create ATAC assay and add it to the object with the same cells in snRNA object from CellRanger
  srt[["ATAC"]] <- CreateChromatinAssay(
    counts = test_mat$Peaks[, colnames(srt)],
    sep = c(":", "-"),
    fragments = fragpath[i],
    annotation = annotation
  )
  rm(list = "test_mat")
  gc()
  
  DefaultAssay(srt)<- "ATAC"
  srt <- NucleosomeSignal(srt)
  srt <- TSSEnrichment(srt)

  
  srt <- FindTopFeatures(srt, min.cutoff = "q0")
  srt <- RunTFIDF(srt)
  srt <- RunSVD(srt)
  srt <- RunUMAP(srt, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
  srt <- FindNeighbors(srt, reduction = "lsi", dims = 2:30, graph.name = c("atac.nn", "atac.snn"))
  srt <- FindClusters(srt, graph.name = "atac.snn")
  cat("==========ATAC assay built =============\n")
  # joint UMAP between CellRanger 
  srt <- FindMultiModalNeighbors(
    object = srt,
    reduction.list = list("pca", "lsi"), 
    dims.list = list(1:20, 2:30),
    modality.weight.name = "RNA.weight",
    verbose = TRUE
  )
  
  # build a joint UMAP visualization
  srt <- RunUMAP(
    object = srt,
    nn.name = "weighted.nn",
    assay = "RNA",
    reduction.name = "weighted.umap",
    verbose = TRUE
  )
  # find clusters based on combined weighted graph of RNA and ATAC
  srt <- FindClusters(srt, graph.name = "wsnn")
  cat("==========Joint assay built =============\n")
  #Grouping variable to use. If set, peaks will be called independently on each group of cells and then combined. 
  peaks <- CallPeaks(srt, macs2.path = "~/.conda/envs/macs3_env/bin/macs3") ## ??? doubt whether to do by group

  # remove peaks on nonstandard chromosomes and in genomic blacklist regions
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
  saveRDS(peaks, paste0(i, "/peaks.Rds"))
  
  # quantify counts in each peak
  macs2_counts <- FeatureMatrix(
    fragments = Fragments(srt),
    features = peaks,
    cells = colnames(srt)
  )
  
  # create a new assay using the MACS2 peak set and add it to the Seurat object
  srt[["peaks"]] <- CreateChromatinAssay(
    counts = macs2_counts,
    fragments = fragpath[[i]],
    annotation = annotation
  )
  
  #processing peaks, dim reduction
  DefaultAssay(srt) <- "peaks"
  srt <- FindTopFeatures(srt, min.cutoff = "q0")
  srt <- RunTFIDF(srt)
  srt <- RunSVD(srt, reduction.name = "lsi.peak", reduction.key = "LSIpeak_")
  
  # DepthCor(srt, assay = "peaks", reduction = "lsi.peak")
  srt <- RunUMAP(srt, reduction = "lsi.peak", dims = 2:30, reduction.name = "umap.peak") %>%
    FindNeighbors(reduction = "lsi.peak", dims = 2:30, graph.name = c("peak.nn", "peak.snn")) %>% 
    FindClusters(graph.name = "peak.snn")
  srt <- RegionStats(srt, genome = BSgenome.Hsapiens.UCSC.hg38)
  cat("==========Peak assay built =============\n")
  srt <- FindMultiModalNeighbors(
    object = srt,
    reduction.list = list("pca", "lsi.peak"), 
    dims.list = list(1:50, 2:40),
    modality.weight.name = "RNA.peak.weight",
    weighted.nn.name = "weighted.peak.nn",
    verbose = TRUE
  )
  
  # build a joint UMAP visualization
  srt <- RunUMAP(
    object = srt,
    nn.name = "weighted.peak.nn",
    assay = "RNA",
    reduction.name = "weighted.peak.umap",
    verbose = TRUE
  )
  srt <- FindClusters(srt, graph.name  = "wknn")
  cat("==========Peak joint assay built =============\n")
  module_anno <- give_module_anno(srt@meta.data[, paste(names(clean_module), "Module")])
  srt$module_group <- module_anno
  srt$simple_group <- ifelse(grepl("&", srt$module_group), "Transition", srt$module_group )
  srt$simple_group[srt$module_group == "NE1&NE2"] <- "NE1&NE2"
  
  cat("\n")
  cat("save in middle steps to in case\n")
  qsave(srt, paste0("~/ATAC/", i, "/", i, "_multi_srt.qs"))
  
  p <- DimPlot(srt, group.by = "simple_group",cols = module_col, pt.size = 0.3)+ ggtitle("Simple Module Group")
  ggsave(paste0(i, "/Module_group.pdf"), width = 4, height = 3, plot = p)
  p1 <- CoveragePlot(
    object = srt,
    region = "FOLH1",
    features = "FOLH1",
    expression.assay = "RNA",
    group.by = "simple_group",
    extend.upstream = 500,
    extend.downstream = 500
  )
  ggsave(paste0(i, "/FOLH1_peak.pdf"), width = 6, height = 4, plot = p1)
  
  cat("===========Plot all UMAPs and save to in case=============\n")
  plotlist <- mapply(plot_fun, srt = list(srt), 
                     reduction = c("umap", "umap.atac", "weighted.umap", "umap.peak","weighted.peak.umap"), 
                     title = c("RNA UMAP", "ATAC UMAP", "Joint UMAP", "Peak UMAP", "Joint Peak UMAP"),
                     SIMPLIFY = F)
  g <- wrap_plots(plotlist, ncol = 1)
  ggsave(paste0(i ,"/Dim_featureplot.pdf"), width = 18, height = 12)
  # qsave(srt, paste0("~/ATAC/", i, "/", i, "_multi_srt.qs"))
  
  
  cat("==========Differential peaks between clusters=============\n")
  DE_peak <- FindAllMarkers(srt, test.use = "LR", latent.vars = "nCount_peaks",
                            only.pos = TRUE,min.pct = 0.05, logfc.threshold = 0.25)
  
  
  DE_peak$closest_gene <- ClosestFeature(srt, DE_peak$gene)$gene_name
  DE_peak$distance <- ClosestFeature(srt, DE_peak$gene)$distance
  saveRDS(DE_peak, paste0("~/ATAC/", i, "/DE_peaks.Rds")) #save 
  
  cat("==========Try find enriched motifs=============\n")
  
  srt <- AddMotifs(srt, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
  meta.feature <- GetAssayData(srt, assay = "peaks", layer = "meta.features")  
  # find peaks open in Pvalb or Sst cells
  open.peaks <- AccessiblePeaks(srt)
  top.da.peak <- unique(DE_peak$gene[DE_peak$p_val_adj < 0.05&DE_peak$avg_log2FC>0.25])
  peaks.matched <- MatchRegionStats(
    meta.feature = meta.feature[open.peaks, ],
    query.feature = meta.feature[top.da.peak, ],
    features.match = "GC.percent",
    n = 50000, 
    verbose = T
  )
  
  sig.motifs <- lapply(unique(DE_peak$cluster), function(clus){
    top.da.peak <- DE_peak$gene[DE_peak$p_val_adj < 0.05&DE_peak$avg_log2FC>0.25&DE_peak$cluster == clus]
    # compare the da.peak and all the open peaks
    enriched.motifs  <- FindMotifs(
      object = srt,
      features = top.da.peak, 
      background = peaks.matched
    )
     sig.motif <- enriched.motifs[enriched.motifs$p.adjust < 0.05,]
     return(sig.motif)
  })
  names(sig.motifs) <- unique(DE_peak$cluster)
  saveRDS(sig.motifs, paste0("~/ATAC/", i, "/sig.motifs.Rds")) #save 
  
  cat("==========Differential peaks between module group=============\n")
  DE_peak <- FindAllMarkers(srt, test.use = "LR", latent.vars = "nCount_peaks", group.by = "simple_group",
                            only.pos = TRUE,min.pct = 0.05, logfc.threshold = 0.25)
  
  
  DE_peak$closest_gene <- ClosestFeature(srt, DE_peak$gene)$gene_name
  DE_peak$distance <- ClosestFeature(srt, DE_peak$gene)$distance
  saveRDS(DE_peak, paste0("~/ATAC/", i, "/module_DE_peaks.Rds")) #save 
  
  
  
  srt <- RunChromVAR(
    object = srt,
    genome = BSgenome.Hsapiens.UCSC.hg38
  )
  
  #change to motif assay
  DefaultAssay(srt) <- "chromvar"
  motif_name <- sapply(row.names(srt), function(x) {name(getMatrixByID(JASPAR2024@db, ID = x))})[rownames(srt)]
  srt[["chromvar"]]@meta.features <- as.data.frame(motif_name)
  
  
  de.motif.activity <- FindAllMarkers(
    object = srt,
    fc.name = "avg_diff", 
    group.by = "simple_group"
  )
  
  
  motifs_names <- sapply(row.names(de.motif.activity), function(x) {name(getMatrixByID(JASPAR2024@db, ID = x))})
  de.motif.activity$motif.name <- motifs_names
  saveRDS(de.motif.activity, paste0("~/ATAC/", i, "/DE_motifs_module.Rds")) #save 
  
  qsave(srt, paste0("~/ATAC/", i, "/", i, "_multi_srt.qs"))
  print(i)











