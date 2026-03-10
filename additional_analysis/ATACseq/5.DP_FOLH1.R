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
library(data.table)
library(pals)
setwd("~/ATAC")
source("functions.R")
paths <- system("realpath ~/ATAC/*/*multi_srt.qs", intern = T)
paths <- setNames(paths, sapply(strsplit(paths, split  = "/"), '[', 5))
i = names(paths)[arg]

JASPAR2024 <- JASPAR2024()
pfm <- getMatrixSet(
  x = JASPAR2024@db,
  opts = list(collection = "CORE", species ='Homo sapiens')
)

net <- dorothea_hs %>% dplyr::filter(confidence != "D") %>% dplyr::select(tf, target)

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

cat("Start making peaks for", paths[arg], "\n===========================================\n")

srt <- qread(paths[arg])
DefaultAssay(srt) <- "peaks"
srt$FOLH1_exp <- as.numeric(srt@assays$RNA@data["FOLH1", ])

if (min(srt$FOLH1_exp) ==0){
  if (sum(srt$FOLH1_exp) > 50){
    srt$FOLH1_group <- ifelse(srt$FOLH1_exp >median(srt$FOLH1_exp) , "high", "low")
    srt$FOLH1_group <- ifelse(srt$FOLH1_exp ==0 , "neg", srt$FOLH1_group)
    cat("There are FOLH1 negative cells\n")
  }else{
    srt$FOLH1_group <- ifelse(srt$FOLH1_exp >median(srt$FOLH1_exp) , "high", "low")
    cat("There are less than 50 FOLH1 negative cells\n")
  }
}else{
  srt$FOLH1_group <- ifelse(srt$FOLH1_exp >median(srt$FOLH1_exp) , "high", "low")
  cat("There are no FOLH1 negative cells\n")
}
Idents(srt) <- srt$FOLH1_group

DE_peak <- FindAllMarkers(srt, test.use = "LR", latent.vars = "nCount_peaks",
                          only.pos = TRUE,min.pct = 0.05, logfc.threshold = 0.25)
if (nrow(DE_peak) < 1){
  cat("No differiantial peaks", i)
}else  {
  DE_peak$closest_gene <- ClosestFeature(srt, DE_peak$gene)$gene_name
  DE_peak$distance <- ClosestFeature(srt, DE_peak$gene)$distance
  saveRDS(DE_peak, paste0("~/ATAC/", i, "/FOLH1_DE_peaks.Rds")) #save 
  
  
  meta.feature <- GetAssayData(srt, assay = "peaks", layer = "meta.features")  
  # find peaks open in Pvalb or Sst cells
  open.peaks <- AccessiblePeaks(srt)
  top.da.peak <- unique(DE_peak$gene[DE_peak$p_val_adj < 0.05])
  if (length(top.da.peak) >0){
     peaks.matched <- MatchRegionStats(
    meta.feature = meta.feature[open.peaks, ],
    query.feature = meta.feature[top.da.peak, ],
    features.match = "GC.percent",
    n = 50000, 
    verbose = T
  )
  
  sig.motifs <- lapply(as.character(unique(DE_peak$cluster)), function(clus){
    top.da.peak <- DE_peak$gene[DE_peak$p_val_adj < 0.05&DE_peak$cluster == clus]
    # compare the da.peak and all the open peaks
    enriched.motifs  <- FindMotifs(
      object = srt,
      features = top.da.peak, 
      background = peaks.matched
    )
    sig.motif <- enriched.motifs[enriched.motifs$p.adjust < 0.05,]
    return(sig.motif)
  })
  names(sig.motifs) <- as.character(unique(DE_peak$cluster))
  saveRDS(sig.motifs, paste0("~/ATAC/", i, "/FOLH1_sig.motifs.Rds")) #save 
  } else {
    cat("No significant DE peaks in ", i , "\n")
  }
 
}



