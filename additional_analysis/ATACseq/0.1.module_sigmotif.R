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

clean_module <- readRDS("~/CASCADEpaper/paper/Fig5_archetype/clean_module.Rds")
names(clean_module) <- c("AR", "Inflammation", "NE1", "NE2", "Cycling", "Glycolysis")

net <- dorothea_hs %>% dplyr::filter(confidence != "D") %>% dplyr::select(tf, target)

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

paths <- system("realpath ~/ATAC/*/*_multi_srt.qs", intern = T)
paths_split <- unlist(lapply(strsplit(paths, split = "/"), '[', 5))
paths <- setNames(paths, paths_split)

sample <- paths_split[arg]
srt <- qread(paths[arg])
setwd(sample)
DefaultAssay(srt) <- "peaks"


cells <- apply(srt@meta.data[, paste(names(clean_module),"Module")], 2, function(x){
  bins <- cut(x, breaks = 10) 
  return(colnames(srt)[bins %in% levels(bins)[7:10]])
})

DE_peak <- lapply(cells, function(x){
  if (length(x) >20 & length(x) < ncol(srt)){
    srt$group <- ifelse(colnames(srt) %in% x , "high", "other")
    #srt$group[!colnames(srt)%in%unlist(cells)] <- "background"
    marker <- try(FindMarkers(srt, test.use = "LR", latent.vars = "nCount_peaks", group.by = "group", only.pos = TRUE,
                              ident.1 = "high", ident.2 = "other", logfc.threshold = 0.25))
    if (nrow(marker) >0){
      marker$gene <- rownames(marker)
      marker$closest_gene <- ClosestFeature(srt, marker$gene)$gene_name
      marker$distance <- ClosestFeature(srt, marker$gene)$distance
    }
    
    return(marker)}
})
DE_peak <- rbindlist(DE_peak, fill = T, idcol = "cluster")
saveRDS(DE_peak, "module_DE_peak_binary.Rds")

srt <- AddMotifs(srt, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
meta.feature <- GetAssayData(srt, assay = "peaks", layer = "meta.features")  
# find peaks open in Pvalb or Sst cells
open.peaks <- AccessiblePeaks(srt)
top.da.peak <- unique(DE_peak$gene[DE_peak$p_val_adj < 0.05])
if (length(top.da.peak)>0){
  peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top.da.peak, ],
  features.match = "GC.percent",
  n = 50000, 
  verbose = T
)

sig.motifs <- lapply(as.character(unique(DE_peak$cluster)), function(clus){
  top.da.peak <- DE_peak$gene[DE_peak$p_val_adj < 0.05&DE_peak$cluster == clus]
  if (length(top.da.peak) > 0){
  # compare the da.peak and all the open peaks
  enriched.motifs  <- FindMotifs(
    object = srt,
    features = top.da.peak, 
    background = peaks.matched
  )
  sig.motif <- enriched.motifs[enriched.motifs$p.adjust < 0.05,]
  return(sig.motif)
  }else{
    cat("No significant peaks for", clus, "\n")
    return(data.frame())
  }
})
names(sig.motifs) <- as.character(unique(DE_peak$cluster))
sig.motifs <- rbindlist(sig.motifs, idcol = "group", fill = T)
saveRDS(sig.motifs, paste0("module_sig.motifs_binary.Rds")) #save 
}else{
  cat("No significant peaks")
}




# later analysis
paths <- system("realpath ~/ATAC/*/module_sig.motifs_binary.Rds", intern = T)
paths_split <- unlist(lapply(strsplit(paths, split = "/"), '[', 5))
paths <- setNames(paths, paths_split)
motif <- lapply(paths, readRDS)
motifdf <- rbindlist(motif, idcol = "sample", fill = T)
motifdf <- motifdf[motifdf$p.adjust <0.05,]
countdf <- as.data.frame(table(motifdf$group, motifdf$motif.name))
countable <- as.data.frame.matrix(table(motifdf$group, motifdf$motif.name))
countable <- countable[, apply(countable, 2, max)>2]
countable <- countable[c("AR Module", "Inflammation Module", "NE1 Module", "NE2 Module", "Glycolysis Module", "Cycling Module"), ]



avgF <- motifdf %>% group_by( group, motif.name) %>% 
  summarise(avg.fold.enrichment = mean(fold.enrichment)) %>% 
  dplyr::filter(motif.name %in% colnames(countable))
avgmat <- reshape2::dcast(avgF, group~motif.name, value.var = "avg.fold.enrichment", fill = 0) %>% column_to_rownames("group")
avgmat <- t(scale(t(avgmat)))
pdf("avgFC_module_motif.pdf", width = 70, height = 6)

ht <- Heatmap(avgmat,  cluster_rows = F, 
        name = "scaled average\nfold enrichment")
draw(ht)
write.table(avgmat[, column_order(ht)], "module_motif_avgFC.csv", quote = F, sep = "\t")
dev.off()


pdf("n_sample_module_motif.pdf", width = 70, height = 6)
ht2 <- Heatmap(countable, name = "n of sample",
        cluster_rows = F)
draw(ht2)
write.table(countable[, column_order(ht)], "module_motif_countable.csv", quote = F, sep = "\t")
dev.off()

rec_count <- apply(countable, 1, function(x) colnames(countable)[which(x >6)])

rec_avg <- apply(avgmat, 1, function(x) colnames(avgmat)[which(x >0)])

clean_motif <- list()
for (i in names(rec_count)){
  clean_motif[[i]] <- intersect(rec_count[[i]], rec_avg[[i]])  
}

saveRDS(clean_motif, "clean_motif.Rds")
shared_motif <- lapply(clean_motif, function(x) lapply(clean_motif, function(y) intersect(y, x)))

table(unlist(clean_motif))
