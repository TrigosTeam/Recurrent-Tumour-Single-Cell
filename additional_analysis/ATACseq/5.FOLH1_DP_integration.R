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
setwd("~/ATAC/integration")
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

combined_obj <- qread("integrated_srt.qs")


cat("==========Differential peaks between module group=============\n")
DefaultAssay(combined_obj) <- "ruv3"
Annotation(combined_obj[["ruv3"]]) <- annotation
Annotation(combined_obj[["peaks"]]) <- annotation
combined_obj$FOLH1_exp <- as.numeric(combined_obj@assays$RNA@data["FOLH1", ])

if (min(combined_obj$FOLH1_exp) ==0){
  if (sum(combined_obj$FOLH1_exp==0) > 50){
    combined_obj$FOLH1_group <- ifelse(combined_obj$FOLH1_exp >mean(combined_obj$FOLH1_exp) , "high", "low")
    combined_obj$FOLH1_group <- ifelse(combined_obj$FOLH1_exp ==0 , "neg", combined_obj$FOLH1_group)
    cat("There are FOLH1 negative cells\n")
  }else{
    combined_obj$FOLH1_group <- ifelse(combined_obj$FOLH1_exp >mean(combined_obj$FOLH1_exp) , "high", "low")
    cat("There are less than 50 FOLH1 negative cells\n")
  }
}else{
  combined_obj$FOLH1_group <- ifelse(combined_obj$FOLH1_exp >mean(combined_obj$FOLH1_exp) , "high", "low")
  cat("There are no FOLH1 negative cells\n")
}
Idents(combined_obj) <- combined_obj$FOLH1_group

DE_peak <- FindAllMarkers(combined_obj, test.use = "LR", latent.vars = "nCount_peaks", group.by = "FOLH1_group",
                          only.pos = TRUE,min.pct = 0.05, logfc.threshold = 0.25)


DE_peak$closest_gene <- ClosestFeature(combined_obj, DE_peak$gene)$gene_name
DE_peak$distance <- ClosestFeature(combined_obj, DE_peak$gene)$distance

saveRDS(DE_peak, "FOLH1_DE_peaks.Rds") #save
DE_peak <- readRDS("~/ATAC/integration/FOLH1_DE_peaks.Rds")
cat("==========Try find enriched motifs=============\n")

combined_obj <- AddMotifs(combined_obj, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)

meta.feature <- GetAssayData(combined_obj, assay = "ruv3", layer = "meta.features")  

combined_obj@assays$ruv3@counts <- combined_obj@assays$ruv3@data
open.peaks <- AccessiblePeaks(combined_obj, assay = "ruv3")
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
      object = combined_obj,
      features = top.da.peak, 
      background = peaks.matched
    )
    sig.motif <- enriched.motifs[enriched.motifs$p.adjust < 0.05,]
    return(sig.motif)
  })
  names(sig.motifs) <- as.character(unique(DE_peak$cluster))
  saveRDS(sig.motifs, paste0("FOLH1_sig.motifs2.Rds")) #save 
} else {
  cat("No significant DE peaks\n")
}

write.table(rbindlist(sig.motifs, idcol = "FOLH1_group"), "significant_FOLH1_motif2.csv", row.names = F, quote = F)


#change to motif assay
combined_obj <- RunChromVAR(
  object = combined_obj,
  genome = BSgenome.Hsapiens.UCSC.hg38
)
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
saveRDS(de.motif.activity, paste0("DE_motifs_activity_FOLH1.Rds")) #save


cat("All finished!")