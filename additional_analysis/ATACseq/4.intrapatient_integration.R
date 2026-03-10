library(ruvIIInb)
library(Seurat)
library(Signac)
library(qs)
library(DelayedArray)
library(dplyr)
library(harmony)
library(ggplot2)
library(pals)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(data.table)
library(JASPAR2024)
library(TFBSTools)
library(GenomicRanges)

source("~/ATAC/functions.R")
setwd("~/ATAC/intrapatient")
# merge all srt to one object and subset to one patients -------
paths <- system("realpath ~/ATAC/*/*multi_srt.qs", intern = T)
paths <- setNames(paths, sapply(strsplit(paths, split  = "/"), '[', 5))
patients <- unique(substr(names(paths), 1, 6))

args <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(args[1])

patient <- patients[arg]
ind <- grep(patient, names(paths))
cat("============start analysis on", patient, "samples===========\n")

if (!dir.exists(patient)) dir.create(patient, recursive = TRUE)
setwd(paste0("~/ATAC/intrapatient/", patient))

# Define overlap peak ranges + control peaks for later normalization---------
peaks <- system("realpath ~/ATAC/*/peaks.Rds", intern = T,)
peaks <- lapply(peaks, readRDS)
peaks <- peaks[ind]

combined.peaks <- GenomicRanges::reduce(x = as(peaks, "GRangesList"))
combined.peaks <- combined.peaks[[1]]
saveRDS(combined.peaks, "combined.peaks.Rds")

ctl <- readRDS("~/ATAC/atac_control_peaks.Rds")
control_shared <- intersect(combined.peaks, ctl)
hits <- findOverlaps(control_shared, combined.peaks)
ref_width <- width(combined.peaks[to(hits)])
control <- width(control_shared)
ctl <- combined.peaks[to(hits)[which(control/ref_width > 0.5)]]
ctl <- paste(ctl@seqnames, ctl@ranges, sep = "-")
saveRDS(ctl, "ctl.Rds")

all_srt <-  lapply(paths[ind], qread)

all_srt <- lapply(all_srt, function(srt){
  DefaultAssay(srt) <- "peaks"
  
  called_peak <- as.matrix(GetAssayData(srt, assay = "peaks", layer = "counts"))
  gr <- StringToGRanges(rownames(srt), sep = c("-", "-"))
  shared_peaks <- intersect(combined.peaks, gr)
  hits <- findOverlaps(shared_peaks, gr)
  gr_width <- width(gr[to(hits)])
  shared_width <- width(shared_peaks)
  right_shared_peaks <-  shared_peaks[which(shared_width/gr_width > 0.4)]
  
  
  counts <- matrix(0, nrow = length(combined.peaks), ncol = ncol(srt), dimnames = list(paste(combined.peaks@seqnames, combined.peaks@ranges, sep = "-"), colnames(srt)))
  counts[findOverlaps(right_shared_peaks, combined.peaks, type = "within", select = "first"), ] <- called_peak[findOverlaps(right_shared_peaks, gr, select = "first"),]
  
  assay<- CreateChromatinAssay(
    counts = counts,
    fragments = Fragments(srt)
  )
  
  srt <- CreateSeuratObject(counts = assay, assay = "peaks")
  return(srt)
})

for(i in seq(length(ind))){
  all_srt[[i]]$sample <- names(paths)[ind[i]]
}


combined_obj <- merge(all_srt[[1]], all_srt[2:length(all_srt)], merge.data = F, project = patient, add.cell.ids = names(paths)[ind])
combined_obj <- combined_obj[rowSums(combined_obj) != 0, combined_obj$nFeature_peaks>100]
combined_obj <- FindTopFeatures(combined_obj, min.cutoff = "q0") %>%
  RunTFIDF() %>%
  RunSVD(features = rownames(combined_obj)[!rownames(combined_obj) %in% ctl] )
combined_obj <- RunHarmony(combined_obj, group.by.vars = "sample", reduction.use = "lsi", dims.use = 2:30,project.dim = FALSE)
combined_obj <-  RunUMAP(combined_obj, dims =1:20,reduction = "harmony", reduction.name = "umap.harmony", reduction.key = "harmonyUMAP_") %>%
  FindNeighbors(reduction = "harmony", dims = 1:20, graph.name = c("harmony.nn", "harmony.snn")) %>%
  FindClusters(graph.name = "harmony.snn")

combined_obj <-  RunUMAP(combined_obj, dims =2:30,reduction = "lsi", reduction.name = "umap.peak", reduction.key = "peakUMAP_") %>%
  FindNeighbors(reduction = "lsi", dims = 2:30, graph.name = c("peak.nn", "peak.snn")) %>%
  FindClusters(graph.name = "peak.snn", resolution = 1.5)

sites<- combined_obj$sample
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
combined_obj$site <- sites

combined_obj$patient <- gsub("00", "", substr(combined_obj$sample, 1, 6))
combined_obj$group <- ifelse(combined_obj$patient%in% c("CA0046", "CA0090"), "NE",ifelse(combined_obj$patient%in% c("CA0058", "CA0027"), "Mixed", "AD"))
combined_obj$sample_id <- paste(combined_obj$site, sub(".*_(\\d+)$", "\\1", combined_obj$sample))
combined_obj$sample_id <- paste(combined_obj$patient, combined_obj$sample_id)

if (!dir.exists("pdf")) dir.create("pdf", recursive = TRUE)
f1 <- DimPlot(combined_obj, reduction = "umap.peak", label = T) + NoLegend()
ggsave("pdf/umap.peak.bycluster.pdf", plot = coneraxes(f1, print = F), width = 5.5, height = 5.5 )
f2 <- DimPlot(combined_obj, group.by = "patient", reduction = "umap.peak", label = T)
ggsave("pdf/umap.peak.bypatient.pdf", plot = coneraxes(f2, print = F), width = 5.5, height = 5.5 )
f3 <- DimPlot(combined_obj, group.by = "sample_id", reduction = "umap.peak", cols = "polychrome", label = T)
ggsave("pdf/umap.peak.bysample.pdf", plot = coneraxes(f3, print = F), width = 5.5, height = 5.5 )
f4 <- DimPlot(combined_obj, group.by = "sample_id", reduction = "umap.peak", cols = "polychrome", split.by = "patient", label = T)+ NoLegend()
ggsave("pdf/umap.peak.bysample.splitpatient.pdf", plot = coneraxes(f4, print = F), width = 18, height = 5.5 )
p1 <- DimPlot(combined_obj, group.by = "patient", reduction = "umap.harmony", label = T)
ggsave("pdf/umap.harmony.bypatient.pdf", plot = coneraxes(p1, print = F), width = 5.5, height = 5.5 )

p2 <- DimPlot(combined_obj, group.by = "harmony.snn_res.0.8", reduction = "umap.harmony", label = T) + NoLegend()
ggsave("pdf/umap.harmony.bycluster.pdf", plot = coneraxes(p2, print = F), width = 5.5, height = 5.5 )

p3 <- DimPlot(combined_obj, reduction = "umap.harmony", group.by = "sample_id",cols = "polychrome")
ggsave("pdf/umap.harmony.bysample.pdf", plot = coneraxes(p3, print = F), width = 5.5, height = 5.5 )

qsave(combined_obj, "combined_obj.qs")

samples <- combined_obj$sample
temp <- structure(seq(length(unique(samples))), names = unique(samples))
batch <- temp[samples]


M <- matrix(0,ncol(combined_obj),length(unique(combined_obj$patient)))
colnames(M) <- unique(combined_obj$patient)

for(p in colnames(M)){
  M[which(combined_obj$patient==p),p] <- 1
}

ctl <- intersect(ctl, rownames(combined_obj))
ruv3nb_out<-try(ruvIIInb::fastruvIII.nb(Y= DelayedArray(as.array(GetAssayData(combined_obj, assay = "peaks", layer = "counts"))), # count matrix with genes as rows and cells as columns
                                        M=M, #Replicate matrix constructed as above
                                        ctl=ctl, #A vector denoting control genes
                                        k=2, # dimension of unwanted variation factors
                                        ncores = 2,
                                        batch = batch
))

#saveRDS(ruv3nb_out, "~/2024_06/tumor_only_ruvout.Rds")
print("tumor only ruv finished")
qsave(ruv3nb_out, "ruv3nb_out.qs")

logPAC <- ruv3nb_out@assays@data$logPAC


rownames(logPAC) <- names(ruv3nb_out@rowRanges)
colnames(logPAC) <- rownames(ruv3nb_out@colData)

qsave(logPAC, "logPAC.qs")
rm( list = "ruv3nb_out")
gc()


combined_obj[["ruv3"]]<- CreateChromatinAssay(data = logPAC, fragments = Fragments(combined_obj))
DefaultAssay(combined_obj) <-"ruv3"
combined_obj <- FindTopFeatures(combined_obj, min.cutoff = "q0") %>%
  RunSVD(features = rownames(combined_obj)[!rownames(combined_obj) %in% ctl], reduction.name= "lsi.ruv3")
combined_obj <-  RunUMAP(combined_obj, dims =2:30,reduction = "lsi.ruv3", reduction.name = "umap.peak.ruv3", reduction.key = "ruv3peakUMAP_") %>%
  FindNeighbors(reduction = "lsi.ruv3", dims = 2:30, graph.name = c("ruv3peak.nn", "ruv3peak.snn")) %>%
  FindClusters(graph.name = "ruv3peak.snn", resolution = 1.5)

qsave(combined_obj, "ruv_srt.qs")

f1 <- DimPlot(combined_obj, reduction = "umap.peak.ruv3", label = T) + NoLegend()
ggsave("pdf/umap.peak.ruv3.bycluster.pdf", plot = coneraxes(f1, print = F), width = 6, height = 5.5 )
f2 <- DimPlot(combined_obj, group.by = "patient", reduction = "umap.peak.ruv3", label = T, raster = F)
ggsave("pdf/umap.peak.ruv3.bypatient.pdf", plot = coneraxes(f2, print = F), width = 6, height = 5.5 )
f3 <- DimPlot(combined_obj, group.by = "sample_id", reduction = "umap.peak.ruv3", label = T, raster = F, cols = as.character(brewer.set3(6))) + labs(title = "ATACseq")
f3 <- coneraxes(f3, print = T)
ggsave("pdf/umap.peak.ruv3.bysample.pdf", plot = f3, width = 5, height = 5.5 )


cat("Complete RUNVIIIbn \n")

# run later Signac steps
## get individual sample ids and cluster based on TF only res = 0.5
clean_module <- readRDS("~/CASCADEpaper/paper/Fig5_archetype/clean_module.Rds")
names(clean_module) <- c("AR", "Inflammation", "NE1", "NE2", "Cycling", "Glycolysis")
# plan("multicore", workers = 4)
# options(future.globals.maxSize = 8000 * 1024^2)

module_col <- setNames( c("dodgerblue3", brewer.dark2(6)[2:6]),names(clean_module))
module_col <- c(module_col, setNames(c("grey90", "grey30", "red"), c("Background","Transition", "NE1&NE2")))
# annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
annotation <- readRDS("~/ATAC/annotation.Rds")
JASPAR2024 <- JASPAR2024()
pfm <- getMatrixSet(
  x = JASPAR2024@db,
  opts = list(collection = "CORE", species ='Homo sapiens')
)

combined_obj <- qread("ruv_srt.qs")
RNAmat <- readRDS("~/integration/2024_06/tumor_only_logPAC.Rds")
combined_obj[["RNA"]] <- CreateAssayObject(data = RNAmat[, intersect(colnames(combined_obj), colnames(RNAmat))])
DefaultAssay(combined_obj) <- "RNA"
combined_obj[["RNA"]]@scale.data <- as.matrix(combined_obj[["RNA"]]@data)
combined_obj <- FindVariableFeatures(combined_obj)
combined_obj <- RunPCA(combined_obj,npc = 30)  



combined_obj <- FindMultiModalNeighbors(
  object = combined_obj,
  reduction.list = list("pca", "lsi.ruv3"), 
  dims.list = list(1:20, 2:30),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
combined_obj <- RunUMAP(
  object = combined_obj,
  nn.name = "weighted.nn",
  assay = "RNA",
  reduction.name = "weighted.umap",
  verbose = TRUE
)
# find clusters based on combined weighted graph of RNA and ATAC
combined_obj <- FindClusters(combined_obj, graph.name = "wsnn")

f2 <- DimPlot(combined_obj, group.by = "sample_id", reduction = "weighted.umap", label = T, raster = F, cols = as.character(brewer.set3(6)))
f2 <-  coneraxes(f2, print = T, labelx = "Weighted UMAP1", labely = "Weighted UMAP2")
ggsave("pdf/RNAweighted.umap.bysample.pdf", plot = f2, width = 6, height = 5.5 )

cat("==========Joint assay built =============\n")


combined_obj <- AddModuleScore(combined_obj, features = clean_module, assay = "RNA")
colnames(combined_obj@meta.data)[grep("Cluster",colnames(combined_obj@meta.data))] <- paste(names(clean_module), "Module")

module_anno <- give_module_anno(combined_obj@meta.data[, paste(names(clean_module), "Module")])
combined_obj$module_group <- module_anno
combined_obj$simple_group <- ifelse(grepl("&", combined_obj$module_group), "Transition", combined_obj$module_group )
combined_obj$simple_group[combined_obj$module_group == "NE1&NE2"] <- "NE1&NE2"
Idents(combined_obj) <- combined_obj$simple_group
combined_obj <- SortIdents(combined_obj)
p <- DimPlot(combined_obj, group.by = "simple_group",cols = module_col, pt.size = 0.3, reduction = "weighted.umap", raster = F)+ ggtitle("Simple Module Group")
p <-  coneraxes(p, print = T, labelx = "Weighted UMAP1", labely = "Weighted UMAP2")
ggsave(paste0("pdf/integrated_module_group.pdf"), width = 5, height = 5, plot = p)


DefaultAssay(combined_obj) <- "peaks"
Annotation(combined_obj[["peaks"]]) <- annotation
combined_obj <- RegionStats(combined_obj, genome = BSgenome.Hsapiens.UCSC.hg38)

#load("~/CASCADE/Analysis/signatures_analysis/objects/all_gene_sets.rda")
#genelist <- lapply(all_gene_sets, function(x) unique(unlist(x)))


# link peaks to genes
combined_obj <- LinkPeaks(
  object = combined_obj,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = unlist(clean_module, use.name = F)
)
qsave(combined_obj, "integrated_srt.qs")
p1 <- CoveragePlot(
  object = combined_obj,
  
  region = "FOLH1",
  features = "FOLH1",
  expression.assay = "RNA",
  group.by = "simple_group",
  extend.upstream = 500,
  extend.downstream = 500
)
ggsave("pdf/integrated_module_group_FOLH1_peak.pdf", width = 5, height = 5.5, plot = p1)

# links <- Links(combined_obj)
# combined_obj <- AddModuleScore(combined_obj, assay = "ruv3", features = list(links$peak[2:7]))
# combined_obj$FOLH1_peak <- combined_obj$Cluster1
# f1 <-FeaturePlot(combined_obj, features = "FOLH1_peak", reduction = "weighted.umap", raster = F)
# f2 <- FeaturePlot(combined_obj, features = "FOLH1", reduction = "weighted.umap", raster = F)
cat("==========Differential peaks between module group=============\n")
DefaultAssay(combined_obj) <- "ruv3"
Annotation(combined_obj[["ruv3"]]) <- annotation
DE_peak <- FindAllMarkers(combined_obj, test.use = "LR", latent.vars = "nCount_peaks", group.by = "simple_group",
                          only.pos = TRUE,min.pct = 0.05, logfc.threshold = 0.25)


DE_peak$closest_gene <- ClosestFeature(combined_obj, DE_peak$gene)$gene_name
DE_peak$distance <- ClosestFeature(combined_obj, DE_peak$gene)$distance

saveRDS(DE_peak, "module_DE_peaks.Rds") #save 

gene_activity <- GeneActivity(combined_obj)
saveRDS(gene_activity, "gene_activity.Rds")
combined_obj[["gene_activity"]] <- CreateAssayObject(counts = gene_activity)

# combined_obj <- NormalizeData(
#   object = combined_obj,
#   assay = "gene_activity",
#   normalization.method = "LogNormalize",
#   scale.factor = median(combined_obj$nCount_gene_activity)
# )

cat("==========Try find enriched motifs=============\n")

combined_obj <- AddMotifs(combined_obj, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
qsave(combined_obj, "integrated_srt.qs")
meta.feature <- GetAssayData(combined_obj, assay = "peaks", layer = "meta.features")  
# find peaks open in Pvalb or Sst cells
open.peaks <- AccessiblePeaks(combined_obj, assay = "peaks")
top.da.peak <- unique(DE_peak$gene[DE_peak$p_val_adj < 0.05&DE_peak$avg_log2FC>0.25])
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top.da.peak, ],
  features.match = "GC.percent",
  n = 50000, 
  verbose = T
)

sig.motifs <- lapply(unique(as.character(DE_peak$cluster)), function(clus){
  top.da.peak <- DE_peak$gene[DE_peak$p_val_adj < 0.05&DE_peak$avg_log2FC>0.25&DE_peak$cluster == clus]
  # compare the da.peak and all the open peaks
  enriched.motifs  <- FindMotifs(
    object = combined_obj,
    features = top.da.peak, 
    background = peaks.matched
  )
  sig.motif <- enriched.motifs[enriched.motifs$p.adjust < 0.05,]
  return(sig.motif)
})
names(sig.motifs) <- unique(DE_peak$cluster)
saveRDS(sig.motifs, paste0("sig.motifs.Rds")) #save 





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
saveRDS(de.motif.activity, paste0("DE_motifs_module.Rds")) #save

qsave(combined_obj, "integrated_srt.qs")
cat("All completed!!\n")
