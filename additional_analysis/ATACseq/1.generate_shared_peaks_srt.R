library(GenomicRanges)
library(Seurat)
library(Signac)
library(qs)

setwd("~/ATAC")
paths <- system("realpath ~/ATAC/*/*multi_srt.qs", intern = T)
paths <- setNames(paths, sapply(strsplit(paths, split  = "/"), '[', 5))
# Define overlap peak ranges + control peaks for later normalization---------
peaks <- system("realpath ~/ATAC/*/peaks.Rds", intern = T)
peaks <- lapply(peaks, readRDS)
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


combined.peaks <- readRDS("combined.peaks.Rds")

args <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(args[1])
i <- names(paths)[arg]

cat("Start making peaks for", paths[arg], "\n===========================================\n")

srt <- qread(paths[arg])
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
qsave(srt,  paste0("~/ATAC/", i, "/", i, "_merge_srt.qs"))