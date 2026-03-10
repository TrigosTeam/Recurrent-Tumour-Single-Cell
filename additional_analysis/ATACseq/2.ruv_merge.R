library(ruvIIInb)
library(Seurat)
library(Signac)
library(qs)
library(DelayedArray)
library(dplyr)
library(harmony)
library(ggplot2)
source("~/ATAC/functions.R")
setwd("~/ATAC")

ctl <- readRDS("~/ATAC/ctl.Rds")

# merge all srt to one object and run harmony -------
paths <- system("realpath ~/ATAC/*/*merge_srt.qs", intern = T)
paths <- setNames(paths, sapply(strsplit(paths, split  = "/"), '[', 5))# all_srt <- lapply(paths, qread)
all_srt <-  lapply(paths, qread)
for(i in 1:24){
  all_srt[[i]]$sample <- names(paths)[i]
}

combined_obj <- merge(all_srt[[1]], all_srt[2:length(all_srt)], merge.data = F, project = "tumor_only", add.cell.ids = names(paths))
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
f1 <- DimPlot(combined_obj, reduction = "umap.peak", label = T) + NoLegend()
ggsave("integration/pdf/umap.peak.bycluster.pdf", plot = coneraxes(f1, print = F), width = 5.5, height = 5.5 )
f2 <- DimPlot(combined_obj, group.by = "patient", reduction = "umap.peak", label = T)
ggsave("integration/pdf/umap.peak.bypatient.pdf", plot = coneraxes(f2, print = F), width = 5.5, height = 5.5 )
f3 <- DimPlot(combined_obj, group.by = "sample_id", reduction = "umap.peak", cols = "polychrome", label = T)
ggsave("integration/pdf/umap.peak.bysample.pdf", plot = coneraxes(f3, print = F), width = 5.5, height = 5.5 )
f4 <- DimPlot(combined_obj, group.by = "sample_id", reduction = "umap.peak", cols = "polychrome", split.by = "patient", label = T)+ NoLegend()
ggsave("integration/pdf/umap.peak.bysample.splitpatient.pdf", plot = coneraxes(f4, print = F), width = 18, height = 5.5 )
p1 <- DimPlot(combined_obj, group.by = "patient", reduction = "umap.harmony", label = T)
ggsave("integration/pdf/umap.harmony.bypatient.pdf", plot = coneraxes(p1, print = F), width = 5.5, height = 5.5 )

p2 <- DimPlot(combined_obj, group.by = "harmony.snn_res.0.8", reduction = "umap.harmony", label = T) + NoLegend()
ggsave("integration/pdf/umap.harmony.bycluster.pdf", plot = coneraxes(p2, print = F), width = 5.5, height = 5.5 )

p3 <- DimPlot(combined_obj, reduction = "umap.harmony", group.by = "sample_id",cols = "polychrome")
ggsave("integration/pdf/umap.harmony.bysample.pdf", plot = coneraxes(p3, print = F), width = 5.5, height = 5.5 )

qsave(combined_obj, "integration/combined_obj.qs")

# start ruvnb iii in R 4.2.0
combined_obj <- qread("integration/combined_obj.qs")
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

#saveRDS(ruv3nb_out, "~/integration/2024_06/tumor_only_ruvout.Rds")
print("tumor only ruv finished")
qsave(ruv3nb_out, "integration/ruv3nb_out.qs")

logPAC <- ruv3nb_out@assays@data$logPAC


rownames(logPAC) <- names(ruv3nb_out@rowRanges)
colnames(logPAC) <- rownames(ruv3nb_out@colData)

qsave(logPAC, "integration/logPAC.qs")
rm( list = "ruv3nb_out")
gc()


combined_obj[["ruv3"]]<- CreateChromatinAssay(data = logPAC, fragments = Fragments(combined_obj))
DefaultAssay(combined_obj) <-"ruv3"
combined_obj <- FindTopFeatures(combined_obj, min.cutoff = "q0") %>%
  RunSVD(features = rownames(combined_obj)[!rownames(combined_obj) %in% ctl], reduction.name= "lsi.ruv3")
combined_obj <-  RunUMAP(combined_obj, dims =2:30,reduction = "lsi.ruv3", reduction.name = "umap.peak.ruv3", reduction.key = "ruv3peakUMAP_") %>%
  FindNeighbors(reduction = "lsi.ruv3", dims = 2:30, graph.name = c("ruv3peak.nn", "ruv3peak.snn")) %>%
  FindClusters(graph.name = "ruv3peak.snn", resolution = 1.5)

qsave(combined_obj, "integration/ruv_srt.qs")

f1 <- DimPlot(combined_obj, reduction = "umap.peak.ruv3", label = T) + NoLegend()
ggsave("integration/pdf/umap.peak.ruv3.bycluster.pdf", plot = coneraxes(f1, print = F), width = 6, height = 5.5 )
f2 <- DimPlot(combined_obj, group.by = "patient", reduction = "umap.peak.ruv3", label = T, raster = F)
ggsave("integration/pdf/umap.peak.ruv3.bypatient.pdf", plot = coneraxes(f2, print = F), width = 6, height = 5.5 )
f3 <- DimPlot(combined_obj, group.by = "sample_id", reduction = "umap.peak.ruv3", cols = "polychrome", label = F, raster = F)+NoLegend() + labs(title = "ATACseq")
f3 <- coneraxes(f3, print = T)
ggsave("integration/pdf/umap.peak.ruv3.bysample.pdf", plot = f3, width = 5, height = 5.5 )





