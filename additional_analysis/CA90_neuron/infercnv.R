setwd("~/snRNA/integration/Neuron_infercnv")
library(Seurat)
library(infercnv)

srt <- readRDS("~/snRNA/integration/perpatient/CA0090_ruvsrt.Rds")
normalsrt <- readRDS("~/snRNA/tumour_normal/final_normalsrt.Rds")

merge_srt <- merge(srt, normalsrt)

counts <- GetAssayData(merge_srt, layer = "counts", assay = "RNA")
anno <- rbind(data.frame(anno = srt$sample), data.frame(anno =normalsrt$cell_anno))
CA90_neu <- colnames(normalsrt)[normalsrt$cell_anno == "Neurons"& normalsrt$patient == "CA0090"]
anno[CA90_neu, "anno"] <- "CA90 Neurons"

infercnvobject = CreateInfercnvObject(raw_counts_matrix=counts,
                                      annotations_file=anno,
                                      delim="\t",
                                      gene_order_file= "multi_g.txt",
                                      ref_group_names= levels(normalsrt$cell_anno))
if (!file.exists("infercnv")){
  dir.create("infercnv")
}

infercnvobject = infercnv::run(infercnvobject,
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir="infercnv",
                               analysis_mode = "cells",
                               cluster_by_groups= T,
                               denoise=F,
                               HMM=F, 
                               output_format = "pdf",
                               num_threads = 1)
cat("Done")

Neuron_feature <- c( "GPC5", "NRXN1", "PCDH9", #astrocyte #G6PAM endothelial marker
                     "RBFOX1", "CADM2", "CSMD1", "KCNIP4")#, #exitatory neu
#"NRXN3", "CNTNAP2", "ROBO2") 

DimPlot(srt)
