library(patchwork)
library(tidyverse)
library(MAST)
library(pals)
library(Seurat)
#include new sample
srt <- readRDS("~/integration/ruvIII/integrated_ruv/tumor_only_ruvsrt.Rds")
clean_module <- readRDS("~/regulon/new_regulon/archetype/2024_06/clean_module.Rds")
srt <- AddModuleScore(srt, features = clean_module, name = "Module")
saveRDS(srt@meta.data, "~/regulon/new_regulon/archetype/2024_06/5_integrated_final_module_meta_old.Rds")


pdf("~/regulon/new_regulon/archetype/2024_06/5_integrated_plot_old.pdf", width = 11, height = 4)
FeaturePlot(srt, c("AR", "ASCL1"),  raster = F)

for (i in paste0("Module", 1:6)){
  p <- FeaturePlot(srt, i, cols =  rev(brewer.rdylbu(10)),raster = F)+FeaturePlot(srt, i,raster = F)
  print(p+plot_layout(ncol = 2))
  # p <-VlnPlot(srt, group.by = "sample", split.by = "patient", features = i, cols = brewer.pal(9,"Set3"))+
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  # print(p)
}
dev.off()