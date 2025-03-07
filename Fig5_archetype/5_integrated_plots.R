library(Seurat)
library(patchwork)
library(tidyverse)
library(MAST)
library(pals)
#include new sample
srt <- readRDS("~/integration/2024_06/tumor_only_ruvsrt.Rds")
clean_module <- readRDS("~/regulon/new_regulon/archetype/2024_06/clean_module.Rds")
srt <- AddModuleScore(srt, features = clean_module, name = "Module")
saveRDS(srt@meta.data, "~/regulon/new_regulon/archetype/2024_06/5_integrated_final_module_meta.Rds")
source("~/CASCADEpaper/paper/functions.R")

pdf("~/regulon/new_regulon/archetype/2024_06/5_integrated_plot.pdf", width = 8, height = 4)
f <-  FeaturePlot(srt, c("AR", "ASCL1"),  raster = F)

for (i in paste0("Module", 1:6)){
  p <- FeaturePlot(srt, i, cols =  rev(brewer.rdylbu(10)),raster = F)+FeaturePlot(srt, i,raster = F)
  print(p+plot_layout(ncol = 2))
  # p <-VlnPlot(srt, group.by = "sample", split.by = "patient", features = i, cols = brewer.pal(9,"Set3"))+
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  # print(p)
}
dev.off()

cells <- apply(srt@meta.data[, grep("Module", colnames(srt@meta.data))], 2, function(x){
  bins <- cut(x, breaks = 10)
  return(colnames(srt)[bins %in% levels(bins)[7:10]])
})


markerlist <- lapply(cells, function(x){
  srt$group <- ifelse(colnames(srt) %in% x, "high", "other")
  Idents(srt) <- srt$group
  return(FindMarkers(srt, ident.1 = "high", ident.2 = "other", test.use ="MAST"))
})
saveRDS(markerlist, "~/regulon/new_regulon/archetype/2024_06/5_integrated_markerlist.Rds")

# paper plot
meta <- readRDS("~/CASCADEpaper/paper/Fig5_archetype/5_integrated_final_module_meta.Rds")
srt <- AddMetaData(srt, meta)
source("~/CASCADEpaper/paper/functions.R")
f1 <- coneraxes(FeaturePlot(srt, "Module1", cols =  rev(brewer.rdylbu(10)),raster = F)+ggtitle("AR")+NoLegend()) 
f2 <- rmaxes(FeaturePlot(srt, "Module2", cols =  rev(brewer.rdylbu(10)),raster = F)+ggtitle("Inflammation")+NoLegend()) 
f3 <- rmaxes(FeaturePlot(srt, "Module3", cols =  rev(brewer.rdylbu(10)),raster = F)+ggtitle("NE1")+NoLegend()) 
f4 <- rmaxes(FeaturePlot(srt, "Module4", cols =  rev(brewer.rdylbu(10)),raster = F)+ggtitle("NE2")+NoLegend()) 
f5 <- rmaxes(FeaturePlot(srt, "Module5", cols =  rev(brewer.rdylbu(10)),raster = F)+ggtitle("Cycling")+NoLegend()) 
f6 <- rmaxes(FeaturePlot(srt, "Module6", cols =  rev(brewer.rdylbu(10)),raster = F)+ggtitle("Hypoxia")+NoLegend()) 

pdf("pdf/integrated_expression.pdf", width = 20, height = 4)
print(f1+f2+f3+f4+f5+f6+plot_layout(nrow = 1, ncol = 6))
dev.off()
