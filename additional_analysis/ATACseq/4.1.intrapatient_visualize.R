library(ggplot2)
library(patchwork)
library(qs)
library(Seurat)
setwd("~/ATAC/intrapatient")
patients <- list.files(".")
source("~/ATAC/functions.R")
source("~/CASCADEpaper/paper/cols.R")
clean_module <- readRDS("~/CASCADEpaper/paper/Fig5_archetype/clean_module.Rds")
names(clean_module) <- c("AR", "Inflammation", "NE1", "NE2", "Cycling", "Glycolysis")

load("/trigos_team/CASCADE/Analysis/signatures_analysis/objects/all_gene_sets.rda")
genelist <- lapply(all_gene_sets, function(x) unique(unlist(x))) # union of different signatures

tang_2022 <- all_gene_sets$tang_2022
names(tang_2022) <- paste0("tang_2022_", names(tang_2022))

temp <- lapply(all_gene_sets, function(x){
  names(which(table(unlist(x)) >= length(x)/2))
})

temp <- temp[!names(temp)%in% c("CRPC", "cancer", "hillock", "metastasis", "club", "invasion", "tang_2022")]
genelist <- c(temp, tang_2022)
names(genelist) <- paste0(names(genelist), "_signature")

# signature ploting
for (p in patients){
  srt <- qread(paste0(p, "/integrated_srt.qs"))
  f <- coneraxes(DimPlot(srt, raster = F, reduction = "umap.peak.ruv3",group.by = "site",label = T)+scale_color_manual(values = site_cols)+ggtitle(paste(p, "ATACseq"))+NoLegend(), front.size =16, print = F)
  ggsave(paste0(p, "/pdf/umap.peak.ruv3.bysite.pdf"),width =2.4, height = 4, plot = f)
  
  DefaultAssay(srt) <- "gene_activity"
  temp <- genelist[sapply(genelist, function(x) any(x %in% rownames(srt)))]
  srt <- AddModuleScore(srt, features = temp)
  colnames(srt@meta.data)[grep("Cluster", colnames(srt@meta.data))] <- names(temp)
  
  if(all(c("AR_signature", "NE_signature","cell_cycle_signature" ,"GI_signature_signature" ,"tang_2022_WNT_markers_signature","tang_2022_stem_cell_markers_signature")%in% colnames(srt@meta.data))){
    f1 <- rmaxes(FeaturePlot(srt, raster = F, reduction = "umap.peak.ruv3",features =  c("AR_signature"))+ggtitle("AR"), front.size =14) 
    f2 <- coneraxes(FeaturePlot(srt, raster = F, reduction = "umap.peak.ruv3",features = "NE_signature")+ggtitle("NE"), front.size =14, labelx = "ruv3.peak UMAP1", labely = "ruv3.peak UMAP2")
    f3<- rmaxes(FeaturePlot(srt, raster = F, reduction = "umap.peak.ruv3",features = "cell_cycle_signature")+ggtitle("cell cycle"), front.size =14)
    f4 <- rmaxes(FeaturePlot(srt, raster = F, reduction = "umap.peak.ruv3",features = "GI_signature_signature")+ggtitle("GI subtype"), front.size =14)
    f5 <- rmaxes(FeaturePlot(srt, raster = F, reduction = "umap.peak.ruv3",features = "tang_2022_WNT_markers_signature")+ggtitle("Tang 2022 WNT"), front.size =14)
    f6 <- rmaxes(FeaturePlot(srt, raster = F, reduction = "umap.peak.ruv3",features = "tang_2022_stem_cell_markers_signature")+ggtitle("Tang 2022 stem-cell-like"), front.size =14)
    f <- (f1 +f2+f3+f4+f5+f6)  + plot_layout(ncol = 3, nrow = 2, byrow = F)&scale_color_gradient2(mid = "grey80", high = "darkred", low = "blue", midpoint = 0)
  }else{
    srt@meta.data[, setdiff(c("AR_signature", "NE_signature","cell_cycle_signature" ,"GI_signature_signature" ,"tang_2022_WNT_markers_signature","tang_2022_stem_cell_markers_signature"), grep("_signature", colnames(srt@meta.data), value = T))] <- 0
    f1 <- rmaxes(FeaturePlot(srt, raster = F, reduction = "umap.peak.ruv3",features =  c("AR_signature"))+ggtitle("AR"), front.size =14) 
    f2 <- coneraxes(FeaturePlot(srt, raster = F, reduction = "umap.peak.ruv3",features = "NE_signature")+ggtitle("NE"), front.size =14)
    f3<- rmaxes(FeaturePlot(srt, raster = F, reduction = "umap.peak.ruv3",features = "cell_cycle_signature")+ggtitle("cell cycle"), front.size =14)
    f4 <- rmaxes(FeaturePlot(srt, raster = F, reduction = "umap.peak.ruv3",features = "GI_signature_signature")+ggtitle("GI subtype"), front.size =14)
    f5 <- rmaxes(FeaturePlot(srt, raster = F, reduction = "umap.peak.ruv3",features = "tang_2022_WNT_markers_signature")+ggtitle("Tang 2022 WNT"), front.size =14)
    f6 <- rmaxes(FeaturePlot(srt, raster = F, reduction = "umap.peak.ruv3",features = "tang_2022_stem_cell_markers_signature")+ggtitle("Tang 2022 stem-cell-like"), front.size =14)
    f <- (f1 +f2+f3+f4+f5+f6)  + plot_layout(ncol = 3, nrow = 2, byrow = F)&scale_color_gradient2(mid = "grey80", high = "darkred", low = "blue", midpoint = 0)
  
  }
  ggsave(paste0(p, "/pdf/Fig2signature.pdf"),width = 6.5, height = 4, plot = f)
  
  srt <- AddModuleScore(srt,clean_module)
  colnames(srt@meta.data)[grep("Cluster", colnames(srt@meta.data))] <- paste(names(clean_module), "Module")
  p1<- rmaxes(FeaturePlot(srt, raster = F, reduction = "umap.peak.ruv3",features =  "AR Module"), front.size =14) 
  p2<- coneraxes(FeaturePlot(srt, raster = F, reduction = "umap.peak.ruv3",features =  "Inflammation Module"), front.size =14, labelx = "ruv3.peak UMAP1", labely = "ruv3.peak UMAP2") 
  p3<- rmaxes(FeaturePlot(srt, raster = F, reduction = "umap.peak.ruv3",features =  "NE1 Module"), front.size =14) 
  p4<- rmaxes(FeaturePlot(srt, raster = F, reduction = "umap.peak.ruv3",features =  "NE2 Module"), front.size =14) 
  p5<- rmaxes(FeaturePlot(srt, raster = F, reduction = "umap.peak.ruv3",features =  "Cycling Module"), front.size =14) 
  p6<- rmaxes(FeaturePlot(srt, raster = F, reduction = "umap.peak.ruv3",features =  "Glycolysis Module" ), front.size =14) 
  g <- p1 +p2+p3+p4+p5+p6 + plot_layout(ncol = 3, nrow = 2, byrow = F)&scale_color_gradient2(mid = "grey80", high = "darkred", low = "blue", midpoint = 0)
  ggsave(paste0(p, "/pdf/Modulesignature.pdf"),width = 9.5, height = 4, plot = g)
  
  saveRDS(srt@meta.data, paste0(p,"/meta.Rds"))
  cat(p, "\n")
}
