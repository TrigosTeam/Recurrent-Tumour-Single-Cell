library(Seurat)
library(Signac)
library(qs)
library(dplyr)
library(ggplot2)
library(pals)
library(patchwork)
library(ggrastr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(data.table)
library(JASPAR2024)
library(TFBSTools)
source("~/ATAC/functions.R")
setwd("~/ATAC/integration")
annotation <- readRDS("~/ATAC/annotation.Rds")
srt <- qread("integrated_srt.qs")
DefaultAssay(srt) <- "ruv3"
Annotation(srt[["ruv3"]]) <- annotation
gene_activity <- GeneActivity(srt)
saveRDS(gene_activity, "gene_activity.Rds")
srt[["gene_activity"]] <- CreateAssayObject(counts = gene_activity)
srt <- NormalizeData(
  object = srt,
  assay = "gene_activity",
  normalization.method = "LogNormalize",
  scale.factor = median(srt$nCount_gene_activity)
)

DefaultAssay(srt) <- "gene_activity"

clean_module <- readRDS("~/CASCADEpaper/paper/Fig5_archetype/clean_module.Rds")
names(clean_module) <- c("AR", "Inflammation", "NE1", "NE2", "Cycling", "Glycolysis")

srt <- AddModuleScore(srt, features = clean_module, assay = "gene_activity")
colnames(srt@meta.data)[grep("Cluster",colnames(srt@meta.data))] <- paste(names(clean_module), "Module")

f1 <- FeaturePlot(srt, paste(names(clean_module), "Module"), pt.size = 0.1, reduction = "umap.peak.ruv3", raster = F, ncol = 1,alpha = 1)& 
  scale_colour_gradientn(colors = rev(brewer.rdylbu(10)))
f <- rmaxes(f1, print = F)

rasterize_all_patches <- function(patch_obj, dpi = 300) {
  # Loop through every plot in the patchwork object
  for (i in seq_along(patch_obj)) {
    # Check if the item is a ggplot object before trying to rasterize
    if (inherits(patch_obj[[i]], "gg")) {
      patch_obj[[i]] <- rasterise(patch_obj[[i]], layers = "Point", dpi = dpi)
    }
  }
  return(patch_obj)
}
f <- rasterize_all_patches(f, dpi = 400)

p <- DimPlot(srt, group.by = "patient",pt.size = 0.1, reduction = "umap.peak.ruv3", raster = F, label = T)+ ggtitle("Patient")+NoLegend()
p <-  coneraxes(p, print = T, labelx = "ruv3.peak UMAP1", labely =  "ruv3.peak UMAP2")

pdf("pdf/module_gene_activity_umap.peak.ruv3.pdf", width = 18, height = 4)
f
p+f+plot_layout(widths = c(1, 6))
dev.off()
pdf("pdf/2module_gene_activity_umap.peak.ruv3.pdf", width = 3, height = 12)
f
dev.off()
DefaultAssay(srt) <- "RNA"
srt <- AddModuleScore(srt, features = clean_module, assay = "RNA")
colnames(srt@meta.data)[grep("Cluster",colnames(srt@meta.data))] <- paste(names(clean_module), "Module-RNA")

saveRDS(srt@meta.data, "~/ATAC/integration/meta.rds")
plotlist <- list()
for( m in names(clean_module)){
  df <- data.frame(Module = m, RNA =  srt@meta.data[,paste(m, "Module-RNA")] , Gene_activity = srt@meta.data[,paste(m, "Module.2")])
  
  p <- ggplot(df, aes(x = RNA, y = Gene_activity)) +
    rasterise(geom_point(alpha = 0.6, color = "darkblue", size = 0.2), dpi = 400) +  # Points with slight transparency
    geom_smooth(method = "lm", color = "red", fill = "pink") + # Linear regression line
    
    # Automatically add correlation coefficient and p-value
    stat_cor(method = "pearson", # Change to "spearman" if data is non-normal
             label.x.npc = "left", 
             label.y.npc = "top",
             size = 5) +
    
    theme_bw(base_size = 18) +
    labs(title = paste("Correlation of", unique(df$Module)),
         x = "RNA Expression Score",
         y = "Gene Activity Score")
  plotlist[[m]] <- p
}

pdf("pdf/correlation_gene_activity_vs_RNA.pdf", width = 12, height = 7)
wrap_plots(plotlist, ncol = 3, axes = "collect")
dev.off()
