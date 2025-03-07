library(Seurat)

srt <- readRDS("~/integration/ruvIII/integrated_ruv/tumor_only_ruvsrt.Rds")
srt$patient <- substr(srt$sample, 1, 6)
Hawley_signature <- readRDS("~/CASCADEpaper/paper/Fig2/Hawley_signature.Rds")
genelist <- split(Hawley_signature, Hawley_signature$cluster)
genelist <- lapply(genelist, function(x) x[, "gene"])
genelist <- lapply(genelist, function(x) x[x%in% rownames(srt)])
genelist <- genelist[sapply(genelist, length)>0]


srt <- AddModuleScore(srt, features = genelist)
srt@meta.data[, grepl("Cluster", colnames(srt@meta.data))] <- names(genelist)
saveRDS(srt@meta.data, "~/CASCADEpaper/paper/Fig2/Hawley_signature_meta.Rds")

pdf("~/CASCADEpaper/paper/Fig2/Hawley_signature.pdf", width = 6, height = 5)
for (i in names(genelist)){
  f <- FeaturePlot(srt, raster = F, i)
  print(f)
}
dev.off()