library(ggplot2)
library(pals)
library(scales)
library(Seurat)
library(colorRamps)
library(ggpubr)
library(geomtextpath)
library(ggpubr)
library(dplyr)
library(circlize)
library(patchwork)
library(data.table)
library(MAST)
setwd("~/CASCADEpaper/paper/PSMA")

srt <- readRDS("~/integration/2024_06/tumor_only_ruvsrt.Rds")
phenotype_meta <- readRDS("~/CASCADEpaper/paper/PSMA/phenotype_meta.Rds")
srt <- AddMetaData(srt, phenotype_meta)

srt$pheno <- paste0(srt$AR_pheno, srt$FOLH1_pheno)
Idents(srt) <- srt$group
com1 <- combn(grep("ARhigh", unique(srt$pheno), value = T), 2, simplify = F)
deg1 <- lapply(com1, function(x){
  marker <- FindMarkers(srt, ident.1 =  x[1], ident.2 = x[2], test.use = "MAST", 
                        min.cells.feature = 20,min.cells.group = 20, group.by = "pheno", subset.ident = c("AR+/NE-", "AR+/NE+"))
  marker$group1 <- x[1]
  marker$group2 <- x[2]
  return(marker)
})
saveRDS(deg1, "ARhighDEG.Rds")

com2 <- combn(grep("ARlow", unique(srt$pheno), value = T), 2, simplify = F)

deg2 <- lapply(com2, function(x){
  marker <- FindMarkers(srt, ident.1 =  x[1], ident.2 = x[2], test.use = "MAST", 
                        min.cells.feature = 20,min.cells.group = 20, group.by = "pheno", subset.ident = c("AR+/NE-", "AR+/NE+"))
  marker$group1 <- x[1]
  marker$group2 <- x[2]
  return(marker)
})
saveRDS(deg2, "ARlowDEG.Rds")

com3 <- combn(grep("ARneg", unique(srt$pheno), value = T), 2, simplify = F)

deg3 <- lapply(com3, function(x){
  marker <- FindMarkers(srt, ident.1 =  x[1], ident.2 = x[2], test.use = "MAST", 
                        min.cells.feature = 20,min.cells.group = 20, group.by = "pheno", subset.ident = c("AR+/NE-", "AR+/NE+"))
  marker$group1 <- x[1]
  marker$group2 <- x[2]
  return(marker)
})
saveRDS(deg3, "ARnegDEG.Rds")



