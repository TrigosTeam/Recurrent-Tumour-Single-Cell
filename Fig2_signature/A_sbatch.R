library(Seurat)
library(ggplot2)
library(dplyr)
library(data.table)
library(scran)
library(patchwork)
library(pals)
library(scales)
library(ggpubr)
library(stringr)
load("/trigos_team/CASCADE/Analysis/signatures_analysis/objects/all_gene_sets.rda")
source("~/CASCADEpaper/paper/functions.R")
setwd("~/CASCADEpaper/paper/Fig2_signature")
genelist <- lapply(all_gene_sets, function(x) unique(unlist(x))) # union of different signatures

tang_2022 <- all_gene_sets$tang_2022
names(tang_2022) <- paste0("tang_2022_", names(tang_2022))

temp <- lapply(all_gene_sets, function(x){
  names(which(table(unlist(x)) >= length(x)/2))
})

temp <- temp[!names(temp)%in% c("CRPC", "cancer", "hillock", "metastasis", "club", "invasion", "tang_2022")]
genelist <- c(temp, tang_2022)
names(genelist) <- paste0(names(genelist), "_signature")


# tumor normal all --------------------------------------------------------
srt <- readRDS("~/integration/2024_06/tumor_normal_ruvsrt.Rds")
genelist2 <- lapply(genelist, function(x) x[x%in% rownames(srt)])

tb <- as.data.frame.matrix(table(srt$seurat_clusters, srt$Cell.Identity))
temp <- apply(tb, 1, function(x) names(which.max(x)))
srt$cluster_identity <- temp[srt$seurat_clusters]
srt$cluster_identity[srt$cluster_identity %in% c("Epithelial cells", "Neurons") & srt$seurat_clusters != "22"] <- srt$patient[srt$cluster_identity %in% c("Epithelial cells", "Neurons") & srt$seurat_clusters != "22"]


srt <- AddModuleScore(srt, features = genelist2)
colnames(srt@meta.data)[grep("Cluster", colnames(srt@meta.data))] <- names(genelist)

saveRDS(srt@meta.data, "signature_meta_tumor_normal.Rds")
pdf("pdf/signature_tumor_normal.pdf", width = 6, height = 5)
DimPlot(srt, group.by = "cluster_identity", raster = F)
DimPlot(srt, group.by = "cluster_identity", raster = F, label = T) + NoLegend()
for(i in names(genelist)){
f <- FeaturePlot(srt, raster = F, features = i)
print(f)
}
dev.off()

# rm(list = "srt")
# gc()

srt$cell_anno <- gsub("00", "", srt$cell_anno)
srt$cell_anno <- str_wrap(srt$cell_anno, 20)


cols <- c(brewer.set1(14), hue_pal()(9))

g <- DimPlot(srt, group.by = "cell_anno", raster = F)+ 
  scale_color_manual(values = cols, breaks  = unique(srt$cell_anno) )+
  ggtitle("")+NoLegend()
g <- coneraxes(g)
ggsave("pdf/tumor_normal_cellanno.pdf", width = 5, height = 5)


g1 <- DimPlot(srt, group.by = "cell_anno", raster = F, cells = colnames(srt)[srt$cell_anno %in% unique(srt$cell_anno)[1:14]])+ 
  scale_color_manual(values = cols[1:14], breaks  = unique(srt$cell_anno)[1:14], name = "Normal cells")+
  theme(legend.text = element_text(size = 16), legend.title = element_text(size = 18))

leg1 <- as_ggplot(get_legend(g1))
g2 <- DimPlot(srt, group.by = "cell_anno", raster = F, cells = colnames(srt)[!srt$cell_anno %in% unique(srt$cell_anno)[1:14]])+ 
  scale_color_manual(values = cols[15:23], breaks  = unique(srt$cell_anno)[15:23], name = "Tumour cells")+
  theme(legend.position = "bottom",legend.text = element_text(size = 16), legend.title = element_text(size = 18))
leg2 <- as_ggplot(get_legend(g2))

saveRDS(g, "plot_obj/tumor_normal_cellanno.Rds")
saveRDS(leg1, "plot_obj/tumor_normal_cellannoleg1.Rds")
saveRDS(leg2, "plot_obj/tumor_normal_cellannoleg2.Rds")

ggsave("pdf/tumor_normal_leg1.pdf", leg1, width = 2, height = 5)
ggsave("pdf/tumor_normal_leg2.pdf", leg2, width = 7, height = 1)

# tumor only all ----------------------------------------------------------

srt <- readRDS("~/integration/2024_06/tumor_only_ruvsrt.Rds")
genelist2 <- lapply(genelist, function(x) x[x%in% rownames(srt)])

srt <- AddModuleScore(srt, features = genelist2)
colnames(srt@meta.data)[grep("Cluster", colnames(srt@meta.data))] <- names(genelist)

pdf("~/CASCADEpaper/paper/Fig2_signature/signature_tumor_only.pdf", width = 6, height = 5)
DimPlot(srt, group.by = "patient", raster = F)
for(i in names(genelist)){
  f <- FeaturePlot(srt, raster = F, features = i)
  print(f)
}
dev.off()

rm(list = "srt")
gc()
meta <- readRDS("~/CASCADEpaper/paper/Fig2_signature/signature_meta_tumor_only.Rds")
srt <- AddMetaData(srt, meta)

f1 <- rmaxes(FeaturePlot(srt, raster = F, features =  c("AR_signature"))+ggtitle("AR")) 
f2 <- coneraxes(FeaturePlot(srt, raster = F, features = "NE_signature")+ggtitle("NE"))
f3<- rmaxes(FeaturePlot(srt, raster = F, features = "cell_cycle_signature")+ggtitle("cell cycle"))
f4 <- rmaxes(FeaturePlot(srt, raster = F, features = "GI_signature_signature")+ggtitle("GI subtype"))
f5 <- rmaxes(FeaturePlot(srt, raster = F, features = "tang_2022_WNT_markers_signature")+ggtitle("Tang 2022 WNT"))
f6 <- rmaxes(FeaturePlot(srt, raster = F, features = "tang_2022_stem_cell_markers_signature")+ggtitle("Tang 2022 stem-cell-like"))

f1 +f2+f3+f4+f5+f6 + plot_layout(ncol = 3, nrow = 2, byrow = F)
ggsave("pdf/tumor_only_signatures.pdf", width = 10, height = 6.5)


srt$patient <- gsub("00", "", srt$patient)
srt$cell_anno <- srt$patient 
srt$cell_anno[srt$sample == "CA0058_hilar_50"] <- "CA58 hilar"
srt$cell_anno[srt$sample %in% c("CA0058_liver_29", "CA0058_liver_38") ] <- "CA58 livers"

g <- coneraxes(DimPlot(srt, group.by = "cell_anno", raster = F, label = T, label.size = 10)+NoLegend()+ggtitle("Tumour only cells"))

ggsave("pdf/tumor_only_cellanno.pdf", width = 6, height = 6)

g <- FeaturePlot(srt, c("AR", "ASCL1"), blend = T, cols = c("lightgrey","red", "blue" ))
f3 <- coneraxes(g[[3]])
leg <- g[[4]]+theme_minimal(base_size = 13)+ggtitle("Color\nthreshold:0.5")
free(f3)+leg + plot_layout(widths = c(4.5, 1))
ggsave("pdf/AR_ASCL1_blend.pdf", height =  6 , width = 6.5)


# cycling cell proportion -------------------
signature_meta_tumor_only <- readRDS("~/CASCADEpaper/paper/Fig2_signature/signature_meta_tumor_only.Rds")
bins <- cut(signature_meta_tumor_only$cell_cycle_signature, breaks = 100)
signature_meta_tumor_only$cycling<- bins %in% levels(bins)[50:100]

signature_meta_tumor_only$pathology <- ifelse(signature_meta_tumor_only$patient %in% c("CA0027", "CA0058"), "Mixed", ifelse(signature_meta_tumor_only$patient %in% c("CA0046", "CA0090"), "NE", "AD"))
df <- signature_meta_tumor_only %>% group_by(pathology, sample, cycling) %>% summarise(n = n()) %>% mutate(prop = n/sum(n))

cyclings <- c(464, 728)
pathology <- c(78296, 34791)
prop.test(cyclings, pathology, alternative = "less")

df <- filter(df, cycling == TRUE)
ggplot(df, aes(x = pathology, y = prop, fill = pathology)) + 
  geom_boxplot()+
  geom_jitter()+
  theme_bw(base_size = 14)
ggsave("pdf/cycling_proption_boxplot.pdf", width = 5, height = 5)

my_comparisons <- list( c("AD", "NE"), c("AD", "Mixed"), c("Mixed", "NE") )
ggpubr::ggboxplot(df, x = "pathology", y = "prop",
          color = "pathology", palette = "jco",add = "jitter")+ 
  stat_compare_means(comparisons = my_comparisons)

ggsave("pdf/cycling_proption_boxplot.pdf", width = 5, height = 5)

ggpubr::ggboxplot(signature_meta_tumor_only, x = "pathology", y = "cell_cycle_signature",
                  color = "pathology", palette = "jco",add = "jitter")+ 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox")
ggsave("pdf/cycling_expression_boxplot.pdf", width = 5, height = 5)

# ruv per patient signature -----------------------------------------------
paths <- system("realpath ~/integration/2024_06/perpatient/tumor_only/*_ruvsrt.Rds", intern = T)
paths_split <- substr(unlist(lapply(strsplit(paths, split = "/"), function(x) x[8])),1 ,6) 
paths <- setNames(paths, paths_split)



for (p in paths_split){
  perp_srt <- readRDS(paths[p])
  temp <- genelist[sapply(genelist, function(x) any(x %in% rownames(perp_srt)))]
  perp_srt <- AddModuleScore(perp_srt, features = temp)
  colnames(perp_srt@meta.data)[grep("Cluster", colnames(perp_srt@meta.data))] <- names(temp)
  saveRDS(perp_srt@meta.data, paste0("~/CASCADEpaper/paper/Fig2_signature/perpatient_signature_tumor_only/",p, "_meta.Rds"))
  pdf(paste0("~/CASCADEpaper/paper/Fig2_signature/perpatient_signature_tumor_only/",p, "_signature.pdf"), width = 4, height = 3)
  for (i in names(temp)){
  f <- FeaturePlot(perp_srt, raster = F, features = i)
  print(f)
  }
  dev.off()
  print(p)
}
signature_meta_tumor_only <- readRDS("~/CASCADEpaper/paper/Fig2_signature/signature_meta_tumor_only.Rds")
source("~/CASCADEpaper/paper/cols.R")
for(p in paths_split){
  perp_srt <- readRDS(paths[p])
  perp_srt <- AddMetaData(perp_srt, metadata = readRDS(paste0("~/CASCADEpaper/paper/Fig2_signature/perpatient_signature_tumor_only/",p, "_meta.Rds")))
  g <- DimPlot(perp_srt, group.by = "site", cols = site_cols ) + ggtitle(gsub("00", "", p))+
    theme(legend.position = "top",legend.text=element_text(size=15))
  g <- coneraxes(g)
  
  if(all(c("AR_signature", "NE_signature","cell_cycle_signature" ,"GI_signature_signature" ,"tang_2022_WNT_markers_signature","tang_2022_stem_cell_markers_signature")%in% colnames(perp_srt@meta.data))){
    f1 <- rmaxes(FeaturePlot(perp_srt, raster = F, features =  c("AR_signature"))+ggtitle("AR"), front.size =14) 
    f2 <- coneraxes(FeaturePlot(perp_srt, raster = F, features = "NE_signature")+ggtitle("NE"), front.size =14)
    f3<- rmaxes(FeaturePlot(perp_srt, raster = F, features = "cell_cycle_signature")+ggtitle("cell cycle"), front.size =14)
    f4 <- rmaxes(FeaturePlot(perp_srt, raster = F, features = "GI_signature_signature")+ggtitle("GI subtype"), front.size =14)
    f5 <- rmaxes(FeaturePlot(perp_srt, raster = F, features = "tang_2022_WNT_markers_signature")+ggtitle("Tang 2022 WNT"), front.size =14)
    f6 <- rmaxes(FeaturePlot(perp_srt, raster = F, features = "tang_2022_stem_cell_markers_signature")+ggtitle("Tang 2022 stem-cell-like"), front.size =14)
    f <- f1 +f2+f3+f4+f5+f6 + plot_layout(ncol = 3, nrow = 2, byrow = F)
    }else{
    perp_srt@meta.data[, setdiff(c("AR_signature", "NE_signature","cell_cycle_signature" ,"GI_signature_signature" ,"tang_2022_WNT_markers_signature","tang_2022_stem_cell_markers_signature"), grep("_signature", colnames(perp_srt@meta.data), value = T))] <- 0
    f1 <- rmaxes(FeaturePlot(perp_srt, raster = F, features =  c("AR_signature"))+ggtitle("AR"), front.size =14) 
    f2 <- coneraxes(FeaturePlot(perp_srt, raster = F, features = "NE_signature")+ggtitle("NE"), front.size =14)
    f3<- rmaxes(FeaturePlot(perp_srt, raster = F, features = "cell_cycle_signature")+ggtitle("cell cycle"), front.size =14)
    f4 <- rmaxes(FeaturePlot(perp_srt, raster = F, features = "GI_signature_signature")+ggtitle("GI subtype"), front.size =14)
    f5 <- rmaxes(FeaturePlot(perp_srt, raster = F, features = "tang_2022_WNT_markers_signature")+ggtitle("Tang 2022 WNT"), front.size =14)
    f6 <- rmaxes(FeaturePlot(perp_srt, raster = F, features = "tang_2022_stem_cell_markers_signature")+ggtitle("Tang 2022 stem-cell-like"), front.size =14)
    f <- f1 +f2+f3+f4+f5+f6 + plot_layout(ncol = 3, nrow = 2, byrow = F)
  }
  
  
  tot <- (free(g) |f) +plot_layout(widths = c(1, 2.5))
  # saveRDS(tot, paste0("~/CASCADEpaper/paper/Fig2_signature/plot_obj/", p, ".Rds"))
  pdf(paste0("pdf/", p, ".pdf"), width = 9.5, height = 4)
  print(tot)
  dev.off()
  print(p)
}

# per sample signatures ---------------------------------------------------
paths <- system("realpath /trigos_team/CASCADE/Analysis/230113_seurat_intron/*/*/*non_tumour_removed.Rds", intern = T)
paths_split <- unlist(lapply(strsplit(paths, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
paths2 <- system("realpath /trigos_team/CASCADE/Analysis/240603_seurat_intron/CA0027/dura*/*non_tumour_removed_v4.Rds", intern = T)
paths_split2 <- unlist(lapply(strsplit(paths2, split = "/"), function(x) paste(x[6], x[7], sep = "_")))

paths <- c(paths, paths2)
paths_split <- c(paths_split, paths_split2)
paths <- setNames(paths, paths_split)

for(i in paths_split){
  srt <- readRDS(paths[i])

  DefaultAssay(srt) <- "RNA"
  lib.median <- median(srt$nCount_RNA)
  target_pseudocount <- 1
  srt <- srt %>% 
    NormalizeData(scale.factor = lib.median / target_pseudocount) %>% 
    ScaleData(do.scale = F)
  
 set.seed(100)
  poisson_fit2 <- modelGeneVarByPoisson(srt@assays$RNA@counts, dispersion = 0.2, size.factors = target_pseudocount * srt$nCount_RNA / lib.median)
  residuals2 <- poisson_fit2$total - poisson_fit2@metadata$trend(poisson_fit2$mean)
  names(residuals2) <- rownames(srt)
  top_genes2 <- rownames(srt)[order(residuals2, decreasing = TRUE)]
  
  srt <- srt %>% 
    RunPCA(features = top_genes2[1:2000]) %>% 
    FindNeighbors(dims = 1:20) %>% 
    FindClusters(resolution = 1) %>% 
    RunUMAP(dims = 1:20)
  
  temp <- genelist[sapply(genelist, function(x) any(x %in% rownames(srt)))]
  srt <- AddModuleScore(srt, features = temp)
  colnames(srt@meta.data)[grep("Cluster", colnames(srt@meta.data))] <- names(temp)
  saveRDS(srt@meta.data, paste0("~/CASCADEpaper/paper/Fig2_signature/persample_signature/",i, "_meta.Rds"))
  pdf(paste0("~/CASCADEpaper/paper/Fig2_signature/persample_signature/",i, ".pdf"), width = 4, height = 3)
  for (n in names(temp)){
    f <- FeaturePlot(srt, raster = F, features = n)
    print(f)
  }
  dev.off()
  print(i)
  
}

srt <- AddMetaData(srt, meta)
f1 <- rmaxes(FeaturePlot(srt, raster = F, features =  c("AR_signature"))+ggtitle("AR")) 
f2 <- coneraxes(FeaturePlot(srt, raster = F, features = "NE_signature")+ggtitle("NE"))
f3<- rmaxes(FeaturePlot(srt, raster = F, features = "cell_cycle_signature")+ggtitle("cell cycle"))
f4 <- rmaxes(FeaturePlot(srt, raster = F, features = "GI_signature_signature")+ggtitle("GI subtype"))
f5 <- rmaxes(FeaturePlot(srt, raster = F, features = "tang_2022_WNT_markers_signature")+ggtitle("Tang 2022 WNT"))
f6 <- rmaxes(FeaturePlot(srt, raster = F, features = "tang_2022_stem_cell_markers_signature")+ggtitle("Tang 2022 stem-cell-like"))

f1 +f2+f3+f4+f5+f6 + plot_layout(ncol = 3, nrow = 2, byrow = F)
ggsave(paste0("~/CASCADEpaper/paper/Fig2_signature/pdf/",i, "_signatures.pdf"),  width = 9.5, height = 4)


# saveRDS(g2, paste0("~/CASCADEpaper/paper/Fig2_signature/plot_obj/", i, ".Rds"))
# png(paste0("~/CASCADEpaper/paper/Fig2_signature/png/", i, ".png"), width = 700, height = 300)
# print(g2)
# dev.off()
