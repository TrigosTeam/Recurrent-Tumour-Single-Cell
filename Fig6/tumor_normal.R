
# library -----------------------------------------------------------------
library(Seurat)
library(dplyr)
library(harmony)
library(data.table)
library(SingleR)
library(patchwork)
library(ggplot2)
refs <- readRDS("~/integration/tumor_normal_seperation/atlasref.Rds")

tumor_normal <- readRDS("/trigos_team/CASCADE/Analysis/integration/ruvIII/integrated_all_sample/integrated_ruv_tumor_normal/tumor_normal_ruvsrt.Rds")

tumor_normal$Cluster.Identity[tumor_normal$Cluster.Identity == "tumor" & tumor_normal$patient %in% c("CA0034", "CA0035", "CA0043", "CA0076", "CA0083")] <- "Adenocarcinoma"
tumor_normal$Cluster.Identity[tumor_normal$Cluster.Identity == "tumor" & tumor_normal$patient %in% c("CA0027", "CA0058")] <- "Mixing"
tumor_normal$Cluster.Identity[tumor_normal$Cluster.Identity == "tumor" & tumor_normal$patient %in% c("CA0046", "CA0090")] <- "Neuroendocrine"

# DimPlot(tumor_normal, raster = F, group.by = "Cluster.Identity")
# DimPlot(tumor_normal, raster = F, group.by = "seurat_clusters", label = T) + NoLegend()
Idents(tumor_normal) <- tumor_normal$seurat_clusters
normal <- subset(tumor_normal, idents = c("33", "35", "12", "24", "22", "20", "14", "29", "27", "5", "23"))
DimPlot(normal, raster = F)
saveRDS(normal, "~/CASCADEpaper/paper/Fig2/normal_pre_ruv.Rds")




# preliminary infercnv test include tumor the rest should meet CNV cut off -----------------------------------------------
# rm(list = "tumor_normal")
# 
# normal[["RNA"]] <- normal[["ruv3"]]
# DefaultAssay(normal) <- "RNA"
# normal <- DietSeurat(normal, counts = T, data = T, assays = "RNA")
# 
# 
# normal <- FindVariableFeatures(normal, selection.method = "vst", nfeatures = 2000)
# normal <- ScaleData(normal, features = rownames(normal))
# normal <- RunPCA(normal, features = VariableFeatures(object = normal))
# normal <- RunHarmony(normal, "patient", assay.use = "RNA")
# normal <- normal %>% 
#   RunUMAP(reduction = "harmony", dims = 1:20) %>% 
#   FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
#   FindClusters(resolution = 0.8) %>% 
#   identity()
# saveRDS(normal, "~/CASCADEpaper/paper/Fig2/normal_harmony_srt.Rds")
# 
# DimPlot(normal)
# 
# library(SingleR)
# refs <- readRDS("~/integration/tumor_normal_seperation/atlasref.Rds")
# pred <- SingleR(normal[["RNA"]]@data, refs, labels = list(refs[[1]]$label.main,refs[[2]]$label.main))
# normal[["Cell.Identity"]] <- pred[colnames(normal), "labels"]
# ids <- unique(normal$Cell.Identity)
# ids <- gsub("_", " ", fixed = T, ids)
# ids <- gsub("-", " ", fixed = T, ids)
# ids <- gsub("s$", "", ids)
# ids <- setNames(ids, unique(normal$Cell.Identity))
# normal[["Cell.Identity"]] <- ids[normal$Cell.Identity]
# DimPlot(normal, group.by = "Cell.Identity")
# meta <- normal@meta.data
# max_mat <- meta %>% group_by(seurat_clusters, Cell.Identity) %>% summarise(n = n())
# max_mat <- max_mat %>% group_by(seurat_clusters) %>% mutate(max_id = max(n)) %>% mutate(max_cell = Cell.Identity[n == max_id])
# max_mat <- max_mat[, c("seurat_clusters", "max_cell")] %>% distinct()
# max_id <- structure(max_mat$max_cell, names = max_mat$seurat_clusters)
# normal$max_id <- max_id[normal$seurat_clusters]
# 
# DimPlot(normal, group.by = "max_id")
# 
# final_normal_v1 <- readRDS("~/integration/tumor_normal_seperation/scATOMIC/final_normal_v1.Rds")
# 
# anno <- as.data.frame(as.character(normal$seurat_clusters))
# rownames(anno) <- colnames(normal)
# anno_1 <- readRDS("~/integration/tumor_normal_seperation/normal_anno_by_sample_id.Rds")
# rownames(anno_1) <- substr(rownames(anno_1), 2, nchar(rownames(anno_1)))
# ref_cells <- rownames(anno_1)[anno_1[, 1] %in% c("B cell", "T cell", "Plasma B cell", "Macrophage", "Monocyte")]
# ref_cells <- ref_cells[ref_cells %in% colnames(normal)]
# 
# anno[ref_cells, 1] <- "ref"
# colnames(anno) <- NULL
# 
# saveRDS(anno, "normal_harmony_anno.Rds")
# 
# 
# 
# library(phylogram)
# library(dendextend)
# dend <- read.dendrogram("~/CASCADEpaper/paper/Fig2/normal_cells/infercnv/unclustered/infercnv.preliminary.observations_dendrogram.txt")
# dend <- color_branches(dend, h = 100)
# plot(dend)
# 
# tree <- cutree(dend, h = 100)
# rm_cells <- names(tree)[tree == 1]
# normal_cells <- names(tree)[tree ==2 ]
# normal_clean <- subset(normal, cells = colnames(normal)[!colnames(normal) %in%rm_cells])
# normal_clean <- FindVariableFeatures(normal_clean, selection.method = "vst", nfeatures = 2000)
# normal_clean <- ScaleData(normal_clean, features = rownames(normal_clean))
# normal_clean <- RunPCA(normal_clean, features = VariableFeatures(object = normal_clean))
# normal_clean <- RunHarmony(normal_clean, "patient", assay.use = "RNA")
# normal_clean <- normal_clean %>% 
#   RunUMAP(reduction = "harmony", dims = 1:20) %>% 
#   FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
#   FindClusters(resolution = 0.5) %>% 
#   identity()
# 
# normal_clean <- FindClusters(normal_clean, resolution = 1)
# 
# pred <- SingleR(normal_clean[["RNA"]]@data, refs, labels = list(refs[[1]]$label.main,refs[[2]]$label.main))
# normal_clean[["Cell.Identity"]] <- pred[colnames(normal_clean), "labels"]
# ids <- unique(normal_clean$Cell.Identity)
# ids <- gsub("_", " ", fixed = T, ids)
# ids <- gsub("-", " ", fixed = T, ids)
# ids <- gsub("s$", "", ids)
# ids <- setNames(ids, unique(normal_clean$Cell.Identity))
# normal_clean[["Cell.Identity"]] <- ids[normal_clean$Cell.Identity]
# DimPlot(normal_clean, label = T)
# meta <- normal_clean@meta.data
# max_mat <- meta %>% group_by(seurat_clusters, Cell.Identity) %>% summarise(n = n())
# max_mat <- max_mat %>% group_by(seurat_clusters) %>% mutate(max_id = max(n)) %>% mutate(max_cell = Cell.Identity[n == max_id])
# max_mat <- max_mat[, c("seurat_clusters", "max_cell")] %>% distinct()
# max_id <- structure(max_mat$max_cell, names = max_mat$seurat_clusters)
# normal_clean$max_id <- max_id[normal_clean$seurat_clusters]
# DimPlot(normal_clean, group.by = "max_id") 
# 
# sites<- normal_clean$sample
# sites[grep("brain", sites)] <- "brain"
# sites[grep("dura", sites)] <- "dura"
# sites[grep("prostate", sites)] <- "prostate"
# sites[grep("liver", sites)] <- "liver"
# sites[grep("fat", sites)] <- "fat"
# sites[grep("lymph|hilar", sites)] <- "LN"
# sites[grep("rib|vertebra", sites)] <- "bone"
# sites[grep("lung", sites)] <- "lung"
# sites[grep("abdomen", sites)] <- "abdomen"
# sites[grep("bladder", sites)] <- "bladder"
# normal_clean$site <- sites
# DimPlot(normal_clean, group.by = "max_id") +DimPlot(normal_clean, group.by = "site")
# 
# anno <- as.data.frame(as.character(normal_clean$seurat_clusters))
# rownames(anno) <- colnames(normal_clean)
# colnames(anno) <- NULL
# anno[anno[, 1] %in% as.character(c(18, 1, 2, 16, 19, 10, 13, 0, 3, 8, 15, 7)), 1] <- "ref"
# saveRDS(anno, "normal_inferncv_clean_annov1.Rds")
# saveRDS(normal_clean[["RNA"]]@counts, "normal_inferncv_clean_countsv1.Rds")
# saveRDS(normal_clean, "normal_infercnv_cleanv1_srt.Rds")
# 
# 
# dend <- read.dendrogram("~/CASCADEpaper/paper/Fig2/normal_cells/infercnv/unclustered/infercnv.preliminary.observations_dendrogram.txt")
# dend <- color_branches(dend, h = 25)
# plot(dend)
# tree <- cutree(dend, h = 25)
# rm_cells <- names(tree)[tree == 3]
# normal_clean2 <- subset(normal_clean, cells = colnames(normal_clean)[!colnames(normal_clean)%in% rm_cells])
# 
# normal_clean2[["ruv3"]] <- normal_clean2[["RNA"]]
# DefaultAssay(normal_clean2) <- "RNA"
# 
# normal_clean2 <- NormalizeData(normal_clean2)
# normal_clean2 <- FindVariableFeatures(normal_clean2, selection.method = "vst", nfeatures = 2000)
# normal_clean2 <- ScaleData(normal_clean2, features = rownames(normal_clean2))
# normal_clean2 <- RunPCA(normal_clean2, features = VariableFeatures(object = normal_clean2))
# normal_clean2 <- RunHarmony(normal_clean2, "patient", assay.use = "RNA")
# normal_clean2 <- normal_clean2 %>% 
#   RunUMAP(reduction = "harmony", dims = 1:20) %>% 
#   FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
#   FindClusters(resolution = 1) %>% 
#   identity()
# 
# pred <- SingleR(normal_clean2[["RNA"]]@data, refs, labels = list(refs[[1]]$label.main,refs[[2]]$label.main))
# normal_clean2[["Cell.Identity"]] <- pred[colnames(normal_clean2), "labels"]
# ids <- unique(normal_clean2$Cell.Identity)
# ids <- gsub("_", " ", fixed = T, ids)
# ids <- gsub("-", " ", fixed = T, ids)
# ids <- gsub("s$", "", ids)
# ids <- setNames(ids, unique(normal_clean2$Cell.Identity))
# normal_clean2[["Cell.Identity"]] <- ids[normal_clean2$Cell.Identity]
# DimPlot(normal_clean2, label = T)
# meta <- normal_clean2@meta.data
# max_mat <- meta %>% group_by(seurat_clusters, Cell.Identity) %>% summarise(n = n())
# max_mat <- max_mat %>% group_by(seurat_clusters) %>% mutate(max_id = max(n)) %>% mutate(max_cell = Cell.Identity[n == max_id])
# max_mat <- max_mat[, c("seurat_clusters", "max_cell")] %>% distinct()
# max_id <- structure(max_mat$max_cell, names = max_mat$seurat_clusters)
# normal_clean2$max_id <- max_id[normal_clean2$seurat_clusters]
# DimPlot(normal_clean2, group.by = "max_id") 
# DimPlot(normal_clean2, group.by = "Cell.Identity") 
# 
# write.csv(normal_clean2[["RNA"]]@counts, "normal_clean2.csv")
# saveRDS(normal_clean2, "normal_infercnv_cleanv2_srt.Rds")


# previous was done through sbatch next remove tumor and rest poor quality cells----------------------------------------

normal_pre_ruv <- 
normal_clean <- subset(normal_pre_ruv, cells = colnames(normal_clean3)) #normal_clean 3 is just cells remove tumor cell annotation in the integrated data, no need to record
DefaultAssay(normal_clean) <- "ruv3"

normal_clean <- FindVariableFeatures(normal_clean, selection.method = "vst", nfeatures = 2000)
normal_clean <- RunPCA(normal_clean, features = VariableFeatures(object = normal_clean))
normal_clean <- RunHarmony(normal_clean, "patient", assay.use = "ruv3")
normal_clean <- normal_clean %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 1) %>% 
  identity()

DimPlot(normal_clean, label = T)
DimPlot(normal_clean, group.by = "patient")
pred <- SingleR(normal_clean[["RNA"]]@data, refs, labels = list(refs[[1]]$label.main,refs[[2]]$label.main))
normal_clean[["Cell.Identity"]] <- pred[colnames(normal_clean), "labels"]
ids <- unique(normal_clean$Cell.Identity)
ids <- gsub("_", " ", fixed = T, ids)
ids <- gsub("-", " ", fixed = T, ids)
ids <- gsub("s$", "", ids)
ids <- setNames(ids, unique(normal_clean$Cell.Identity))
normal_clean[["Cell.Identity"]] <- ids[normal_clean$Cell.Identity]
meta <- normal_clean@meta.data
max_mat <- meta %>% group_by(seurat_clusters, Cell.Identity) %>% summarise(n = n())
max_mat <- max_mat %>% group_by(seurat_clusters) %>% mutate(max_id = max(n)) %>% mutate(max_cell = Cell.Identity[n == max_id])
max_mat <- max_mat[, c("seurat_clusters", "max_cell")] %>% distinct()
max_id <- structure(max_mat$max_cell, names = max_mat$seurat_clusters)
normal_clean$max_id <- max_id[normal_clean$seurat_clusters]
DimPlot(normal_clean, group.by = "max_id") 
sites<- normal_clean$sample
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
normal_clean$site <- sites
DimPlot(normal_clean, group.by = "site") + DimPlot(normal_clean, group.by = "max_id") 

normal_clean <- subset(normal_clean, cells = WhichCells(normal_pre_ruv, idents = 27, invert = T)) #cluster 27 is mergeing to CA76 tumor

f <- DimPlot(normal_clean)
normal_cells <- CellSelector(f) # exclude all the cells clustered with tumor
normal_clean <- subset(normal_clean, cells = intersect(colnames(normal_clean), normal_cells))

normal_clean <- FindVariableFeatures(normal_clean, selection.method = "vst", nfeatures = 2000)
normal_clean <- RunPCA(normal_clean, features = VariableFeatures(object = normal_clean))
normal_clean <- RunHarmony(normal_clean, "patient", assay.use = "ruv3")
normal_clean <- normal_clean %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 1) %>% 
  identity()

DimPlot(normal_clean, label = T)
DimPlot(normal_clean, group.by = "patient")
pred <- SingleR(normal_clean[["RNA"]]@data, refs, labels = list(refs[[1]]$label.main,refs[[2]]$label.main))
normal_clean[["Cell.Identity"]] <- pred[colnames(normal_clean), "labels"]
ids <- unique(normal_clean$Cell.Identity)
ids <- gsub("_", " ", fixed = T, ids)
ids <- gsub("-", " ", fixed = T, ids)
ids <- gsub("s$", "", ids)
ids <- setNames(ids, unique(normal_clean$Cell.Identity))
normal_clean[["Cell.Identity"]] <- ids[normal_clean$Cell.Identity]
meta <- normal_clean@meta.data
max_mat <- meta %>% group_by(seurat_clusters, Cell.Identity) %>% summarise(n = n())
max_mat <- max_mat %>% group_by(seurat_clusters) %>% mutate(max_id = max(n)) %>% mutate(max_cell = Cell.Identity[n == max_id])
max_mat <- max_mat[, c("seurat_clusters", "max_cell")] %>% distinct()
max_id <- structure(max_mat$max_cell, names = max_mat$seurat_clusters)
normal_clean$max_id <- max_id[normal_clean$seurat_clusters]
DimPlot(normal_clean, group.by = "max_id") 
sites<- normal_clean$sample
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
normal_clean$site <- sites
DimPlot(normal_clean, group.by = "site") + DimPlot(normal_clean, group.by = "max_id") 
FeaturePlot(normal_clean, c("nFeature_RNA"))+
FeaturePlot(normal_clean, c("S.Score", "G2M.Score"))

# FeaturePlot(normal_clean, c("IL7R", "CCR7", "S100A4", "CD8A", "CD14")) # Naive T cells
# FeaturePlot(normal_clean, c("CD19", "CD25", "CD30", "CD20", "CD27"))
DimPlot(normal_clean, cells.highlight = colnames(normal_clean)[normal_clean$nFeature_RNA<200|normal_clean$nFeature_RNA>7000])
rm_cells <- colnames(normal_clean)[normal_clean$nFeature_RNA<200|normal_clean$nFeature_RNA>7000]


# normal_clean2 -----------------------------------------------------------

normal_clean2 <- subset(normal_clean, cells = WhichCells(normal_clean, cells = rm_cells, invert = T))

normal_clean2 <- FindVariableFeatures(normal_clean2, selection.method = "vst", nfeatures = 2000)
normal_clean2 <- ScaleData(normal_clean2, features = rownames(normal_clean2))
normal_clean2 <- RunPCA(normal_clean2, features = VariableFeatures(object = normal_clean2))
normal_clean2 <- RunHarmony(normal_clean2, "patient", assay.use = "ruv3")
normal_clean2 <- normal_clean2 %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 1) %>% 
  identity()

pred <- SingleR(normal_clean2[["ruv3"]]@data, refs, labels = list(refs[[1]]$label.main,refs[[2]]$label.main))
normal_clean2[["Cell.Identity"]] <- pred[colnames(normal_clean2), "labels"]
ids <- unique(normal_clean2$Cell.Identity)
ids <- gsub("_", " ", fixed = T, ids)
ids <- gsub("-", " ", fixed = T, ids)
ids <- gsub("s$", "", ids)
ids <- setNames(ids, unique(normal_clean2$Cell.Identity))
normal_clean2[["Cell.Identity"]] <- ids[normal_clean2$Cell.Identity]

meta <- normal_clean2@meta.data
max_mat <- meta %>% group_by(seurat_clusters, Cell.Identity) %>% summarise(n = n())
max_mat <- max_mat %>% group_by(seurat_clusters) %>% mutate(max_id = max(n)) %>% mutate(max_cell = Cell.Identity[n == max_id])
max_mat <- max_mat[, c("seurat_clusters", "max_cell")] %>% distinct()
max_id <- structure(max_mat$max_cell, names = max_mat$seurat_clusters)
normal_clean2$max_id <- max_id[normal_clean2$seurat_clusters]
DimPlot(normal_clean2, group.by = "max_id") + DimPlot(normal_clean2, label = T)

sites<- normal_clean2$sample
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
normal_clean2$site <- sites

saveRDS(normal_clean2, "normal_finalv1.Rds")
markers <- FindAllMarkers(normal_clean2)
saveRDS(markers, "normal_finalv1_marker.Rds")


sub_marker <- markers %>% filter(avg_log2FC > 0 & p_val_adj < 0.05 )#& abs(pct.1-pct.2)>0.5




# check with Anna's previous version --------------------------------------
load("/trigos_team/CASCADE/Analysis/Normal/srt_normal_RNA_harmony.Rds")
anna_anno <- structure(srt_normal_RNA_harmony$Cluster.Identity, names = colnames(srt_normal_RNA_harmony))
normal_clean2$anna_anno <- anna_anno[colnames(srt_normal_RNA_harmony)]
DimPlot(normal_clean2, group.by = "max_id", label = T) + 
  DimPlot(normal_clean2, group.by = "anna_anno", label = T)
#all the cells are in Anna's previous harmony obj, most epithelial and neuron mixing is removed.



# looking for clean normal cell marker ------------------------------------


cell_type_marker <- fread("rna_single_cell_type.tsv")
marker_ref <- cell_type_marker %>% filter(`Gene name` %in% sub_marker$gene)

top_marker <- marker_ref %>% group_by(`Cell type`) %>% filter(nTPM %in% head(sort(nTPM, decreasing = T)))
test <- lapply(split(top_marker, top_marker$`Cell type`), function(x){
 temp <- apply(x, 1, function(y){
   if(sum(sub_marker$gene %in% y[2])>0){
     
   return(paste(sub_marker$cluster[sub_marker$gene %in% y[2]], collapse = ","))
   }
 })
 if( class(temp) == "list"){
   if(sum(sapply(temp, is.null))>0){
 temp[sapply(temp, is.null)] <- "no"
  }
 }
 
 x$cluster <- unlist(temp)
 return(x)
})
test <- rbindlist(test, fill = T)
View(test)

pdf("normal_marker.pdf", width = 8, height = 8)
for ( i in names(test)){
  f <- FeaturePlot(normal_clean2, features = test[[i]]$`Gene name`, cols = rev(brewer.rdylbu(5))) + plot_annotation(title  = i)
  print(f)
}
dev.off()


library(readxl)
markers <- read_excel("encyclopedia_table.xlsx", sheet = 1)
genes <- strsplit(markers$`Top 10 important genes from the CellTypist model`, split = ",")
genes <- lapply(genes, function(x){
  temp <- gsub("^ ", "", x)
  temp <- temp[temp %in% rownames(normal_clean2)]
})
names(genes) <- markers$`Low-hierarchy cell types`
library(patchwork)
genes <- genes[sapply(genes, length)>0]

pdf("normal_marker_encyclopidia.pdf", width = 8, height = 8)
for (i in names(genes)){
  f <- FeaturePlot(normal_clean2, genes[[i]]) + plot_annotation(title = i)
  print(f)  
}
dev.off()


c("TIGIT","NECTIN2","TGFB1")

FeaturePlot(normal_clean2, c("CD38", "CD49", "CD34", "CD10"))
DimPlot(normal_clean2, group.by = "max_id", label = T)+DimPlot(normal_clean2, label = T)
Erythroid_marker <- c("HBA2", "HBA1","HBB")
Epithelial_cell <- c("WFDC2","EPCAM", "SLPI")
Fibroblast_feature <- c("DCN","LUM", "CFD", "VCAN", "ALDH1A1")
Endothelial_feature <- c("RAMP2","VWF", "ST6GALNAC3", "NOTCH1")

Epithelial_feature <- c("EPCAM","PCA3","KRT18","KRT8", "WFDC2")

Macrophage_feature <- c("C1QA","C1QC","LYZ", "CD86", "CD68", "CD14", "FCGR3A","FCGR2B")
focused_B_cell_feature <- c("MS4A1","BANK1", "TNFRSF13C", "CD74", "IGKC", "IGHA1")
plasma_feature <- c("CD38","MZB1","XBP1")
T_cell_feature <- c("CD3G", 'CD4', "CD8A", "CCL5", "FYN",
                    "CD247", "IL7R") # naive T cell

Hepatocytes_feature <- c("APOC3", "APOA2", "ALB", "SAA1", "ORM1", "HP")
Hepatocytes_feature <- c(Hepatocytes_feature, "TF", "CYP3A4", "CYP2E1", "ASS1", "APOE")

Neuron_feature <- c("GPM6A", "GPC5", "NRXN1", "PCDH9", #astrocyte
                    "RBFOX1", "CADM2", "CSMD1", "KCNIP4", #exitatory neu
                    "NRXN3", "RBFOX1", "CNTNAP2", "CSMD1", "ROBO2") # inhibitory neu
smooth_muscle <- c("TAGLN")


all_marker <- c(Erythroid_marker, Epithelial_cell, Fibroblast_feature, Endothelial_feature, 
                Epithelial_feature, Macrophage_feature, focused_B_cell_feature, 
                plasma_feature, T_cell_feature, Hepatocytes_feature, Neuron_feature, smooth_muscle)

DoHeatmap(normal_clean2, features = all_marker)
DotPlot(normal_clean2, features = unique(all_marker))+ theme(axis.text.x = element_text(angle = 90))





# normal_clean 3: filter further based on marker and sample information -------------------
#### cluster 20 is either lymphoid prenitor, or confounded by somehow cycling states, showing B cell and macrophage markers but all not strong, better remove 
rm_cells <- WhichCells(normal_clean2, idents = "20")
#### cluster 19 is annotated as neuron but don't have strong neuron markers like others and they are not from brain tissue
rm_cells <- c(rm_cells, WhichCells(normal_clean2, idents = "19"))
#### by checking hepatocytes markers, and sample origin, remove non-liver cells in cluster 13
rm_cells <- c(rm_cells, colnames(normal_clean2)[normal_clean2$site != "liver" & normal_clean2$seurat_clusters == '13'])
#### the chondrocytes look very suspicious as well, showing weak fibroblast markers, but we dont have chondrocytes marker so let's temprary leave them

DimPlot(normal_clean2, cells.highlight = colnames(normal_clean2)[!colnames(normal_clean2)%in%rm_cells])

normal_clean3 <- subset(normal_clean2, cells = colnames(normal_clean2)[!colnames(normal_clean2)%in%rm_cells])

normal_clean3 <- FindVariableFeatures(normal_clean3, selection.method = "vst", nfeatures = 2000)
normal_clean3 <- ScaleData(normal_clean3, features = rownames(normal_clean3))
normal_clean3 <- RunPCA(normal_clean3, features = VariableFeatures(object = normal_clean3))
normal_clean3 <- RunHarmony(normal_clean3, "patient", assay.use = "ruv3")
normal_clean3 <- normal_clean3 %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 1) %>% 
  identity()

pred <- SingleR(normal_clean3[["ruv3"]]@data, refs, labels = list(refs[[1]]$label.main,refs[[2]]$label.main))
normal_clean3[["Cell.Identity"]] <- pred[colnames(normal_clean3), "labels"]
ids <- unique(normal_clean3$Cell.Identity)
ids <- gsub("_", " ", fixed = T, ids)
ids <- gsub("-", " ", fixed = T, ids)
ids <- gsub("s$", "", ids)
ids <- setNames(ids, unique(normal_clean3$Cell.Identity))
normal_clean3[["Cell.Identity"]] <- ids[normal_clean3$Cell.Identity]

meta <- normal_clean3@meta.data
max_mat <- meta %>% group_by(seurat_clusters, Cell.Identity) %>% summarise(n = n())
max_mat <- max_mat %>% group_by(seurat_clusters) %>% mutate(max_id = max(n)) %>% mutate(max_cell = Cell.Identity[n == max_id])
max_mat <- max_mat[, c("seurat_clusters", "max_cell")] %>% distinct()
max_id <- structure(max_mat$max_cell, names = max_mat$seurat_clusters)
normal_clean3$max_id <- max_id[normal_clean3$seurat_clusters]
DimPlot(normal_clean3, group.by = "max_id") + DimPlot(normal_clean3, label = T)

sites<- normal_clean3$sample
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
normal_clean3$site <- sites
DotPlot(normal_clean3, features = unique(all_marker))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))


FeaturePlot(normal_clean3, cell_type_marker %>% filter(`Cell type` =="Adipocytes") %>% arrange(nTPM) %>% pull(`Gene name`) %>% tail())
FeaturePlot(normal_clean3, c("LEPTIN", "LPL", "FATP1", "NRG4", "P2RX5", "PAT2"))
normal_clean3$cell_id <- NA

anno2 <- normal_clean2$max_id
normal_clean3$anno2 <- anno2[colnames(normal_clean3)]
unmatched_cells<- names(which(normal_clean3$max_id != normal_clean3$anno2))
head(unmatched_cells)
DimPlot(normal_clean3, cells.highlight =unmatched_cells, sizes.highlight = 0.5)
DimPlot(normal_clean3, group.by = "anno2", label = T)+DimPlot(normal_clean3, group.by = "max_id", label = T) 

normal_clean3$cell_id[normal_clean3$seurat_clusters %in% c("0", "3", "13", "18", "19", "20") ] <- "Fibroblast"
normal_clean3$cell_id[normal_clean3$seurat_clusters %in% c("1", "5", "21") ] <- "Macrophage"
normal_clean3$cell_id[normal_clean3$seurat_clusters %in% c("2", "4") ] <- "T cell"
normal_clean3$cell_id[normal_clean3$seurat_clusters %in% c("24") ] <- "Naive cell" # showing marker of T cell, B cell and Erythroblast
normal_clean3$cell_id[normal_clean3$seurat_clusters %in% c("6", "10") ] <- "B cell"
normal_clean3$cell_id[normal_clean3$seurat_clusters %in% c("7", "8", "17") ] <- "Endothelial cell"

normal_clean3$cell_id[normal_clean3$seurat_clusters %in% c("9") ] <- "Adipocyte"
normal_clean3$cell_id[normal_clean3$seurat_clusters %in% c("11", "12") ] <- "Epithelial cell"
normal_clean3$cell_id[normal_clean3$seurat_clusters %in% c("14") ] <- "Hepatocyte"
normal_clean3$cell_id[normal_clean3$seurat_clusters %in% c("15") ] <- "Plasma cell"
normal_clean3$cell_id[normal_clean3$seurat_clusters %in% c("16") ] <- "Erythroblast"
normal_clean3$cell_id[normal_clean3$seurat_clusters %in% c("22", "23", "25") ] <- "Neuron"

DimPlot(normal_clean3, group.by = "cell_id", label = T) # based on marker expression
saveRDS(normal_clean3, "normal_clean3.Rds")

# current markers work well
# worth to remove cells cluster differently between normal_clean2 and normal_clean3 with very differnt cell annotation.

rm_cells <- unmatched_cells[!unmatched_cells %in% WhichCells(normal_clean3, idents = c("0", "3", "13", "9"))]


f <- DimPlot(normal_clean3, group.by = "cell_id", label = T)
suspicious <- CellSelector(f)
rm_cells <- unique(c(rm_cells, intersect(suspicious, WhichCells(normal_clean3, idents = "5"))))

normal_clean4 <- subset(normal_clean3, cells = colnames(normal_clean3)[!colnames(normal_clean3)%in%rm_cells])

normal_clean4 <- FindVariableFeatures(normal_clean4, selection.method = "vst", nfeatures = 2000)
normal_clean4 <- ScaleData(normal_clean4, features = rownames(normal_clean4))
normal_clean4 <- RunPCA(normal_clean4, features = VariableFeatures(object = normal_clean4))
normal_clean4 <- RunHarmony(normal_clean4, "patient", assay.use = "ruv3")
normal_clean4 <- normal_clean4 %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 1) %>% 
  identity()

pred <- SingleR(normal_clean4[["ruv3"]]@data, refs, labels = list(refs[[1]]$label.main,refs[[2]]$label.main))
normal_clean4[["Cell.Identity"]] <- pred[colnames(normal_clean4), "labels"]
ids <- unique(normal_clean4$Cell.Identity)
ids <- gsub("_", " ", fixed = T, ids)
ids <- gsub("-", " ", fixed = T, ids)
ids <- gsub("s$", "", ids)
ids <- setNames(ids, unique(normal_clean4$Cell.Identity))
normal_clean4[["Cell.Identity"]] <- ids[normal_clean4$Cell.Identity]

meta <- normal_clean4@meta.data
max_mat <- meta %>% group_by(seurat_clusters, Cell.Identity) %>% summarise(n = n())
max_mat <- max_mat %>% group_by(seurat_clusters) %>% mutate(max_id = max(n)) %>% mutate(max_cell = Cell.Identity[n == max_id])
max_mat <- max_mat[, c("seurat_clusters", "max_cell")] %>% distinct()
max_id <- structure(max_mat$max_cell, names = max_mat$seurat_clusters)
normal_clean4$max_id <- max_id[normal_clean4$seurat_clusters]
DimPlot(normal_clean4, group.by = "max_id") + DimPlot(normal_clean4, label = T)

sites<- normal_clean4$sample
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
normal_clean4$site <- sites
DotPlot(normal_clean4, features = unique(all_marker))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))


DimPlot(normal_clean4, group.by = "max_id", label = T) + DimPlot(normal_clean4, group.by = "cell_id", label = T)

DimPlot(normal_clean4, group.by = "max_id", label = T) + DimPlot(normal_clean4, group.by = "site", label = T)


normal_clean4$anno3 <- normal_clean4$cell_id

normal_clean4$cell_id <- NA
normal_clean4$cell_id[normal_clean4$seurat_clusters %in% c("1", "12") ] <- "Fibroblast"
normal_clean4$cell_id[normal_clean4$seurat_clusters %in% c("20", "16", "2") ] <- "Chondrocyte"
normal_clean4$cell_id[normal_clean4$seurat_clusters %in% c("10") ] <- "Adipocyte"

normal_clean4$cell_id[normal_clean4$seurat_clusters %in% c("0", "5", "21") ] <- "Macrophage"
normal_clean4$cell_id[normal_clean4$seurat_clusters %in% c("3", "4") ] <- "T cell"
normal_clean4$cell_id[normal_clean4$seurat_clusters %in% c("22") ] <- "Naive cell" # showing marker of T cell, B cell and Erythroblast
normal_clean4$cell_id[normal_clean4$seurat_clusters %in% c("7", "11") ] <- "B cell"
normal_clean4$cell_id[normal_clean4$seurat_clusters %in% c("13") ] <- "Plasma cell"
normal_clean4$cell_id[normal_clean4$seurat_clusters %in% c("15") ] <- "Erythroblast"

normal_clean4$cell_id[normal_clean4$seurat_clusters %in% c( "8", "9", "17") ] <- "Endothelial cell"
normal_clean4$cell_id[normal_clean4$seurat_clusters %in% c("6") ] <- "Epithelial cell"
normal_clean4$cell_id[normal_clean4$seurat_clusters %in% c("14") ] <- "Hepatocyte"
normal_clean4$cell_id[normal_clean4$seurat_clusters %in% c("19") ] <- "Neuron"
DimPlot(normal_clean4, group.by = "cell_id", label = T) + DimPlot(normal_clean4, group.by = "site", label = T)
normal_clean4$cell_id[is.na(normal_clean4$cell_id)] <- "unclear"
Idents(normal_clean4) <- normal_clean4$cell_id
DotPlot(normal_clean4, features = unique(all_marker))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
saveRDS(normal_clean4, "normal_clean4.Rds")

# normal_clean5: last removal cells ---------------------------------------
unmatched_cells <- colnames(normal_clean4)[normal_clean4$cell_id != normal_clean4$anno3]
unmatched_cells <- unmatched_cells[!unmatched_cells %in% colnames(normal_clean4)[normal_clean4$cell_id == "Chondrocyte"]]
DimPlot(normal_clean4, cells.highlight = unmatched_cells, sizes.highlight = 0.2)
f <- DimPlot(normal_clean4, group.by = "cell_id", label = T)
temp <- CellSelector(f)

rm_cells <- c(colnames(normal_clean4)[normal_clean4$cell_id == "unclear"], unmatched_cells)

normal_clean5 <- subset(normal_clean4, cells = colnames(normal_clean4)[!colnames(normal_clean4)%in%rm_cells])

normal_clean5 <- FindVariableFeatures(normal_clean5, selection.method = "vst", nfeatures = 2000)
normal_clean5 <- ScaleData(normal_clean5, features = rownames(normal_clean5))
normal_clean5 <- RunPCA(normal_clean5, features = VariableFeatures(object = normal_clean5))
normal_clean5 <- RunHarmony(normal_clean5, "patient", assay.use = "ruv3")
normal_clean5 <- normal_clean5 %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 1) %>% 
  identity()

pred <- SingleR(normal_clean5[["ruv3"]]@data, refs, labels = list(refs[[1]]$label.main,refs[[2]]$label.main))
normal_clean5[["Cell.Identity"]] <- pred[colnames(normal_clean5), "labels"]
ids <- unique(normal_clean5$Cell.Identity)
ids <- gsub("_", " ", fixed = T, ids)
ids <- gsub("-", " ", fixed = T, ids)
ids <- gsub("s$", "", ids)
ids <- setNames(ids, unique(normal_clean5$Cell.Identity))
normal_clean5[["Cell.Identity"]] <- ids[normal_clean5$Cell.Identity]

meta <- normal_clean5@meta.data
max_mat <- meta %>% group_by(seurat_clusters, Cell.Identity) %>% summarise(n = n())
max_mat <- max_mat %>% group_by(seurat_clusters) %>% mutate(max_id = max(n)) %>% mutate(max_cell = Cell.Identity[n == max_id])
max_mat <- max_mat[, c("seurat_clusters", "max_cell")] %>% distinct()
max_id <- structure(max_mat$max_cell, names = max_mat$seurat_clusters)
normal_clean5$max_id <- max_id[normal_clean5$seurat_clusters]
DimPlot(normal_clean5, group.by = "max_id") + DimPlot(normal_clean5, label = T)

sites<- normal_clean5$sample
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
normal_clean5$site <- sites
DotPlot(normal_clean5, features = unique(all_marker))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
saveRDS(normal_clean5, "normal_clean5.Rds")

normal_clean5$anno4 <- normal_clean5$cell_id
normal_clean5$cell_id <- NA
normal_clean5$cell_id[normal_clean5$seurat_clusters %in% c("0", "1") ] <- "Fibroblast"
normal_clean5$cell_id[normal_clean5$seurat_clusters %in% c("13") ] <- "Smooth muscle"
normal_clean5$cell_id[normal_clean5$seurat_clusters %in% c("21", "19") ] <- "Chondrocyte"
normal_clean5$cell_id[normal_clean5$seurat_clusters %in% c("11") ] <- "Adipocyte"

normal_clean5$cell_id[normal_clean5$seurat_clusters %in% c("4", "6", "10", "14") ] <- "Macrophage"
normal_clean5$cell_id[normal_clean5$seurat_clusters %in% c("2", "3") ] <- "T cell"
normal_clean5$cell_id[normal_clean5$seurat_clusters %in% c("23") ] <- "Naive cell" # showing marker of T cell, B cell and Erythroblast
normal_clean5$cell_id[normal_clean5$seurat_clusters %in% c("5", "12") ] <- "B cell"
normal_clean5$cell_id[normal_clean5$seurat_clusters %in% c("15") ] <- "Plasma cell"
normal_clean5$cell_id[normal_clean5$seurat_clusters %in% c("17") ] <- "Erythroblast"

normal_clean5$cell_id[normal_clean5$seurat_clusters %in% c( "7", "9", "20") ] <- "Endothelial cell"
normal_clean5$cell_id[normal_clean5$seurat_clusters %in% c("8","18") ] <- "Epithelial cell"
normal_clean5$cell_id[normal_clean5$seurat_clusters %in% c("16") ] <- "Hepatocyte"
normal_clean5$cell_id[normal_clean5$seurat_clusters %in% c("22") ] <- "Neuron"

Idents(normal_clean5) <- normal_clean5$cell_id
DimPlot(normal_clean5, label = T)
#lack effective chondrocytes and adipocyte marker
library(MAST)
markers <- FindAllMarkers(normal_clean5, test.use = "MAST", only.pos = T, return.thresh = 0.05,logfc.threshold = 0.5)
FeaturePlot(normal_clean5, c("THSD4", "NTN4", "EYA2", "CXADR", "AC013652.1", "EHF", "CPA6", "AC104041.1", "FAM160A1", "ELF3", "FAAH2", "HIST1H2AC", "TMC5")

)



Epithelial_cell <- c("WFDC2","EPCAM", "SLPI", "EHF", "TMC5")#"CXADR", developing neuron system

Hepatocytes_feature <- c("APOC3", "APOA2", "ALB", "SAA1", "ORM1", "HP")
Hepatocytes_feature <- c(Hepatocytes_feature, "TF", "CYP3A4", "CYP2E1")#, "ASS1", "APOE")

Neuron_feature <- c( "GPC5", "NRXN1", "PCDH9", #astrocyte #G6PAM endothelial marker
                    "RBFOX1", "CADM2", "CSMD1", "KCNIP4")#, #exitatory neu
                    #"NRXN3", "CNTNAP2", "ROBO2") # inhibitory neu

Fibroblast_feature <- c("DCN","LUM", "CFD", "VCAN", "ALDH1A1", "PDGFRA", "COL3A1")
smooth_muscle <- c("TAGLN", "MYH11","ACTA2", "KCNJ8", "CNN1")
Adipocyte_feature <- c("ABCC9", "RGS5", "GJC1", "ADCY3")# ADCY3, neuronal primary cilia and obesity
Chondrocyte_feature <- c("COL1A1", "KAZN")


Endothelial_feature <- c("RAMP2","VWF", "ST6GALNAC3", "NOTCH1", "FLT1")


Erythroid_marker <- c("HBA2", "HBA1","HBB", "ANK1")
Macrophage_feature <- c("C1QA","C1QC","LYZ", "CD86", "CD68", "CD14", "FCGR3A","FCGR2B")
focused_B_cell_feature <- c("MS4A1","BANK1", "TNFRSF13C", "CD74", "IGKC", "IGHA1")
plasma_feature <- c("CD38","MZB1","XBP1","PIM2")
T_cell_feature <- c("CD3G", 'CD4', "CD8A", "CCL5",# "FYN",
                    "CD247", "IL7R") # naive T cell


all_marker <- list(Hepatocytes_feature, Neuron_feature,Epithelial_cell, Fibroblast_feature, smooth_muscle, Adipocyte_feature,Chondrocyte_feature, Endothelial_feature, 
                   Macrophage_feature, focused_B_cell_feature, 
                   plasma_feature, T_cell_feature, Erythroid_marker)
names(all_marker) <- c("Hepatocyte", "Neuron", "Epithelial", "Fibroblast", "Smooth\nmuscle", "Adipocyte", "Chondrocyte", "Edothelial", "Macrophage", "B cell", "Plasma cell", "T cell", "Erythroblst")
# all_marker <- c(Hepatocytes_feature, Neuron_feature,Epithelial_cell, Fibroblast_feature, smooth_muscle, Adipocyte_feature,Chondrocyte_feature, Endothelial_feature, 
#                  Macrophage_feature, focused_B_cell_feature, 
#                 plasma_feature, T_cell_feature, Erythroid_marker)
saveRDS(all_marker, "normal_markers.Rds")


DotPlot(normal_clean5, features = unique(all_marker), group.by = "cell_id")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

saveRDS(normal_clean5, "normal_final_srt.Rds")
saveRDS(markers, "normal_final_marker.Rds")

DimPlot(normal_clean5)
DefaultAssay(normal_clean5) <- "ruv3"
Idents(normal_clean5, "cell_id")
DotPlot(normal_clean5, feature= unique(all_marker))

DefaultAssay(normal_clean5) <- "RNA"
light_srt <- DietSeurat(normal_clean5, assays = "RNA")
light_sce <- as.SingleCellExperiment(light_srt)
saveRDS(light_sce, "normal_sce.Rds")


### run harmony on original RNA assay
DefaultAssay(normal_final_srt) <- "RNA"
normal_final_srt <- NormalizeData(normal_final_srt)
normal_final_srt <- ScaleData(normal_final_srt)

normal_final_srt <- FindVariableFeatures(normal_final_srt, selection.method = "vst", nfeatures = 2000)
normal_final_srt <- RunPCA(normal_final_srt, features = VariableFeatures(object = normal_final_srt))
normal_final_srt <- RunHarmony(normal_final_srt, "patient", assay.use = "RNA")
normal_final_srt <- normal_final_srt %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 1) %>% 
  identity()
DimPlot(normal_final_srt, group.by = "cell_id", label = T) +
DimPlot(normal_final_srt, group.by = "patient", label = T) + DimPlot(normal_final_srt, label = T)+NoLegend()
DotPlot(normal_final_srt, features = unique(all_marker))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

### see only ruv assay
DefaultAssay(normal_final_srt) <- "ruv3"
normal_final_srt <- FindVariableFeatures(normal_final_srt, selection.method = "vst", nfeatures = 2000)
normal_final_srt <- RunPCA(normal_final_srt, features = VariableFeatures(object = normal_final_srt))
normal_final_srt <- normal_final_srt %>% 
  RunUMAP(reduction = "pca", dims = 1:20) %>% 
  FindNeighbors(reduction = "pca", dims = 1:20) %>% 
  FindClusters(resolution = 1) %>% 
  identity()
DimPlot(normal_final_srt, group.by = "cell_id", label = T) +
  DimPlot(normal_final_srt, group.by = "patient", label = T) + DimPlot(normal_final_srt, label = T)+NoLegend()
DotPlot(normal_final_srt, features = unique(all_marker))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))


