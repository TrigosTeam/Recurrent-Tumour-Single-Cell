library(SeuratObject)
library(Seurat)
library(scran)
library(simspec)
library(dplyr)
library(patchwork)
library(ggplot2)

paths <- system("realpath ~/230113_seurat_intron/*/*/*non_tumour_removed.Rds", intern = T)
paths_split <- unlist(lapply(strsplit(paths, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
paths2 <- system("realpath ~/240603_seurat_intron/CA0027/dura*/*non_tumour_removed_v4.Rds", intern = T)
paths_split2 <- unlist(lapply(strsplit(paths2, split = "/"), function(x) paste(x[6], x[7], sep = "_")))

paths <- c(paths, paths2)
paths_split <- c(paths_split, paths_split2)

intron_tumor <-lapply(paths, readRDS)
names(intron_tumor) <- paths_split
print("read in done")

renamed_tumor <- lapply(names(intron_tumor), function(name){
  temp <- RenameCells(intron_tumor[[name]], add.cell.id = name)
  temp$sample <- name
  DefaultAssay(temp) <- "RNA"
  return(temp)
})
rm(list = "intron_tumor")
gc()
print("rename done")

test <- merge(renamed_tumor[[1]], renamed_tumor[2:length(renamed_tumor)], merge.data = F, project = "tumor_only")

print("merge done")
# test <- Merge_Seurat_List(renamed_tumor, project = "merge_all_tumor", merge.data = F)
rm(list = "renamed_tumor")
gc()
DefaultAssay(test) <- "RNA"
lib.gmean3 <- exp(mean(log(colSums(test[["RNA"]]@counts))))
target_pseudocount = 1
test <- test %>%
  NormalizeData(scale.factor = lib.gmean3 / target_pseudocount) %>%
  ScaleData(do.scale = FALSE)
set.seed(100)
poisson_fit3 <- modelGeneVarByPoisson(test[["RNA"]]@counts, dispersion = 0.1, size.factors = colSums(test[["RNA"]]@counts)/lib.gmean3)
residuals3 <- poisson_fit3$total - poisson_fit3@metadata$trend(poisson_fit3$mean)
names(residuals3) <- rownames(test)
top_genes3 <- rownames(test)[order(residuals3, decreasing = TRUE)]
test <- test %>%
  RunPCA(features = top_genes3[1:2000]) %>%
  FindVariableFeatures(nfeatures = 2000)
test <- cluster_sim_spectrum(object = test, label_tag = "sample")
test <- RunUMAP(test, reduction = "css", dims = 1:ncol(Embeddings(test, "css")))
test <- FindNeighbors(test, reduction = "css", dims = 1:ncol(Embeddings(test, "css")))
test <- FindClusters(test, resolution = 1)
test$patient <- substr(test$sample, 1, 6)
saveRDS(test, "~/Fig2_signature/sc_integration/tumor_only_RNA_simspec_srt.Rds")
# saveRDS(test, "~/integration/objects/tumor_merged_27_34_35_43_46_58_76_83_90_srt.Rds")
# test <- readRDS("~/integration/objects/tumor_merged_27_34_35_43_46_58_76_83_90_srt.Rds")
pdf( "~/Fig2_signature/sc_integration/tumor_only_RNA_simspec_UMAP.pdf", width = 20, height = 10)
f <- DimPlot(test, group.by = "sample", label = T)
print(f)

f <- DimPlot(test, group.by = "patient", label = T)
print(f)
dev.off()
