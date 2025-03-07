paths <- system("realpath /trigos_team/CASCADE/Analysis/230113_seurat_intron/*/*/*non_tumour_removed.Rds", intern = T)
paths_split <- unlist(lapply(strsplit(paths, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
paths2 <- system("realpath /trigos_team/CASCADE/Analysis/240603_seurat_intron/CA0027/dura*/*non_tumour_removed_v4.Rds", intern = T)
paths_split2 <- unlist(lapply(strsplit(paths2, split = "/"), function(x) paste(x[6], x[7], sep = "_")))

paths <- c(paths, paths2)
paths_split <- c(paths_split, paths_split2)
paths <- setNames(paths,paths_split)
meta <- readRDS("~/regulon/new_regulon/archetype/2024_06/5_integrated_final_module_meta.Rds")
srt<- readRDS(paths[1])
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

clean_module <- readRDS("~/regulon/new_regulon/archetype/2024_06/clean_module.Rds")
srt <- AddModuleScore(srt, features  = clean_module)
names(clean_module)
colnames(srt@meta.data)[grep("Cluster", colnames(srt@meta.data))] <- paste(names(clean_module), "module")


arc <- readRDS("~/regulon/new_regulon/archetype/2024_06/2_arc/2_arcCA0027_dura_base_skull_13.Rds")
PCs4arch= t(srt[["pca"]]@cell.embeddings)
meta <- readRDS("~/regulon/new_regulon/archetype/2024_06/2_meta/CA0027_dura_base_skull_13.Rds")
srt <- AddMetaData(srt, meta)

p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = srt$seurat_clusters,
                 data_alpha = 0.2, 
                 colors = c(brewer.accent(13), "red"),
                 text_size = 60, data_size = 6) 