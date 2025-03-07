library(Seurat)
library(RColorBrewer)
library(data.table)
library(dplyr)
library(tidyverse)
library(patchwork)
library(mmand)
library(pals)
clean_module <- readRDS("~/CASCADEpaper/paper/Fig5_archetype/clean_module.Rds")
paths <- system("realpath /trigos_team/CASCADE/Analysis/230113_seurat_intron/*/*/*non_tumour_removed.Rds", intern = T)
paths_split <- unlist(lapply(strsplit(paths, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
paths2 <- system("realpath /trigos_team/CASCADE/Analysis/240603_seurat_intron/CA0027/dura*/*non_tumour_removed_v4.Rds", intern = T)
paths_split2 <- unlist(lapply(strsplit(paths2, split = "/"), function(x) paste(x[6], x[7], sep = "_")))

paths <- c(paths, paths2)
paths_split <- c(paths_split, paths_split2)
paths <- setNames(paths,paths_split)
srt <- readRDS(paths[1])

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

srt <- AddModuleScore(srt, features  = clean_module)
names(clean_module)
colnames(srt@meta.data)[grep("Cluster", colnames(srt@meta.data))] <- paste(names(clean_module), "module")



auc_paths <- system("realpath /trigos_team/CASCADE/Analysis/230303_enrich_regulons/*/*/srt_AUC.Rds", intern = T)
auc_paths_split <- unlist(lapply(strsplit(auc_paths, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
auc_paths <- setNames(auc_paths, auc_paths_split)
auc_srt <- readRDS(auc_paths[1])
auc_srt[["new_embeddings"]] <- CreateDimReducObject(embeddings = srt[["umap"]]@cell.embeddings, key = "RNA_", assay = "AUC")



pdf("regulon test.pdf", width =10, height = 8)
TFs <- lapply(names(clean_module), function(n){
  x <- clean_module[[n]]
  temp <- intersect(x, gsub(".pos", "",rownames(auc_srt)))
  sapply(temp, function(gene){
    f <- FeaturePlot(auc_srt, features = paste0(gene, ".pos"), reduction = "new_embeddings") + 
      FeaturePlot(srt, features = gene) +
      FeaturePlot(srt, features = paste(n, "module")) + plot_annotation(title = n) + plot_layout(ncol = 3)
    g <- FeaturePlot(auc_srt, features = paste0(gene, ".pos"), reduction = "new_embeddings", cols = rev(brewer.rdylbu(10))) + 
      FeaturePlot(srt, features = gene, cols = rev(brewer.rdylbu(10))) +
      FeaturePlot(srt, features = paste(n, "module"), cols = rev(brewer.rdylbu(10))) + plot_annotation(title = n) + plot_layout(ncol = 3)
    print(f/g)
  })
  return(temp)
})
dev.off()




library(clusterProfiler)
library(enrichplot)
library(dorothea)
library(decoupleR)
library(ggplot2)
library(dplyr)
library(ggpubr)

marker <- list.files("~/regulon/new_regulon/archetype/2024_06/4_markers/", pattern = "_marker.Rds", full.names = T) %>% lapply(., readRDS)
DEG <- lapply(marker, lapply, function(x) x$DEG )
cepo <- lapply(marker, lapply, function(x) x$cepo )
names(DEG) <- list.files("~/regulon/new_regulon/archetype/2024_06/4_markers/", pattern = "_marker.Rds") %>% gsub("_marker.Rds", "", .)
names(cepo) <- names(DEG)
markerlist <- DEG

generate_genelist <- function(markerdf){
  temp <- setNames(markerdf$avg_log2FC, rownames(markerdf))
  temp <- sort(temp, decreasing = T)
  return(temp)
}

net <- dorothea_hs %>% filter(confidence != "D") %>% select(tf, target)

lapply(clean_module, function(x) {
  sort(table(net$tf[net$target %in% x]), decreasing = T)
})


#### fisher exact test ----------
netlist <- split(net$target, net$tf)
allgenes <- unlist(clean_module, use.names = F)
fishtest <- function(cat1, cat2, total = allgenes){
  return(fisher.test(table(total %in% cat1, total%in% cat2)))
}

test <- lapply(clean_module, function(m){
  lapply(netlist, function(n){
    if(sum(allgenes %in% n) >0 ){
      f<- fishtest(m, n, total = allgenes)
      
    }else{print("tf not related")}
  })
})

enlist <- lapply(markerlist, function(markers){
  enresult <- lapply(markers[!sapply(markers, is.null)], function(x){
    hall <- try(GSEA(generate_genelist(x), TERM2GENE = net))
  return(hall@result)
  })
  return(rbindlist(enresult, idcol = "module", fill = T))
})
saveRDS(enlist, "regulonenich.Rds")

endf <- rbindlist(regulonenich, idcol = "sample")

temp <- endf %>% filter(NES > 0) %>% select(sample,module, ID) %>% distinct() %>% group_by(module, ID) %>% summarise(n = n()) 
temp <- temp[temp$ID %in% names(which(table(temp$ID) == 1)),]
%>% arrange(desc(n)) %>% slice(1:5)%>% 
  mutate(Description = gsub("_", " ", ID))

NES <- apply(temp,1, function(x) endf %>% filter(module == x[1]&ID == x[2]) %>% summarise(mean = mean(NES)) %>% pull(mean))
temp$avg_NES <- NES
temp <- temp[order(temp$module, temp$n),]
temp$Description_lv <- factor(paste(temp$ID, temp$module), levels = paste(temp$ID, temp$module))


ggplot2::ggplot(temp, aes(n, Description_lv, fill=avg_NES)) + 
  geom_col(orientation='y') + 
  scale_fill_continuous(low="#FEE0D2", high="#DE2D26", guide=guide_colorbar(reverse=TRUE)) + 
  theme_bw() + ylab(NULL)+
  facet_grid(module~., scales = "free", space = "free") + 
  theme(text = element_text(size = 10), legend.position = "bottom") + 
  scale_y_discrete(labels = setNames(temp$Description, temp$Description_lv))
