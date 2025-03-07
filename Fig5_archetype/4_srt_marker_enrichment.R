library(patchwork)
library(dplyr)
library(Seurat)
library(scran)
library(data.table)
library(pals)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(ggnewscale)
library(dorothea)
library(dplyr)
library(tibble)
library(data.table)
library(ggstance)
library(MAST)
library(gridExtra)
library(stringr)
library(ggplot2)
library(tidyverse)
library(Cepo)
library(grid)

setwd("~/regulon/new_regulon/archetype/2024_06")
marker_cluster <- readRDS("~/regulon/new_regulon/archetype/2024_06/clean_module.Rds")

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
  target_pseudocount = 1
  set.seed(100)
  poisson_fit <- modelGeneVarByPoisson(srt@assays$RNA@counts, dispersion = 0.2, size.factors = target_pseudocount * srt$nCount_RNA / lib.median)
  residuals  <- poisson_fit$total - poisson_fit@metadata$trend(poisson_fit$mean)
  names(residuals ) <- rownames(srt )
  top_genes  <- rownames(srt)[order(residuals , decreasing = TRUE)]
  
  srt <- srt %>% 
    RunPCA(features = top_genes[1:2000]) %>% 
    FindNeighbors(dims = 1:20) %>% 
    FindClusters(resolution = 1) %>% 
    RunUMAP(dims = 1:20)
  
  # marker_cluster2 <- lapply(marker_cluster, function(x) x[x %in% rownames(srt)])
  # 
  # print(paste("module", which(sapply(marker_cluster2, length)==0), "is empty"))
  # 
  # marker_cluster2 <- marker_cluster2[sapply(marker_cluster2, length)>0]
  
  srt <- AddModuleScore(srt, marker_cluster, name = "Module")
  saveRDS(srt@meta.data[, grep("Module", colnames(srt@meta.data))], paste0("~/regulon/new_regulon/archetype/2024_06/4_srt_meta/", i, ".Rds"))
  
  pdf(paste0("~/regulon/new_regulon/archetype/2024_06/4_srt_plot/", i, ".pdf"), width = 12, height = 6)
  f <- FeaturePlot(srt, features = paste0("Module",1:7), cols = rev(brewer.rdylbu(10))) + plot_layout(nrow = 2, ncol = 4)
  print(f)
  f <- FeaturePlot(srt, features = paste0("Module",1:7))+ plot_layout(nrow = 2, ncol = 4)
  print(f)
  dev.off()
  
  cells <- apply(srt@meta.data[, grep("Module", colnames(srt@meta.data))], 2, function(x){
    bins <- cut(x, breaks = 10) 
    return(colnames(srt)[bins %in% levels(bins)[7:10]])
  })
  saveRDS(cells, paste0("~/regulon/new_regulon/archetype/2024_06/4_population/", i, "_cell.Rds") )
  
  cepoout <-lapply(cells, function(x){
    if (length(x) >20 & length(x) < ncol(srt)){
      srt$group <- ifelse(colnames(srt) %in% x , "high", "other")
      #srt$group[!colnames(srt)%in%unlist(cells)] <- "background"
      out <- Cepo(srt[["RNA"]]@data, srt$group, useNames = T)
      Idents(srt) <- srt$group
      marker <- try(FindMarkers(srt, ident.1 = "high", ident.2 = "other", test.use = "MAST"))
      print(length(x))
      return(list(cepo = as.data.frame(out$stats), DEG = marker))}
  })
  saveRDS(cepoout, paste0("~/regulon/new_regulon/archetype/2024_06/4_markers/", i, "_marker.Rds") )
  print(i)
}

rm(list = ls())
gc()


Hall <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
KEGG <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")%>%
  dplyr::select(gs_name, gene_symbol)
Reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME")%>%
  dplyr::select(gs_name, gene_symbol)

C5 <- msigdbr(species = "Homo sapiens", category = "C5") %>%
  dplyr::select(gs_name, gene_symbol)
C6 <- msigdbr(species = "Homo sapiens", category = "C6") %>%
  dplyr::select(gs_name, gene_symbol)
C7 <- msigdbr(species = "Homo sapiens", category = "C7") %>%
  dplyr::select(gs_name, gene_symbol)
all <- rbind(Hall, KEGG, Reactome,C5,C6)

C8 <- msigdbr(species = "Homo sapiens", category = "C8") %>%
  dplyr::select(gs_name, gene_symbol)

generate_genelist_cepo <- function(markerdf){
  temp <- setNames(markerdf$high, rownames(markerdf))
  temp <- sort(temp, decreasing = T)
  return(temp)
}


generate_genelist <- function(markerdf){
  temp <- setNames(markerdf$avg_log2FC, rownames(markerdf))
  temp <- sort(temp, decreasing = T)
  return(temp)
}


marker <- list.files("~/regulon/new_regulon/archetype/2024_06/4_markers/", pattern = "_marker.Rds", full.names = T) %>% lapply(., readRDS)
DEG <- lapply(marker, lapply, function(x) x$DEG )
cepo <- lapply(marker, lapply, function(x) x$cepo )
names(DEG) <- list.files("~/regulon/new_regulon/archetype/2024_06/4_markers/", pattern = "_marker.Rds") %>% gsub("_marker.Rds", "", .)
names(cepo) <- names(DEG)
markerlist <- DEG

enlist <- lapply(markerlist, function(markers){
  enresult <- lapply(markers[!sapply(markers, is.null)], function(x){
    hall <- try(GSEA(generate_genelist(x), TERM2GENE = Hall))
    kegg <- try(GSEA(generate_genelist(x), TERM2GENE = KEGG))
    reactome <- try(GSEA(generate_genelist(x), TERM2GENE = Reactome))
    return(list(hall = hall, kegg = kegg, reactome = reactome))
  })
  return(enresult)
})


saveRDS(enlist, "4_enlist.Rds")
cells <- list.files("~/regulon/new_regulon/archetype/2024_06/4_population", pattern = "_cell.Rds", full.names = T) %>% lapply(., readRDS)
names(cells) <- list.files("~/regulon/new_regulon/archetype/2024_06/4_population", pattern = "_cell.Rds") %>% gsub("_cell.Rds", "",.)



enrichment_barplots_new <- function(endf, module, cells){
  
  df <- endf[endf$module== module, ]
  if(nrow(df)>1){
    num <- sum(sapply(cells, function(x) length(x[[module]])>20))
    temp <- as.data.frame(head(sort(table(df$ID), decreasing = T), 5)) %>%
      `colnames<-`(c("Description", "Freq")) %>% 
      mutate(Description = gsub("_", " ", Description)) %>%
      mutate(Description = str_wrap(Description, width = 20))
    
    
    y <- df  %>% 
      mutate(Description = gsub("_", " ", Description)) %>%
      mutate(Description = str_wrap(Description, width = 20)) %>% 
      filter(Description %in% temp$Description)%>%
      group_by(Description)%>% summarise(mean_NES = mean(NES))
    y <- merge(y, temp, by = "Description", all = T)
    y$proportion <- y$Freq/num
    y$Description <- factor(y$Description, levels = y$Description[order(y$proportion, y$mean_NES, decreasing = T)])
    
    f <- ggplot(y, aes(proportion, Description, fill=mean_NES)) + 
      geom_col(orientation='y') + 
      scale_fill_continuous(low='steelblue1', high='tomato2', guide=guide_colorbar(reverse=TRUE)) + 
      theme_minimal() + ylab(NULL) + ggtitle(module)
    print(f)}
}



endf <- lapply(enlist, lapply, lapply, function(x) if(class(x)!= "try-error") x@result)

endf <- lapply(endf, lapply, function(x) rbindlist(x, idcol = "db", fill = T))
endf <- lapply(endf, function(x) rbindlist(x, idcol = "module", fill = T))
endf <- rbindlist(endf, idcol = "sample", fill = T)
endf <- as.data.frame(endf)
endf <- endf[endf$NES > 0,]

pdf("4_enrichment.pdf", width = 4, height = 3)
temp <- sapply(unique(endf$module), function(x) {
  print(x)
  enrichment_barplots_new(endf, x,cells)
})
dev.off()




enlist_cepo <- lapply(cepo, function(markers){
  enresult <- lapply(markers[!sapply(markers, is.null)], function(x){
    hall <- try(GSEA(generate_genelist_cepo(x), TERM2GENE = Hall))
    kegg <- try(GSEA(generate_genelist_cepo(x), TERM2GENE = KEGG))
    reactome <- try(GSEA(generate_genelist_cepo(x), TERM2GENE = Reactome))
    return(list(hall = hall, kegg = kegg, reactome = reactome))
  })
  return(enresult)
})
saveRDS(enlist_cepo, "4_enlist_cepo.Rds")


endf <- lapply(enlist_cepo, lapply, lapply, function(x) if(class(x)!= "try-error") x@result)

endf <- lapply(endf, lapply, function(x) rbindlist(x, idcol = "db", fill = T))
endf <- lapply(endf, function(x) rbindlist(x, idcol = "module", fill = T))
endf <- rbindlist(endf, idcol = "sample", fill = T)
endf <- as.data.frame(endf)
endf <- endf[endf$NES > 0,]

pdf("4_enrichment_cepo.pdf", width = 4, height = 3)
temp <- sapply(unique(endf$module), function(x) enrichment_barplots_new(endf, x,cells))
dev.off()

# metas <- list.files("~/regulon/new_regulon/archetype/2024_02/sup/newtral/srt_meta", full.names = T) %>% lapply(. , readRDS)
# names(metas) <- gsub(".Rds", "", list.files("~/regulon/new_regulon/archetype/2024_02/sup/newtral/srt_meta", full.names = F))
# 
# 
# 
# head(metas[[1]])
# metas <- rbindlist(metas,idcol  = "sample")
# head(metas)
# metas <- as.data.frame(metas)
# par(mfrow = c(3, 3))
# plot(density(metas$Module1))
# abline(v = 0.7)
# plot(density(metas$Module2))
# abline(v = 0.75)
# plot(density(metas$Module3))
# abline(v =1.4)
# plot(density(metas$Module4))
# abline(v =0.55)
# plot(density(metas$Module5))
# abline(v =0.5)
# plot(density(metas$Module6))
# abline(v =0.6)
# plot(density(metas$Module7))
# abline(v =0.4)
# plot(density(metas$Module8))
# abline(v =0.4)
# 
# cutoff <- setNames(c(0.7, 0.75, 1.4, 0.55, 0.5, 0.6, 0.4, 0.4), paste0("Module", 1:8))
# 
# sapply(paste0("Module", 1:8), function(i){
#   table(data.frame(sample = metas$sample, cell = metas[,i] >=  as.numeric(cutoff[i])))[,2]
# })
# 


# combine plots -----------------------------------------------------------
library(stringr)
enlist <- readRDS("~/CASCADEpaper/paper/Fig5_archetype/4_enlist.Rds")
endf <- lapply(enlist, lapply, lapply, function(x) if(class(x)!= "try-error") x@result)

endf <- lapply(endf, lapply, function(x) rbindlist(x, idcol = "db", fill = T))
endf <- lapply(endf, function(x) rbindlist(x, idcol = "module", fill = T))
endf <- rbindlist(endf, idcol = "sample", fill = T)
endf <- as.data.frame(endf)
endf <- endf[endf$NES > 0,]
temp <- endf %>% select(sample,module, ID) %>% distinct() %>% group_by(module, ID) %>% summarise(n = n()) %>% arrange(desc(n)) %>% dplyr::slice(1:5)%>% 
  mutate(Description = gsub("_", " ", ID)) %>%
  mutate(Description = str_wrap(Description, width = 30))

NES <- apply(temp,1, function(x) endf %>% filter(module == x[1]&ID == x[2]) %>% summarise(mean = mean(NES)) %>% pull(mean))
temp$avg_NES <- NES
temp <- temp[order(temp$module, temp$n),]
temp$Description_lv <- factor(paste(temp$ID, temp$module), levels = paste(temp$ID, temp$module))


ggplot2::ggplot(temp, aes(n, Description_lv, fill=avg_NES)) + 
    geom_col(orientation='y') + 
    scale_fill_continuous(low="#FEE0D2", high="#DE2D26", guide=guide_colorbar(reverse=TRUE)) + 
    theme_classic(base_size = 16) + ylab(NULL)+
  facet_grid(module~., scales = "free", space = "free") + 
  theme(legend.position = "top", legend.key.width = unit(2,"cm")) + 
  scale_y_discrete(labels = setNames(temp$Description, temp$Description_lv))+
  labs(fill = "mean NES", x = "n")
ggsave("pdf/combined_enrich.pdf", width = 8, height = 12)

endf <- rbind(mastdf, cepodf)
temp2 <- endf %>% select(sample,module, ID) %>% distinct() %>% group_by(module, ID) %>% summarise(n = n()) %>% arrange(desc(n)) %>% slice(1:10)
temp <- endf %>% select(sample,module, ID) %>% distinct() %>% filter(ID %in% temp2$ID)
temp <- as.data.frame.matrix(table(temp$module, temp$ID))
temp <- melt(t(temp))
NES <- apply(temp,1, function(x) endf %>% filter(module == x[2]&ID == x[1]) %>% summarise(mean = mean(NES)) %>% pull(mean))
temp$avg_NES <- NES
temp <- temp %>% mutate(Description = gsub("_", " ", Var1)) %>%
  mutate(Description = str_wrap(Description, width = 20))
ggplot2::ggplot(temp , aes(value, Description, fill=avg_NES)) + 
  geom_col(orientation='y') + 
  scale_fill_continuous(low="#FEE0D2", high="#DE2D26", guide=guide_colorbar(reverse=TRUE)) + 
  theme_bw() + ylab(NULL)+
  facet_grid(.~Var2, scales = "free", space = "free") + 
  theme(text = element_text(size = 10), legend.position = "bottom") 
