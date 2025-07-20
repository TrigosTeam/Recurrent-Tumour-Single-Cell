library(Seurat)
library(scran)
library(scater)
library(data.table)
library(dplyr)
library(tibble)
library(patchwork)
library(tidyverse)
library(MAST)
library(pals)
library(circlize)
library(ComplexHeatmap)
source("~/functions.R")
library(Hmisc)

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

paths <- system("realpath ~/230113_seurat_intron/*/*/*non_tumour_removed.Rds", intern = T)
paths_split <- unlist(lapply(strsplit(paths, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
paths2 <- system("realpath ~/240603_seurat_intron/CA0027/dura*/*non_tumour_removed_v4.Rds", intern = T)
paths_split2 <- unlist(lapply(strsplit(paths2, split = "/"), function(x) paste(x[6], x[7], sep = "_")))

paths <- c(paths, paths2)
paths_split <- c(paths_split, paths_split2)
paths <- setNames(paths,paths_split)
meta <- readRDS("~/Fig5_archetype/5_integrated_final_module_meta.Rds")
subclones <- readRDS("~/Fig3_genomics_normal_cells/ATAClone/subclones.Rds")

sub_path <- system("realpath ~/240628_ATAClone/*/*/leiden_cluster.Rds", intern = T)

gf = rtracklayer::import("~/reference_annotation/multiome/genes.gtf")
gf_genes <- gf[gf$type =="gene"]


DEGs <- list()
DE_region <- list()
for (p in  sub_path){
  sample = str_split(p, "/")[[1]][7]
  patient = str_split(p, "/")[[1]][6]
  i = paste(patient, sample,sep = "_")
  subclone <- subclones[[paste(patient, sample, sep = '_')]]
  srt<- readRDS(paths[i])
  
  srt <- subset(srt, cells = names(subclone))
  DefaultAssay(srt) <- "RNA"
  lib.median <- median(srt$nCount_RNA)
  target_pseudocount <- 1
  srt <- srt %>%
    NormalizeData(scale.factor = lib.median / target_pseudocount) %>%
    ScaleData(do.scale = F)

  
  # 
  # srt <- srt %>% 
  #   RunPCA(features = top_genes2[1:2000]) %>% 
  #   FindNeighbors(dims = 1:20) %>% 
  #   FindClusters(resolution = 1) %>% 
  #   RunUMAP(dims = 1:20)
  # 
  
  
  
  # joint_cn <- readRDS(paste("~/240628_ATAClone", patient, sample, "joint_cn.Rds", sep = "/"))
  # subcloneall <- readRDS(p)
  srt$subclone <- subclone[colnames(srt)]
  Idents(srt) <- srt$subclone
  if (length(unique(as.character(subclone)))>1){
    ncomb <- combn(unique(as.character(subclone)), 2, simplify = F)
    
    markerdf <- lapply(ncomb, function(comb){
      df <- FindMarkers(srt,ident.1 = comb[1], ident.2 = comb[2],  test.use = "MAST")
      if (nrow(df) > 1){
        return(data.frame(df, comp = paste(comb, collapse = "|")) %>% rownames_to_column("gene"))
      }else(
        return(data.frame(comp = paste(comb, collapse = "|")))
      )
      
    })
    markerdfs <- rbindlist(markerdf, fill = T)
    DEGs[[i]] <- markerdfs
    
    single_cn <- readRDS(paste("~/240628_ATAClone", patient, sample, "single_cn.Rds", sep = "/"))
    single_cn <- single_cn[, colnames(srt)]
    
    srt[["ATAClone"]] <-CreateAssayObject(data = single_cn) 
    regiondf <- lapply(ncomb, function(comb){
      df <- FindMarkers(srt,ident.1 = comb[1], ident.2 = comb[2],  test.use = "MAST", assay = "ATAClone")
      if (nrow(df) > 1){
        return(data.frame(df, comp = paste(comb, collapse = "|")) %>% rownames_to_column("gene"))
      }else(
        return(data.frame( comp = paste(comb, collapse = "|")))
      )
    })
    regiondfs <- rbindlist(regiondf, fill = T)
    DE_region[[i]][["region"]] <- regiondfs
    
    if(any(ncol(markerdfs)<2, ncol(regiondfs) <2)){
      print("no DEG/ no DE regions")
    }else{
      chrom <- sapply(strsplit(unique(regiondfs$gene), split = ".", fixed = T), "[", 1)
      ranges<- sapply(strsplit(unique(regiondfs$gene), split = ":", fixed = T), "[", 2)
      start <- sapply(strsplit(ranges, split = "-", fixed = T), "[", 1)
      end <- sapply(strsplit(ranges, split = "-", fixed = T), "[", 2)
    
      bingr <- as(data.frame(chrom = chrom, start = start, end = end), "GRanges")
      hit <- findOverlaps(bingr, gf_genes)
      ind <- split( to(hit), from(hit))
    
      genes <- sapply(ind, function(x) intersect(gf_genes$gene_name[x], rownames(srt)))
      names(genes) <- unique(regiondfs$gene)[as.numeric(names(genes))]
      DE_region[[i]][["genes"]] <- genes
      
      prop <- sapply(unique(markerdfs$comp), function(cl){
        deg <- markerdfs$gene[markerdfs$comp == cl&markerdfs$p_val_adj<0.05]
        gFC <- setNames(markerdfs$avg_log2FC[markerdfs$comp == cl&markerdfs$p_val_adj<0.05], deg)
        region <- regiondfs$gene[regiondfs$comp == cl&regiondfs$p_val_adj <0.05]
        region_genes <- unique(unlist(genes[region]))
        rFC <- setNames(regiondfs$avg_log2FC[regiondfs$comp == cl&regiondfs$p_val_adj <0.05], region)
        rFC <- lapply(region, function(r) setNames( rep(rFC[r], length(genes[[r]])), as.character(genes[[r]])))
        rFC <- unlist(rFC)
        ingenes <- intersect(deg, region_genes)
        return(data.frame(prop = sum(sign(gFC[ingenes]) == sign(rFC[ingenes]))/length(deg), cor =  cor(gFC[ingenes], rFC[ingenes])))
      })
      DE_region[[i]][["prop"]] <- prop 
      
    }
    
  }

  genes <- genes[sapply(genes, length)>0]


  single_cn_sub <- single_cn[names(genes), names(subclone)]

  srt <- AddModuleScore(srt, genes)
  exp_mat <- srt@meta.data[names(subclone), grep("Cluster", colnames(srt@meta.data))]
  colnames(exp_mat) <- names(genes)
  corlist[[i]] <- list(CN = single_cn_sub, exp = exp_mat)
  print(i)
} 







saveRDS(DEGs, "~/Fig3_genomics_normal_cells/ATAClone/DEGs.Rds")
saveRDS(DE_region, "~/Fig3_genomics_normal_cells/ATAClone/DE_region.Rds")




test <- lapply(DE_region, '[[', "prop")

test <- lapply(1:22, function(ind){
  markerdfs <- rbindlist(DEGs[[ind]])
  regiondfs <- DE_region[[ind]]$region
  genes <- DE_region[[ind]]$genes
  prop <- sapply(unique(markerdfs$comp), function(cl){
    deg <- markerdfs$gene[markerdfs$comp == cl&markerdfs$p_val_adj<0.05]
    gFC <- setNames(markerdfs$avg_log2FC[markerdfs$comp == cl&markerdfs$p_val_adj<0.05], deg)
    region <- regiondfs$gene[regiondfs$comp == cl&regiondfs$p_val_adj <0.05]
    region_genes <- unique(unlist(genes[region]))
    rFC <- setNames(regiondfs$avg_log2FC[regiondfs$comp == cl&regiondfs$p_val_adj <0.05], region)
    rFC <- lapply(region, function(r) setNames( rep(rFC[r], length(genes[[r]])), as.character(genes[[r]])))
    rFC <- unlist(rFC)
    ingenes <- intersect(deg, region_genes)
    if(length(ingenes) > 3){
      
      cortest <- cor.test(gFC[ingenes], rFC[ingenes], method = "spearman")
      return(data.frame(prop = sum(sign(gFC[ingenes]) == sign(rFC[ingenes]))/length(deg), cor =  cortest$estimate, corP = cortest$p.value))
    }else if(length(ingenes) == 0 ){
      return(data.frame(prop = NA, cor =  NA, corP = NA))
    }else{
      return(data.frame(prop = sum(sign(gFC[ingenes]) == sign(rFC[ingenes]))/length(deg), cor =  NA, corP = NA))
    }
  })
})
names(test) <- names(DEGs)

hist(unlist(sapply(test[sapply(test, class)!= "logical"], function(x) x["cor",])), breaks = 100)
range(unlist(sapply(test[sapply(test, class)!= "logical"], function(x) x["cor",]))[unlist(sapply(test[sapply(test, class)!= "logical"], function(x) x["corP",]))<0.05], na.rm = T)
# [1] -0.1902691  0.8780730

summary(unlist(sapply(test[sapply(test, class)!= "logical"], function(x) x["cor",])))
summary(unlist(sapply(test[sapply(test, class)!= "logical"], function(x) x["prop",])))
hist(unlist(sapply(test[sapply(test, class)!= "logical"], function(x) x["prop",])), breaks = 100)
range(unlist(sapply(test[sapply(test, class)!= "logical"], function(x) x["prop",])), na.rm = T)
# [1] 0.0000000 0.8653846

lapply(DE_region, lapply, function(y) lapply(y[["prop"]], function(x) data.frame(avg_prop = mean(x["prop",], n = ncol(x)))))

test3 <- lapply(DE_region, function(x) x[["region"]]%>% filter(p_val_adj < 0.05) %>% group_by(comp) %>% summarise(n = n()))
test2 <- lapply(test, function(x) data.frame(prop = as.numeric(x["prop", ]), ncomp = ncol(x), comp = colnames(x)))

temp <- lapply(names(test), function(x) merge(test2[[x]], test3[[x]], by = "comp"))
names(temp) <- names(test)
temp <- rbindlist(temp, idcol = "sample")
head(temp)

plot(temp$prop, temp$n)
cor.test(temp$prop, temp$n, method = "spearman", alternative = "greater")
saveRDS(temp, "~/Fig3_genomics_normal_cells/ATAClone/prop_vs_nregion.Rds")


library(patchwork)
library(dplyr)
library(Seurat)
library(data.table)
library(pals)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(ggnewscale)
library(dplyr)
library(tibble)
library(data.table)
library(ggstance)
library(MAST)
library(gridExtra)
library(stringr)
library(ggplot2)
library(tidyverse)


Hall <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
KEGG <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")%>%
  dplyr::select(gs_name, gene_symbol)
Reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME")%>%
  dplyr::select(gs_name, gene_symbol)

C2 <- msigdbr(species = "Homo sapiens", category = "C4") %>%
  dplyr::select(gs_name, gene_symbol) # curated gene sets
C4 <- msigdbr(species = "Homo sapiens", category = "C4") %>%
  dplyr::select(gs_name, gene_symbol) # computational gene sets
C5 <- msigdbr(species = "Homo sapiens", category = "C5") %>%
  dplyr::select(gs_name, gene_symbol) # GO term ontology
C6 <- msigdbr(species = "Homo sapiens", category = "C6") %>%
  dplyr::select(gs_name, gene_symbol) #oncology
C7 <- msigdbr(species = "Homo sapiens", category = "C7") %>%
  dplyr::select(gs_name, gene_symbol) # immunologic signature
C8 <- msigdbr(species = "Homo sapiens", category = "C8") %>%
  dplyr::select(gs_name, gene_symbol) #cell type


DEGs <- readRDS("~/Fig3_genomics_normal_cells/ATAClone/DEGs.Rds")
DEGs <- lapply(DEGs, function(x) {
  temp <- rbindlist(x, fill = T)
  return(split(temp, temp$comp))
  })
generate_genelist <- function(markerdf){
  df <- markerdf[markerdf$p_val_adj < 0.05,]
  temp <- setNames(markerdf$avg_log2FC, markerdf$gene)
  temp <- sort(temp, decreasing = T)
  return(temp)
}


genelists <- lapply(DEGs, lapply, generate_genelist)
enrichlist <- lapply(genelists, lapply, function(x) try(GSEA(geneList = x, TERM2GENE = Hall)))
endf <- lapply(enrichlist, function(x){
  temp <- x[sapply(x, class) == "gseaResult"]
  temp <- rbindlist(lapply(temp, function(y) y@result), fill = T, idcol = "comp")
})
endf <- rbindlist(endf, idcol = "sample", fill = T)
saveRDS(endf, "~/Fig3_genomics_normal_cells/ATAClone/endf.Rds")
#### DEGs proportion
# DEGs <- readRDS("~/Fig3_genomics_normal_cells/ATAClone/DEGs.Rds")
# 
# for (p in  sub_path){
#   sample = str_split(p, "/")[[1]][7]
#   patient = str_split(p, "/")[[1]][6]





#   joint_cn <- readRDS(paste("~/240628_ATAClone", patient, sample, "joint_cn.Rds", sep = "/"))
#   subclone <- subclones[[paste(patient, sample, sep = "_")]]
#   subclone <- split(names(subclone), subclone)
#   test <- sapply(subclone, function(cell) rowMeans(joint_cn[, cell]))
#   
#   chrom <- sapply(strsplit(rownames(single_cn), split = ".", fixed = T), "[", 1)
#   ranges<- sapply(strsplit(rownames(single_cn), split = ":", fixed = T), "[", 2)
#   start <- sapply(strsplit(ranges, split = "-", fixed = T), "[", 1)
#   end <- sapply(strsplit(ranges, split = "-", fixed = T), "[", 2)
# 
#   bingr <- as(data.frame(chrom = chrom, start = start, end = end), "GRanges")
#   hit <- findOverlaps(bingr, gf_genes)
#   ind <- split( to(hit), from(hit))
  