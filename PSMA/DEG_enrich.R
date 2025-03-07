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
library(scran)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(ggnewscale)
library(dorothea)
library(tibble)
library(ggstance)
library(gridExtra)
library(stringr)
library(grid)

setwd("~/CASCADEpaper/paper/PSMA")


phenotype_meta <- readRDS("~/CASCADEpaper/paper/PSMA/phenotype_meta.Rds")
samples <- unique(phenotype_meta$sample[phenotype_meta$group %in% c("AR+/NE-", "AR+/NE+")])
phenotype_meta <- split(phenotype_meta, phenotype_meta$sample)



paths <- system("realpath /trigos_team/CASCADE/Analysis/230113_seurat_intron/*/*/*non_tumour_removed.Rds", intern = T)
paths_split <- unlist(lapply(strsplit(paths, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
paths2 <- system("realpath /trigos_team/CASCADE/Analysis/240603_seurat_intron/CA0027/dura*/*non_tumour_removed_v4.Rds", intern = T)
paths_split2 <- unlist(lapply(strsplit(paths2, split = "/"), function(x) paste(x[6], x[7], sep = "_")))

paths <- c(paths, paths2)
paths_split <- c(paths_split, paths_split2)

paths <- setNames(paths, paths_split)

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
C8 <- msigdbr(species = "Homo sapiens", category = "C8") %>%
  dplyr::select(gs_name, gene_symbol)

enrichment_barplots <- function(geneList, geneSets, groupname){
  enresult <- try(GSEA(geneList, TERM2GENE = geneSets))
  if (class(enresult)!= "try-error"){
    if(nrow(enresult@result)>0 ){
      y <- arrange(enresult, abs(NES)) %>% 
        group_by(sign(NES)) %>% 
        mutate(Description = gsub("_", " ", Description)) %>%
        mutate(Description = str_wrap(Description, width = 20)) %>% 
        dplyr::slice(1:5)
      f <- ggplot(y, aes(NES, fct_reorder(Description, NES), fill=p.adjust), showCategory=10) + 
        geom_col(orientation='y') + 
        scale_fill_continuous(low='blue', high='red', guide=guide_colorbar(reverse=TRUE)) + 
        theme_minimal() + ylab(NULL) + ggtitle(groupname)
      print(f)
      return(enresult)}
    else{
      print("no enrichment")
    }
  }else{print("not enough gene mapped")}
}

generate_genelist <- function(markerdf){
  temp <- setNames(markerdf$avg_log2FC, markerdf$gene)
  temp <- sort(temp, decreasing = T)
  return(temp)
}

degenrich_list <- list()
for (i in samples){
  srt <- readRDS(paths[i])
  meta <- phenotype_meta[[i]]
  rownames(meta) <-gsub(paste0(unique(meta$sample), "_"), "", rownames(meta))
  srt <- AddMetaData(srt, metadata = meta)
  srt$pheno <- factor(paste0(srt$AR_pheno, srt$FOLH1_pheno))
  
  com <- lapply(c("ARhigh", "ARlow", "ARneg"), function(x){
    if(length(unique(grep(x, srt$pheno, value = T)))>1){
      return(combn(grep(x, levels(srt$pheno), value = T), 2, simplify = F))
    }
  })
  com <- unlist(com, recursive = F)
  names(com) <- sapply(com, function(x) paste(x, collapse = " vs "))
  
  deg <- lapply(com, function(x){
    marker <- try(FindMarkers(srt, ident.1 =  x[1], ident.2 = x[2], test.use = "MAST", 
                          min.cells.feature = 20,min.cells.group = 20, group.by = "pheno"))
    if(class(marker)!= "try-error"){
    marker$group1 <- x[1]
    marker$group2 <- x[2]
    marker <- marker %>% filter(p_val_adj < 0.05) %>% rownames_to_column("gene")
    return(marker)} else{return(NULL)}
  })
  
  enlist <- lapply(deg, function(markerdf){
    if(!is.null(markerdf)){
      enresult <- lapply(list(Hall,Reactome, KEGG, C5, C6), function(geneSets){
        try(GSEA(generate_genelist(markerdf), TERM2GENE = geneSets))
      })
      names(enresult) <- c("Hallmaker", "Reactome", "KEGG", "GOterm", "Oncology")
      enresult <- lapply(enresult[sapply(enresult, class)!= "try-error"], function(x) x@result)
      return(rbindlist(enresult, idcol = "geneSet", fill = T))
    }
  })
  endf <- rbindlist(enlist, idcol = "compare", fill = T)
  endf$group <- NA
  endf$group[grep("ARhigh", endf$compare, fixed = T)] <- "ARhigh"
  endf$group[grep("ARlow", endf$compare, fixed = T)] <- "ARlow"
  endf$group[grep("ARneg", endf$compare, fixed = T)] <- "ARneg"
  degenrich_list[[i]][["DEG"]] <- deg
  degenrich_list[[i]][["enrich"]] <- endf
  print(i)
}

saveRDS(degenrich_list, "degenrich_list.Rds")



enrichlist <- lapply(degenrich_list, function(x) x[["enrich"]])
enrichlist <- rbindlist(enrichlist, idcol = "sample", fill = T)
enrichlist$sign <- ifelse(sign(enrichlist$NES)>0,  "increase", "decrease")

rm <- as.data.frame.matrix(table(enrichlist$ID, enrichlist$sign))
rm <- names(which(apply(rm, 1, function(x) all(x>0))))

enrichlist %>%  filter(!ID %in% rm) %>% 
  group_by(sign, group,ID) %>% 
  summarise(avgNES= mean(NES), n = n()) -> temp

temp <- as.data.frame(table(temp$sign, temp$ID))
highlightname <- gsub("_", " ", as.character(temp$Var2[temp$Freq==3]))

enrichlist %>%  filter(!ID %in% rm) %>% 
  group_by(sign, group,ID) %>% 
  summarise(avgNES= mean(NES), n = n()) %>% 
  group_by(sign, group) %>% arrange(desc(n)) %>% 
  slice(1:5) %>%
  mutate(n = sign(avgNES)*n) %>% 
  mutate(ID = str_wrap(gsub("_", " ", ID), 50)) %>% 
  mutate(name = paste(group, ID)) %>% 
  group_by(group) %>% 
  arrange(n) %>% 
  mutate(name = factor(name, levels = name)) -> test
test$color <- ifelse(test$ID %in% highlightname&test$sign=="increase", "red", 
                     ifelse(test$ID %in% highlightname&test$sign=="decrease", "blue", "black"))

g <- ggplot(test, aes(x = n, y = name, fill= avgNES)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  facet_grid(group~. ,space = "free", scales = "free")+
  scale_fill_gradient2(high = "tomato3", mid = "white", low = "steelblue3")+
  labs( x = "number of samples with enrichment", y = "") + 
  scale_y_discrete(labels = test$ID, breaks = test$name) +
  scale_x_continuous(breaks = c(-8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4), labels = abs(c(-8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4)))+
  

pdf("enrichment.pdf", height = 6, width = 9)
print(g)
dev.off()


