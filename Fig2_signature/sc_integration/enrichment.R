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
  df <- markerdf[markerdf$p_val_adj < 0.05,]
  temp <- setNames(markerdf$avg_log2FC, markerdf$gene)
  temp <- sort(temp, decreasing = T)
  return(temp)
}

inter_patient_DEG <- readRDS("~/Fig2_signature/sc_integration/DEG/inter_patient_DEG.Rds")

all <- rbind(Hall, KEGG, Reactome, C6, C8 )
# genelists <- lapply(inter_patient_DEG, generate_genelist)
# enlists <- lapply(genelists, function(x) try(GSEA(geneList = x, TERM2GENE = all)))
# saveRDS(enlists, "~/Fig2_signature/sc_integration/DEG/enlist_interpatient.Rds")

cna_compare <- function(degs,cna1, cna2, test = c("greater", "less")){
  if(is.null(degs)){
    return(0)
  }
  temp <- sapply(degs, function(gene){
    CN1 <- mean(cna1$SCN[cna1$gene == gene])
    CN2 <- mean(cna2$SCN[cna2$gene == gene])
    if(test == "greater"){
      return(CN1 > CN2)}
    else if(test == "less"){
      return(CN1 < CN2)
    }
  }, USE.NAMES = F)
  return(temp)
}

geneCNAlist <- readRDS("~/Fig3_genomics_normal_cells/copy_number/geneCNAlist.Rds")
names(geneCNAlist)[8] <- "CA0035_paraaortic_lymph_node_1"

patients <- substr(names(geneCNAlist),1,6)
geneCNApp <- split(geneCNAlist, patients)
cna_aff <- lapply(names(inter_patient_DEG), function(p){
  p1 <- strsplit(p, "-")[[1]][1]
  p2 <- strsplit(p, "-")[[1]][2]
  DEG <- inter_patient_DEG[[p]]
  cna1 <- rbindlist(geneCNApp[[p1]], idcol = "sample")
  cna2 <- rbindlist(geneCNApp[[p2]], idcol = "sample")
  sign_pos<-  DEG$avg_log2FC > 0
  sign_neg <- DEG$avg_log2FC < 0
  degs <- rownames(DEG)
  result_pos <- cna_compare(degs, cna1, cna2, test = "greater")
  result_neg <- cna_compare(degs, cna1, cna2, test = "less")
  
  return(DEG[(sign_pos&result_pos)|(sign_neg&result_neg), ])
  })
names(cna_aff) <- names(inter_patient_DEG)
cna_free <- lapply(names(inter_patient_DEG), function(p){
  p1 <- strsplit(p, "-")[[1]][1]
  p2 <- strsplit(p, "-")[[1]][2]
  DEG <- inter_patient_DEG[[p]]
  cna1 <- rbindlist(geneCNApp[[p1]], idcol = "sample")
  cna2 <- rbindlist(geneCNApp[[p2]], idcol = "sample")
  sign_pos<-  DEG$avg_log2FC > 0
  sign_neg <- DEG$avg_log2FC < 0
  degs <- rownames(DEG)
  result_pos <- cna_compare(degs, cna1, cna2, test = "greater")
  result_neg <- cna_compare(degs, cna1, cna2, test = "less")
  
  return(DEG[!((sign_pos&result_pos)|(sign_neg&result_neg)), ])
})
names(cna_free) <- names(inter_patient_DEG)
  
genelists <- lapply(cna_aff, generate_genelist)
enlists <- lapply(genelists, function(x) try(GSEA(geneList = x, TERM2GENE = all)))
saveRDS(enlists, "~/Fig2_signature/sc_integration/DEG/enlist_interpatient_cna_aff.Rds")
genelists <- lapply(cna_free, generate_genelist)
enlists <- lapply(genelists, function(x) try(GSEA(geneList = x, TERM2GENE = all)))
saveRDS(enlists, "~/Fig2_signature/sc_integration/DEG/enlist_interpatient_cna_free.Rds")

# for(i in list.files("~/Fig2_signature/sc_integration/DEG/inter_lesion")){
#   degs <- readRDS(paste0("~/Fig2_signature/sc_integration/DEG/inter_lesion/", i))
#   genelists <- lapply(degs, generate_genelist)
#   enlists <- lapply(genelists, function(x) try(GSEA(geneList = x, TERM2GENE = all)))
#   saveRDS(enlists, paste0("~/Fig2_signature/sc_integration/DEG/inter_lesion/enrich_", i))
#   print(i)
# }

for(i in list.files("~/Fig2_signature/sc_integration/DEG/inter_lesion", pattern ="^CA")){
  DEGs <- readRDS(paste0("~/Fig2_signature/sc_integration/DEG/inter_lesion/", i))
  cna_aff <- lapply(names(DEGs), function(p){
    p1 <- strsplit(p, "-")[[1]][1]
    p2 <- strsplit(p, "-")[[1]][2]
    DEG <- DEGs[[p]]
    cna1 <- geneCNAlist[[p1]]
    cna2 <- geneCNAlist[[p2]]
    sign_pos<-  DEG$avg_log2FC > 0
    sign_neg <- DEG$avg_log2FC < 0
    degs <- rownames(DEG)
    result_pos <- cna_compare(degs, cna1, cna2, test = "greater")
    result_neg <- cna_compare(degs, cna1, cna2, test = "less")
    
    return(DEG[(sign_pos&result_pos)|(sign_neg&result_neg), ])
  })
  names(cna_aff) <- names(DEGs)
  cna_free <- lapply(names(DEGs), function(p){
    p1 <- strsplit(p, "-")[[1]][1]
    p2 <- strsplit(p, "-")[[1]][2]
    DEG <- DEGs[[p]]
    cna1 <- geneCNAlist[[p1]]
    cna2 <- geneCNAlist[[p2]]
    sign_pos<-  DEG$avg_log2FC > 0
    sign_neg <- DEG$avg_log2FC < 0
    degs <- rownames(DEG)
    result_pos <- cna_compare(degs, cna1, cna2, test = "greater")
    result_neg <- cna_compare(degs, cna1, cna2, test = "less")
    
    return(DEG[!((sign_pos&result_pos)|(sign_neg&result_neg)), ])
  })
  names(cna_free) <- names(DEGs)
  genelists <- lapply(cna_aff, generate_genelist)
  enlists <- lapply(genelists, function(x) try(GSEA(geneList = x, TERM2GENE = all)))
  saveRDS(enlists, paste0("~/Fig2_signature/sc_integration/DEG/inter_lesion/enrich_cna_aff", i))
  genelists <- lapply(cna_free, generate_genelist)
  enlists <- lapply(genelists, function(x) try(GSEA(geneList = x, TERM2GENE = all)))
  saveRDS(enlists, paste0("~/Fig2_signature/sc_integration/DEG/inter_lesion/enrich_cna_free", i))
  print(i)
}

# downstream -------------------------------------------------------------
# 
# enlist_interpatient <- readRDS("~/Fig2_signature/sc_integration/DEG/enlist_interpatient.Rds")
# enlist <- lapply(enlist_interpatient, function(x) x@result)
# endf <- rbindlist(enlist, idcol = "comp", fill = T)
# View(endf)
# 
# inter_patient_DEG <- readRDS("~/Fig2_signature/sc_integration/DEG/inter_patient_DEG.Rds")
# inter_patient_DEG <- lapply(inter_patient_DEG, function(x) x[x$p_val_adj<0.05,])
# 
# ndeg <- sapply(inter_patient_DEG, nrow)
# patients <- gsub(".Rds","",list.files("~/Fig2_signature/sc_integration/DEG/inter_lesion", pattern = "^CA"))
# mat <- matrix(0, nrow = 9, ncol = 9, dimnames = list(patients, patients))
# for(i in names(ndeg)){
#   p1 <- strsplit(i, "-")[[1]][1]
#   p2 <- strsplit(i, "-")[[1]][2]
#   mat[p1, p2] <- ndeg[[i]]
#   mat[p2, p1] <- ndeg[[i]]
# }
# 
# dmat <- as.dist(mat)
# hc <- hclust(dmat)
# plot(hc)
# library(ComplexHeatmap)
# 
# Heatmap(mat, name = "No.DEG")

enrich_list <- list()
for( p in names(geneCNApp)){
  DEG <- readRDS(paste0("~/Fig2_signature/sc_integration/DEG/inter_lesion/",p, ".Rds"))
  DEG <- DEG[DEG$cluster %in% names(geneCNAlist),]
  DEG <- DEG %>% filter(p_val_adj < 0.05)
  comp <- comps[[p]]
  enrich <- lapply(unique(DEG$cluster), function(i){
    subcomp <- split(rownames(comp[[i]]), comp[[i]]$id)
    cna_aff <- unique(unlist(subcomp[c("pos >1","neg <1")], use.names = F))
    cna_aff_list <- generate_genelist(DEG %>% filter(cluster == i & gene %in% cna_aff))
    cna_affenrich <- try(GSEA(geneList = cna_aff_list, TERM2GENE = all))
    cna_free <- unique(unlist(subcomp[c("pos <1","neg >1")]))
    cna_free_list <- generate_genelist(DEG %>% filter(cluster == i & gene %in% cna_free))
    cna_freeenrich <- try(GSEA(geneList = cna_free_list, TERM2GENE = all))
    print(i)
    return(rbindlist(list(cna_aff = cna_affenrich@result, cna_free = cna_freeenrich@result), idcol = "CNA", fill = T))
  })
  names(enrich) <- unique(DEG$cluster)
  enrich_list[[p]] <- rbindlist(enrich, idcol = "sample", fill = T)
  print(p)
}
saveRDS(enrich_list, "~/Fig2_signature/sc_integration/DEG/inter_lesion/enrich_list.Rds")



### inter-organ DEGs --------

inter_organ_DEG <- readRDS("~/Fig2_signature/sc_integration/DEG/inter_organ_DEG.Rds")
markerlist <- split(inter_organ_DEG, inter_organ_DEG$cluster)
genelists <- lapply(markerlist, generate_genelist)
enrich_list <- lapply(genelists, function(x) {
  try(GSEA(geneList = x, TERM2GENE = all))
  })

enrichdf <- rbindlist( lapply(enrich_list, function(x) x@result), idcol = "organ", fill = T)
saveRDS(enrichdf, "~/Fig2_signature/sc_integration/DEG/inter-organ_enrichment.Rds")
