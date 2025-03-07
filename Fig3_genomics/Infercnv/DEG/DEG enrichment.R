library(patchwork)
library(dplyr)
library(Seurat)
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

paths <- system("realpath /trigos_team/CASCADE/Analysis/230113_seurat_intron/*/*/*non_tumour_removed.Rds", intern = T)
paths_split <- unlist(lapply(strsplit(paths, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
paths <- setNames(paths, paths_split)
subclone_anno <- readRDS("~/CASCADEpaper/paper/Fig3/Infercnv/subclone_anno.Rds")
sig_region <- readRDS("~/CASCADEpaper/paper/Fig3/Infercnv/region_test/sig_region.Rds")
multi_g <- readRDS("~/genomic/exon_only_reference.Rds")

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



for(i in names(subclone_anno)){
  srt <- readRDS(paths[i])
  if(sum(is.na(names(subclone_anno[[i]])))>0){
    subclone_anno[[i]][is.na(names(subclone_anno[[i]]))] <- "filtered"
    names(subclone_anno[[i]])[is.na(names(subclone_anno[[i]]))] <- setdiff(colnames(srt), names(subclone_anno[[i]]))
  }
  srt$subclone <- subclone_anno[[i]][colnames(srt)]
  DefaultAssay(srt) <- "RNA"
  Idents(srt) <- srt$subclone
  comb <- combn(names(which(table(subclone_anno[[i]]) > 5)), 2, simplify = F)
  DEG <- lapply(comb, function(x){
    temp <- FindMarkers(srt, ident.1 = x[1], ident.2 = x[2], test.use = "MAST")
    temp$gene <- rownames(temp)
    temp$comp <- paste(x[1], "vs", x[2])
    temp$ident1 <- x[1]
    temp$ident2 <- x[2]
    return(temp)
  })
  DEG <- rbindlist(DEG,fill = T)
  DEG <- DEG %>% filter(p_val_adj<=0.05)
  if (nrow(DEG) >0){
    regions <- sig_region[[i]]
    #check position of DEG
    genes <- apply(regions,1, function(regions){
      chrom <- as.numeric(regions["chr"])
      start <- as.numeric(regions["start"])
      end <- as.numeric(regions["end"])
      return(unique(multi_g$gene.name[multi_g$chrom==chrom&multi_g$start>=start&multi_g$end<=end]))
    })
    
    regions <- as.data.frame(regions)
    regions$genes <- sapply(genes, paste, collapse = ";")
    judge <- apply(DEG, 1,  function(x){
      gene <- x["gene"]
      ident1 <- x["ident1"]
      ident2 <- x["ident2"]
      temp <- which(sapply(genes, function(y) any(y == gene)))
      if(length(temp) >0){
        sbdf <- regions[temp, !colnames(regions)%in%c("chr", "start", "end", "comp", "genes")]
        
        if (x["avg_log2FC"] < 0){
          test <- sbdf[,paste(ident1,ident1, sep = ".")] < sbdf[,paste(ident2,ident2, sep = ".")] 
        }else{
          test <- sbdf[,paste(ident1,ident1, sep = ".")] > sbdf[,paste(ident2,ident2, sep = ".")]
        }
        temp2 <- as.data.frame(x)
        temp2["test",] <- test
        return(cbind(data.table::transpose(temp2) %>% `colnames<-`(rownames(temp2)), regions[temp, ]))}
      else{
        temp2 <- as.data.frame(x)
        temp2["test",] <- "not in"
        return(data.table::transpose(temp2) %>% `colnames<-`(rownames(temp2)))}
    })
    judge <- rbindlist(judge, fill = TRUE)
    saveRDS(judge, paste0("~/CASCADEpaper/paper/Fig3/Infercnv/DEG/", i, ".Rds"))
    # functional enrichment
    markers <- split(DEG, DEG$comp)
    geneLists <- lapply(markers, generate_genelist)
    enres <- list()
    pdf(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/enrichment/", i, ".pdf"), width = 5, height = 4)
    for(n in names(geneLists)){
      
      enres[[n]] <- mapply(enrichment_barplots, list(geneLists[[n]]), list(Hall, KEGG, C5, C6), as.list(paste("subclone",n, c("Hallmark", "KEGG", "GO term", "Oncology"))), SIMPLIFY = F)
    }
    
    dev.off()
    saveRDS(enres,paste0("~/CASCADEpaper/paper/Fig3/Infercnv/enrichment/", i, ".Rds" ))
    print(i)
  }else{
    print(i)
    print("no DEG")
  }
  
}



# later analysis ----------------------------------------------------------
DEGs <- list.files("~/CASCADEpaper/paper/Fig3/Infercnv/DEG", full.names = T, pattern = ".Rds") %>% lapply(., readRDS)
names(DEGs) <- gsub(".Rds", "", list.files("~/CASCADEpaper/paper/Fig3/Infercnv/DEG", full.names = F, pattern = ".Rds"))

DEGs2 <- lapply(DEGs, function(x) x[x$p_val_adj <=0.05,])
# #CA0046_paraaortic_lymph_node_15, CA0046_prostate_17, CA0046_lung_7,CA0083_paraaortic_lymph_node_12
# 
# 
# 
tb <- lapply(DEGs[sapply(DEGs, nrow)>0], function(x) as.data.frame(table(x$test, x$comp)))
tb <- lapply(tb, function(x){
  temp <- x %>% group_by(Var2) %>% mutate(total = sum(Freq)) %>% mutate(prop = Freq/total)
  temp <- dcast(temp, Var2 ~ Var1, value.var = "Freq")
  temp <- as.data.frame(temp) %>% column_to_rownames(var = "Var2")
  
  temp$total <- rowSums(temp)
  temp$`within%` <- 1 - temp$`not in`/temp$total
  #return(temp)
  if ("TRUE" %in% colnames(temp)){
  temp$agree <- temp$`TRUE`/temp$total
  }else{temp$agree <- 0}
  temp <- temp %>% rownames_to_column("comp")
  return(temp)
} )


tb <- rbindlist(tb, idcol = "sample", fill = T)
# tb <- dcast(tb, sample ~ Var1, value.var = "Freq")
# tb <- as.data.frame(tb)
# tb[, 2:4] <- apply(tb[, 2:4], 2, function(x){
#   x[is.na(x)] <- 0
#   return(x)})
# 
# tb$total <- tb$`FALSE`+ tb$`not in` +tb$`TRUE`
# tb$`within %` <- (tb$`FALSE`+tb$`TRUE`)/tb$total *100
# tb$`direction agree` <- tb$`TRUE`/tb$total * 100
# tb$`direction agree within` <- tb$`direction agree`/tb$`within %`*100
# region_num <- sapply(sig_region, nrow)
# tb$`within % adjusted by region number` <-   tb$`within %`/ region_num[tb$sample]
# summary(tb$`within %`)
# summary(tb$`within % adjusted by region number`)
saveRDS(tb, "~/CASCADEpaper/paper/Fig3/Infercnv/DEG_table.Rds")
write.table(tb, "~/CASCADEpaper/paper/Fig3/Infercnv/DEG_table.csv", quote = F, row.names = F, col.names = T)

# enrichment check --------------------------------------------------------
enrich <- list.files("~/CASCADEpaper/paper/Fig3/Infercnv/enrichment", full.names = T, pattern = ".Rds") %>% lapply(., readRDS)
names(enrich) <- gsub(".Rds", "", list.files("~/CASCADEpaper/paper/Fig3/Infercnv/enrichment", full.names = F, pattern = ".Rds"))

temp <- lapply(enrich, lapply, lapply, function(x){
  if(class(x) == "gseaResult"){
    return(x@result)
  }
})

temp <- lapply(temp, lapply, function(x) {
  names(x) <- c("Hallmark", "KEGG", "GO term", "Oncology")
  return(rbindlist(x, idcol = "db", fill = T))
})

enrichment <- lapply(temp, function(x) rbindlist(x, fill = T, idcol  = "subclone"))
enrichment <- rbindlist(enrichment, idcol = "sample", fill = T)
test <- apply(enrichment, 1, function(x){
  sample <- x["sample"]
    temp <- unlist(strsplit(x["core_enrichment"], split = "/"))
    subdf <- DEGs[[sample]][, 1:10]
    
    test <-subdf %>% filter(gene %in% temp & comp == x["subclone"]) %>% pull(test) %>% table %>% as.matrix()
    core_genes <- subdf %>% filter(gene %in% temp & comp == x["subclone"] & test == "TRUE") %>% pull(gene)
    test <- as.data.frame(t(test))
    test$sum <- sum(test)
    test$`% DEG`  <- test$sum/ nrow(subdf %>% filter(comp == x["subclone"]))
    test$core_genes <- paste(core_genes, collapse = "/")
return(test)
})
test <- rbindlist(test, fill = T)
test2 <- as.matrix(test[, c(1:4,6)])
test2[is.na(test2)] <- 0
test2 <- as.data.frame(test2)
test <- cbind (test2, test[, "core_genes"])
test$`% within` <- (1 - test$`not in`/test$sum)
test$`% TURE within` <- test$`TRUE` / test$sum 
enrichment <- cbind(enrichment, test)
summary(test$`% within`)

summarise_enrichment <- enrichment %>% group_by(sample, subclone) %>% summarize( n = n(), within= mean(`FALSE` + `TRUE`), agree = mean(`TRUE`))

saveRDS(enrichment, "~/CASCADEpaper/paper/Fig3/Infercnv/enrichment table.Rds")
write.table(enrichment, "~/CASCADEpaper/paper/Fig3/Infercnv/enrichment_table.csv", quote = F, row.names = F, col.names = T, sep = ",")
write.table(summarise_enrichment, "~/CASCADEpaper/paper/Fig3/Infercnv/summarise_enrichment_table.csv", quote = F, row.names = F, col.names = T, sep = ",")


tb$ID <- paste(tb$sample, tb$comp)
summarise_enrichment$ID <- paste(summarise_enrichment$sample, summarise_enrichment$subclone)


tb2 <- merge(tb, summarise_enrichment, by = "ID", all.x = T)
write.table(tb2, "~/CASCADEpaper/paper/Fig3/Infercnv/DEG_enrichment_table.csv", quote = F, row.names = F, col.names = T, sep = ",")


db1 <- as.data.frame(table(enrichment$ID, enrichment$sample)) %>% mutate(ID = paste(Var1, Var2))
db2 <- as.data.frame(table(Tenrich$ID, Tenrich$sample))%>% mutate(ID = paste(Var1, Var2))
db_tb <- merge(db1, db2, by = "ID", all = T)
db_tb$diff <- db_tb$Freq.x-db_tb$Freq.y
