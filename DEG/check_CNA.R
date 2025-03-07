library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)

geneCNAlist <- readRDS("~/CASCADEpaper/paper/Fig3_genomics/copy_number/geneCNAlist.Rds")
names(geneCNAlist)[8] <- "CA0035_paraaortic_lymph_node_1"
geneCNAlist <- geneCNAlist[-9]
names(geneCNAlist)[7] <- "CA0035_bladder_2"


purplereport <- readRDS("~/CASCADEpaper/paper/Fig3_genomics/purplereport.Rds")

ploidy <- lapply(purplereport, function(x) x$ploidy)
ploidy <- ploidy[-9]
names(ploidy) <- names(geneCNAlist)

geneCNAlist <- lapply(seq_along(geneCNAlist), function(x){
  temp <- geneCNAlist[[x]]
  temp$adj_CN <- temp$SCN/ploidy[[x]]
  return(temp)
})
names(geneCNAlist) <- names(ploidy)

patients <- substr(names(geneCNAlist),1,6)
geneCNApp <- split(geneCNAlist, patients)
geneCNApp <- geneCNApp[c("CA0027","CA0034","CA0035", "CA0043", "CA0046", "CA0058", "CA0076" , "CA0083","CA0090")]


# inter patient -----------------------------------------------------------
clonal_CNApp <- lapply(geneCNApp, function(x){
  temp <- rbindlist(x, idcol = "sample")
  temp <- split(temp, temp$gene)
  genelist <- sapply(temp, function(y){
    if(nrow(y) == length(x) & var(y$CN) == 0){
      return(unique(y$CN))
    }
  })
  return(unlist(genelist))
})

clonalgenes <- sapply(clonal_CNApp, names)

all_CNApp <-lapply(geneCNApp, function(x){
  temp <- rbindlist(x, idcol = "sample")
  temp <- split(temp, temp$gene)
  genelist <- sapply(temp, function(y){
    if(nrow(y) == length(x)){
      return(mean(y$CN))
    }
  })
  return(unlist(genelist))
})

inter_patient_DEG <- readRDS("~/integration/2024_06/DEG/inter_patient_DEG.Rds")
DEG <- inter_patient_DEG %>% filter(p_val_adj < 0.05)

pdf("interpatient_comp.pdf", width = 8, height = 4)
comp <- lapply(names(geneCNApp), function(i){
  DEGp <- DEG %>% filter(cluster == i )
  genes <- DEGp$gene
  CNp <- sapply(genes, function(g) mean(sapply(geneCNApp[[i]], function(y) y$CN[y$gene == g])))
  CN.na <- CNp[is.na(CNp)]
  CNp <- CNp[!is.na(CNp)]
  CNs <- lapply(genes, function(g) sapply(geneCNApp[setdiff(names(geneCNApp),i)], sapply, function(y) y$CN[y$gene == g]))
  CNs <- sapply(CNs, function(y) mean(unlist(y), na.rm = T))
  CNs <- setNames(CNs, genes)
  CN_ratio <- CNp/CNs[names(CNp)]
  df <- data.frame(patient = i , NOdeg = length(genes), degnoCN = length(CN.na), CN_ratio = CN_ratio, avglogFC = DEGp$avg_log2FC[match(names(CNp), DEGp$gene)])
  df$sign <- ifelse(sign(df$avglogFC) == 1, "pos", "neg")
  df$CN <- ifelse(df$CN_ratio >1, ">1", "<1")
  df$cor <- round(cor(df$CN_ratio, df$avglogFC), 3)
  df$id <- paste(df$sign, df$CN)
  df$clonal <- ifelse(names(CNp) %in% clonalgenes[[i]], "clonal", "subclonal")
  df$CN_clonal <- paste(df$CN, df$clonal)
  sumtb <-  df %>% group_by(sign, CN) %>% summarise(n = n()) %>% mutate(group_prop = round(n/sum(n),2), id = paste(sign, CN), pos =cumsum(group_prop) - (0.5 * group_prop))
  sumtb$prop <- round(sumtb$n/sum(sumtb$n), 2)
  sumtb$CN_ratio = ifelse(sumtb$CN == "<1", 0.5, 1.5)
  sumtb$avglogFC = ifelse(sumtb$sign == "neg", 2.5, -2.5)
  sumtb$sign <- factor(sumtb$sign)
  sumtb$CN <- factor(sumtb$CN, levels = c(">1", "<1"))
  f <- ggplot(df, aes(x = CN_ratio, y = avglogFC)) + geom_point(aes(color = CN_clonal)) + 
    geom_text(data = sumtb, aes(label = prop))+
    scale_color_manual(values = c( "tomato","royalblue","darkred", "navy" ), breaks = c(">1 subclonal", "<1 subclonal",">1 clonal", "<1 clonal" ))+ggtitle(paste(i, "correlation =", round(cor(df$CN_ratio, df$avglogFC), 3))) + theme_bw()

  p <- ggplot(sumtb, aes(x = sign)) + geom_bar(aes(fill= CN, y = group_prop), stat = "identity") + 
    geom_text(aes(label = n, y = pos),position = position_stack(vjust = 0.5))+
    scale_fill_manual(values = c( "royalblue","tomato"), breaks = c("<1", ">1"))+ theme_bw()
  print(f+p)
  return(df)
})
dev.off()
saveRDS(comp, "inter_patient_comp.Rds")



# inter lesion ------------------------------------------------------------
pdf("interlesion_comp.pdf", width = 8, height = 4)
comps <- list()
for( p in names(geneCNApp)){
  DEG <- readRDS(paste0("~/integration/2024_06/DEG/inter_lesion/",p, ".Rds"))
  DEG <- DEG[DEG$cluster %in% names(geneCNAlist),]
  DEG <- DEG %>% filter(p_val_adj < 0.05)
  comp <- lapply(unique(DEG$cluster), function(i){
    DEGp <- DEG %>% filter(cluster == i )
    genes <- DEGp$gene
    CNp <- unlist(sapply(genes, function(g) geneCNAlist[[i]]$CN[geneCNAlist[[i]]$gene == g]))
    CN.na <- CNp[is.na(CNp)]
    CNp <- CNp[!is.na(CNp)]
    CNs <- lapply(genes, function(g) sapply(geneCNApp[[substr(i, 1, 6)]][setdiff(names(geneCNApp[[substr(i, 1, 6)]]),i)],function(y) y$CN[y$gene == g]))
    CNs <- sapply(CNs, function(y) mean(unlist(y), na.rm = T))
    CNs <- setNames(CNs, genes)
    CN_ratio <- CNp/CNs[names(CNp)]
    df <- data.frame(patient = i , NOdeg = length(genes), degnoCN = length(CN.na), CN_ratio = CN_ratio, avglogFC = DEGp$avg_log2FC[match(names(CNp), DEGp$gene)])
    df$sign <- ifelse(sign(df$avglogFC) == 1, "pos", "neg")
    df$CN <- ifelse(df$CN_ratio >1, ">1", "<1")
    df$cor <- round(cor(df$CN_ratio, df$avglogFC), 3)
    df$id <- paste(df$sign, df$CN)

    sumtb <-  df %>% group_by(sign, CN) %>% summarise(n = n()) %>% mutate(group_prop = round(n/sum(n),2), id = paste(sign, CN), pos =cumsum(group_prop) - (0.5 * group_prop))
    sumtb$prop <- round(sumtb$n/sum(sumtb$n), 2)
    sumtb$CN_ratio = ifelse(sumtb$CN == "<1", 0.5, 1.5)
    sumtb$avglogFC = ifelse(sumtb$sign == "neg", 2.5, -2.5)
    sumtb$sign <- factor(sumtb$sign)
    sumtb$CN <- factor(sumtb$CN, levels = c(">1", "<1"))
    f <- ggplot(df, aes(x = CN_ratio, y = avglogFC)) + geom_point(aes(color = CN)) + 
      geom_text(data = sumtb, aes(label = prop))+
      scale_color_manual(values = c( "royalblue","tomato"), breaks = c("<1", ">1"))+ggtitle(paste(i, "correlation =", round(cor(df$CN_ratio, df$avglogFC), 3))) + theme_bw()
    
    g <- ggplot(sumtb, aes(x = sign)) + geom_bar(aes(fill= CN, y = group_prop), stat = "identity") + 
      geom_text(aes(label = n, y = pos),position = position_stack(vjust = 0.5))+
      scale_fill_manual(values = c( "royalblue","tomato"), breaks = c("<1", ">1"))+ theme_bw()
    print(f+g)
    return(df)
  })
  names(comp) <- unique(DEG$cluster)
  comps[[p]] <- comp
  print(p)
}
dev.off()

saveRDS(comps, "inter_interlesion_comp.Rds")

# enrichment of CNA agreed DEGs
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

C4 <- msigdbr(species = "Homo sapiens", category = "C4") %>%
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
      return(enresult@result)}
    else{
      print("no enrichment")
    }
  }else{print("not enough gene mapped")}
}

generate_genelist <- function(markerdf){
  markerdf <- markerdf %>% filter(id %in% c("pos >1", "neg <1"))
  temp <- setNames(markerdf$avglogFC, rownames(markerdf))
  temp <- sort(temp, decreasing = T)
  return(temp)
}


comp <- readRDS("~/integration/2024_06/DEG/inter_interlesion_comp.Rds")

genelists <- lapply(comp, lapply, generate_genelist)
pdf("enrichment of CN influenced DEGs of C6 Oncology.pdf")
enlist <- lapply(genelists, lapply, function(x) enrichment_barplots(x, C6, "oncology"))
dev.off()

pdf("enrichment of CN influenced DEGs of C4 computional gene sets.pdf")
enlist <- lapply(genelists, lapply, function(x) enrichment_barplots(x, C4, "C4"))
dev.off()


###### correlation between CN and expression ---------
cnlist <- readRDS("~/CASCADEpaper/paper/Fig3_genomics/copy_number/cnlist.Rds")
tumor_only_meta <- readRDS("~/CASCADEpaper/paper/Fig5_archetype/5_integrated_final_module_meta.Rds")

names(cnlist)[8:9] <- c("CA0035_bladder_2","CA0035_paraaortic_lymph_node_1" )


tumor_only_logPAC <- readRDS("~/integration/2024_06/tumor_only_logPAC.Rds")
library(GenomicRanges)
library(ggplot2)
gf = rtracklayer::import("/trigos_team/CASCADE/Analysis/reference_annotation/multiome/genes.gtf")
multi_g <- gf[gf$type =="gene"]

library(RColorBrewer)
n <- 24
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

purplereport <- readRDS("~/CASCADEpaper/paper/Fig3_genomics/purplereport.Rds")

ploidy <- lapply(purplereport, function(x) x$ploidy)
names(ploidy) <- names(cnlist)

pdf("sample_CN_expression_regression.pdf", width = 5, height = 5)
 tiled_cn_exp <- lapply(unique(tumor_only_meta$sample), function(i){
  cn <- as(as.data.frame(cnlist[[i]]), "GRanges")
  seql <- split(end(cn), seqnames(cn))
  seql <- sapply(seql, max)
  tg <- tileGenome(seql, tilewidth = 50000000, cut.last.tile.in.chrom = T)
  ind <- findOverlaps(cn, tg)
  ind <- split(from(ind), to(ind))
  agg_cn <- sapply(ind, function(x) mean(cn$copyNumber[x]))
  agg_cn <- agg_cn/ploidy[[i]]*2
  tg$cn <- agg_cn
  log_cn <- log(agg_cn + 1)
  tg$log_cn <- log_cn
  
  hits <- findOverlaps(multi_g, tg, select = "all", minoverlap = 1)
  ind <- split(from(hits), to(hits))
  gene_tg <- sapply(ind, function(x) unique(multi_g$gene_name[x]))
  gene_avg <- rowMeans(as.matrix(tumor_only_logPAC[, grep(i, colnames(tumor_only_logPAC))]))
  agg_gene <- sapply(gene_tg, function(x) mean(gene_avg[intersect(x, names(gene_avg))]))
  tg$agg_gene <- NA
  tg$agg_gene[as.numeric(names(agg_gene))] <- agg_gene
  
  df <- as.data.frame(tg)
  lms <- lm(agg_gene~cn, df[df$cn<=5,])
  jt <- DescTools::JonckheereTerpstraTest(df[df$cn<=5,"cn"], g = df[df$cn<=5,"agg_gene"], alternative = "increasing", nperm = 5000)
  cor <- cor.test(df$cn, df$agg_gene, alternative ="greater", method = "spearman")
  tg$intercept <- round(lms$coefficients["(Intercept)"], 3)
  tg$slope <- round(lms$coefficients["cn"], 3)
  tg$sample <- i
  tg$patient <- substr(i, 1, 6)
  tg$jt.p <- round(jt$p.value, 3)
  tg$rho <- round(cor$estimate, 3)
  f <- ggplot(df, aes(x = cn, y = agg_gene)) + geom_point(aes(color = seqnames))+ 
    scale_color_manual(values = col_vector)+
    stat_smooth(method = "lm", 
                formula = y ~ x, 
                geom = "smooth") + 
    xlim(0, 5)+
    annotate("text", x = 4, y = 0.075, label = paste("rho:", round(cor$estimate, 3), "p:", round(cor$p.value, 3)))+
    ggtitle(i, subtitle = paste("intercept:", round(lms$coefficients["(Intercept)"], 3), 
                                "slope:", round(lms$coefficients["cn"], 3),
                                "JT pval:", round(jt$p.value, 3))) + theme_bw()
  print(f)
  print(i)
  return(tg)
 })
dev.off()
saveRDS(tiled_cn_exp, "tiled_cn_exp.Rds")


cordf <- lapply(unique(tumor_only_meta$sample), function(i){
  cn <- as(as.data.frame(cnlist[[i]]), "GRanges")
  seql <- split(end(cn), seqnames(cn))
  seql <- sapply(seql, max)
  tg <- tileGenome(seql, tilewidth = 50000000, cut.last.tile.in.chrom = T)
  ind <- findOverlaps(cn, tg)
  ind <- split(from(ind), to(ind))
  agg_cn <- sapply(ind, function(x) mean(cn$copyNumber[x]))
  agg_cn <- agg_cn/ploidy[[i]]*2
  tg$cn <- agg_cn
  log_cn <- log(agg_cn + 1)
  tg$log_cn <- log_cn
  
  hits <- findOverlaps(multi_g, tg, select = "all", minoverlap = 1)
  ind <- split(from(hits), to(hits))
  gene_tg <- sapply(ind, function(x) unique(multi_g$gene_name[x]))
  gene_avg <- rowMeans(as.matrix(tumor_only_logPAC[, grep(i, colnames(tumor_only_logPAC))]))
  agg_gene <- sapply(gene_tg, function(x) mean(gene_avg[intersect(x, names(gene_avg))]))
  tg$agg_gene <- NA
  tg$agg_gene[as.numeric(names(agg_gene))] <- agg_gene
  
  df <- as.data.frame(tg)
  cor <- cor.test(df$cn, df$agg_gene, alternative ="greater", method = "spearman")
  return(data.frame(rho = cor$estimate, p = cor$p.value, sample = i))
})

df <- rbindlist(cordf)
df$patient <- substr(df$sample, 1, 6)
df$patient <- gsub("00", "",df$patient)
df$adj_p <- p.adjust(df$p, method = "bonferroni")
df$p_judge <- ifelse(df$adj_p < 0.05, "significant", "non significant")

xlab <- setNames(paste(tumor_only_meta$site, sapply(strsplit(tumor_only_meta$sample, split = "_"), function(x) tail(x, 1))), tumor_only_meta$sample)
xlab <- xlab[unique(tumor_only_meta$sample)]
saveRDS(df, "expression_CN_correlation.Rds")
df$patient <- gsub("CA", "", df$patient)

ggplot(df, aes(x = sample, y = rho, fill= p_judge)) + geom_bar(stat = "identity") + 
  theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        strip.background = element_rect(fill = "white"), 
        panel.spacing = unit(0.1, "cm"), 
        legend.direction = "horizontal", 
        legend.position = "top", 
        legend.margin = margin(0, 0, -10,0)) +
  facet_grid(.~patient, space = "free", scales = "free")+
  scale_fill_manual(values = c("red",  "lightseagreen"), name = "p-value")+
  scale_x_discrete(labels = xlab)+
  labs(x = "Sample")
ggsave("expression_CN_correlation_bar.pdf", width = 7.3, height = 4)

