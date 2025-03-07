require(ggplot2)

library(reshape2)
library(data.table)
library(scales)
library(ggpubr)
library(dplyr)
library(phylogram)
library(dendextend)
library(patchwork)
library(ComplexHeatmap)
library(pals)
library(circlize)

# PURPLE visual with gene segmentation-------------
genelist <- readRDS("~/CASCADEpaper/paper/Fig3_genomics/copy_number/geneCNAlist.Rds")
cnlist <- readRDS("~/CASCADEpaper/paper/Fig3_genomics/copy_number/cnlist.Rds")

paths <- system("realpath /trigos_team/CASCADE/Sequencing_updated_reference/CA*/*/purple/*purity.tsv", intern = T)
paths_split <- data.table::transpose(strsplit(paths, "/"))
ids <- paths_split[[5]]
sites <- paths_split[[6]]
purplereport <- lapply(paths, read_tsv)
names(purplereport) <- paste(ids, sites, sep= "_")

### use gene based CN -> too refined, not clear
sapply(seq_along(genelist), function(x){
  ploidy <- purplereport[[x]]$ploidy
  cndf <- genelist[[x]]
  cndf <- cndf[order(cndf$chrom, cndf$g_start),]
  cndf$chrom <- gsub("chr", "", cndf$chrom)
  cndf$chrom <- factor(cndf$chrom, levels = c(as.character(1:22), "X"))
  mat <- matrix(cndf$SCN, nrow = 1)
  mat <- mat/ploidy*2
  col_fun <- colorRamp2(c(0, 2, 4), c("steelblue3", "white", "red3"))
  ht  <-  Heatmap(mat, border = T, cluster_columns = F, column_split = cndf$chrom, col = col_fun, name = "normalized\nCN", cluster_column_slices = F,
                  heatmap_legend_param = list(labels_gp = gpar(fontsize = 12), title_gp = gpar(fontsize = 16)),
                  column_title_gp = gpar(fontsize= 12))
  pdf(paste0("~/CASCADEpaper/paper/Fig3_genomics/copy_number/purple_fig/", names(genelist)[x], ".pdf"), width = 10 , height= 2.7)
  draw(ht, column_title = "WGS PURPLE output",column_title_gp = gpar(fontsize = 16, fontface = "bold"))
  dev.off()
})

### use original segments (too many bins from some chromosomes)
sapply(seq_along(genelist), function(x){
  ploidy <- purplereport[[x]]$ploidy
  cndf <- cnlist[[x]]
  cndf <- cndf[order(cndf$chromosome, cndf$start),]
  cndf$chromosome <- gsub("chr", "", cndf$chromosome)
  cndf$chromosome <- factor(cndf$chromosome, levels = c(as.character(1:22), "X", "Y"))
  mat <- matrix(cndf$copyNumber, nrow = 1)
  mat <- mat/ploidy*2
  col_fun <- colorRamp2(c(0, 2, 4), c("steelblue3", "white", "red3"))
  ht  <-  Heatmap(mat, border = T, cluster_columns = F, column_split = cndf$chromosome, col = col_fun, name = "normalized\nCN", cluster_column_slices = F,
                  heatmap_legend_param = list(labels_gp = gpar(fontsize = 12), title_gp = gpar(fontsize = 16)),
                  column_title_gp = gpar(fontsize= 12))
  pdf(paste0("~/CASCADEpaper/paper/Fig3_genomics/copy_number/purple_fig/", names(genelist)[x], ".pdf"), width = 10 , height= 2.7)
  draw(ht, column_title = "WGS PURPLE output",column_title_gp = gpar(fontsize = 16, fontface = "bold"))
  dev.off()
})



# PURPLE viual use ATAClone bins size (in paper)-------------
library(GenomicRanges)
library(IRanges)
cndf <- cnlist$CA0027_dura_base_skull_13
gr <- as(cndf, "GRanges")

joint_cn <- readRDS("/trigos_team/CASCADE/Analysis/240628_ATAClone/CA0027/dura_base_skull_13/joint_cn.Rds")
chrom <- sapply(strsplit(rownames(joint_cn), split = ".", fixed = T), "[", 1)
ranges<- sapply(strsplit(rownames(joint_cn), split = ":", fixed = T), "[", 2)
start <- sapply(strsplit(ranges, split = "-", fixed = T), "[", 1)
end <- sapply(strsplit(ranges, split = "-", fixed = T), "[", 2)

bingr <- as(data.frame(chrom = chrom, start = start, end = end), "GRanges")

cols <- factor(gsub("chr", "", chrom), levels = c(as.character(1:22),"X", "Y"))
col_fun <- colorRamp2(c(0, 2, 4), c("steelblue3", "white", "red3"))
sapply(seq_along(genelist), function(x){
ploidy <- purplereport[[x]]$ploidy
cndf <- cnlist[[x]]
gr <- as(cndf, "GRanges")
hit <- findOverlaps(bingr, gr)
ind <- split( to(hit), from(hit))
cns <- sapply(ind, function(x) mean(cndf$copyNumber[x]))
mat <- matrix(cns, nrow = 1)
mat <- mat/ploidy*2
ht  <-  Heatmap(mat, border = T, cluster_columns = F, column_split = cols, col = col_fun, name = "CN", cluster_column_slices = F,
                heatmap_legend_param = list(labels_gp = gpar(fontsize = 13), title_gp = gpar(fontsize = 14)),
                column_title_gp = gpar(fontsize= 13))
pdf(paste0("purple_fig/", names(genelist)[x], ".pdf"), width = 7.5 , height= 2.5)
draw(ht, column_title = "WGS output",column_title_gp = gpar(fontsize = 16, fontface = "bold"))
dev.off()

})

# subdf <- lapply(names(genelist), function(sample) {
#   x <- genelist[[sample]]
#   temp <- x[, c("gene", "SCN")]
#   colnames(temp) <- c("gene", sample)
#   return(temp)
#   })
# jointdf <- Reduce(function(x, y) merge(x, y, all = T, by= "gene"), subdf)
# 
# jointdf <- jointdf[order(jointdf$chrom.x, jointdf$g_start.x),]
# jointmat <- t(as.matrix(jointdf[, grep("CA00", colnames(jointdf))]))
# 
# Heatmap(jointmat, border = T, cluster_columns = F,show_column_names = F,  column_split = jointdf$chrom.x, col = col_fun)
# 

