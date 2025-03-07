require(ggplot2)
library(infercnv)
library(reshape2)
library(data.table)
library(scales)
library(ggpubr)
library(dplyr)
library(phylogram)
library(dendextend)
library(patchwork)
genelist <- readRDS("~/CASCADEpaper/paper/Fig3/copy_number/geneCNAlist.Rds")
genelist <- lapply(genelist, function(x){
  x$CN <- x$SCN
  x <- x[, -c(7, 11, 12)]
})

genelist <- genelist[-9]
names(genelist)[7] <- "CA0035_bladder_2" 
names(genelist)[8] <- "CA0035_paraaortic_lymph_node_1"


for (i in names(genelist)){
  df <- genelist[[i]]
  df_list <- list()
  for (n in unique(df$chrom)){
    subdf <- data.table(subset(df, df$chrom == n), key = "g_start")
    subdf$bin_width <- subdf$g_end-subdf$g_start
    
    #delete duplicated start at the same position with smaller bin width
    x <- subdf$g_start[duplicated(subdf$g_start)]
    if(length(x)>0){
      temp <- aggregate(bin_width~g_start, subdf[subdf$g_start%in%x,], min)
      subdf <-  subdf[-which(subdf$g_start%in%x&subdf$bin_width%in%temp$bin_width),]
    }
    
    l <- nrow(subdf)
    y <- 0
    while (y < l){
      l <- nrow(subdf)
      subdf <- subdf[order(subdf$g_start), ]
      last_end <- subdf$g_end[1:(nrow(subdf)-1)]
      subdf$last_end <- c(NA, last_end)
      next_end <- subdf$g_end[2:nrow(subdf)]
      subdf$next_end <- c(next_end,NA)
      ind3 <- which(subdf$next_end<=subdf$g_end)
      if (length(ind3)>0){
        subdf <- subdf[-(ind3+1), ]  
      }
      y <- nrow(subdf)
    }
    #select genes that overlop part of them, use the earliest start and lastest end
    y <- 0
    while (y < l){
      l <- nrow(subdf)
      subdf <- subdf[order(subdf$g_start), ]
      subdf$last_end[2:nrow(subdf)] <- subdf$g_end[1:(nrow(subdf)-1)]
      next_end <- subdf$g_end[2:nrow(subdf)]
      subdf$next_end <- c(next_end,NA)
      ind <- which(subdf$g_start<subdf$last_end&subdf$g_end>subdf$last_end)
      if(length(ind)>0){
        ind2 <- ind-1
        temp <- data.frame(gene = subdf$gene[ind2], chrom = subdf$chrom[ind2], g_start = subdf$g_start[ind2], g_end = subdf$g_end[ind], 
                           f_start = NA, f_end = NA, CN = subdf$CN[ind], p_gene = NA, region = NA, 
                           bin_width = subdf$g_end[ind]-subdf$g_start[ind2], last_end = NA, next_end= NA)
        subdf <- subdf[-unique(c(ind, ind2)),]
        subdf <- rbind(subdf, temp)  
      }
      y <- nrow(subdf)
    }
    
    subdf$seg <- seq(1, nrow(subdf))
    subdf$samples <- "sample"
    df_list[[n]] <- subdf[order(subdf$g_start), ]
    print(n)
  }
  genelist[[i]]<- do.call(rbind.data.frame, df_list)
  print(i)
}

for (i in names(genelist)[names(genelist)%in%list.files("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/clustered")][c(7, 8, 27:32)]){
  #for (i in grep("CA0058", names(genelist)[names(genelist)%in%list.files("~/CN/pooling_normal")], value = T)){
  print(i)
  temp <- readRDS(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/clustered/", i, "/preliminary.infercnv_obj"))
  anno <- readRDS(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/anno/", i, ".Rds"))
  colnames(anno) <- c("ID")
  anno$Cell <- rownames(anno)
  anno_t <- anno[nchar(anno$ID)<3, ]
  mat <- temp@expr.data
  mat <- mat[, intersect(anno_t$Cell, colnames(mat))]
  
  
  cell_ord_uncluster <- read.dendrogram(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/", i, "/infercnv.preliminary.observations_dendrogram.txt"))
  cell_ord <-list()
  hc <- temp@tumor_subclusters$hc
  
  cell_ord <- lapply(unique(anno_t$ID), function(n){
    test <- data.frame(cell = hc[[n]]$labels, order =hc[[n]]$order, cluster = n)
    return(test)
  })
  names(cell_ord) <- unique(anno_t$ID)
  
  cell_ord <- rbindlist(cell_ord)
  cell_ord <- data.table(cell_ord, key = "cell")
  gene_ord<- temp@gene_order
  gene_ord$width <- gene_ord$stop-gene_ord$start
  gene_ord$gene <- rownames(gene_ord)
  gene_ord <- data.table(gene_ord, key = "gene")
  df_list <- list()
  for (n in unique(gene_ord$chr)){
    subdf <- data.table(subset(gene_ord, gene_ord$chr == n), key = "start")
    x <- subdf$start[duplicated(subdf$start)]
    if(length(x)>0){
      temp <- aggregate(width~start, subdf[subdf$start%in%x,], min)
      subdf <-  subdf[-which(subdf$start%in%x&subdf$width%in%temp$width),]
    }
    
    l <- nrow(subdf)
    y <- 0
    while (y < l){
      l <- nrow(subdf)
      subdf <- subdf[order(subdf$start), ]
      last_end <- subdf$stop[1:(nrow(subdf)-1)]
      subdf$last_end <- c(NA, last_end)
      next_end <- subdf$stop[2:nrow(subdf)]
      subdf$next_end <- c(next_end,NA)
      ind3 <- which(subdf$next_end<=subdf$stop)
      if (length(ind3>0)){
        subdf <- subdf[-(ind3+1), ]  
      }
      y <- nrow(subdf)
    }
    
    y <- 0
    while (y < l){
      l <- nrow(subdf)
      subdf <- subdf[order(subdf$start), ]
      subdf$last_end[2:nrow(subdf)] <- subdf$stop[1:(nrow(subdf)-1)]
      next_end <- subdf$stop[2:nrow(subdf)]
      subdf$next_end <- c(next_end,NA)
      ind <- which(subdf$start<subdf$last_end&subdf$stop>subdf$last_end)
      if(length(ind)>0){
        ind2 <- ind-1
        temp <- data.frame(gene = subdf$gene[ind2], chr = subdf$chr[ind2], start = subdf$start[ind2], stop = subdf$stop[ind], 
                           width = subdf$stop[ind]-subdf$start[ind2], last_end = NA, next_end = NA)
        subdf <- subdf[-unique(c(ind, ind2)),]
        subdf <- rbind(subdf, temp)  
      }
      y <- nrow(subdf)
    }
    df_list[[n]] <- subdf[order(subdf$start)]
  }
  gene_ord<- rbindlist(df_list)
  gene_ord_f <- gene_ord[,c("gene",'chr', "width", "start")]
  gene_ord_f <- gene_ord_f[gene_ord_f$gene%in%genelist[[i]]$gene, ]
  mat <- mat[rownames(mat)%in%gene_ord_f$gene,]
  
  
  
  df <- reshape2::melt(mat)
  colnames(df) <- c("gene", "cell", "value")
  df <- data.table(df, key = "gene")
  df <- merge(df, gene_ord_f, by="gene", all = T)
  df <- merge(df, cell_ord, by  = "cell",all = T)
  df$cell <- factor(df$cell, levels = unique(df$cell[order(df$cluster, df$order)]))
  df$gene <- factor(df$gene, levels = unique(df$gene[order(df$start)]))
  
  
  
  p<- ggplot(data = df, aes(x = gene, y = cell, fill = value, width = log10(width)))+ 
    geom_tile()+
    # scale_fill_gradientn(na.value = "#D3D3D3",
    # colors=c("#00008B","blue","white","#DC143C", "black"),
    # values=rescale(c(min(df$value),quantile(df$value[df$value<1], 0.1),1,quantile(df$value[df$value>1], 0.9),max(df$value[!is.na(df$value)]))),
    # limits=c(min(df$value[!is.na(df$value)]), max(df$value[!is.na(df$value)])))+
    scale_fill_gradientn(na.value = "transparent",
                         colors=c("#00009D","white","#9D0002"),
                         values=rescale(c(0.85, 1, 1.15)),
                         limits=c(0.85, 1.15), 
                         breaks = c(0.85, 0.95, 1.05, 1.15))+
    facet_grid(cluster~chr, scales = "free")+
    theme(strip.text.x = element_text(size = 10), axis.title = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(), 
          axis.text.y =element_blank(), 
          panel.spacing = unit(0.1, "lines"))+ ggtitle(paste("Infercnv plots:", i, "grouped by clusters of tumor only srt object"))

  
  df$cell <- factor(df$cell, levels = labels(cell_ord_uncluster))
  p2<- ggplot(data = df, aes(x = gene, y = cell, fill = value, width = log10(width)))+ 
    geom_tile()+
    scale_fill_gradientn(na.value = "transparent",
                         colors=c("#00009D","white","#9D0002"),
                         values=rescale(c(0.85, 1, 1.15)),
                         limits=c(0.85, 1.15), 
                         breaks = c(0.85, 0.95, 1.05, 1.15))+
    facet_grid(.~chr, scales = "free_x")+
    theme(strip.text.x = element_text(size = 10), axis.title = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(), 
          axis.text.y =element_blank(), 
          panel.spacing = unit(0.1, "lines"))+ ggtitle(paste("Infercnv plots:",i, "without grouping"))

  
  genelist[[i]] <- genelist[[i]][genelist[[i]]$gene%in%gene_ord_f$gene,]
  # ind <- max(genelist[[i]]$CN[!is.na(genelist[[i]]$CN)])
  # ind2 <- ind-2
  # colors <- c("blue","#0000FF4D", "#FFFFFF", tail(colorspace::diverge_hsv(ind2*2+1), ind2))
  #adjustcolor("blue", alpha.f = 0.3)
  #show_col()
  patient <- substr(i, 1, 6)
  sample <- sub(paste0(patient, "_"), "", i, fixed = T)
  if( length(list.files(paste0("/trigos_team/CASCADE/Sequencing_updated_reference/", patient,"/",sample,"/purple"),pattern =  "purity.tsv", full.names = T) ) == 0) { midpoint = 2}else{
    midpoint <-list.files(paste0("/trigos_team/CASCADE/Sequencing_updated_reference/", patient,"/",sample,"/purple"),pattern =  ".purple.purity.tsv", full.names = T) %>% read.table(sep = "\t", header = T) %>% dplyr::select(ploidy) %>% as.numeric
  }
  g <- ggplot(data =genelist[[i]],aes(x = seg, y = samples, fill = CN, width = log10(bin_width)))+
    geom_tile()+
    facet_grid(.~chrom, scale = "free_x")+
    theme(strip.text.x = element_text(size = 10), axis.title = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(), 
          axis.text.y =element_blank(), 
          panel.spacing = unit(0.1, "lines"))+
    scale_fill_gradientn(na.value = "transparent",
                         colors=c("#00009D","white","#9D0002"),
                         values=rescale(c(min(genelist[[i]]$CN),midpoint, max(genelist[[i]]$CN))),
                         limits=c(min(genelist[[i]]$CN), max(genelist[[i]]$CN)), 
                         breaks = c(min(genelist[[i]]$CN),midpoint, max(genelist[[i]]$CN)))+ ggtitle(paste("PURPLE plot:",i))
  
  #scale_fill_manual(values = colors, na.value = "#D3D3D3", breaks = as.character(seq(0,ind)),name = "CN")+ggtitle(paste(i, "gene CN"))
  

  pdf(paste0("~/CASCADEpaper/paper/Fig3/copy_number/svg/", i, ".png"), width = 12, height = 10)
  # print(g)
  print(p + p2 + g + plot_layout(nrow = 3))
  dev.off()
}