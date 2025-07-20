rmaxes <- function(g, front.size = 18, print = T){
  if(class(g)[1] == "patchwork"){
    f1 <- g
    for (i in 1:length(g)){
      f1[[i]] <- rmaxes(f1[[i]], print = F, front.size = front.size)
    }
    # f1 <- lapply(f1, rmaxes, print = F)
    f1 <- f1+plot_layout(axes = "collect")
    if(print){
      print(f1)
    }
    return(f1)
  }else{
  temp <- g + theme(axis.title = element_blank(), text = element_text(size = front.size),
                    plot.margin = margin(t = 0, l = -0.3, r = 0, b = -0.3, "cm"),
                    legend.key.size = unit(0.2,"inches"))+
    scale_x_continuous(breaks = NULL, limits = c(min(g$data[, 1]-1), max(g$data[, 1])))+
    scale_y_continuous(breaks = NULL, limits = c(min(g$data[, 2]-1), max(g$data[, 2])))
  if(print){
    print(temp)
  }
  return(temp)
  }
}

coneraxes <- function(p, labelx = "UMAP1", labely ="UMAP2", colorbar = TRUE, front.size = 18, print = T){
  if(colorbar){
    arr <- list(x = min(p$data[, 1])-1, y = min(p$data[, 2])-1, 
                x_len = diff(range(p$data[, 1]))/5, y_len = diff(range(p$data[, 1]))/5)
    f <- p + annotate("segment", 
                      x = arr$x, xend = arr$x + c(arr$x_len, 0), 
                      y = arr$y, yend = arr$y + c(0, arr$y_len), 
                      arrow = arrow(type = "closed", length = unit(3, 'pt')))+
      annotate(geom = "text", x = arr$x, y = arr$y-0.3, label = labelx, hjust = "left", vjust = "top", size = front.size/.pt)+
      annotate(geom = "text", x = arr$x-0.3, y = arr$y, label = labely, angle = 90, hjust = "left", vjust = "bottom", size = front.size/.pt)+
      theme(axis.title = element_blank(), text = element_text(size = front.size),
            plot.margin = margin(t = 0, l = 0.1, r = 0, b = 0.1, "cm"),
            legend.key.size = unit(0.2,"inches"), 
            legend.margin = margin(l = -10, t = -10))+
      scale_x_continuous(breaks = NULL) +
      scale_y_continuous(breaks = NULL)+
      coord_cartesian(clip = 'off')
    if(print){
    print(f)
    }
    return(f)
  }
}


#### normal cell plot functions ----------
source("~/CASCADEpaper/paper/cols.R")
fplot <- function(srt,features, title = NULL, reduction = c("css_umap", "harmony_umap")[2]){
  DefaultAssay(srt) <- "ruv3"
  if(class(features) == "list"){
    sapply(names(features), function(x){
      FeaturePlot(srt, features = features[[x]], reduction = reduction) + plot_annotation(title = x) + plot_layout(guides = "collect", axes = "collect")
    })
  }else{
    features <- features[features %in% rownames(srt)]
    fs <- split(features, ceiling(seq_along(features)/6))
    lapply(fs, function(x) FeaturePlot(srt, features = x, reduction = reduction) + plot_annotation(title = title) + plot_layout(nrow = 2, ncol = 3, guides = "collect", axes = "collect"))
    
  }
  
}

barplot_bytissue <- function(endo3){
  df <- endo3@meta.data %>% group_by(subtype, site) %>% summarise(n = n()) %>% mutate (freq = n/sum(n))
  df$subtype <- factor(df$subtype, levels = sort(unique(endo3$subtype)))
  df$xlab <- str_wrap(df$subtype, 10)
  df2 <- df %>% group_by(subtype) %>% mutate(n = sum(n))
  
  df$labs <- paste0(signif(df$freq*100, 2))
  df$labs[df$freq < 0.03] <- NA
  
  ggplot(df, aes(x = xlab, y = freq, fill = site)) + geom_bar(stat = "identity")  +  theme_classic2(base_size = 26)+
    geom_text(aes(label = labs),
              position = position_stack(vjust = 0.5),        # Position text in the middle of the bars
              color = "black")+
    theme(axis.text.x = element_text(), 
          axis.title= element_blank()) + 
    scale_fill_manual(values =site_cols)+
    geom_text(data=df2,aes(x=xlab,y=1,label=n),hjust=0, angle = 90)+
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1.2), position = "right")+
    scale_x_discrete(guide = guide_axis(angle = -30), position = "top", labels)#guide_axis(n.dodge = 2)
}

barplot_bytissue_byrow <- function(endo3){
  df <- endo3@meta.data %>% group_by(site, subtype) %>% summarise(n = n()) %>% mutate (freq = n/sum(n))
  df$subtype <- factor(df$subtype, levels = sort(unique(endo3$subtype)))
  df2 <- df %>% group_by(site) %>% mutate(n = sum(n))
  
  df$labs <- paste0(signif(df$freq*100, 2))
  df$labs[df$freq < 0.03] <- NA
  
  ggplot(df, aes(x = site, y = freq, fill = subtype)) + geom_bar(stat = "identity")  +  theme_bw(base_size = 18)+
    geom_text(aes(label = labs),
              position = position_stack(vjust = 0.5),        # Position text in the middle of the bars
              color = "white")+
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), 
          axis.title.x = element_blank(), 
          plot.margin = margin(-10, -10, -10, 5)) + 
    scale_fill_manual(values =ggsci::pal_observable()(nlevels(df$subtype)))+theme(legend.position = "none")+
    geom_text(data=df2,aes(x=site,y=1,label=n),hjust=0, angle = 90)+
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1.2))+
    labs(y = "Cell Proportion")
  
}

barplot_bypatient <- function(endo3){
  df <- endo3@meta.data %>% group_by(patient, subtype) %>% summarise(n = n()) %>% mutate (freq = n/sum(n))
  df$subtype <- factor(df$subtype, levels = sort(unique(endo3$subtype)))
  df$patient <- gsub("00", "", df$patient)
  df2 <- df %>% group_by(patient) %>% mutate(n = sum(n))

  df$labs <- paste0(signif(df$freq*100, 2))
  df$labs[df$freq < 0.03] <- NA
  
  ggplot(df, aes(x = patient, y = freq, fill = subtype)) + geom_bar(stat = "identity")  +  theme_bw(base_size = 18)+
    geom_text(aes(label = labs),
              position = position_stack(vjust = 0.5), 
              color = "white")+
    scale_fill_manual(values =ggsci::pal_observable()(nlevels(df$subtype)))+
    theme(legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
          axis.title.x = element_blank(), 
          plot.margin = margin(-10, r = 5, -10, -10)) +
    geom_text(data=df2,aes(x=patient,y=1,label=n),hjust=0, angle = 90)+
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1.2))+
    labs(y = "Cell Proportion")
}

# barplot_bypatient <- function(endo3){
#   df <- endo3@meta.data %>% group_by(patient, subtype) %>% summarise(n = n()) %>% mutate (freq = n/sum(n))
#   df$subtype <- factor(df$subtype, levels = sort(unique(endo3$subtype)))
#   df$patient <- gsub("00", "", df$patient)
#   df2 <- df %>% group_by(patient) %>% mutate(n = sum(n))
#   
#   df$labs <- paste0(signif(df$freq*100, 2), "%")
#   df$labs[df$freq < 0.03] <- NA
#   
#   ggplot(df, aes(x = patient, y = freq, fill = subtype)) + geom_bar(stat = "identity")  +  theme_bw(base_size = 18)+
#     geom_text(aes(label = labs),
#               position = position_stack(vjust = 0.5),        # Position text in the middle of the bars
#               color = "white",                               # Set text color
#               size = 4)+
#     scale_fill_manual(values =ggsci::pal_observable()(nlevels(df$subtype)))+
#     theme(legend.position = "none",axis.text.x = element_text(angle = 90)) +
#     geom_text(data=df2,aes(x=patient,y=1,label=n),hjust=0, angle = 90)+
#     scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1.2))
# }

barplot_bypathology <- function(endo3){
  endo3$pathology <-ifelse(endo3$patient %in% c("CA0090", "CA0046"), "NE", ifelse(endo3$patient %in% c("CA0027", "CA0058"), "Mixed", "AD"))
  df <- endo3@meta.data %>% group_by(pathology, subtype) %>% summarise(n = n()) %>% mutate (freq = n/sum(n))
  df$subtype <- factor(df$subtype, levels = sort(unique(endo3$subtype)))
  df2 <- df %>% group_by(pathology) %>% mutate(n = sum(n))
  
  df$labs <- paste0(signif(df$freq*100, 2), "%")
  df$labs[df$freq < 0.03] <- NA
  
  ggplot(df, aes(x = pathology, y = freq, fill = subtype)) + geom_bar(stat = "identity")  +  theme_bw(base_size = 18)+
    geom_text(aes(label = labs),
              position = position_stack(vjust = 0.5),        # Position text in the middle of the bars
              color = "white",                               # Set text color
              size = 4)+
    scale_fill_manual(values =ggsci::pal_observable()(nlevels(df$subtype)))+
    geom_text(data=df2,aes(x=pathology,y=1,label=n),hjust=0, angle = 90)+
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1.2))+
    theme(legend.position = "none")
}


vlnplot <- function(srt, markers, group){
  exp <- srt[["ruv3"]]@data[intersect(markers, rownames(srt[["ruv3"]])), ]
  df <- reshape2::melt(as.matrix(exp))
  colnames(df) <- c("gene", "cell", "ruv_exp")
  df$patient <- srt$patient[match(df$cell, colnames(srt))]
  df$sample <- srt$sample[match(df$cell, colnames(srt))]
  df$site <- srt$site[match(df$cell, colnames(srt))]
  df$cluster_marker  <- srt$cluster_marker[match(df$gene, srt$marker)]
  df$group <- factor(as.character(group[as.character(df$gene)]), levels = rev(group[levels(df$gene)]))
  lab <- df %>% group_by(sample, site) %>% summarise(n = n())
  lab <- setNames(lab$site, lab$sample)
  # f <- ggplot(df, aes(x = patient, fill = patient, y = ruv_exp))+
  #   geom_violin(position=position_dodge(),draw_quantiles=c(0.5))+
  #   facet_grid(cluster_marker+gene~., switch = "y") + theme_bw(base_size = 18) + 
  #   scale_fill_manual(values = brewer.set1(9))+
  #   scale_y_continuous(position = "right") + 
  #   theme(strip.background =element_rect(fill="white"), strip.text = element_text(colour = 'black', size = 7))
  # # print(f)
  # f <- ggplot(df, aes(x = sample, fill = patient, y = ruv_exp))+
  #   geom_violin(position=position_dodge(),draw_quantiles=c(0.5))+
  #   scale_x_discrete(labels = lab)+
  #   facet_grid(cluster_marker+gene~patient, switch = "y", scales = "free_x", space = "free_x") + theme_bw(base_size = 18) + 
  #   scale_fill_manual(values = brewer.set1(9))+
  #   scale_y_continuous(position = "right") + 
  #   theme(strip.background =element_rect(fill="white"), strip.text = element_text(colour = 'black', size = 7),
  #         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  #print(f)
  df$facet <- paste(df$group, df$gene, sep = "\n")
  f <- ggplot(df, aes(x = site, fill = site, y = ruv_exp))+
    geom_violin(position=position_dodge(),draw_quantiles=c(0.5))+
    facet_grid(facet~., switch = "y", ) + theme_bw(base_size = 26) + 
    scale_fill_manual(values = brewer.set3(9))+
    scale_y_continuous(position = "right", breaks = c(0, round(max(df$ruv_exp)/2), round(max(df$ruv_exp)))) + 
    scale_x_discrete(guide = guide_axis(angle = 30))+
    theme(strip.background =element_rect(fill="white"), 
          strip.text = element_text(colour = 'black'), 
          strip.text.y.left = element_text(angle = 0),
          legend.position = "none", 
          axis.title = element_blank(), 
          panel.spacing = unit(0.2, "cm"))
  print(f)
  return(f)
}

plot_cluster <- function(Tcells){
  f <- DimPlot(Tcells, reduction = "harmony_umap", label = T, group.by = "harmony_cluster")
  print(f)
  Tcells <- FindClusters(Tcells, resolution = 0.5, graph.name = "harmony_snn")
  DefaultAssay(Tcells) <- "ruv3"
  markers <- FindAllMarkers(Tcells, test.use = "MAST" , only.pos = T)
  f <- DimPlot(Tcells, reduction = "harmony_umap", label = T)
  print(f)
  df <- markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% arrange(desc(avg_log2FC)) %>% dplyr::slice(1:3)
  print(df, n = 50)
  marker <- df %>% group_by(cluster) %>% arrange(desc(avg_log2FC)) %>% dplyr::slice(1)
  Tcells$cluster_marker <- paste(Tcells$seurat_clusters, marker$gene[match(Tcells$seurat_clusters, marker$cluster)], sep = ":")
  Tcells$marker <- marker$gene[match(Tcells$seurat_clusters, marker$cluster)]
  f <- DimPlot(Tcells, reduction = "harmony_umap", label = T, group.by = "cluster_marker", repel = T)
  print(f)
  return(Tcells)
}


plot_azimuth <- function(azimuth_pred, immunes, reduction = c("css_umap", "harmony_umap")[2]){
  predictions <- read.delim(azimuth_pred, row.names = 1)
  immunes <- AddMetaData(
    object = immunes,
    metadata = predictions)
  df <- immunes@meta.data
  df <- df %>% group_by(harmony_cluster, predicted.celltype.l2) %>% summarise(n = n()) %>% mutate(prop = n/sum(n))
  df2 <- df%>% group_by(harmony_cluster) %>% filter(prop >0.4) %>% mutate(anno = paste(predicted.celltype.l2, collapse = "/"))
  anno <- setNames(df2$anno, df2$harmony_cluster)
  immunes$azimuth_cluster <- as.character(anno[as.character(immunes$harmony_cluster)])
  f <- DimPlot(immunes, reduction = reduction, group.by = "predicted.celltype.l2")
  print(f)
  f <- DimPlot(immunes, reduction = reduction, group.by = "azimuth_cluster")
  print(f)
  return(immunes)
}

preprocess_srt <- function(normal2){
  DefaultAssay(normal2) <- "RNA"
  lib.median <- median(normal2$nCount_RNA)
  target_pseudocount <- 1
  normal2 <- normal2 %>% 
    NormalizeData(scale.factor = lib.median / target_pseudocount) %>% 
    ScaleData(do.scale = F)
  
  set.seed(100)
  poisson_fit2 <- modelGeneVarByPoisson(normal2@assays$RNA@counts, dispersion = 0.2, size.factors = target_pseudocount * normal2$nCount_RNA / lib.median)
  residuals2 <- poisson_fit2$total - poisson_fit2@metadata$trend(poisson_fit2$mean)
  names(residuals2) <- rownames(normal2)
  top_genes2 <- rownames(normal2)[order(residuals2, decreasing = TRUE)]
  normal2 <- normal2 %>% 
    NormalizeData(scale.factor = lib.median / target_pseudocount) %>% 
    ScaleData(do.scale = F)
  
  normal2 <- normal2 %>% 
    RunPCA(features = top_genes2[1:2000]) %>% 
    FindNeighbors(dims = 1:20) %>% 
    FindClusters(resolution = 1) %>% 
    RunUMAP(dims = 1:20)
  normal2$RNA_cluster <- normal2$seurat_clusters
  
  normal2 <- cluster_sim_spectrum(object = normal2, label_tag = "sample")
  normal2 <- FindNeighbors(normal2, reduction = "css", dims = 1:ncol(Embeddings(normal2, "css")), graph.name = c("css_nn", "css_snn" ))
  normal2 <- FindClusters(normal2, resolution = 1, graph.name = "css_snn")
  normal2 <- RunUMAP(normal2, reduction = "css", dims = 1:ncol(Embeddings(normal2, "css")), reduction.name = "css_umap")
  normal2$css_cluster <- normal2$seurat_clusters
  
  normal2 <- RunHarmony(normal2, "sample") %>% 
    RunUMAP(reduction = "harmony", dims = 1:30,reduction.name = "harmony_umap") %>% 
    FindNeighbors(reduction = "harmony", dims = 1:30, graph.name = c("harmony_nn","harmony_snn")) %>% 
    FindClusters(resolution = 1, graph.name = "harmony_snn")
  normal2$harmony_cluster <- normal2$seurat_clusters
  Idents(normal2) <- normal2$harmony_cluster
  
  return(normal2)
}

#subc <- readRDS(paste0("~/normal_cell_annotation/subtype/", celltype, ".Rds"))
plot_subtype <- function(subc, celltype, markerlist){
  p1 <- coneraxes(DimPlot(subc, reduction = "harmony_umap", group.by = "subtype", cols = pal_observable()(6)), front.size = 16)+ theme(legend.position = "top", plot.title = element_blank())
  p2 <- DotPlot(subc, features = markerlist, group.by = "subtype")+theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
                                                                           axis.title = element_blank(), 
                                                                           # axis.text.y = element_text(angle = -30, hjust = 1, vjust = 1), 
                                                                           legend.position = "none", 
                                                                         panel.spacing = unit(0.1, "cm"))
  
  free(p1)+p2 +plot_layout(widths = c(1, 1.5))
  ggsave(paste0("subtypepdf/", celltype , "1.pdf"), width = 5.6, height = 3)
  
  p3 <- barplot_bypatient(subc)
  p4 <- barplot_bytissue_byrow(subc)
  p3 + p4+ plot_layout(axes = "collect")
  ggsave(paste0("subtypepdf/", celltype , "2.pdf"), width = 6, height = 3)
  
  g1 <- rmaxes(DimPlot(subc, reduction = "harmony_umap", group.by = "site")+scale_color_manual(values = site_cols)+theme(plot.title = element_blank(), legend.position = "bottom"), front.size = 26)
  g2 <- barplot_bytissue(subc)+NoLegend()
  markers <- readRDS(paste0("~/normal_cell_annotation/subtype/",celltype,"_tissuemarker.Rds"))
  features <- markers %>% group_by(cluster)%>% arrange(desc(pct.1-pct.2)) %>% arrange(desc(avg_log2FC)) %>% dplyr::slice(1)
  features$cluster <- as.character(features$cluster)
  features$cluster <- factor(features$cluster, levels =sort(as.character(features$cluster)) ) 
  group <- setNames(as.character(features$cluster), as.character(features$gene))
  g3 <- vlnplot(subc, features$gene[order(features$cluster)], group = group)
  
  g1+inset_element(g2, left = 0.6, top = 1, right = 1, bottom = 0)
  # g1+free(g2)+plot_layout(widths = c(1.5, 1))&theme( plot.title = element_blank(),legend.margin = margin(t = -10))
  ggsave(paste0("subtypepdf/", celltype , "3.pdf"), width = 6, height = 5)
  ggsave(paste0("subtypepdf/", celltype , "4.pdf"), g3, width = 6, height = 7)
}

##### ggplot theme ----------
phenotype_meta <- readRDS("~/Fig6_PSMA/phenotype_meta.Rds")
xlab <- setNames(paste(phenotype_meta$site, sapply(strsplit(phenotype_meta$sample, split = "_"), function(x) tail(x, 1))), phenotype_meta$sample)
