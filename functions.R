require(Seurat)
require(ggplot2)
require(patchwork)
require(pals)
require(ggrastr)
require(ggpubr)


rmaxes <- function(g, front.size = 18, print = TRUE) {
  
  remove_axes <- theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks   = element_blank(),
    axis.line = element_blank(),
    plot.margin  = margin(0, 0, 0, 0, "cm"),
    text = element_text(size = front.size),
    legend.key.size = unit(0.2, "inches")
  )
  
  # If it's a patchwork object:
  if (inherits(g, "patchwork")) {
    # Apply theme to every plot inside patchwork
    g2 <- g & remove_axes   # <-- patchwork syntax to apply theme to all subplots
    if (print) print(g2)
    return(g2)
  }
  
  # If it's a regular ggplot:
  g2 <- g + remove_axes
  g2 <-  rasterise(g2, layers = "Point", dpi = 300)
  if (print) print(g2)
  return(g2)
}


coneraxes <- function(p, which_plot = 1, labelx = "UMAP1", labely = "UMAP2", 
                      colorbar = TRUE, front.size = 18, print = T) {
  
  # 1. Determine if input is a Patchwork object or a single ggplot
  is_patchwork <- inherits(p, "patchwork")
  
  # 2. Select the target plot to annotate
  if (is_patchwork) {
    # Extract the specific plot from the patchwork list
    target_p <- p[[which_plot]]
  } else {
    target_p <- p
  }
  
  # 3. Apply logic to the target plot
  if (colorbar) {
    # Ensure we are using the data from the specific subplot
    # Note: Assumes columns 1 and 2 are the coordinates (common in Seurat/DimPlot)
    x_min <- min(target_p$data[, 1], na.rm = TRUE)
    y_min <- min(target_p$data[, 2], na.rm = TRUE)
    x_range <- diff(range(target_p$data[, 1], na.rm = TRUE))
    
    # Calculate arrow dimensions
    arr <- list(
      x = x_min - 1, 
      y = y_min - 1, 
      x_len = x_range / 5, 
      y_len = x_range / 5
    )
    
    # Handle the 'rmaxes' function: 
    # Use it if it exists in your environment, otherwise manually remove axes
    if (exists("rmaxes")) {
      f <- rmaxes(target_p, print = FALSE)
    } else {
      f <- target_p + theme_void() + theme(legend.position = "none") 
    }
    
    # Add annotations to the subplot
    f <- f + 
      annotate("segment", 
               x = arr$x, xend = arr$x + c(arr$x_len, 0), 
               y = arr$y, yend = arr$y + c(0, arr$y_len), 
               arrow = arrow(type = "closed", length = unit(3, 'pt'))) + 
      
      # X-axis Label
      annotate(geom = "text", 
               x = arr$x, 
               y = arr$y - 0.2, 
               label = labelx, 
               hjust = "left", 
               vjust = "top", 
               size = 15/.pt) + 
      
      # Y-axis Label
      annotate(geom = "text", 
               x = arr$x - 0.2, 
               y = arr$y, 
               label = labely, 
               angle = 90, 
               hjust = "left", 
               vjust = "bottom", 
               size = 15/.pt) + 
      
      # Apply Theme settings
      theme(axis.title = element_blank(), 
            text = element_text(size = front.size),
            plot.margin = margin(t = 0, l = 0.1, r = 0, b = 0.1, "cm"),
            legend.key.size = unit(0.2, "inches"),
            legend.margin = margin(l = -10, t = -10)) + 
      scale_x_continuous(breaks = NULL) + 
      scale_y_continuous(breaks = NULL) + 
      coord_cartesian(clip = 'off')
  } else {
    # If colorbar is FALSE, just return original (or formatted) plot
    f <- target_p
  }
  p <- rmaxes(p, print = F)
  # 4. Return based on object type
  if (is_patchwork) {
    # Put the modified plot back into the patchwork object
    p[[which_plot]] <- f
    if (print) print(p)
    return(p)
  } else {
    if (print) print(f)
    return(f)
  }
}

preprocess_srt <- function(srt){
  require(Seurat)
  require(scran)
  DefaultAssay(srt) <- "RNA" 
  lib.median <- median(srt$nCount_RNA)
  target_pseudocount = 1
  srt <- srt %>%
    NormalizeData(scale.factor = lib.median / target_pseudocount) %>%
    ScaleData(do.scale = FALSE)
  set.seed(100)
  poisson_fit <- modelGeneVarByPoisson(GetAssayData(srt, slot = "counts"), dispersion = 0.2, size.factors = target_pseudocount * srt$nCount_RNA / lib.median)
  residuals <- poisson_fit$total - poisson_fit@metadata$trend(poisson_fit$mean)
  names(residuals) <- rownames(srt)
  top_genes <- rownames(srt)[order(residuals, decreasing = TRUE)]
  
  srt <- srt %>%
    RunPCA(features = top_genes[1:2000]) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 1) %>%
    RunUMAP(dims = 1:20)
  return(srt)
}

give_module_anno <- function(meta){
  df <- lapply(seq(from = 45, to = 70, by =1), function(p){
    group<- apply(meta[, 1:6], 2, function(x){
      bins <- cut(x, breaks = 100)
      group <- bins %in% levels(bins)[p:100]
    })
    ind_high <- apply(group, 1, function(x) sum(x))
    return(as.data.frame(table(ind_high)) %>% mutate(cutoff = p))
  })
  df <- rbindlist(df)
  df$cutoff <- factor(df$cutoff, levels = df %>% filter(ind_high ==1) %>% arrange(Freq) %>% pull(cutoff))
  df %>% filter(ind_high == "1") -> subdf
  
  cutoff <- as.numeric(as.character(subdf$cutoff[which.max(subdf$Freq)]))
  
  
  group<-  apply(meta[, 1:6], 2, function(x){
    bins <- cut(x, breaks = 100)
    group <- bins %in% levels(bins)[cutoff:100]
  })
  ind_high <- apply(group, 1, function(x) if(sum(x) > 1) data.frame( group = paste(c("AR","Inflammation", "NE1","NE2", "Cycling","Glycolysis")[x], collapse = "&"), freq = sum(x)) else if(sum(x)==1) data.frame(group = c("AR","Inflammation", "NE1","NE2", "Cycling","Glycolysis")[x], freq = 1) else data.frame(group = "background", freq = 0))
  ind_high <- rbindlist(ind_high)
  meta$oldgroup <- ind_high$group
  ind <- grepl("AR&", ind_high$group)&ind_high$freq == 2
  ind_high$group[ind] <- gsub("AR&", "",ind_high$group[ind])
  ind_high$freq[ind] <- ind_high$freq[ind]-1
  cat("cutoff above", cutoff, "percentile\n")
  return(ind_high$group)
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

#subc <- readRDS(paste0("~/CASCADEpaper/paper/normal_cells_202406/subtype/", celltype, ".Rds"))
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
  markers <- readRDS(paste0("~/CASCADEpaper/paper/normal_cells_202406/subtype/",celltype,"_tissuemarker.Rds"))
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
phenotype_meta <- readRDS("~/CASCADEpaper/paper/PSMA/phenotype_meta.Rds")
xlab <- setNames(paste(phenotype_meta$site, sapply(strsplit(phenotype_meta$sample, split = "_"), function(x) tail(x, 1))), phenotype_meta$sample)
