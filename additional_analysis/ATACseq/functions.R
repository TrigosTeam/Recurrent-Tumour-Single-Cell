require(patchwork)
require(ggplot2)
require(ggpubr)
rmaxes <- function(g, front.size = 18, print = F) {
  
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

plot_fun <- function(srt, reduction, title){
  f1 <- DimPlot(srt, reduction = reduction, label = T) + NoLegend() + ggtitle(title)
  f1 <- coneraxes(f1, print = F, labelx = paste0(srt@reductions[[reduction]]@key, "1"), labely = paste0(srt@reductions[[reduction]]@key, "2"))
  f2 <- FeaturePlot(srt, reduction = reduction,  paste(names(clean_module), "Module"), ncol = 6, raster = F)
  f2 <- rmaxes(f2, print = F)
  p <- f1+f2 +plot_layout(widths = c(1, 6))
  return(p)
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
    RunPCA(features = top_genes2[1:2000]) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 1) %>%
    RunUMAP(dims = 1:20)
}

give_module_anno <- function(meta){
  df <- lapply(seq(from = 45, to = 70, by =1), function(p){
    group<- apply(meta, 2, function(x){
      bins <- cut(x, breaks = 100)
      group <- bins %in% levels(bins)[p:100]
    })
    ind_high <- apply(group[, 2:6], 1, function(x) sum(x))
    return(as.data.frame(table(ind_high)) %>% mutate(cutoff = p))
  })
  df <- rbindlist(df)
  df$cutoff <- factor(df$cutoff, levels = df %>% dplyr::filter(ind_high ==1) %>% arrange(Freq) %>% pull(cutoff))
  df %>% dplyr::filter(ind_high == "1") -> subdf
  
  cutoff <- as.numeric(as.character(subdf$cutoff[which.max(subdf$Freq)]))
  
  group<-  apply(meta, 2, function(x){
    bins <- cut(x, breaks = 100)
    group <- bins %in% levels(bins)[cutoff:100]
  })
  ind_high <- apply(group, 1, function(x) if(sum(x) > 1) data.frame( group = paste(names(clean_module)[x], collapse = "&"), freq = sum(x)) else if(sum(x)==1) data.frame(group = names(clean_module)[x], freq = 1) else data.frame(group = "Background", freq = 0))
  ind_high <- rbindlist(ind_high)
  ind <- grepl("AR&", ind_high$group)&ind_high$freq == 2
  ind_high$group[ind] <- gsub("AR&", "",ind_high$group[ind])
  ind_high$freq[ind] <- ind_high$freq[ind]-1
  return(ind_high$group)
}



  
  
  
  
  
  
  
  
