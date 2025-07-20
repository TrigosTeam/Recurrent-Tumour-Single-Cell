library(Seurat)
library(scran)
library(scater)
library(data.table)
library(dplyr)
library(patchwork)
library(tidyverse)
library(MAST)
library(pals)
library(circlize)
library(ComplexHeatmap)
source("~/functions.R")
# subclone visual UMAP, AR/ASCL1 expression, heatmap---------
paths <- system("realpath ~/230113_seurat_intron/*/*/*non_tumour_removed.Rds", intern = T)
paths_split <- unlist(lapply(strsplit(paths, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
paths2 <- system("realpath ~/240603_seurat_intron/CA0027/dura*/*non_tumour_removed_v4.Rds", intern = T)
paths_split2 <- unlist(lapply(strsplit(paths2, split = "/"), function(x) paste(x[6], x[7], sep = "_")))

paths <- c(paths, paths2)
paths_split <- c(paths_split, paths_split2)
paths <- setNames(paths,paths_split)
meta <- readRDS("~/Fig5_archetype/5_integrated_final_module_meta.Rds")

sub_path <- system("realpath ~/240628_ATAClone/*/*/leiden_cluster.Rds", intern = T)
subclones <- list()
for (p in  sub_path){
  sample = str_split(p, "/")[[1]][7]
  patient = str_split(p, "/")[[1]][6]
  i = paste(patient, sample,sep = "_")
  
  srt<- readRDS(paths[i])
  
  DefaultAssay(srt) <- "RNA"
  lib.median <- median(srt$nCount_RNA)
  target_pseudocount <- 1
  srt <- srt %>% 
    NormalizeData(scale.factor = lib.median / target_pseudocount) %>% 
    ScaleData(do.scale = F)
  
  set.seed(100)
  poisson_fit2 <- modelGeneVarByPoisson(srt@assays$RNA@counts, dispersion = 0.2, size.factors = target_pseudocount * srt$nCount_RNA / lib.median)
  residuals2 <- poisson_fit2$total - poisson_fit2@metadata$trend(poisson_fit2$mean)
  names(residuals2) <- rownames(srt)
  top_genes2 <- rownames(srt)[order(residuals2, decreasing = TRUE)]
  
  srt <- srt %>% 
    RunPCA(features = top_genes2[1:2000]) %>% 
    FindNeighbors(dims = 1:20) %>% 
    FindClusters(resolution = 1) %>% 
    RunUMAP(dims = 1:20)
  
  joint_cn <- readRDS(paste("~/240628_ATAClone", patient, sample, "joint_cn.Rds", sep = "/"))
  subclone <- readRDS(p)
  srt$subclone <- subclone[colnames(srt)]
  Idents(srt) <- srt$subclone
  ind <- setdiff(levels(srt$subclone), names(which.max(sapply(levels(srt$subclone), function(x) sum(joint_cn[, colnames(srt)[which(srt$subclone == x)]] %in% c(2, 1))))))
  
  g <- DimPlot(srt, group.by = "subclone", cells =colnames(srt)[srt$subclone %in% ind])
  g <- coneraxes(g)
  g2 <- VlnPlot(srt, features = c("AR", "ASCL1"), group.by = "subclone", idents = ind)
  f <- free(g)|g2
  ggsave(paste0("sublcone_fig/", i,".pdf"), f, width = 9 , height= 3)
  # saveRDS(f, paste0("~/Fig3/ATAClone/plot_obj/", i, "_srt.Rds"))
  
  
  chrom <- sapply(strsplit(rownames(joint_cn), split = ".", fixed = T), "[", 1)
  chrom <- factor(gsub("chr", "", chrom), levels = c(as.character(1:22), "X", "Y"))
  
  cnmat <- joint_cn[, colnames(srt)[srt$subclone %in% ind]]
  rols <- as.character(subclone[colnames(srt)[srt$subclone %in% ind]])
  color_fun <- circlize::colorRamp2(breaks = c(0,2,4), colors = c("blue", "white", "red"))
  
  ht <- Heatmap(t(cnmat), cluster_columns = F, cluster_rows = F, column_split = chrom, col = color_fun, show_row_names = F, show_column_names = F, 
          row_split = rols, border =  T, name = "CN", 
          right_annotation =  rowAnnotation(subclone = rols, col = list(subclone = setNames(ggsci::pal_futurama()(length(ind)), unique(rols))),
                                            show_annotation_name=F, annotation_legend_param = list(
                                              labels_gp = gpar(fontsize = 13),
                                            title_gp = gpar(fontsize = 14)
                                            )),
          row_title_gp = gpar(fontsize = 13),
          column_title_gp = gpar(fontsize= 13),
          heatmap_legend_param = list(labels_gp = gpar(fontsize = 13), title_gp = gpar(fontsize = 14)))
  
  pdf(paste0("sublcone_fig/", i,"_hm.pdf"),width = 9 , height= 2.5)
  draw(ht, column_title = "ATAClone output",column_title_gp = gpar(fontsize = 16, fontface = "bold"),merge_legend = TRUE)
  dev.off()
  
  subclones[[i]] <- subclone[subclone %in% ind & names(subclone) %in% colnames(srt)]
  print(i)
}
  
saveRDS(subclones, "subclones.Rds")


# saveRDS(ht, paste0("~/Fig3/ATAClone/plot_obj/", i, "_heatmap.Rds"))





# per-patient joint cn profile---------
subclones <- readRDS("~/Fig3_genomics_normal_cells/ATAClone/subclones.Rds")
samples <- list()
for (i in names(subclones)){
  names(subclones[[i]]) <- paste(i, names(subclones[[i]]), sep = "_")
  samples[[i]] <- setNames(rep(i, length(subclones[[i]])), names(subclones[[i]]))
}

gf = rtracklayer::import("~/reference_annotation/multiome/genes.gtf")
gf_genes <- gf[gf$type =="gene"]


genelist <- lapply(all_gene_sets, function(x) unique(unlist(x))) # union of different signatures

tang_2022 <- all_gene_sets$tang_2022
names(tang_2022) <- paste0("tang_2022_", names(tang_2022))

temp <- lapply(all_gene_sets, function(x){
  names(which(table(unlist(x)) >= length(x)/2))
})

temp <- temp[!names(temp)%in% c("CRPC", "cancer", "hillock", "metastasis", "club", "invasion", "tang_2022")]
genelist <- c(temp, tang_2022)
names(genelist) <- paste0(names(genelist), "_signature")

gf_sub <- gf_genes[gf_genes$gene_name %in% c(genelist$AR_signature, genelist$NE_signature)]
phenotype_meta <- readRDS("~/PSMA/phenotype_meta.Rds")
xlab <- setNames(paste( phenotype_meta$site, sapply(strsplit(phenotype_meta$sample, split = "_"), function(x) tail(x, 1))), phenotype_meta$sample)


for (p in unique(substr(names(subclones), 1, 6))[c(1,6)]){
  paths <- system(paste0("realpath ~/240628_ATAClone/",p,"/*/joint_cn.Rds"), intern = T)
  joint_cn <- lapply(paths, readRDS)
  sample <- sapply(str_split(paths, "/"), "[", 7)
  
  names(joint_cn) <- paste(p, sample, sep = "_")
  sub_sublcone <- subclones[names(joint_cn)]
  names(sub_sublcone) <- NULL
  sub_sublcone <- unlist(sub_sublcone)
  sub_samples <- samples[names(joint_cn)]
  names(sub_samples) <- NULL
  sub_samples <- unlist(sub_samples)
  
  for (i in names(joint_cn)){
    colnames(joint_cn[[i]]) <- paste(i, colnames(joint_cn[[i]]), sep = "_")
  }
  
  test <- Reduce(cbind, joint_cn)
  
  
  cnmat <- t(test)
  cnmat <- cnmat[names(sub_sublcone), ]
  rols <- sub_sublcone[rownames(cnmat)]
  rols2 <- xlab[sub_samples[rownames(cnmat)]]
  
  color_fun <- circlize::colorRamp2(breaks = c(0,2,4), colors = c("blue", "white", "red"))
  chrom <- sapply(strsplit(colnames(cnmat), split = ".", fixed = T), "[", 1)
  chrom <- gsub("chr", "",chrom)
  cmat <- split(cnmat,  paste(rols2, rols, sep = "-"))
  cs <- table(paste(rols2, rols, sep = "-"))
  cmat <- lapply(names(cmat), function(x) matrix(cmat[[x]], nrow= cs[x])[1, ])
  names(cmat) <- names(cs)
  cdf <- as.data.frame(cmat)
  cdf <- na.omit(cdf)
  colnames(cdf) <- names(cmat)
  hc <- hclust(dist(t(cdf)),"centroid")
  order <- gsub(".", " ", fixed= T, labels(as.dendrogram(hc)))
  rol <- factor(paste(rols2, rols, sep = "-"), levels = order)
  
  ht <- Heatmap(cnmat,col = color_fun, border =  T, name = "CN est", 
                cluster_columns = F, cluster_rows = F,cluster_row_slices = FALSE, 
                column_split = factor(chrom, levels = c(1:22, "X", "Y")),  row_split = rol,  row_gap = unit(0, "mm"),
                show_row_names = F, show_column_names = F, row_title_rot = 0,
                right_annotation =  rowAnnotation( Sample = rols2,
                                                  col = list(Sample = setNames(pals::brewer.set3(9)[1:n_distinct(rols2)], unique(rols2))),
                                                  show_annotation_name=F, 
                                                  annotation_legend_param = list(
                                                    labels_gp = gpar(fontsize = 13),
                                                    title_gp = gpar(fontsize = 16)
                                                  )),
                row_title_gp = gpar(fontsize = 13),
                column_title_gp = gpar(fontsize= 13),
                heatmap_legend_param = list(labels_gp = gpar(fontsize = 13), title_gp = gpar(fontsize = 16)))
  
  pdf(paste0("perpatient_subclone_fig/", p,".pdf"),width = 10 , height= 6)
  #plot(hc)
  draw(ht, column_title = gsub("00", "", p),column_title_gp = gpar(fontsize = 13, fontface = "bold"))
  dev.off()
  
  srt <- readRDS(paste0("~/Fig2_signature/sc_integration/perpatient/tumor_only/", p, "_ruvsrt.Rds"))
  names(rol) <- rownames(cnmat)
  srt$subclone <- NA
  srt$subclone[names(rol)] <- as.character(rol)
  cols <- randcol(n_distinct(rol))
  coneraxes(DimPlot(srt, group.by = "subclone", cells = intersect(colnames(srt), rownames(cnmat)), 
                    cols = cols)) +
    theme(legend.position = "bottom")+
    guides(color = guide_legend(ncol = 3,override.aes = list(size = 2)))
  ggsave(paste0("perpatient_subclone_fig/", p,"_UMAP.pdf"),width = 4 , height= 6)
  f1 <- rmaxes(FeaturePlot(srt, features = "AR", cells = intersect(colnames(srt), rownames(cnmat)))) 
  f2 <- coneraxes(FeaturePlot(srt, features = "ASCL1", cells = intersect(colnames(srt), rownames(cnmat))))
  print(f1/f2)
  ggsave(paste0("perpatient_subclone_fig/", p,"_Feature.pdf"),width = 4 , height= 6)
  
  VlnPlot(srt,group.by = "subclone", features = c("AR", "ASCL1"), ncol = 1, cols = cols)
  ggsave(paste0("perpatient_subclone_fig/", p,"_Vln.pdf"),width = 4 , height= 6)
  # hit <- findOverlaps(bingr, gf_sub)
  # ind <- split( to(hit), from(hit))
  # cols <- sapply(ind, function(x) paste(gf_sub$gene_name[x], collapse = ","))
  # chrom <- sapply(strsplit(colnames(cnmat[, unique(from(hit))]), split = ".", fixed = T), "[", 1)
  # chrom <- factor(chrom, levels = paste0("chr", c(1:22, "X")))
  # text <- lapply(ind, function(x) data.frame(text = gf_sub$gene_name[x], fontsize = 6))
  # names(text) <- chrom
  # 
  # 
  # ht2 <- Heatmap(cnmat[, unique(from(hit))],col = color_fun, border =  T, name = "CN est", 
  #                cluster_columns = F, cluster_rows = F,cluster_row_slices = FALSE, 
  #                column_split = chrom,  row_split = rol,  row_gap = unit(0, "mm"),
  #                show_row_names = F, show_column_names = F, row_title_rot = 0,
  #                right_annotation =  rowAnnotation(subclone = rols, sample = rols2,
  #                                                  col = list(subclone = setNames(ggsci::pal_futurama()(length(levels(rols))), levels(rols)),
  #                                                             sample = setNames(pals::brewer.set3(9)[1:length(joint_cn)], names(joint_cn))),
  #                                                  show_annotation_name=F, 
  #                                                  annotation_legend_param = list(
  #                                                    labels_gp = gpar(fontsize = 12),
  #                                                    title_gp = gpar(fontsize = 16)
  #                                                  )),
  #                top_annotation =  columnAnnotation(textbox = anno_textbox(chrom, text, by = "anno_block")),
  #                row_title_gp = gpar(fontsize = 12),
  #                column_title_gp = gpar(fontsize= 10),
  #                heatmap_legend_param = list(labels_gp = gpar(fontsize = 12), title_gp = gpar(fontsize = 16)))
  # 
  # pdf(paste0(" perpatient_subclone_fig/", p,"_NE.pdf"),width = 10 , height= 6)
  # draw(ht2, column_title = p,column_title_gp = gpar(fontsize = 12, fontface = "bold"))
  # dev.off()
  print(p)
}



