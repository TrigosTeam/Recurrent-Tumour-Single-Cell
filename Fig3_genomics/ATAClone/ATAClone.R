
  
  ## Set-up
  
  ### Libraries
  
  #### Code
  

library("ggplot2")
library("igraph")


### ATAClone functions

#### Code


source("/trigos_team/CASCADE/Analysis/240628_ATAClone/scripts/make_peaks_altius.R")
source("/trigos_team/CASCADE/Analysis/240628_ATAClone/scripts/ataclone_binning.R")
source("/trigos_team/CASCADE/Analysis/240628_ATAClone/scripts/ataclone_cn.R")


### Sample information

#### Parameters


id <- "CA0027"
site <- "bladder_2"
override_purple_path <- "/trigos_team/CASCADE/Sequencing_updated_reference/CA0035/bladder_base_2/purple/CA0035_bladder_base_2.purple.cnv.somatic.tsv"
override_frag_path <- NULL
override_srt_path <- NULL
override_cr_arc_barcode_stats_path <- NULL
override_cr_atac_barcode_stats_path <- NULL
```


### Sample paths


if (is.null(override_purple_path)){
  purple_path <- paste0("/trigos_team/CASCADE/Sequencing_updated_reference/", id, "/", site, "/purple/", id, "_", site, ".purple.cnv.somatic.tsv")
} else {
  purple_path <- override_purple_path
}

if (is.null(override_frag_path)){
  frag_path <- paste0("/trigos_team/CASCADE/SingleCell/", id, "/", site, "/multiome/", id, "_", site, "_multi/outs/atac_fragments.tsv.gz")
} else {
  frag_path <- override_frag_path
}

if (is.null(override_srt_path)){
  srt_path <- paste0("/trigos_team/CASCADE/Analysis/230113_seurat_intron/", id, "/", site, "/", id, "_", site, "_srt_new_filtering_doublets_removed.Rds")
} else {
  srt_path <- override_srt_path
}

if (is.null(override_cr_arc_barcode_stats_path)){
  cr_arc_barcode_stats_path <- paste0("/trigos_team/CASCADE/SingleCell/", id, "/", site, "/multiome/", id, "_", site, "_multi/outs/per_barcode_metrics.csv")
} else {
  cr_arc_barcode_stats_path <- override_cr_arc_barcode_stats_path
}

if (is.null(override_cr_atac_barcode_stats_path)){
  cr_atac_barcode_stats_path <- paste0("/trigos_team/CASCADE/SingleCell/", id, "/", site, "/ATACseq/", id, "_", site, "_atac/outs/singlecell.csv")
} else {
  cr_atac_barcode_stats_path <- override_cr_atac_barcode_stats_path
}


## Pre-processing

### Load data

#### Code


atac_fragments <- rtracklayer::import(frag_path, format = "bed")
srt <- readRDS(srt_path)
bin_width <- 10000000
atac_peaks <- get_peaks_altius("/trigos_team/CASCADE/Analysis/240628_ATAClone/files/dat_bin_FDR01_hg38.mtx.gz?download=1", "/trigos_team/CASCADE/Analysis/240628_ATAClone/files/DHS_Index_and_Vocabulary_metadata.tsv?download=1", "/trigos_team/CASCADE/Analysis/240628_ATAClone/files/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz?download=1", 15)

barcodes <- colnames(srt)

atac_fragments_filtered <- atac_fragments[atac_fragments@elementMetadata$name %in% barcodes]
atac_fragments_filtered_peak <- get_peak_overlaps(atac_fragments_filtered, atac_peaks, "PIC")

chr_arm_bins <- get_chr_arm_bins(bin_width)
all_counts <- get_bin_counts(atac_fragments_filtered, chr_arm_bins, barcodes)
peak_counts <- get_bin_counts(atac_fragments_filtered_peak, chr_arm_bins, barcodes)

saveRDS(all_counts, "all_counts.Rds")
saveRDS(peak_counts, "peak_peaks.Rds")

min_peak_count_cell <- 300
min_peak_frac_cell <- 0.05

peak_frac_cell <- colSums(peak_counts) / colSums(all_counts)

if (file.exists(cr_arc_barcode_stats_path) & file.exists(cr_atac_barcode_stats_path)){
  cr_atac_barcode_stats <- get_atac_barcode_stats(cr_arc_barcode_stats_path, cr_atac_barcode_stats_path, colnames(peak_counts))
} else {
  cr_atac_barcode_stats <- NULL
}

if (is.null(cr_atac_barcode_stats)){
  is_cell <- rep(T, ncol(peak_counts))
} else {
  is_cell <- cr_atac_barcode_stats$is__cell_barcode == 1
}

if (is.null(cr_atac_barcode_stats)){
  is_multiplet <- rep(T, ncol(peak_counts))
} else {
  is_multiplet <- cr_atac_barcode_stats$excluded_reason == 3  
}
filtering_df <- data.frame(peak_frac_cell, log10.peak.count = log10(colSums(peak_counts)), is_cell, is_multiplet)

ggplot(filtering_df, aes(x = peak_frac_cell, y = log10.peak.count, color = is_cell)) + geom_point() + geom_vline(xintercept = min_peak_frac_cell) + geom_hline(yintercept = log10(min_peak_count_cell))

ggplot(filtering_df, aes(x = peak_frac_cell, y = log10.peak.count, color = is_multiplet)) + geom_point() + geom_vline(xintercept = min_peak_frac_cell) + geom_hline(yintercept = log10(min_peak_count_cell))

peak_count_overdispersion <- 0.025
sim_peak_count_overdispersion <- 0.025

keep_idx <- which(colSums(peak_counts) > min_peak_count_cell & peak_frac_cell > min_peak_frac_cell)
peak_counts_filtered <- peak_counts[,keep_idx]
all_counts_filtered <- all_counts[,keep_idx]

peak_counts_filtered_norm <- transformGamPoi::acosh_transform(peak_counts_filtered, overdispersion = peak_count_overdispersion, size_factors = F)
peak_counts_filtered_norm_pca <- prcomp(t(peak_counts_filtered_norm), center = T, scale. = F)
peak_counts_filtered_norm_cor <- t(peak_counts_filtered_norm_pca$x[,2:ncol(peak_counts_filtered_norm_pca$x)] %*% t(peak_counts_filtered_norm_pca$rotation[,2:ncol(peak_counts_filtered_norm_pca$x)])) + peak_counts_filtered_norm_pca$center

#use to check mean-variance relationship. Also use for thresholding PCA
peak_counts_sim <- simulate_counts(peak_counts_filtered, sim_peak_count_overdispersion)
peak_counts_sim_norm <- transformGamPoi::acosh_transform(peak_counts_sim, overdispersion = peak_count_overdispersion, size_factors = F)
peak_counts_sim_norm_pca <- prcomp(t(peak_counts_sim_norm), center = T, scale. = F)
peak_counts_sim_norm_cor <- t(peak_counts_sim_norm_pca$x[,2:ncol(peak_counts_sim_norm_pca$x)] %*% t(peak_counts_sim_norm_pca$rotation[,2:ncol(peak_counts_sim_norm_pca$x)])) + peak_counts_sim_norm_pca$center

#use as alternative for thresholding PCA
peak_counts_filtered_norm_cor_permuted <- permute_rows(peak_counts_filtered_norm_cor)
peak_counts_filtered_norm_cor_permuted_pca <- prcomp(t(peak_counts_filtered_norm_cor_permuted), center = T, scale. = F)

loess_df <- data.frame(sim.mean = rowMeans(peak_counts_sim_norm_cor), sim.var = rowVars(peak_counts_sim_norm_cor))
loess_sim <- loess(sim.var ~ sim.mean, loess_df, span = 0.2)
plt_df <- data.frame(mean = rowMeans(peak_counts_filtered_norm_cor), var = rowVars(peak_counts_filtered_norm_cor))

ggplot(loess_df, aes(x = sim.mean, y = sim.var)) + geom_point() + stat_function(fun = function(x){predict(loess_sim, x)})
ggplot(plt_df, aes(x = mean, y = var)) + geom_point() + stat_function(fun = function(x){predict(loess_sim, x)})

knn_k <- 11

keep_pcs_sim <- choose_pcs(peak_counts_filtered_norm_pca$sdev, peak_counts_sim_norm_pca$sdev[2:length(peak_counts_sim_norm_pca$sdev)], tolerance = 1)
keep_pcs_permute <- choose_pcs(peak_counts_filtered_norm_pca$sdev, peak_counts_filtered_norm_cor_permuted_pca$sdev, tolerance = 1)

use_pcs <- 2:max(keep_pcs_sim, keep_pcs_permute)

peak_counts_filtered_norm_cor_denoised <- t(peak_counts_filtered_norm_pca$x[,use_pcs] %*% t(peak_counts_filtered_norm_pca$rotation[,use_pcs])) + peak_counts_filtered_norm_pca$center

peak_counts_knn <- scran::buildKNNGraph(t(peak_counts_filtered_norm_pca$x[,use_pcs]), k = knn_k, directed = F)

#leidenAlg implements a method to identify a 'significant' resolution parameter for a given graph. It turns out a similar result can be obtained by scanning resolutions and finding the one where internal edges plateaus?
#See: https://www.nature.com/articles/srep02930 and https://leidenalg.readthedocs.io/en/stable/advanced.html
try_k <- seq(knn_k - 20, knn_k + 20, 1)
try_k <- try_k[try_k > 0]
leiden_resolutions <- try_k / (gorder(peak_counts_knn) - 1)
leiden_sweep <- sweep_leiden_resolution(peak_counts_knn, "CPM", leiden_resolutions, 20, 2)
plot(leiden_resolutions, sapply(leiden_sweep, mean), ylab = "mean.internal.edges", main = paste0("k = ", knn_k))
abline(v = knn_k / (gorder(peak_counts_knn) - 1))

keep_pcs_sim <- choose_pcs(peak_counts_filtered_norm_pca$sdev, peak_counts_sim_norm_pca$sdev[2:length(peak_counts_sim_norm_pca$sdev)], tolerance = 1)
keep_pcs_permute <- choose_pcs(peak_counts_filtered_norm_pca$sdev, peak_counts_filtered_norm_cor_permuted_pca$sdev, tolerance = 1)

use_pcs <- 2:max(keep_pcs_sim, keep_pcs_permute)

peak_counts_knn <- scran::buildKNNGraph(t(peak_counts_filtered_norm_pca$x[,use_pcs]), k = knn_k, directed = F)

#resolution is interpretable as a community density threshold - discovered communities will be at least this density
#appropriate resolution is related to overall density of graph (i.e. 2|E|/(n(n - 1)) or mean.degree/(n - 1). See: https://github.com/igraph/rigraph/issues/529
#because k influences density of graph, such that mode of (degree / (n - 1)) is always k / (n - 1), I choose (k - 1) / (n - 1) as the threshold
leiden_resolution <- (knn_k - 1) / (gorder(peak_counts_knn) - 1)

#this is trying to recreate find_partition_with_rep from leidenAlg package. It's actually still very stochastic - maybe parameters are not optimised?
peak_counts_leiden <- setNames(choose_best_leiden_partition(peak_counts_knn, "CPM", leiden_resolution, 50, 20), colnames(peak_counts_filtered))

peak_counts_umap <- umap::umap(peak_counts_filtered_norm_pca$x[,use_pcs])
colnames(peak_counts_umap$layout) <- c("UMAP_1", "UMAP_2")

plt_df <- data.frame(peak_counts_umap$layout, cluster = peak_counts_leiden, log10.peak.counts = log10(colSums(peak_counts_filtered)), droplet.peak.frac = colSums(peak_counts_filtered) / colSums(all_counts_filtered))
ggplot(plt_df, aes(x = UMAP_1, y= UMAP_2, color = cluster)) + geom_point()
ggplot(plt_df, aes(x = UMAP_1, y= UMAP_2, color = log10.peak.counts)) + geom_point() + scale_color_gradient(low = "blue", high = "yellow")
ggplot(plt_df, aes(x = UMAP_1, y= UMAP_2, color = droplet.peak.frac)) + geom_point() + scale_color_gradient(low = "blue", high = "yellow")

plt_df2 <- data.frame(srt@reductions$umap@cell.embeddings, cluster.atac = peak_counts_leiden[colnames(srt)], cluster.rna = srt$seurat_clusters)
ggplot(plt_df2, aes(x = UMAP_1, y= UMAP_2, color = cluster.rna)) + geom_point()
ggplot(plt_df2, aes(x = UMAP_1, y= UMAP_2, color = cluster.atac)) + geom_point()

saveRDS(peak_counts_leiden, "leiden_cluster.Rds")

ref.cluster <- 4

single_cn_estimates <- get_copy_number(peak_counts_filtered, peak_counts_leiden, ref.cluster = ref.cluster, joint = F)
joint_cn_estimates <- get_copy_number(peak_counts_filtered, peak_counts_leiden, ref.cluster = ref.cluster, joint = T)
#new_orders <- get_new_order(joint_cn_estimates, peak_counts_filtered_norm_cor, peak_counts_leiden)
new_feature_names <- get_new_feature_names(rownames(joint_cn_estimates))
feature_factor <- get_feature_factor(new_feature_names)
new_cell_order <- get_new_cell_order(joint_cn_estimates, peak_counts_filtered_norm_cor_denoised, peak_counts_leiden)
cell_factor <- get_cell_factor(joint_cn_estimates[,new_cell_order], peak_counts_leiden[new_cell_order])
chr_arm <- ComplexHeatmap::HeatmapAnnotation(chr_arm = get_chr_arms(new_feature_names), col = list(chr_arm = c(p = "black", q = "grey")), show_legend = F)

color_fun <- circlize::colorRamp2(breaks = c(0,2,4), colors = c("blue", "white", "red"))
ComplexHeatmap::Heatmap(t(single_cn_estimates[new_feature_names,new_cell_order]), cluster_rows = F, cluster_columns = F, column_split = feature_factor, row_split = cell_factor, row_labels = rep('', ncol(joint_cn_estimates)), column_labels = rep('', nrow(joint_cn_estimates)), border = T, col = color_fun, top_annotation = chr_arm)
ComplexHeatmap::Heatmap(t(joint_cn_estimates[new_feature_names,new_cell_order]), cluster_rows = F, cluster_columns = F, column_split = feature_factor, row_split = cell_factor, row_labels = rep('', ncol(joint_cn_estimates)), column_labels = rep('', nrow(joint_cn_estimates)), border = T, col = color_fun, top_annotation = chr_arm)
#need to denoise!! - also affects hclust for reordering cells
ComplexHeatmap::Heatmap(scale(t(peak_counts_filtered_norm_cor_denoised[new_feature_names,new_cell_order]), scale = F), cluster_rows = F, cluster_columns = F, column_split = feature_factor, row_split = cell_factor, row_labels = rep('', ncol(joint_cn_estimates)), column_labels = rep('', nrow(joint_cn_estimates)), border = T, top_annotation = chr_arm)

saveRDS(single_cn_estimates, "single_cn.Rds")
saveRDS(joint_cn_estimates, "joint_cn.Rds")
saveRDS(peak_counts_filtered_norm_cor_denoised, "denoised_norm.Rds")


subclones <- split(names(leiden_cluster), leiden_cluster)


submat <- sapply(subclones, function(x){
  rowMeans(joint_cn[, x])
})

submat <- na.omit(submat)

ref.cluster <- names(which.min(apply(submat, 2, var)))
tumormat <- submat[, setdiff(colnames(submat), ref.cluster)]

vars_df <- apply(tumormat, 1, var)
summary(vars_df, breaks = 100, xlim = c(0, 0.4))
hist(vars_df)

ComplexHeatmap::Heatmap(t(joint_cn_estimates[new_feature_names,]), cluster_rows = F, cluster_columns = F, column_split = feature_factor, row_split = leiden_cluster, row_labels = rep('', ncol(joint_cn_estimates)), column_labels = rep('', nrow(joint_cn_estimates)), border = T, col = color_fun, top_annotation = chr_arm)
