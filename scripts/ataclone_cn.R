get_copy_number <- function(count_matrix, clusters, ref.cluster, joint, ploidy.correction.method = "prop"){
  require("MatrixGenerics")
  #removed bias correction - need to investigate whether it improves calling. Scaling to mean in presence of zeroes may be a problem.
  cor.list <- list()
  for (i in levels(clusters)){
    cor.list[[i]] <- rowMeans(count_matrix[,as.integer(clusters) == i])
  }
  cn.mat <- count_matrix
  if (ploidy.correction.method == "prop"){
    for (i in seq_along(cor.list)){
      cor.list[[i]] <- cor.list[[i]] / sum(cor.list[[i]])
    }
    if (!joint){
      cn.mat <- t(t(count_matrix) / colSums(count_matrix))
    }
  }
  if (joint){
    for (i in levels(clusters)){
      cn.mat[,as.integer(clusters) == i] <- rep(2 * as.numeric(cor.list[[i]]) / as.numeric(cor.list[[ref.cluster]]), sum(as.integer(clusters) == i))
      cn.mat[grep("^chrX|^chrY", rownames(cn.mat)),as.integer(clusters) == i] <- 0.5 * cn.mat[grep("^chrX|^chrY", rownames(cn.mat)),as.integer(clusters) == i]
    }
  } else {
    for (i in levels(clusters)){
      cn.mat[,as.integer(clusters) == i] <- 2 * cn.mat[,as.integer(clusters) == i] / as.numeric(cor.list[[ref.cluster]])
      cn.mat[grep("^chrX|^chrY", rownames(cn.mat)),as.integer(clusters) == i] <- 0.5 * cn.mat[grep("^chrX|^chrY", rownames(cn.mat)),as.integer(clusters) == i]
    }
  }
  cn.mat[is.infinite(cn.mat)] <- NA
  cn.mat
}

get_new_order <- function(cn.mat, count_matrix_cor, clusters){
  #sort within clusteres based on corrected count matrix
  hclust.list <- list()
  for (i in levels(clusters)){
    cell.dist <- dist(t(count_matrix_cor[,as.integer(clusters) == i]))
    hclust.list[[i]] <- colnames(cn.mat)[as.integer(clusters) == i][hclust(cell.dist, method = "ward.D2")$order]
  }
  #perform hierarchical clustering on per-cluster copy-number estimates
  #get first cell for each cluster (they all have the same estimate)
  cluster.idx <- setNames(integer(length(levels(clusters))), levels(clusters))
  for (cluster in levels(clusters)){
    cluster.idx[cluster] <- which(as.integer(clusters) == cluster)[1]
  }
  cluster.dist <- dist(t(cn.mat[,cluster.idx]), method = "manhattan")
  cluster.hclust <- hclust(cluster.dist, method = "complete")
  cluster.order <- names(cluster.idx)[cluster.hclust$order]
  cluster.order <- c(ref.cluster, cluster.order[cluster.order != ref.cluster])
  new.cluster.levels <- factor(clusters[unlist(hclust.list[cluster.order])], levels = cluster.order)
  
  chr.num <- as.integer(gsub("chr", "", sapply(strsplit(rownames(cn.mat), "\\."), `[`, 1)))
  new.rownames <- c(rownames(cn.mat)[order(chr.num)[!is.na(chr.num)]], rownames(cn.mat)[is.na(chr.num)])
  
  new.row.levels <- gsub("\\.", "", gsub("chr", "", sapply(strsplit(new.rownames, ":"), `[`, 1)))
  new.row.levels <- factor(new.row.levels, unique(new.row.levels))
  
  list(new.rownames, unlist(hclust.list), new.cluster.levels, new.row.levels)
}

get_new_feature_names <- function(feature.names){
  chr.num <- as.integer(gsub("chr", "", sapply(strsplit(feature.names, "\\."), `[`, 1)))
  new.rownames <- c(feature.names[order(chr.num)[!is.na(chr.num)]], feature.names[is.na(chr.num)])
  
  new.rownames
}

get_feature_factor <- function(feature.names){
  #chr.num <- as.integer(gsub("chr", "", sapply(strsplit(feature.names, "\\."), `[`, 1)))
  
  #new.row.levels <- gsub("\\.", "", gsub("chr", "", sapply(strsplit(feature.names, ":"), `[`, 1)))
  #new.row.levels <- factor(new.row.levels, unique(new.row.levels))
  #new.row.levels
  new.row.levels <- gsub("chr", "", sapply(strsplit(feature.names, "\\."), `[`, 1))
  factor(new.row.levels, unique(new.row.levels))
}

get_new_cell_order <- function(cn.mat, count_matrix_cor, clusters){
  #sort within clusteres based on corrected count matrix
  hclust.list <- list()
  for (i in levels(clusters)){
    cell.dist <- dist(t(count_matrix_cor[,as.integer(clusters) == i]))
    hclust.list[[i]] <- colnames(cn.mat)[as.integer(clusters) == i][hclust(cell.dist, method = "ward.D2")$order]
  }
  unlist(hclust.list)
}

get_cell_factor <- function(cn.mat, clusters){
  #perform hierarchical clustering on per-cluster copy-number estimates
  #get first cell for each cluster (they all have the same estimate)
  cluster.idx <- setNames(integer(length(levels(clusters))), levels(clusters))
  for (cluster in levels(clusters)){
    cluster.idx[cluster] <- which(as.integer(clusters) == cluster)[1]
  }
  cluster.dist <- dist(t(cn.mat[,cluster.idx]), method = "manhattan")
  cluster.hclust <- hclust(cluster.dist, method = "complete")
  cluster.order <- names(cluster.idx)[cluster.hclust$order]
  cluster.order <- c(ref.cluster, cluster.order[cluster.order != ref.cluster])
  new.cluster.levels <- factor(clusters, levels = as.integer(cluster.order))
  new.cluster.levels
}

get_chr_arms <- function(bin.names){
  sapply(strsplit(sapply(strsplit(bin.names, ":"), `[`, 1), "\\."), `[`, 2)
}

reorder_matrix <- function(){}

#cell.list <- as.list(as.data.frame(peak_counts_filtered[new_feature_names,]))
#cell.arm.count <- do.call(cbind, lapply(cell.list, get_arm_counts, feature_factor, get_chr_arms(new_feature_names)))
get_arm_counts <- function(x, chr.factor, arm.factor){
  sapply(split(x, list(chr.factor, arm.factor)), sum)
}

get_purple_cn <- function(purple_path, bin.names){
  purple.cn <- read.delim(purple_path, header = T)
  purple.cn.granges <- makeGRangesFromDataFrame(purple.cn)
  ataclone.granges <- GRanges(gsub("\\..", "", bin.names))
  cn.overlaps <- findOverlaps(ataclone.granges, purple.cn.granges)
  cn.overlaps.split <- split(cn.overlaps@to, cn.overlaps@from)
  cn.overlaps.widths <- list()
  for (i in seq_along(cn.overlaps.split)){
    cn.overlaps.widths[[i]] <- numeric(length = length(cn.overlaps.split[[i]]))
    for (j in seq_along(cn.overlaps.widths[[i]])){
      cn.overlaps.widths[[i]][j] <- width(GenomicRanges::intersect(ataclone.granges[i], purple.cn.granges[cn.overlaps.split[[i]][j]]))
    }
  }
  ataclone.purple.estimates <- numeric(length(cn.overlaps.split))
  for (i in seq_along(ataclone.purple.estimates)){
    ataclone.purple.estimates[i] <- sum(purple.cn[cn.overlaps.split[[i]],"copyNumber"] * cn.overlaps.widths[[i]] / sum(cn.overlaps.widths[[i]]))
  }
  ataclone.purple.estimates
}
