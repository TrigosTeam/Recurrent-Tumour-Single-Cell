library(Seurat)
library(dplyr)
library(cluster)
library(reshape2)

paths <- system("realpath ~/integration/2024_06/perpatient/tumor_only/*_ruvsrt.Rds", intern = T)
paths_split <- substr(unlist(lapply(strsplit(paths, split = "/"), function(x) x[8])),1 ,6) 
paths <- setNames(paths, paths_split)
dism_snp <- readRDS("~/CASCADEpaper/paper/Fig3_genomics/dism_snp.Rds")

sh_list <- list()
dist_list <- list()
for (p in paths_split){
  perp_srt <- readRDS(paths[p])
  perp_srt$labs <- paste(gsub("00", "", perp_srt$patient), perp_srt$site, sapply(strsplit(perp_srt$sample, split = "_"), function(x) tail(x, 1)) )
  pca <- Embeddings(perp_srt, "pca")
  eu_dist <- as.matrix(dist(pca[, 1:20]))
  sh <- as.data.frame(silhouette(x = as.integer(factor(perp_srt$labs)), dmatrix = eu_dist))
  cluster <- factor(perp_srt$labs)
  patient <- gsub("00", "", p)
  labs <- setNames(levels(cluster),unique(as.integer(cluster)))
  sh$cluster <- labs[as.character(sh$cluster)]
  sh$neighbor <- labs[as.character(sh$neighbor)]
  sh <- reshape2::dcast(sh, value.var = "sil_width", cluster~neighbor, fun.aggregate = mean) 
  rownames(sh) <- sh$cluster
  sh_list[[p]] <- sh
  
  if(p == "CA0035"){
    colnames(dism_snp[[patient]]) <- c("CA35 bladder 2" ,"CA35 LN 1")
    rownames(dism_snp[[patient]]) <- c("CA35 bladder 2" ,"CA35 LN 1")
  }
  xy <- t(combn(colnames(dism_snp[[patient]]), 2))
  snpdf <- data.frame(xy, dist = dism_snp[[patient]][xy]) %>% mutate(obj = paste(X1, X2)) %>% select(c("dist", "obj")) %>% `colnames<-`(c("SNP_dis", "obj"))
  shdf <-  data.frame(xy, dist = sh[xy]) %>% mutate(obj = paste(X1, X2)) %>% select(c("dist", "obj")) %>% `colnames<-`(c("avg_sil_width", "obj"))
  
  df <- merge(snpdf, shdf, by = "obj", all.x = T, all.y = F)
  dist_list[[p]] <- df
  print(p)
}

saveRDS(sh_list, "~/CASCADEpaper/paper/Fig3_genomics/sh_list.Rds")
saveRDS(dist_list, "~/CASCADEpaper/paper/Fig3_genomics/dist_list.Rds")

###### plot--------
dist_list <- readRDS("~/CASCADEpaper/paper/Fig3_genomics/dist_list.Rds")
distdf <- rbindlist(dist_list, idcol = "patient")

ggplot(distdf, aes(x = log10(SNP_dis), y = as.numeric(avg_sil_width), color = gsub("00", "", patient))) + geom_point() +
  theme_classic(base_size = 18)+
  scale_color_manual(values = patient_cols)+
  facet_wrap(.~patient, scales = "free") +
  labs(y = "average silhouette distance", x = "log10(SNP distance)", color = "patient")
ggsave("silhouetteVSdistSNP.pdf", width  = 6, height = 6)
