library(Seurat)
library(dplyr)
library(scales)
library(circlize)
library(dendextend)
library(heatmaply)
library(htmlwidgets)
library(ComplexHeatmap)

make_count_mat <- function(df,name_var, count_var){
  df <- as.data.frame(df)
  mat <- matrix(0, nrow = length(unique(df[, name_var])),
                ncol = length(unique(df[, name_var])), 
                dimnames = list(unique(df[, name_var]), unique(df[, name_var])))
  for (i in rownames(mat)){
    for(n in colnames(mat)){
      var1 <- df[df[,name_var] ==i, count_var]
      var2 <- df[df[,name_var] ==n, count_var]
      mat[i, n ] <-length(intersect(var1, var2))/min(length(var1),length(var2))
    }
  }
  return(mat)
}


# enriced genes df --------------------------------------------------------

enrichlist <- list.files("~/regulon/new_regulon/archetype/2024_06/2_enrichment", full.names = T) %>% lapply(., readRDS)
names(enrichlist) <- gsub(".Rds", "", list.files("~/regulon/new_regulon/archetype/2024_06/2_enrichment"))


enriched_genes <- lapply(enrichlist, function(x) x[["lab"]][["enriched_genes"]])
enriched_genes <- rbindlist(enriched_genes, idcol = "sample")
enriched_genes$arch_group <- paste(enriched_genes$sample, enriched_genes$arch_name, sep = "_")
enriched_genes$genes <- as.character(enriched_genes$genes)


# DEG markers -------------------------------------------------------------

markerlist <- list.files("~/regulon/new_regulon/archetype/2024_06/2_marker", full.names = T) %>% lapply(., readRDS)
names(markerlist) <- gsub(".Rds", "", list.files("~/regulon/new_regulon/archetype/2024_06/2_marker"))
markers <- rbindlist(markerlist, idcol = "sample")
markers$arch_group <- paste(markers$sample, markers$cluster, sep = "_")


arch_group_gene <- sapply(unique(markers$arch_group), function(arch){
  markergene <- markers %>% filter(arch_group == arch) %>% 
    filter(p_val_adj < 0.05 & avg_log2FC  >0) %>%  pull(gene)
  enrichgene <- enriched_genes %>% filter(arch_group == arch) %>%
    filter(p < 0.05 & mean_diff>0) %>% pull(genes)
  return(intersect(enrichgene, markergene))
})


library(dendextend)

derive_groupgene <- function(arch_group_gene, mat, k, cutoff, show_name = F){
  ht <- draw(Heatmap(mat))
  draw(Heatmap(mat, column_split = k, row_split = k, show_row_names = show_name, show_column_names = show_name, name = "prop"))
  tree <- cutree(column_dend(ht), k = k)
  group <- split(names(tree), tree)
  groupgene <- lapply(group, function(x){
    names(which(table(unlist(arch_group_gene[x]))>=(length(x)*cutoff)))
  })
  groupgene <- groupgene[sapply(groupgene, length)>0]
  return(groupgene)
}

test_groupgene <- function(groupgene){
mat <- sapply(groupgene, function(x){
  sapply(groupgene, function(y){
    length(intersect(x, y))/min(length(x), length(y))
  })
})
draw(Heatmap(mat, col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), name = "prop"))
return(mat)
}

arch_group_gene <- arch_group_gene[sapply(arch_group_gene, length)>10]
mat <- sapply(arch_group_gene, function(x){
  sapply(arch_group_gene, function(y){
    return( length(intersect(x, y ))/min(length(x), length(y)))
  })
})
plot(density(mat))
abline(v = 0.25)
Heatmap(mat)
groupgene <- derive_groupgene(mat, k = 9, 0.25)

## remove 3 as it comes from only one archtype (CA27_dura_13)
## 4, 6
## 1 take all, setdiff(8 and 1 )
# 5 take all 

groupgene[[4]] <- unique(unlist(groupgene[c(4,6)]))
groupgene <- groupgene[-c(3, 6)]
rmgenes <- setdiff(names(which(table(unlist(groupgene))>1)), groupgene[[7]])


srt <- AddModuleScore(srt, groupgene, name = "Group")
FeaturePlot(srt, paste0("Group", 1:7), cols = rev(RdYlBu(10)))
saveRDS(groupgene, "module_marker.Rds")


# second round split by patholgoy -----------------------------------------



NE_mat <- mat[grep("CA0046|CA0090|CA0058_liver", colnames(mat)),grep("CA0046|CA0090|CA0058_liver", colnames(mat))]
plot(density(NE_mat))
abline(v = 0.35)
NE_groupgene <- derive_groupgene(arch_group_gene,NE_mat, k = 8, 0.35)
names(NE_groupgene) <- paste0("NE", seq(length(NE_groupgene)))


AD_mat <- mat[grep("CA0034|CA0035|CA0043|CA0076|CA0083", colnames(mat)),grep("CA0034|CA0035|CA0043|CA0076|CA0083", colnames(mat))]
plot(density(AD_mat))
abline(v = 0.3)
Heatmap(AD_mat)
AD_groupgene <- derive_groupgene(arch_group_gene,AD_mat, k = 5,  0.3)
names(AD_groupgene) <- paste0("AD", seq(length(AD_groupgene)))

mix_mat <-mat[grep("CA0027|CA0058_hilar", colnames(mat)),grep("CA0027|CA0058_hilar", colnames(mat))]
plot(density(mix_mat))
abline(v = 0.45 )
Heatmap(mix_mat)
mix_groupgene <- derive_groupgene(arch_group_gene,mix_mat, k  = 5, 0.45)
names(mix_groupgene) <- paste0("mix", seq(length(mix_groupgene)))

testmat <- test_groupgene(c(NE_groupgene, AD_groupgene, mix_groupgene))
plot(density(testmat))
abline(v = 0.3)
Heatmap(ifelse(testmat>0.25, 1, 0))
Heatmap(testmat)

allgroup <- c(NE_groupgene, AD_groupgene, mix_groupgene)
allgroupgene <- derive_groupgene(allgroup, testmat, k=6, 0.3, show_name = T)


unique(unlist(allgroup[c("AD2", "NE6", "NE2")]))
unique(unlist(allgroup[c("mix1", "NE1")]))
unique(unlist(allgroup[c("mix2", "NE3")]))




srt <- AddModuleScore(srt, groupgene, name = "Module")
FeaturePlot(srt, paste0("Module", 1:9))
srt <- AddModuleScore(srt, NE_groupgene, name = "NE")
FeaturePlot(srt, paste0("NE", 1:10))
srt <- AddModuleScore(srt, AD_groupgene, name = "AD")
srt <- AddModuleScore(srt, mix_groupgene, name = "mix")
FeaturePlot(srt, c(paste0("AD", 1:3),paste0("mix", 1:4)))

srt <- AddModuleScore(srt, allgroupgene, name = "all")
FeaturePlot(srt, paste0("all", 1:6))
s

CA0090 <- AddModuleScore(CA0090, allgroupgene, name = "all")
FeaturePlot(CA0090, paste0("all", 1:6))


temp <- sapply(allgroupgene, function(x){
  sapply(allgroupgene, function(y){
    length(intersect(x, y))/min(length(x), length(y))
  })
})
Heatmap(temp)

View(table(unlist(allgroupgene)))
intersect(allgroupgene[[4]], allgroupgene[[3]])

purifing <- function(allgroupgene, ind){
  setdiff(allgroupgene[[ind]], unlist(allgroupgene[setdiff(seq(length(allgroupgene)),ind)]))
}

AR <- purifing(allgroupgene, 5)
inflammation <- purifing(allgroupgene, 6)
NE1 <- purifing(allgroupgene, 3)
NE2 <- purifing(allgroupgene, 4)
cycling <- purifing(allgroupgene, 1)
hypoxia <- purifing(allgroupgene, 2)

clean_module <-  list(AR = AR, inflammation = inflammation, NE1 = NE1, NE2 = NE2, cycling = cycling, hypoxia = hypoxia )
temp <- sapply(clean_module, function(x){
  sapply(clean_module, function(y){
    length(intersect(x, y))/min(length(x), length(y))
  })
})
Heatmap(temp)
srt<- AddModuleScore(srt, clean_module, name = "clean")
FeaturePlot(srt, paste0("clean", 1:6))
saveRDS(clean_module, "clean_module")


# seconde round clean after enrichment result seen ------------------------


module_marker <- readRDS("~/regulon/new_regulon/archetype/2024_06/module_marker.Rds")

# clean module 3 and 4, 7, 3 is androgen respond, 4, 7 are NE
AR_related <- setdiff(module_marker[[3]], c(module_marker[[4]], module_marker[[7]]))
inflammation <- setdiff(module_marker[[2]], c(module_marker[[1]]))
NE_transition <- setdiff(module_marker[[4]], module_marker[[1]])
