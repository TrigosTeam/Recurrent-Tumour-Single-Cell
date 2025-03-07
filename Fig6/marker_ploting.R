library(Seurat)
library(dplyr)
library(harmony)
library(data.table)
library(SingleR)
library(patchwork)
library(readxl)

normal_clean2 <- readRDS("~/CASCADEpaper/paper/Fig6/normal_clean2.Rds")

markers <- read_excel("~/CASCADEpaper/paper/Fig6/encyclopedia_table.xlsx", sheet = 1)
genes <- strsplit(markers$`Top 10 important genes from the CellTypist model`, split = ",")
genes <- lapply(genes, function(x){
  temp <- gsub("^ ", "", x)
  temp <- temp[temp %in% rownames(normal_clean2)]
})
names(genes) <- markers$`Low-hierarchy cell types`

genes <- genes[sapply(genes, length)>0]

pdf("~/CASCADEpaper/paper/Fig6/normal_marker_encyclopidia.pdf", width = 8, height = 10)
for (i in names(genes)){
  f <- FeaturePlot(normal_clean2, genes[[i]]) + plot_annotation(title = i)+ plot_layout(nrow = 4, ncol = 3)
  print(f)  
  print(i)
}
dev.off()

library(patchwork)
library(pals)
library(data.table)
setwd("~/CASCADEpaper/paper/Fig2/normal_cells")

normal_clean2 <- readRDS("~/CASCADEpaper/paper/Fig2/normal_cells/normal_finalv1.Rds")
markers <- readRDS("~/CASCADEpaper/paper/Fig2/normal_cells/normal_finalv1_marker.Rds")
sub_marker <- markers %>% filter(avg_log2FC > 0 & p_val_adj < 0.05 )#& abs(pct.1-pct.2)>0.5
cell_type_marker <- fread("rna_single_cell_type.tsv")


marker_ref <- cell_type_marker %>% filter(`Gene name` %in% sub_marker$gene)

top_marker <- marker_ref %>% group_by(`Cell type`) %>% filter(nTPM %in% head(sort(nTPM, decreasing = T)))
test <- lapply(split(top_marker, top_marker$`Cell type`), function(x){
  temp <- apply(x, 1, function(y){
    if(sum(sub_marker$gene %in% y[2])>0){
      
      return(paste(sub_marker$cluster[sub_marker$gene %in% y[2]], collapse = ","))
    }
  })
  if( class(temp) == "list"){
    if(sum(sapply(temp, is.null))>0){
      temp[sapply(temp, is.null)] <- "no"
    }
  }
  
  x$cluster <- unlist(temp)
  return(x)
})


pdf("~/CASCADEpaper/paper/Fig2/normal_cells/normal_marker.pdf", width = 8, height = 8)
for ( i in names(test)){
  f <- FeaturePlot(normal_clean2, features = test[[i]]$`Gene name`, cols = rev(brewer.rdylbu(5))) + plot_annotation(title  = i)
  print(f)
  f <- FeaturePlot(normal_clean2, features = test[[i]]$`Gene name`) + plot_annotation(title  = i)
  print(f)
}
dev.off()
