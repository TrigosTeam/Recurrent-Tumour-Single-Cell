library(Seurat)
library(ggplot2)
library(CellChat)
library(patchwork)
library(data.table)
library(dplyr)
library(circlize)
library(ComplexHeatmap)
library(pals)
library(ggpubr)
source("~/snRNA/functions.R")
setwd("~/snRNA/cellchat")

# calculate the cell number across cell types
final_normalsrt <- readRDS("~/snRNA/tumour_normal/final_normalsrt.Rds")
normal_col <- setNames(brewer.set1(14), levels(final_normalsrt$cell_anno))
final_normalsrt$group <- ifelse( final_normalsrt$patient %in% c("CA0090", "CA0046")|grepl("CA0058_liver", final_normalsrt$sample), "NE+", "AR+")
df <- final_normalsrt@meta.data %>% group_by(group, sample, cell_anno) %>% summarise(n = n())
my_comparisons <- list(c("AR+", "NE+"))
ggplot(df, aes(x = group, y = n, fill = cell_anno)) +
  geom_boxplot() +
  facet_wrap(.~cell_anno, scales = "free_y")+
  theme_classic()+
  scale_fill_manual(values = normal_col)+
  theme(legend.position = "none") +
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     label.x = 1.5,
                     label.y.npc = 0.8)+
  labs(y = "Number of cell (n)", x = NULL) 
ggsave("celltype number in different group.pdf", height = 4, width = 6)


# enriched pathway matrix -------
## with general cell type annotaions 
### Tumour as receiver
mats <- lapply(list.files("enrichedpathway", full.names = T, pattern = "[0-9+].Rds"), readRDS)
names(mats) <- gsub(".Rds", "", list.files("enrichedpathway", pattern = "[0-9+].Rds"))
mats <- lapply(mats, function(df) df %>% filter(source!=target) %>% filter(target == "Tumour"))
all_cols <- unique(unlist(lapply(mats, function(x) unique(x$source))))
all_rows <- unique(unlist(lapply(mats, function(x) unique(x$pathway_name))))
mat <- generate_feature_mat(mats, all_cols, all_rows, formula = "pathway_name~source")

meta <- readRDS("~/snRNA/integration/tumour_only_meta.Rds")

module_score <- meta[, c( paste(names(clean_module), "Module"))]
module_score <- split(module_score, meta$sample)
avg_exp <- lapply(module_score, colMeans)
anno_col_fun = colorRamp2(c(-0.5, 1), c( "white", "darkblue"))

groups <- ifelse( substr(names(mats),1,6) %in% c("CA0090", "CA0046"), "NE", 
                  ifelse( substr(names(mats),1,6) %in% c("CA0027", "CA0058"), "Mixed", "AD"))

ha2 = HeatmapAnnotation(patient = gsub("00", "", substr(names(mats),1,6)), 
                        site = sites, 
                        group = groups, 
                        AR_Module = as.numeric(sapply(avg_exp[colnames(opmat)], '[',"AR Module")), 
                        Inflammation_Module = as.numeric(sapply(avg_exp[colnames(opmat)], '[',"Inflammation Module")), 
                        NE1_Module = as.numeric(sapply(avg_exp[colnames(opmat)], '[',"NE1 Module")), 
                        NE2_Module = as.numeric(sapply(avg_exp[colnames(opmat)], '[',"NE2 Module")), 
                        Cycling_Modle  = as.numeric(sapply(avg_exp[colnames(opmat)], '[',"Cycling Module")), 
                        Glycolysis_Module = as.numeric(sapply(avg_exp[colnames(opmat)], '[',"Glycolysis Module")), 
                        col = list(patient = setNames(brewer.set1(9), unique(gsub("00", "", substr(names(mats),1,6)))),
                                   site = setNames(brewer.set3(9), unique(sites)), 
                                   group = setNames(c("indianred", "seagreen","steelblue2"), c("AD", "Mixed", "NE")),
                                   AR_Module = colorRamp2(c(-0.5, 1), c( "white", "dodgerblue3")),
                                   Inflammation_Module = colorRamp2(c(0, 0.6), c( "white", "#D95F02")),
                                   NE1_Module = colorRamp2(c(-0.2, 0.6), c( "white", "#7570B3")), 
                                   NE2_Module = colorRamp2(c(0, 0.6), c( "white", "#E7298A")),
                                   Cycling_Modle  =colorRamp2(c(-0.05, 0.15), c( "white", "#66A61E")),
                                   Glycolysis_Module =colorRamp2(c(-0.4, 0.4), c( "white", "#E6AB02"))))
ha = HeatmapAnnotation(patient = substr(names(mats),1,6), 
                       site = sites, 
                       group = groups, 
                       col = list(patient = setNames(brewer.set1(9), unique(substr(names(mats),1,6))),
                                  site = setNames(brewer.set3(9), unique(sites)), 
                                  group = setNames(c("indianred", "seagreen","steelblue2"), c("AD", "Mixed", "NE"))))

htmat <- t(mat[, colMeans(mat)>mean(colMeans(mat)[colMeans(mat)>0]) ])
celltypes= sapply(strsplit(rownames(htmat), split = ":"), '[',1)

ra2 = rowAnnotation(celltype = celltypes[ind], 
                    annotation_name_side = "top", 
                    col = list(celltype =normal_col))
rns <- sapply(strsplit(rownames(htmat), split = ":"), "[", 2)[ind]
ht2 <- Heatmap(ifelse(htmat>0, "enriched", "")[ind,],  
               col = c("enriched" = "red"),
               name = "cellchat prob",
               show_column_names = F,
               show_row_names = T, 
               row_labels = rns,
               column_order = order(as.numeric(sapply(avg_exp[colnames(htmat)], '[',"AR Module"))),
               top_annotation = ha2, 
               left_annotation = ra2)
pdf("heatmap/tumour_as_receiver_pathway_bycelltype.pdf", width = 12, height = 12)
draw(ht2, column_title = "Tumour as receiver", column_title_gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()


htmat <- t(mat[, colMeans(mat)>mean(colMeans(mat)[colMeans(mat)>0]) ])
celltypes= sapply(strsplit(rownames(htmat), split = "-"), '[',1)
ra = rowAnnotation(celltype = celltypes, 
                   annotation_name_side = "top", 
                   col = list(celltype = setNames(brewer.accent(n_distinct(celltypes)), unique(celltypes))))
ht <- Heatmap(htmat, 
        col = colorRamp2(breaks = c(0, 0.1, max(mat)), colors = c("steelblue3", "white", "tomato")),
        name = "cellchat prob",
        show_column_names = F,
        show_row_names = F,
        top_annotation = ha, 
        left_annotation = ra)
pdf("heatmap/tumour_as_receiver_pathway.pdf", width = 6, height = 6)
draw(ht, column_title = "Tumour as receiver", column_title_gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()

g1 <- pcaplot(feature_mat = mat)+ plot_annotation(title = "Tumour as receiver")+plot_layout(ncol = 2)
ggsave("pca/tumour_as_receiver_pathway.pdf", width = 6, height = 6)

### Tumour as sender
mats <- lapply(list.files("enrichedpathway", full.names = T, pattern = "[0-9+].Rds"), readRDS)
names(mats) <- gsub(".Rds", "", list.files("enrichedpathway", pattern = "[0-9+].Rds"))
mats <- lapply(mats, function(df) df %>% filter(source!=target) %>% filter(source == "Tumour"))
all_cols <- unique(unlist(lapply(mats, function(x) unique(x$target))))
all_rows <- unique(unlist(lapply(mats, function(x) unique(x$pathway_name))))

mat <- generate_feature_mat(mats, all_cols, all_rows, formula = "pathway_name~target")
htmat <- t(mat[, colMeans(mat)>mean(colMeans(mat)[colMeans(mat)>0]) ])
celltypes= sapply(strsplit(rownames(htmat), split = "-"), '[',1)
ra = rowAnnotation(celltype = celltypes, 
                   annotation_name_side = "top", 
                   col = list(celltype = setNames(brewer.accent(n_distinct(celltypes)), unique(celltypes))))
ht <- Heatmap(htmat, 
              col = colorRamp2(breaks = c(0, 0.1, max(mat)), colors = c("steelblue3", "white", "tomato")),
              name = "cellchat prob",
              show_column_names = F,
              show_row_names = F,
              top_annotation = ha, 
              left_annotation = ra)
pdf("heatmap/tumour_as_sender_pathway.pdf", width = 6, height = 6)
draw(ht, column_title = "Tumour as sender", column_title_gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()
g2 <- pcaplot(feature_mat = mat)+ plot_annotation(title = "Tumour as sender")+plot_layout(ncol = 2)
ggsave("pca/tumour_as_sender_pathway.pdf", width = 6, height = 6)

# enriched LR pair -------
## with general cell type annotaions 
### Tumour as receiver




mats <- lapply(list.files("enrichedLR", full.names = T, pattern = "[0-9+].Rds"), readRDS)
names(mats) <- gsub(".Rds", "", list.files("enrichedLR", pattern = "[0-9+].Rds"))
mats <- lapply(mats, function(df) df %>% filter(source!=target) %>% filter(target == "Tumour"))
all_cols <- unique(unlist(lapply(mats, function(x) unique(as.character(x$source)))))
all_rows <- unique(unlist(lapply(mats, function(x) unique(as.character(x$interaction_name)))))



module_score <- meta[, c( paste(names(clean_module), "Module"))]
module_score <- split(module_score, meta$sample)
avg_exp <- lapply(module_score, colMeans)
anno_col_fun = colorRamp2(c(-0.5, 1), c( "white", "darkblue"))

mat <- generate_feature_mat(mats, all_cols, all_rows, formula = "interaction_name~source")
htmat <- t(mat[, colMeans(mat)>mean(colMeans(mat)[colMeans(mat)>0]) ])
celltypes= sapply(strsplit(rownames(htmat), split = ":"), '[',1)
opmat <- ifelse(htmat>0, "enriched", "")

ha2 = HeatmapAnnotation(patient = patients, 
                        site = sites, 
                        group = groups, 
                        AR_Module = as.numeric(sapply(avg_exp[colnames(opmat)], '[',"AR Module")), 
                        Inflammation_Module = as.numeric(sapply(avg_exp[colnames(opmat)], '[',"Inflammation Module")), 
                        NE1_Module = as.numeric(sapply(avg_exp[colnames(opmat)], '[',"NE1 Module")), 
                        NE2_Module = as.numeric(sapply(avg_exp[colnames(opmat)], '[',"NE2 Module")), 
                        Cycling_Modle  = as.numeric(sapply(avg_exp[colnames(opmat)], '[',"Cycling Module")), 
                        Glycolysis_Module = as.numeric(sapply(avg_exp[colnames(opmat)], '[',"Glycolysis Module")), 
                        col = list(patient = setNames(brewer.set1(9), unique(patients)),
                                   site = setNames(brewer.set3(9), unique(sites)), 
                                   group = setNames(c("indianred", "seagreen","steelblue2"), c("AD", "Mixed", "NE")),
                                   AR_Module = colorRamp2(c(-0.5, 1), c( "white", "dodgerblue3")),
                                   Inflammation_Module = colorRamp2(c(0, 0.6), c( "white", "#D95F02")),
                                   NE1_Module = colorRamp2(c(-0.2, 0.6), c( "white", "#7570B3")), 
                                   NE2_Module = colorRamp2(c(0, 0.6), c( "white", "#E7298A")),
                                   Cycling_Modle  =colorRamp2(c(-0.05, 0.15), c( "white", "#66A61E")),
                                   Glycolysis_Module =colorRamp2(c(-0.4, 0.4), c( "white", "#E6AB02"))))
alter_fun = list(
  background = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#CCCCCC", col = NA)),
  # red rectangles
  enriched = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "red", col = NA))
)

ind <- rowSums(htmat>0)>5

ra2 = rowAnnotation(celltype = celltypes[ind], 
                   annotation_name_side = "top", 
                   col = list(celltype =normal_col))
rns <- sapply(strsplit(rownames(htmat), split = ":"), "[", 2)[ind]
ht <- Heatmap(htmat[ind, ],
              show_column_names = F,
              clustering_distance_rows = "binary",
              clustering_method_rows = "ward.D2",
              column_order = order(avg_expression$FOLH1[match(colnames(htmat), avg_expression$sample)]))

row_ord <- ifelse(htmat>0 , 1, 0)
row_neg_score <- rowSums(row_ord[ind, order(as.numeric(sapply(avg_exp[colnames(opmat)], '[',"AR Module")))][, 1:13])
row_pos_score <- rowSums(row_ord[ind, order(as.numeric(sapply(avg_exp[colnames(opmat)], '[',"AR Module")))][, 14:33])
row_group <- ifelse(row_pos_score>0&row_neg_score==0,"AR", ifelse(row_pos_score==0&row_neg_score>0, "NE", "Mixed"))

op <- oncoPrint(opmat[ind, ], alter_fun = alter_fun, col = c(enriched = "red"),
                column_order = order(as.numeric(sapply(avg_exp[colnames(htmat)], '[',"AR Module"))),
                top_annotation = ha2, 
                left_annotation = ra2,
                right_annotation = NULL,
                row_split = row_group,
                heatmap_legend_param = list(title = "CellChat\nResults", direction = "horizontal"))

pdf("heatmap/tumour_as_receiver_LR_bycelltype.pdf", width = 12, height = 12)
draw(op, merge_legend = T)
# draw(ht2, column_title = "Tumour as receiver", column_title_gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()

g1 <- pcaplot(feature_mat = mat)+ plot_annotation(title = "Tumour as receiver")+plot_layout(ncol = 2)
ggsave("pca/tumour_as_receiver_LR.pdf", width = 6, height = 6)

all_LR <- lapply(mats, function(x) unique(as.character(x$interaction_name)))
rmat <- matrix(0, ncol =33, nrow = 270, dimnames = list(unique(unlist(all_LR)), names(mats)))
for(i in names(mats)){
  rmat[all_LR[[i]], i] <- 1
}
row_ord <- ifelse(rmat>0 , 1, 0)
row_neg_score <- rowSums(row_ord[rowSums(rmat>0)>5, order(as.numeric(sapply(avg_exp[colnames(opmat)], '[',"AR Module")))][, 1:13])
row_pos_score <- rowSums(row_ord[rowSums(rmat>0)>5, order(as.numeric(sapply(avg_exp[colnames(opmat)], '[',"AR Module")))][, 14:33])
row_group <- ifelse(row_pos_score>0&row_neg_score==0,"AR", ifelse(row_pos_score==0&row_neg_score>0, "NE", "Mixed"))

ht2 <- Heatmap(ifelse(rmat>0, "enriched", "")[rowSums(rmat>0)>5,],  
               col = c("enriched" = "red"),
               name = "cellchat prob",
               show_column_names = F,
               show_row_names = T, 
               row_split = row_group,
               column_order = order(as.numeric(sapply(avg_exp[colnames(rmat)], '[',"AR Module"))),
               top_annotation = ha2)
pdf("heatmap/tumour_as_receiver_LR.pdf", width = 12, height = 12)
draw(ht2, column_title = "Tumour as receiver", column_title_gp = gpar(fontsize = 14, fontface = "bold"), merge_legend = T)
dev.off()


ha3 = HeatmapAnnotation(patient = patients, 
                        site = sites, 
                        group = groups, 
                        FOLH1 = avg_expression$FOLH1[match(colnames(htmat), avg_expression$sample)], 
                        AR = avg_expression$AR[match(colnames(htmat), avg_expression$sample)], 
                        col = list(patient = setNames(brewer.set1(9), unique(patients)),
                                   site = setNames(brewer.set3(9), unique(sites)), 
                                   group = setNames(c("indianred", "seagreen","steelblue2"), c("AD", "Mixed", "NE")),
                                   AR = colorRamp2(c(0, 3), c( "white", "dodgerblue3")),
                                   FOLH1 = colorRamp2(c(0, 3), c( "white", "#E7298A"))))

op2 <- oncoPrint(opmat[ind, ], alter_fun = alter_fun, 
                row_labels = rns,
                column_order = order(avg_expression$FOLH1[match(colnames(htmat), avg_expression$sample)]),
                top_annotation = ha3, 
                left_annotation = ra2,
                right_annotation = NULL,
                row_order = row_order(ht),
                heatmap_legend_param = list(title = "CellChat\nResults", direction = "horizontal"))

pdf("heatmap/tumour_as_receiver_LR_bycelltype_FOLH1rank.pdf", width = 12, height = 12)
draw(op2, merge_legend = T)
dev.off()

#### focus on enrichemnt of receptors on tumour
all_receptor<- lapply(mats, function(x) unique(as.character(x$receptor)))
avg_prob <- lapply(mats, function(x) {
  df <- x %>% group_by(receptor) %>% summarise(mean = mean(prob)) 
  return(setNames(df$mean, df$receptor))
})

max_prob <- lapply(mats, function(x) {
  df <- x %>% group_by(receptor) %>% summarise(mean = max(prob)) 
  return(setNames(df$mean, df$receptor))
})

rmat <- matrix(0, ncol =33, nrow = 101, dimnames = list(unique(unlist(all_receptor)), names(mats)))
for(i in names(mats)){
  rmat[names(max_prob[[i]]), i] <- avg_prob[[i]]
}
opmat <- ifelse(rmat>0 , "enriched", "")

ha = HeatmapAnnotation(patient = substr(names(mats),1,6), 
                       site = sites, 
                       group = groups, 
                       col = list(patient = setNames(brewer.set1(9), unique(substr(names(mats),1,6))),
                                  site = setNames(brewer.set3(9), unique(sites)), 
                                  group = setNames(c("indianred", "seagreen","steelblue2"), c("AD", "Mixed", "NE"))))
ht <- Heatmap(opmat[rowSums(rmat)>3, ], col = c(enriched = "red"),
              top_annotation = ha,
              show_column_names = F,
              #column_order = order(as.numeric(sapply(avg_exp[colnames(opmat)], '[',"AR Module"))),
              #column_order = order(colSums(rmat), decreasing = T),
              heatmap_legend_param = list(title = "CellChat\nResults", direction = "horizontal"))

op <- oncoPrint(opmat[rowSums(rmat)>3, ], alter_fun = alter_fun, col = c(enriched = "red"),
          top_annotation = ha,
          column_order = column_order(ht), 
          #column_order = order(colSums(rmat), decreasing = T),
          heatmap_legend_param = list(title = "CellChat\nResults", direction = "horizontal"))
pdf("heatmap/tumour_as_receiver_receptorenriched.pdf", width = 7, height = 10)
draw(op,merge_legend = T)
dev.off()

# op2 <- oncoPrint(opmat[rowSums(rmat)>3, ], alter_fun = alter_fun, col = c(enriched = "red"),
#                 top_annotation = ha,
#                 column_split = sites,
#                 column_title = NULL,
#                 #column_order = order(colSums(rmat), decreasing = T),
#                 heatmap_legend_param = list(title = "CellChat\nResults", direction = "horizontal"))
# pdf("heatmap/tumour_as_receiver_receptorenriched_groupbysite.pdf", width = 7, height = 10)
# draw(op2,merge_legend = T)
# dev.off()
# 
# op3 <- oncoPrint(opmat[rowSums(rmat)>3, ], alter_fun = alter_fun, col = c(enriched = "red"),
#                  top_annotation = ha,
#                  column_split =groups,
#                  column_title = NULL,
#                  row_order = row_order(ht), 
#                  #column_order = order(colSums(rmat), decreasing = T),
#                  heatmap_legend_param = list(title = "CellChat\nResults", direction = "horizontal"))
# pdf("heatmap/tumour_as_receiver_receptorenriched_groupbydiseasegroup.pdf", width = 7, height = 10)
# draw(op3,merge_legend = T)
# dev.off()

ind <- rowSums(rmat>0)>5
row_ord <- ifelse(rmat>0 , 1, 0)
ht <- Heatmap(row_ord[ind, ],
              show_column_names = F,
              clustering_distance_rows = "binary",
              clustering_method_rows = "ward.D2",
              column_order =  order(avg_expression$FOLH1[match(colnames(opmat), avg_expression$sample)]))

row_neg_score <- rowSums(row_ord[ind, order(as.numeric(sapply(avg_exp[colnames(opmat)], '[',"AR Module")))][, 1:13])
row_pos_score <- rowSums(row_ord[ind, order(as.numeric(sapply(avg_exp[colnames(opmat)], '[',"AR Module")))][, 14:33])

row_group <- ifelse(row_pos_score>0&row_neg_score==0,"AR", ifelse(row_pos_score==0&row_neg_score>0, "NE", "Mixed"))
op <- oncoPrint(opmat[ind, ], alter_fun = alter_fun, col = c(enriched = "red"),
                top_annotation = ha2,
                row_split = row_group,
                right_annotation = NULL,
                column_order = order(as.numeric(sapply(avg_exp[colnames(opmat)], '[',"AR Module"))),
                heatmap_legend_param = list(title = "CellChat\nResults"))

pdf("heatmap/tumour_as_receiver_receptorenriched_bymodulescore.pdf", width = 10, height = 10)
draw(op, column_title = "Receptors when tumour as receiver", column_title_gp = gpar(fontsize = 14, fontface = "bold"),merge_legend = T)
dev.off()

op2 <- oncoPrint(opmat[ind, ], alter_fun = alter_fun, 
                 column_order = order(avg_expression$FOLH1[match(colnames(opmat), avg_expression$sample)]),
                 top_annotation = ha3, 
                 right_annotation = NULL,
                 row_order = row_order(ht),
                 heatmap_legend_param = list(title = "CellChat\nResults", direction = "horizontal"))

pdf("heatmap/tumour_as_receiver_receptorenriched_bymodulescore_FOLH1rank.pdf", width = 12, height = 12)
draw(op2, merge_legend = T)
dev.off()


#### focus on enrichemnt of ligand from normal cells
all_ligand<- lapply(mats, function(x) unique(as.character(x$ligand)))
rmat <- matrix(0, ncol =33, nrow = 134, dimnames = list(unique(unlist(all_ligand)), names(mats)))
for(i in names(mats)){
  rmat[all_ligand[[i]], i] <- 1
}
opmat <- ifelse(rmat ==1 , "enriched", "")

row_ord <- ifelse(rmat>0 , 1, 0)
row_neg_score <- rowSums(row_ord[ind, order(as.numeric(sapply(avg_exp[colnames(opmat)], '[',"AR Module")))][, 1:13])
row_pos_score <- rowSums(row_ord[ind, order(as.numeric(sapply(avg_exp[colnames(opmat)], '[',"AR Module")))][, 14:33])

row_group <- ifelse(row_pos_score>0&row_neg_score==0,"AR", ifelse(row_pos_score==0&row_neg_score>0, "NE", "Mixed"))
op <- oncoPrint(opmat[ind, ], alter_fun = alter_fun, col = c(enriched = "red"),
                top_annotation = ha2,
                row_split = row_group,
                right_annotation = NULL,
                column_order = order(as.numeric(sapply(avg_exp[colnames(opmat)], '[',"AR Module"))),
                heatmap_legend_param = list(title = "CellChat\nResults"))
pdf("heatmap/tumour_as_receiver_ligandenriched_bymodulescore.pdf", width = 10, height = 10)
draw(op, column_title = "Ligands when tumour as receiver", column_title_gp = gpar(fontsize = 14, fontface = "bold"),merge_legend = T)
dev.off()






### Tumour as sender
mats <- lapply(list.files("enrichedLR", full.names = T, pattern = "[0-9+].Rds"), readRDS)
names(mats) <- gsub(".Rds", "", list.files("enrichedLR", pattern = "[0-9+].Rds"))
mats <- lapply(mats, function(df) df %>% filter(source!=target) %>% filter(source == "Tumour"))
all_cols <- unique(unlist(lapply(mats, function(x) unique(as.character(x$target)))))
all_rows <- unique(unlist(lapply(mats, function(x) unique(as.character(x$interaction_name)))))
mat <- generate_feature_mat(mats, all_cols, all_rows, formula = "interaction_name~target")
htmat <- t(mat[, colMeans(mat)>mean(colMeans(mat)[colMeans(mat)>0]) ])
celltypes= sapply(strsplit(rownames(htmat), split = "-"), '[',1)
ra = rowAnnotation(celltype = celltypes, 
                   annotation_name_side = "top", 
                   col = list(celltype = setNames(brewer.accent(n_distinct(celltypes)), unique(celltypes))))
ht <- Heatmap(htmat, 
              col = colorRamp2(breaks = c(0, 0.1, max(mat)), colors = c("steelblue3", "white", "tomato")),
              name = "cellchat prob",
              show_column_names = F,
              show_row_names = F,
              top_annotation = ha, 
              left_annotation = ra)
pdf("heatmap/tumour_as_sender_LR.pdf", width = 6, height = 6)
draw(ht, column_title = "Tumour as sender", column_title_gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()
g2 <- pcaplot(feature_mat = mat)+ plot_annotation(title = "Tumour as sender")+plot_layout(ncol = 2)
ggsave("pca/tumour_as_sender_LR.pdf", width = 6, height = 6)


all_ligand<- lapply(mats, function(x) unique(as.character(x$ligand)))
rmat <- matrix(0, ncol =33, nrow = 81, dimnames = list(unique(unlist(all_ligand)), names(mats)))
for(i in names(mats)){
  rmat[all_ligand[[i]], i] <- 1
}
opmat <- ifelse(rmat ==1 , "enriched", "")
alter_fun = list(
  background = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#CCCCCC", col = NA)),
  # red rectangles
  enriched = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "red", col = NA))
)
op <- oncoPrint(opmat, alter_fun = alter_fun, col = c(enriched = "red"),
                top_annotation = ha,
                #column_order = order(colSums(rmat), decreasing = T),
                heatmap_legend_param = list(title = "CellChat\nResults", direction = "horizontal"))
pdf("heatmap/tumour_as_sender_ligandenriched.pdf", width = 7, height = 10)
draw(op,merge_legend = T)
dev.off()