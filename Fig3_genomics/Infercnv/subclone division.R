library(infercnv)
library(dendextend)
library(phylogram)
library(dendextend)
library(dplyr)
library(Seurat)
library(tidyverse)

sample <- list.files("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered")
paths <- system("realpath /trigos_team/CASCADE/Analysis/230113_seurat_intron/*/*/*non_tumour_removed.Rds", intern = T)
paths_split <- unlist(lapply(strsplit(paths, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
paths <- setNames(paths, paths_split)
subclone_anno <- list()


# CA0027_dura_base_skull_13 -----------------------------------------------


# CA0027_dura
i = sample[1]
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
infercnv_obj <- readRDS( "20_HMM_pred.repr_intensitiesHMMi6.hmm_mode-samples.Pnorm_0.5.infercnv_obj" ))
infercnv_obj@tumor_subclusters <- NULL
plot_cnv(infercnv_obj = infercnv_obj,
         cluster_by_groups = T, 
         out_dir = "bycluster")

plot_per_group(infercnv_obj = infercnv_obj, on_references = F, write_expr_matrix = F, out_dir = "bycluster/")

dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)
tree <- cutree(dend, h  = 7)
table(tree)
tree[tree %in% c(4, 5)] <- 2
table(tree)
names(tree)[tree ==3][names(tree)[tree == 3] %in% colnames(srt)[infercnv_obj@observation_grouped_cell_indices$`7`]]
tree[names(tree)[tree ==3][names(tree)[tree == 3] %in% colnames(srt)[infercnv_obj@observation_grouped_cell_indices$`7`]]] <- 4

tree <- tree[colnames(srt)]
subclone_anno[[i]] <- tree
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)

# CA0027_paraaortic_lymph_node_1 -----------------------------------------------
i = sample[2]
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)
tree <- cutree(dend, h  = 12)

tree <- tree[colnames(srt)]
subclone_anno[[i]] <- tree
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)

# CA0027_prostate_9 -----------------------------------------------
i = sample[3]
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)
tree <- cutree(dend, k = 3)

tree <- tree[colnames(srt)]
subclone_anno[[i]] <- tree
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)

#"CA0034_liver_left_11" -----------------------------------------------
i = sample[4]
print(i)
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)

tree <- cutree(dend, k = 2)
tree <- tree[colnames(srt)]
subclone_anno[[i]] <- tree
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)


#CA0034_liver_right_8 -----------------------------------------------
i = sample[5]
print(i)
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)

tree <- cutree(dend, k = 3)
tree <- tree[colnames(srt)]
subclone_anno[[i]] <- tree
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)


#CA0034_paraaortic_lymph_node_2 -----------------------------------------------
i = sample[6]
print(i)
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)

tree <- cutree(dend, k = 3)
tree <- tree[colnames(srt)]
subclone_anno[[i]] <- tree
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)

#CA0035_bladder_2 -----------------------------------------------
i = sample[7]
print(i)
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)
d1=color_branches(dend,k=4, col = c(1:4))
plot(d1)

tree <- cutree(dend, k = 4)
tree <- tree[colnames(srt)]
table(tree)
tree[tree == 4] <- 2
subclone_anno[[i]] <- tree
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)


#CA0035_paraaortic_lymph_node_1 -----------------------------------------------
i = sample[8]
print(i)
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)
d1=color_branches(dend,k=7, col = c(1:7))
plot(d1)

tree <- cutree(dend, k = 7)
tree[tree == 2] <- 1
tree[tree %in% 4:7] <- 2

tree <- tree[colnames(srt)]
table(tree)
subclone_anno[[i]] <- tree
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)

#CA0043_liver_12 -----------------------------------------------
i = sample[9]
print(i)
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)
d1=color_branches(dend,k=4, col = c(1:4))
plot(d1)

tree <- cutree(dend, k = 2)

tree <- tree[colnames(srt)]
table(tree)
subclone_anno[[i]] <- tree
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)


#CA0043_liver_7 -----------------------------------------------
i = sample[10]
print(i)
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")

tree <- cutree(dend, k = 2)
tree <- tree[colnames(srt)]
table(tree)
subclone_anno[[i]] <- tree
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)


#CA0043_portal_lymph_node_4 -----------------------------------------------
i = sample[11]
print(i)
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")

tree <- cutree(dend, k = 4)
tree <- tree[colnames(srt)]
table(tree)
subclone_anno[[i]] <- tree
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)

#CA0046_hilar_lymph_node_5 -----------------------------------------------
i = sample[12]
print(i)
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")

tree <- cutree(dend, k = 3)
tree <- tree[colnames(srt)]
tree[tree == 2] <- 1
tree[tree == 3] <- 2
table(tree)
subclone_anno[[i]] <- tree
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)


#CA0046_liver_12 -----------------------------------------------
i = sample[13]
print(i)
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")

tree <- cutree(dend, k = 2)
tree <- tree[colnames(srt)]

table(tree)
subclone_anno[[i]] <- tree
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)


#CA0046_lung_7 -----------------------------------------------
i = sample[14]
print(i)
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
d1=color_branches(dend,k=3, col = c(1:3))
plot(d1)

tree <- cutree(dend, k = 3)
tree <- tree[colnames(srt)]

table(tree)
subclone_anno[[i]] <- tree
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)

#CA0046_paraaortic_lymph_node_15 -----------------------------------------------
i = sample[15]
print(i)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)
d1=color_branches(dend,h = 7, col = c(1:7))
plot(d1)

tree <- cutree(dend, k = 6)
tree <- tree[colnames(srt)]
table(tree)
tree[tree == 2 ] <- 1
tree[tree %in% 4:5] <- 2
tree[tree == 6] <- 3

subclone_anno[[i]] <- tree
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)

#CA0046_pelvic_lymph_node_19 -----------------------------------------------
i = sample[16]
print(i)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)

tree <- cutree(dend, k = 2)
tree <- tree[colnames(srt)]
subclone_anno[[i]] <- tree
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)

#CA0046_prostate_17 -----------------------------------------------
i = sample[17]
print(i)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)

tree <- cutree(dend, k = 2)
tree <- tree[colnames(srt)]
subclone_anno[[i]] <- tree
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)

#CA0058_hilar_50 -----------------------------------------------
i = sample[18]
print(i)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)

tree <- cutree(dend, k = 4)
tree <- tree[colnames(srt)]
subclone_anno[[i]] <- tree
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)



# CA0058_liver_29 -----------------------------------------------
i = sample[19]
print(i)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)
d1=color_branches(dend,k = 4, col = c(1:7))
plot(d1)

tree <- cutree(dend, k = 4)
tree <- tree[colnames(srt)]
tree[tree == 4] <- 3

subclone_anno[[i]] <- tree
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)


# CA0058_liver_38 -----------------------------------------------
i = sample[20]
print(i)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)
d1=color_branches(dend,h = 5.5, col = c(1:10))
plot(d1)

tree <- cutree(dend, h = 5.5)
tree <- tree[colnames(srt)]
table(tree)
tree[tree %in% 1:3] <- 1
tree[tree %in% 4:5] <- 2
tree[tree == 6] <- 3 
tree[tree == 9] <- 4
tree[tree %in% 7:8] <- 5
table(tree)

subclone_anno[[i]] <- tree
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)

# CA0076_left_rib_17-----------------------------------------------
i = sample[21]
print(i)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)


tree <- cutree(dend, k = 4)
tree <- tree[colnames(srt)]
table(tree)

subclone_anno[[i]] <- tree
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)



# CA0076_liver_right_40-----------------------------------------------
i = sample[22]
print(i)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)


tree <- cutree(dend, k = 5)
tree <- tree[colnames(srt)]
tree[tree == 5] <- 3
table(tree)

subclone_anno[[i]] <- tree
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)



# CA0076_vertebra_25-----------------------------------------------
i = sample[23]
print(i)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)


tree <- cutree(dend, k = 3)
tree <- tree[colnames(srt)]
table(tree)

subclone_anno[[i]] <- tree
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)

# CA0083_hilar_lymph_node_47 -----------------------------------------------
i = sample[24]
print(i)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)
d1=color_branches(dend,h = 5, col = c(1:10))
plot(d1)

tree <- cutree(dend, h = 5)
tree <- tree[colnames(srt)]
table(tree)

tree[tree %in% c(2:5, 7, 9)] <- 1
tree[tree %in% c(6, 8)] <- 2



subclone_anno[[i]] <- tree
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)

# CA0083_liver_49 -----------------------------------------------
i = sample[25]
print(i)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)
d1=color_branches(dend,h = 5, col = c(1:10))
plot(d1)

tree <- cutree(dend, k =3)
tree <- tree[colnames(srt)]
table(tree)

tree[tree == 3] <- 1

subclone_anno[[i]] <- tree
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)

# CA0083_lung_55-----------------------------------------------
i = sample[26]
print(i)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)

### no subclone finally 



# CA0083_paraaortic_lymph_node_12-----------------------------------------------
i = sample[27]
print(i)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)
d1=color_branches(dend,k = 6, col = c(1:10))
plot(d1)

tree <- cutree(dend, k =6)
tree <- tree[colnames(srt)]
table(tree)

tree[tree != 3] <- 1
tree[tree == 3] <- 2

subclone_anno[[i]] <- tree
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)

# CA0083_perinephric_fat_19-----------------------------------------------
i = sample[28]
print(i)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)

#no subclone
#cluster 6 maybe normal??????


# CA0090_abdomen_13-----------------------------------------------
i = sample[29]
print(i)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)

tree <- cutree(dend, k = 3)
table(tree)
tree <- tree[colnames(srt)]
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
# 2nd branch is possible junk or normal cells
## no subclone
## 

# CA0090_brain_2-----------------------------------------------
i = sample[30]
print(i)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)

tree <- cutree(dend, k = 2)
table(tree)
tree <- tree[colnames(srt)]
srt$subclone <- tree
subclone_anno[[i]] <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)

# CA0090_liver_43-----------------------------------------------
i = sample[31]
print(i)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)

tree <- cutree(dend, k = 3)
table(tree)
tree <- tree[colnames(srt)]
srt$subclone <- tree
subclone_anno[[i]] <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)
# CA0090_paraaortic_lymph_node_39-----------------------------------------------
i = sample[32]
print(i)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)

tree <- cutree(dend, k = 2)
table(tree)
tree <- tree[colnames(srt)]
srt$subclone <- tree
subclone_anno[[i]] <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)

# CA0090_porta_hepatis_lymph_node_52-----------------------------------------------
i = sample[33]
print(i)
setwd(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/",i))
srt <- readRDS(paths[i])
DimPlot(srt, label = T)
dend <- read.dendrogram("infercnv.observations_dendrogram.txt")
plot(dend)
d1=color_branches(dend,k = 13, col = polychrome(13))
plot(d1)

tree <- cutree(dend, k = 13)
table(tree)
tree <- tree[colnames(srt)]
nots <- names(tree[tree == 10])

tree <- cutree(dend, k = 4)
table(tree)
tree <- tree[colnames(srt)]
tree[tree ==2] <- 3
tree[nots] <- 2

subclone_anno[[i]] <- tree
srt$subclone <- tree
DimPlot(srt, label = T) + DimPlot(srt, group.by = "subclone",label = T)
ggsave(filename = paste0( i, ".png"),path = "~/CASCADEpaper/paper/Fig3/Infercnv/srt_plot/",device = "png",  width = 8, height = 4)

saveRDS(subclone_anno, "~/CASCADEpaper/paper/Fig3/Infercnv/subclone_anno.Rds")
