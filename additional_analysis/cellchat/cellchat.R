args <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(args[1])

#https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html#part-v-save-the-cellchat-object
library(Seurat)
library(ggplot2)
library(CellChat)
library(patchwork)
library(data.table)
library(dplyr)
library(pals)
source("~/snRNA/functions.R")

setwd("~/snRNA/cellchat")
normal_meta <- readRDS("~/snRNA/tumour_normal/normal_meta.Rds")
normal_cols <- setNames(c(polychrome(14), "lightgrey", "grey95"), c(levels(normal_meta$cell_anno), "Tumour", "Ambient cells")) 
archetype_anno <- readRDS("~/snRNA/archetype_anno.Rds")
module_cols <- setNames( c("dodgerblue3", brewer.dark2(6)[2:6]),c("AR","Inflammation", "NE1","NE2", "Cycling","Glycolysis"  ))


paths <- list.files("~/snRNA/tumour_normal", full.names = T, pattern = "doublets_removed")
samples <- gsub("_srt_new_filtering_doublets_removed.Rds", "", paths)
samples <- gsub("/home/sweng/snRNA/tumour_normal/", "", samples)
paths <- setNames(paths, samples)

sample <- samples[arg]
cat("Processing sample:", sample, "\n")

srt <- readRDS(paths[sample])

srt$annotation <- NA
srt$annotation <- as.character(normal_meta$cell_anno)[match( paste(sample, colnames(srt), sep = "_"), rownames(normal_meta))]
srt$annotation[is.na(srt$annotation)] <- "Tumour"

anno <- give_module_anno(archetype_anno[[sample]][, 1:6])
cols <- mixture_colors(anno, module_cols)
srt$archetype_anno <- anno[match(colnames(srt),rownames(archetype_anno[[sample]]))]
srt$archetype_anno[is.na(srt$archetype_anno)] <- srt$annotation[is.na(srt$archetype_anno)]
srt$archetype_anno[srt$archetype_anno == "Tumour"] <- "Ambient cells"
srt$annotation[srt$archetype_anno == "Ambient cells"] <- "Ambient cells"

srt <- subset(srt, cells = colnames(srt)[srt$archetype_anno!= "Ambient cells"])

srt <- preprocess_srt(srt)

srt$archetype_anno <- factor(srt$archetype_anno, levels = c(names(cols), names(normal_cols)))

g <- DimPlot(srt, group.by = "annotation", cols = normal_cols)+ggtitle("Cell Type")
ggsave(paste0("Dimplot/", sample, ".png"), plot = g, width = 6, height = 3, dpi = 350)
f <- DimPlot(srt, group.by = "archetype_anno", cols = c(cols, normal_cols))+ggtitle("Full annotation")
f <- f + scale_color_manual(values = c(cols, normal_cols), breaks = c(names(module_cols), names(normal_cols)))
ggsave(paste0("Dimplot/", sample, "_fullanno.png"), plot = f, width = 8, height = 3, dpi = 350)


cat("Finished preprocessing and plotting\n")
cat("Starting CellChat analysis with cell type annotation\n")
# ---- cellchat run with cell type annotation ----
cellchat <- createCellChat(object = srt, group.by = "annotation", assay = "RNA")
cellchat@DB <-CellChatDB.human
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat,do.fast = FALSE)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
#cellchat <- filterCommunication(cellchat, min.cells = names(which.min(table(srt$annotation[srt$annotation!= "Ambient cells"]))))
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
# df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
df_pathway <- subsetCommunication(cellchat, slot.name = "netP")
saveRDS(df_pathway, paste0("enrichedpathway/", sample, ".Rds"))
df_LR <- subsetCommunication(cellchat, slot.name = "net")
saveRDS(df_LR, paste0("enrichedLR/", sample, ".Rds"))
saveRDS(cellchat, paste0("cellchat_obj/", sample, ".Rds"))

# groupSize <- as.numeric(table(cellchat@idents))
# h <- ceiling(groupSize/5)
# pdf(paste0("circleplot/", sample, ".pdf"), width = 15, height = 9)
# mat <- cellchat@net$weight
# par(mfrow = c(3,5), xpd=TRUE)
# for (i in which(rownames(mat) %in% c(names(module_cols), names(normal_cols)))) {
#   mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#   mat2[i, ] <- mat[i, ]
#   netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
# }
# dev.off()


srt$archetype_anno = droplevels(srt$archetype_anno, exclude = setdiff(levels(srt$archetype_anno),unique(srt$archetype_anno)))
cat("Starting CellChat analysis with archetype annotation\n")
# ----cellchat run with archetype annotation----
cellchat <- createCellChat(object = srt, group.by = "archetype_anno", assay = "RNA")
cellchat@DB <-CellChatDB.human
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat,do.fast = FALSE)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
#cellchat <- filterCommunication(cellchat, min.cells = names(which.min(table(srt$annotation[srt$annotation!= "Ambient cells"]))))
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
# df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
df_pathway <- subsetCommunication(cellchat, slot.name = "netP")
saveRDS(df_pathway, paste0("enrichedpathway/", sample, "_archetype_anno.Rds"))
df_LR <- subsetCommunication(cellchat, slot.name = "net")
saveRDS(df_LR, paste0("enrichedLR/", sample, "_archetype_anno.Rds"))
saveRDS(cellchat, paste0("cellchat_obj/", sample, "_archetype_anno.Rds"))

# groupSize <- as.numeric(table(cellchat@idents))
# h <- ceiling(groupSize/5)
# pdf(paste0("circleplot/", sample, "_archetype_anno.pdf"), width = 15, height = 9)
# mat <- cellchat@net$weight
# par(mfrow = c(3,5), xpd=TRUE)
# for (i in which(rownames(mat) %in% c(names(module_cols), names(normal_cols)))) {
#   mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#   mat2[i, ] <- mat[i, ]
#   netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
# }
# dev.off()


