library(Seurat)
library(patchwork)
library(tidyverse)
library(pals)
library(ggridges)

#include new sample
srt <- readRDS("~/integration/2024_06/tumor_only_ruvsrt.Rds")
clean_module <- readRDS("~/CASCADEpaper/paper/Fig5_archetype/clean_module.Rds")
names(clean_module) <- c("AR", "Inflammation", "NE1", "NE2", "Cycling", "Glycolysis")

meta <- readRDS("~/CASCADEpaper/paper/Fig5_archetype/5_integrated_final_module_meta.Rds")
cells <- apply(meta[, grep("Module", colnames(meta))], 2, function(x){
  bins <- cut(x, breaks = 10)
  group <- bins %in% levels(bins)[5:10]
  return(setNames(group, rownames(meta)))
})
colnames(cells) <- paste0(colnames(cells), "_group") 
meta <- cbind(meta, cells)
meta$group <- "low expression"
ind_high <- apply(meta[, paste0("Module", 1:6, "_group")], 1, function(x) if(sum(x) > 1) paste(names(clean_module)[x], collapse = "&"))

for(i in paste0("Module", 1:6, "_group")){
  meta$group[meta[, i]] <- paste(i, "high")
}
meta$group[!sapply(ind_high, is.null)] <- unlist(ind_high)
meta$group <- ifelse(grepl("&", meta$group), "Transition", meta$group)



meta$group <- factor(meta$group, levels = 
c(paste0("Module", 1:6, "_group", " high"),  
unique(sort(grep("&", meta$group, value = T))), 
"low expression"), labels = c(names(clean_module),unique(sort(grep("&", meta$group, value = T))), "Background"))


p1 <- ggplot(meta, aes(y= group, x  = nFeature_RNA))+
  geom_density_ridges() +
  geom_jitter(height = 0.1, size = 0.5, alpha = 0.5, aes(color = patient))+#geom_boxplot( width = 0.2, outliers = F)+
  theme_classic(base_size = 15)+
  scale_color_manual(values = brewer.set1(9))+
  labs(x = "Unique Gene Counts", y = "Module group")

p2 <- ggplot(meta, aes(y= group, x  = nCount_RNA))+
  geom_density_ridges() +
  geom_jitter(height = 0.1, size = 0.5, alpha = 0.5, aes(color = patient))+#geom_boxplot( width = 0.2, outliers = F)+
  theme_classic(base_size = 15)+
  scale_color_manual(values = brewer.set1(9))+
  labs(x = "Total Gene Counts", y = "Module group")

ggsave("ridge_plot.pdf", plot = p1+p2+plot_layout(axes = "collect", guides = "collect"), width = 17, height = 15)
ggsave("ridge_plot.png", plot = p1+p2+plot_layout(axes = "collect", guides = "collect"), width = 17, height = 15, dpi = 350)


