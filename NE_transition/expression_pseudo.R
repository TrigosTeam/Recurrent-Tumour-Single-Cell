library(Seurat)
library(scran)
library(scater)
library(data.table)
library(dplyr)
library(patchwork)
library(tidyverse)
library(MAST)
library(pals)
library(TSCAN)
library(SingleCellExperiment)
library(circlize)
source("~/CASCADEpaper/paper/functions.R", echo=F)
setwd("~/CASCADEpaper/paper/NE_transition")
paths <- system("realpath /trigos_team/CASCADE/Analysis/230113_seurat_intron/*/*/*non_tumour_removed.Rds", intern = T)
paths_split <- unlist(lapply(strsplit(paths, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
paths2 <- system("realpath /trigos_team/CASCADE/Analysis/240603_seurat_intron/CA0027/dura*/*non_tumour_removed_v4.Rds", intern = T)
paths_split2 <- unlist(lapply(strsplit(paths2, split = "/"), function(x) paste(x[6], x[7], sep = "_")))

paths <- c(paths, paths2)
paths_split <- c(paths_split, paths_split2)
paths <- setNames(paths,paths_split)
meta <- readRDS("~/CASCADEpaper/paper/Fig5_archetype/5_integrated_final_module_meta.Rds")
srt<- readRDS(paths[1])
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

clean_module <- readRDS("~/CASCADEpaper/paper/Fig5_archetype/clean_module.Rds")
srt <- AddModuleScore(srt, features  = clean_module)
names(clean_module)
colnames(srt@meta.data)[grep("Cluster", colnames(srt@meta.data))] <- paste(names(clean_module), "module")

high_cells <-  apply(srt@meta.data[, grep("module", colnames(srt@meta.data))], 2, function(x){
  bins <- cut(x, breaks = 100)
  group <- bins %in% levels(bins)[60:100]
  return(colnames(srt)[group])
})

group<-  apply(srt@meta.data[, grep("module", colnames(srt@meta.data))], 2, function(x){
  bins <- cut(x, breaks = 100)
  group <- bins %in% levels(bins)[60:100]
})
ind_high <- apply(group, 1, function(x) if(sum(x) > 1) data.frame( group = paste(names(clean_module)[x], collapse = "&"), freq = sum(x)) else if(sum(x)==1) data.frame(group = names(clean_module)[x], freq = 1) else data.frame(group = "background", freq = 0))
ind_high <- rbindlist(ind_high)
ind <- grepl("AR&", ind_high$group)&ind_high$freq == 2
ind_high$group[ind] <- gsub("AR&", "",ind_high$group[ind])
ind_high$freq[ind] <- ind_high$freq[ind]-1

srt$group <- factor(ind_high$group, levels = ind_high %>% arrange(freq) %>% pull(group) %>% unique())
color_high <- setNames( c("dodgerblue3", brewer.dark2(6)[2:6]),names(clean_module))
color_high <- c(color_high, setNames("lightgrey", "background"))
mixed <- sapply(unique(ind_high$group), function(x){
  if(grepl("&", x)){
    modules <- strsplit(x, split = "&")[[1]]
    if (length(modules) <5){
      colors <- c(color_high[modules], rep("white", 4-length(modules)))
      c1 <- colorRamp2(c(0, 1), colors[1:2])(0.5)
      c2 <- colorRamp2(c(0, 1), colors[3:4])(0.5)
      final_color <-  colorRamp2(c(0, 1), c(c1, c2))(0.5)
      return(final_color)} 
  }else{
    return(as.character(color_high[x]))
  }
})
srt$group2 <- as.character(srt$group)
srt$group2[srt$group2 == "NE1&NE2"] <- "NE1-NE2"
srt$group2[grep("&", srt$group2)] <- "transition"
srt$group2 <- factor(srt$group2, levels = c("AR" ,"inflammation","NE1","NE2" ,"NE1-NE2","cycling","hypoxia","background","transition"  ))



line.data <- readRDS("~/CASCADEpaper/paper/NE_transition/line.data.Rds")
TSCAN.pseudo <- readRDS("~/CASCADEpaper/paper/NE_transition/TSCAN.pseudo.Rds")
srt$pseudo <- TSCAN.pseudo

srt$group_ord <- as.integer(srt$group2)
setNames("slategrey", "transition")
colors <- c(color_high, setNames(mixed["NE1&NE2"], "NE1-NE2"),setNames("slategrey", "transition"))
colors <- colors[levels(srt$group2)]

f1 <- coneraxes(DimPlot(srt, group.by = "group2")+ geom_line(data=line.data, mapping=aes(x=UMAP_1, y=UMAP_2, group=edge), color = "red")+
            ggtitle("Archetype Modules"))+scale_color_manual(values=colors, 
                                                                  label = c("AR", "Inflammation", "NE1", "NE2", "NE1-NE2", "Cycling","Glycolysis", "Background", "Transition"))+theme(legend.position = "bottom")

p1 <- rmaxes(FeaturePlot(srt, "NE1 module", cols = rev(brewer.rdylbu(10)))+NoLegend()+ geom_line(data=line.data, mapping=aes(x=UMAP_1, y=UMAP_2, group=edge), color = "red")+ggtitle("NE1"))
p2 <- rmaxes(FeaturePlot(srt, "NE2 module", cols = rev(brewer.rdylbu(10)))+NoLegend()+ geom_line(data=line.data, mapping=aes(x=UMAP_1, y=UMAP_2, group=edge), color = "red")+ggtitle("NE2"))
p3 <- coneraxes(FeaturePlot(srt, "pseudo")+ geom_line(data=line.data, mapping=aes(x=UMAP_1, y=UMAP_2, group=edge), color = "red")+ggtitle("pseudotime"))

p1+p2+p3+plot_layout(ncol = 1, nrow = 3)
ggsave("pseudoNE12.pdf", width = 4, height = 9)


srt$ASCL1 <- srt@assays$RNA@data["ASCL1", ]

ggplot(srt@meta.data, aes(x = group2, y = ASCL1)) + geom_boxplot()

f2 <- ggplot(srt@meta.data, aes(x=Cluster3, y = Cluster4,color = ASCL1)) + geom_point() + scale_color_distiller(palette = "Spectral", breaks= c(0, 1.5, 2.5))+
  theme_classic(base_size = 18)+
  labs(x = "NE1 Module Score", y = "NE2 Module Score")+
  theme(legend.position = c(0.8, 0.2), 
        legend.direction = "horizontal",
        legend.background = element_blank(), 
        legend.title.position = "top",
        plot.background = element_rect(fill = "transparent", color = NA)) #c(0.6, 0.2), legend.direction = "horizontal"
leg <- as_ggplot(get_legend(f1))
f3 <- coneraxes(FeaturePlot(srt, "ASCL1"))+theme(legend.position = "none")
(f1+f3+f2)/leg+plot_layout( heights = c(4, 1))
ggsave("ASCL1 expression.pdf", width = 9, height = 5)

### signature trend in NE only cells -------
NEsrt <- subset(srt, subset = group%in% c("AR", "NE2", "NE1&NE2"))

df <- lapply(clean_module$NE2, function(x){
  temp <- data.frame(pseudo = srt$pseudo, expression = as.numeric(srt[["RNA"]][x,]), marker =x)
  temp$cor <- cor(temp$pseudo, temp$expression, method = "spearman")
  return(temp)
})
df <- rbindlist(df)

ggplot(df, aes(x = pseudo, y = expression, color = marker)) + geom_point()+
geom_smooth(method='lm', formula= y~x)

genelist <- split(clean_module$NE1,  ceiling(seq_along(clean_module$NE1)/6))

pdf("expression_pseudoNE1.pdf", width = 6, height = 4)
df <- lapply(genelist, function(genes){
  df <- lapply(genes, function(x){
    temp <- data.frame(pseudo = srt$pseudo, expression = as.numeric(srt[["RNA"]][x,]), marker =x)
    temp$cor <- cor(temp$pseudo, temp$expression, method = "spearman")
    return(temp)
  })
  df <- rbindlist(df)
  df$marker <- paste(df$marker, round(df$cor, 3), sep = ":")
  g <- ggplot(df, aes(x = pseudo, y = expression)) + geom_point(aes(color = marker))+
    geom_smooth(method="loess")+facet_wrap(.~marker,ncol = 3) + theme_bw() + NoLegend()
  print(g)
  return(df)
})
dev.off()
NE1df <- df

genelist <- split(clean_module$NE2,  ceiling(seq_along(clean_module$NE2)/6))

df <- rbindlist(df)
subdf <- df %>% select(marker, cor) %>% distinct()


HEPACAM2
CACNA2D1
TNIK
HDAC9
PLCB4
WASF3
PVT1
NFIB

genes <- clean_module$NE2
df <- lapply(genes, function(x){
  temp <- data.frame(pseudo = NEsrt$pseudo, expression = as.numeric(NEsrt[["RNA"]][x,]), marker =x, 
                     gene = x, group = NEsrt$group)
  temp$cor <- cor(temp$pseudo, temp$expression, method = "spearman")
  return(temp)
})

df <- rbindlist(df)
df$marker <- paste(df$marker, round(df$cor, 3), sep = ":")
g <- ggplot(df, aes(x = pseudo, y = expression)) + geom_point(aes(color = group))+
  geom_smooth(method="loess")+facet_wrap(.~marker, ncol = 3) + theme_bw() + 
  scale_color_manual(values = mixed[levels(NEsrt$group)])
print(g)


library(ggpubr)

my_comparisons <- list( c("AR", "NE2"), c("NE2", "NE1&NE2"), c("AR", "NE1&NE2") )
ggboxplot(df, x = "group", y = "expression",color = "group", facet.by = "marker")+ 
  stat_compare_means(group = "marker", comparisons = my_comparisons, label = "p.signif") # Add pairwise comparisons p-value

compdf <- compare_means(expression ~ group, data = df, 
              group.by = "gene")
compdf <- as.data.frame(compdf)


median_df <- df %>% group_by(gene, group) %>% summarise(med = median(expression))
median_df <- as.data.frame(median_df)

increases <- apply(compdf, 1, function(x){
  temp <- median_df$med[median_df$gene==x["gene"]&median_df$group==x["group1"]]
  temp2 <- median_df$med[median_df$gene==x["gene"]&median_df$group==x["group2"]]
  return(ifelse(temp<temp2, "increase", ifelse(temp==temp2, "equal", "decrease")))
})
compdf$trend <- increases
compdf$p.signif <- ifelse(compdf$p.signif=="ns", "ns", "sig")
compdf$type <- paste(compdf$p.signif, compdf$trend)

gene_trend <- sapply(split(compdf, compdf$gene), function(x) paste(x$type, collapse = ";"))
gene_trend[grepl("^sig decrease", gene_trend)&grepl("sig increase$", gene_trend)]

pdf("expression_pseudoNE2_NEsrt.pdf", width = 4, height = 8)
  df <- lapply(clean_module$NE2, function(x){
    temp <- data.frame(pseudo = NEsrt$pseudo, expression = as.numeric(NEsrt[["RNA"]][x,]), marker =x, gene = x, group = NEsrt$group)
    temp$cor <- cor(temp$pseudo, temp$expression, method = "spearman")
    
  temp$marker <- paste(temp$marker, round(temp$cor, 3), sep = ":")
  g <- ggplot(temp, aes(x = pseudo, y = expression)) + geom_point(aes(color = group))+
    geom_smooth(method="loess")+ggtitle(unique(temp$marker)) + theme_bw() + scale_color_manual(values = mixed[levels(NEsrt$group)]) + NoLegend()

  p <- ggboxplot(temp, x = "group", y = "expression",color = "group", palette = mixed[c("AR", "NE2", "NE1&NE2")])+ 
    stat_compare_means(group = "marker", comparisons = my_comparisons, label = "p.signif")
  print(g/p)
  return(temp)
})
dev.off()

genes <- c( "CACNA2D1","NFIB","HEPACAM2", "RAD51B","HDAC9", "FHIT", "PRIM2", "PLCL2")

plot_psuedo <- function(x, df){
  temp <- df %>% filter(gene == x)
  g <- ggplot(temp, aes( y = expression,x = pseudo)) + 
    geom_point(aes(color = group), size = 0.2)+
    geom_smooth(method="loess")+
    ggtitle(x) + 
    theme_bw(base_size = 16) + 
    theme(legend.position = "none", plot.title=element_text(size = 18),
          plot.margin = margin(t = 0.1, l = 0, r = 0.1, b = 0, "cm"), 
          axis.title.y = element_blank(),axis.title.x = element_blank())+
    scale_color_manual(values = mixed[levels(NEsrt$group)])+
    scale_x_continuous(breaks = c(0, 25, 50))
  p <- ggboxplot(temp, x = "group", y = "expression",color = "group", 
                 palette = mixed[c("AR", "NE2", "NE1&NE2")],  outlier.size = 0.1, ggtheme = theme_bw())+
    stat_compare_means(group = "marker", comparisons = my_comparisons, label = "p.signif", method = "wilcox.test",
                       bracket.size = 0.1, label.y.npc = "bottom", hide.ns = F, vjust =1)+
    theme(text = element_text(size = 12), legend.position = "none",
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          plot.title=element_text(size = 18),
          # axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = margin(t = 0, l = 0, r = 0.1, b = 0, "cm")) + 
    scale_y_continuous(breaks = NULL)
  return(g+p+plot_layout(axes = "collect"))
  # return(p)
}

f <- lapply(genes, function(x){
  return(plot_psuedo(x, df))
})
p <- wrap_plots(f, ncol = 2, byrow = T, guides = "collect")
p&plot_layout(axes = "collect")
ggsave("expression_pseudoNE2_NEsrt_selected.pdf", width = 20, height = 20, units = "cm")
ggsave("expression_pseudoNE2_NEsrt_selected_boxplot.pdf", width = 9, height = 8, units = "cm")


df <- rbind(data.frame(exp = srt$Cluster3, pseudotime = srt$pseudo, module = "NE1"),
                  data.frame(exp = srt$Cluster4, pseudotime = srt$pseudo, module = "NE2"))

ggplot(df, aes(x = pseudotime, y = exp, color = module))+geom_point(alpha = 0.4)+
  scale_color_manual(values = color_high)+theme_bw(base_size = 9)+
  theme(legend.position="top",
        legend.justification="right",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))


ggsave("~/CASCADEpaper/paper/NE_transition/NE12pseudo.pdf", width = 5, height = 5, unit = "cm")

