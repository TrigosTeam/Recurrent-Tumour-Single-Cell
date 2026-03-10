source("~/CASCADEpaper/paper/functions.R")
source("~/CASCADEpaper/paper/normal_cells_202406/normal_markers.R")

srt <- readRDS("~/integration/2024_06/perpatient/tumor_normal/CA0090_ruvsrt.Rds")
final_normal_anno <- readRDS("~/CASCADEpaper/paper/normal_cells_202406/subtype/final_normal_anno.Rds")

DimPlot(srt)

srt$cell_anno <- "Tumours"
srt$cell_anno[intersect(colnames(srt), names(final_normal_anno))] <- as.character(final_normal_anno[intersect(colnames(srt), names(final_normal_anno))]) 
srt$cell_anno <- ifelse(srt$cell_anno %in% c("Neurons", "Tumours"), srt$cell_anno, "Other normal cells")
Idents(srt) <-  srt$cell_anno

cols <- c("Tumours" = "grey40", "Neurons" = "red", "Other normal cells" = "steelblue")
p1 = DimPlot(srt, group.by = "cell_anno") + scale_color_manual(values = cols) + ggtitle("CA90 cells integration")

p1 = coneraxes(p1)
p1+theme(legend.position = "bottom")

f = VlnPlot(srt,c("ASCL1", Neuron_feature) )


exp = GetAssayData(srt, layer = "data")[c("ASCL1","SYT1", "RBFOX3", "SNAP25"), ]
data = data.frame(group = srt$cell_anno, t(exp)) 
data = reshape2::melt(data = data)
ggplot(data, aes(x = group, y = value, fill = group))+
  geom_violin()+
  scale_fill_manual(values = cols)+
  geom_jitter(height = 0.1, size = 0.1)+
  facet_grid(.~variable)+
  theme_classic(base_size = 15)+
  theme(legend.position = "none", axis.text.x =  element_text(angle = 45, hjust = 1, vjust = 1))+
  labs(y = "Expression", x = "")
  