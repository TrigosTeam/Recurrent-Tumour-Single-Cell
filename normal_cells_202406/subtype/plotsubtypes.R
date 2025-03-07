######Tcell B----------
Tcells <- readRDS("~/CASCADEpaper/paper/normal_cells_202406/subtype/Tcells.Rds")
f1 <- coneraxes( DimPlot(Tcells, group.by = "subtype", reduction = "harmony_umap", cols = ggsci::pal_observable()(3)), front.size = 16)+ theme(legend.position = "top", plot.title = element_blank())
f1
tcell_marker <- list(Treg = c("IL2RA",  "RTKN2","IKZF2"), 
                     Effector = c("CD8A", "CD8B"),
                     Helper = c("IL4R", "CD4"),
                     'Mem/naive' = c("IL7R", "BACH2", "LEF1"))
f2 <- DotPlot(Tcells, group.by = "subtype", features = tcell_marker)+
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), 
                                                                            axis.title = element_blank(), 
                                                                            axis.text.y = element_text(angle = -30, hjust = 1, vjust = 1), 
                                                                            legend.position = "none")
free(f1)+f2+plot_layout(widths = c(1, 1.5))
ggsave("subtypepdf/Tcell1.pdf", width = 5.6, height = 3)

p3 <- barplot_bypatient(Tcells)
p4 <- barplot_bytissue_byrow(Tcells)
p3 + p4+ plot_layout(axes = "collect")
ggsave("subtypepdf/Tcell2.pdf", width = 6, height = 3)
f3 <- coneraxes(DimPlot(Tcells, reduction = "harmony_umap", label = T, group.by = "cluster_marker", repel = T)+theme(plot.title = element_blank()), front.size = 15) 
ggsave("subtypepdf/Tcell3.pdf", width = 3, height = 3)
f4 <- rmaxes( FeaturePlot(Tcells, unique(Tcells$marker)[1:7], reduction = "harmony_umap", cols = rev(brewer.rdylbu(10))))
ggsave("subtypepdf/Tcell4.pdf", width = 3, height = 3)

#####Bcell ----------
Bcells <- readRDS("~/CASCADEpaper/paper/normal_cells_202406/subtype/Bcells.Rds")
Bcells <- plot_cluster(Bcells)

markerlist <- list(Plasma = c("CD38", "XBP1", "IGHGP"), 
                   `Mem/naive` = c("BANK1", "CD83", "MARCH1"), 
                   Erythroid = c("RHCE", "HBA2", "ANK1"))
p1 <- coneraxes(DimPlot(Bcells, reduction = "harmony_umap", group.by = "subtype", cols = pal_observable()(6)), front.size = 16)+ theme(legend.position = "top", plot.title = element_blank())
p2 <- DotPlot(Bcells, features = markerlist, group.by = "subtype")+theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
                                                                         axis.title = element_blank(), 
                                                                         # axis.text.y = element_text(angle = -30, hjust = 1, vjust = 1), 
                                                                         legend.position = "none")

free(p1)+p2 +plot_layout(widths = c(1, 1.5))
ggsave("subtypepdf/Bcell1.pdf", width = 5.6, height = 3)

p3 <- barplot_bypatient(Bcells)
p4 <- barplot_bytissue_byrow(Bcells)
p3 + p4+ plot_layout(axes = "collect")
ggsave("subtypepdf/Bcell2.pdf", width = 6, height = 3)

g1 <- rmaxes(DimPlot(Bcells, reduction = "harmony_umap", group.by = "site")+scale_color_manual(values = site_cols)+theme(plot.title = element_blank(), legend.position = "bottom"), front.size = 26)
g2 <- barplot_bytissue(Bcells)+NoLegend()
markers <- readRDS("~/CASCADEpaper/paper/normal_cells_202406/subtype/Bcells_tissuemarker.Rds")
features <- markers %>% group_by(cluster)%>% arrange(desc(pct.1-pct.2)) %>% arrange(desc(avg_log2FC)) %>% dplyr::slice(1)
features$cluster <- as.character(features$cluster)
features$cluster <- factor(features$cluster, levels =sort(as.character(features$cluster)) ) 
group <- setNames(as.character(features$cluster), as.character(features$gene))
g3 <- vlnplot(Bcells, features$gene[order(features$cluster)], group = group)

g1+inset_element(g2, left = 0.6, top = 1, right = 1, bottom = 0)
# g1+free(g2)+plot_layout(widths = c(1.5, 1))&theme( plot.title = element_blank(),legend.margin = margin(t = -10))
ggsave("subtypepdf/Bcell3.pdf", width = 6, height = 5)

g3
ggsave("subtypepdf/Bcell4.pdf", width = 6, height = 7)

##### Macrophage ----------
Macro <- readRDS("~/CASCADEpaper/paper/normal_cells_202406/subtype/Macro.Rds")
markerlist <- list( `angio\nTAM` = c("VCAN", "FCN1","MARCO","OLR1"), 
                    TRM = c("F13A1", "LYVE1", "CD163L1"), 
                    TAM = c("SPP1", "CD83", "CD109"))
plot_subtype(Macro, "Macro", markerlist = markerlist)

##### Endothelial --------
Endo <- readRDS("~/CASCADEpaper/paper/normal_cells_202406/subtype/Endo.Rds")
markerlist <- list(`tip\nlike `= c("SLC45A4", "IGFBP3"),#, "ESM1"
                   arterial = c("LTBP4", "FBLN5"), #"GJA5", , "SOX17"
                   vein = c( "VCAM1","TSHZ2"), #,"MCTP1"
                   lymphatic = c("PKHD1L1", "STON2"), #, "PDGFC"
                   `liver\nsinusoidal` = c("ATRNL1", "OIT3")) #, "FGF23", "DNASE1L3"
plot_subtype(Endo, "Endo", markerlist = markerlist)

##### Epithelial --------
Epi <- readRDS("~/CASCADEpaper/paper/normal_cells_202406/subtype/Epi.Rds")

markerlist <- list(basal = c("KRT5", "KRT15", "VAV3", "TP63"), #, "TG"
                   luminal = c( "CPA6", "ALDH1A2"), #"ROBO1",
                   ciliated = c("ADCY2", "PTPRN2"), #"NR4A1", 
                   aveoliar = c("ABCA3", "SFTPB"))# "DOCK2", 
plot_subtype(Epi, "Epi", markerlist = markerlist)

##### Fibroblsat -------
Fibro <- readRDS("~/CASCADEpaper/paper/normal_cells_202406/subtype/Fibro.Rds")
markerlist <- list(adipose =c("LSAMP", "TRPC4"),#, "DACH1"
                   chondrocyte = c("IBSP", "SATB2"), #,"CD36"
                   myofibroblast = c("COL3A1","VCAN"),# "COL6A3",
                   pricyte= c("ADAMTS9","EBF1"),#"DLC1",
                   `smooth\nmuscle` = c("ACTA2", "MYH11")#, "CNN1"
)
plot_subtype(Fibro, "Fibro", markerlist = markerlist)

