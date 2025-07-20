library(dplyr)

meta <- readRDS("~/Fig2_signature/sc_integration/tumor_only_meta_expression.Rds")
sites <- meta$sample
sites[grep("brain", sites)] <- "brain"
sites[grep("dura", sites)] <- "dura"
sites[grep("prostate", sites)] <- "prostate"
sites[grep("liver", sites)] <- "liver"
sites[grep("fat", sites)] <- "fat"
sites[grep("lymph|hilar", sites)] <- "LN"
sites[grep("rib|vertebra", sites)] <- "bone"
sites[grep("lung", sites)] <- "lung"
sites[grep("abdomen", sites)] <- "abdomen"
sites[grep("bladder", sites)] <- "bladder"
meta$site <- sites
meta$pathology <- ifelse(meta$patient %in% c("CA0090", "CA0046"), "NE", ifelse(meta$patient %in% c("CA0027", "CA0058"), "Mixed", "AD"))
meta$patient <- gsub("00", "", meta$patient)
g1 <- ggplot(meta, aes(x = site, y = AR, fill = patient)) + geom_boxplot()+
  geom_hline(yintercept = 0.3*(max(meta$AR) - min(meta$AR))) +
  facet_grid(.~pathology+patient, scales = "free", space = "free")+
  scale_fill_manual(values = brewer.set1(9))+
  theme_bw()+
  theme(axis.text.x = element_blank(), legend.position = "right") + labs(x = "")

g2 <- ggplot(meta, aes(x = site, y = FOLH1, fill = patient)) + geom_boxplot()+
  geom_hline(yintercept = 0.3*(max(meta$FOLH1) - min(meta$FOLH1))) +
  facet_grid(.~pathology+patient, scales = "free", space = "free")+
  scale_fill_manual(values = brewer.set1(9))+
  theme_bw()+
  theme(axis.text.x = element_blank(), legend.position = "right")+labs(x = "")

g3 <- ggplot(meta, aes(x = site, y = ASCL1, fill = patient)) + geom_boxplot()+
  geom_hline(yintercept = 0.3*(max(meta$ASCL1) - min(meta$ASCL1))) +
  facet_grid(.~pathology+patient, scales = "free", space = "free")+
  scale_fill_manual(values = brewer.set1(9))+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90), legend.position = "right")

meta2 <- melt(meta, measure.vars = c("AR", "ASCL1", "FOLH1"), variable.name = "marker", value.name = "expression")

pdf("AR_ASCL_FOLH_expression.pdf", width = 10, height = 7)
g2+g1+g3 + plot_layout(nrow = 3, guides = "collect", axis_titles = "collect", axes = "collect")
ggplot(meta2, aes(x = site, y = expression, fill = patient)) + geom_boxplot()+
  geom_hline(yintercept = 0.3*(max(meta$ASCL1) - min(meta$ASCL1))) +
  facet_grid(marker~pathology+patient, space = "free", scales = "free_x")+
  scale_fill_manual(values = brewer.set1(9))+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90), legend.position = "right")
dev.off()
  
melt(meta, )
#phenotype generation is in ~/Fig2_signature/sc_integration/AR_PSMA_NE_phenotype.R
meta <- readRDS("~/Fig2_signature/sc_integration/tumor_only_phenotype_meta.Rds")
sites <- meta$sample
sites[grep("brain", sites)] <- "brain"
sites[grep("dura", sites)] <- "dura"
sites[grep("prostate", sites)] <- "prostate"
sites[grep("liver", sites)] <- "liver"
sites[grep("fat", sites)] <- "fat"
sites[grep("lymph|hilar", sites)] <- "LN"
sites[grep("rib|vertebra", sites)] <- "bone"
sites[grep("lung", sites)] <- "lung"
sites[grep("abdomen", sites)] <- "abdomen"
sites[grep("bladder", sites)] <- "bladder"
meta$site <- sites
meta %>% group_by(patient, site, adeno_pheno) %>% summarise(n = n()) %>% mutate(Freq= n/sum(n)) -> phenodf

lv <- unique(phenodf$adeno_pheno)[c(1,2,4, 7, 3, 5, 6, 8)]
phenodf$phenotype <- factor(phenodf$adeno_pheno, levels = lv )
cols <- setNames(brewer.set3(8), levels(phenodf$phenotype))
phenodf$pathology <- ifelse(phenodf$patient %in% c("CA0090", "CA0046"), "NE", ifelse(phenodf$patient %in% c("CA0027", "CA0058"), "Mixed", "AD"))
phenodf$patient <- gsub("00","", phenodf$patient)




pdf("PSMA_heterogeneity_phenotype.pdf", width = 10, height = 5)
ggplot(phenodf, aes(x = site, y = Freq, fill = phenotype)) + geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(hjust =1, angle = 90, vjust = 0.5)) +
  facet_grid(.~pathology+patient, space = "free_x", scales = "free_x") + 
  scale_fill_manual(values = cols , labels = c("AR high/FOLH1 high/ASCL1 high","AR high/FOLH1 high/ASCL1 low", 
                                              "AR low/FOLH1 high/ASCL1 high", "AR low/FOLH1 low/ASCL1 low" ,"AR high","FOLH1 high","ASCL1 high","cycling" )) 

dev.off()


# check CN of FOLH1
geneCNAlist <- readRDS("~/Fig3_genomics_normal_cells/copy_number/geneCNAlist.Rds")
CNs <- rbindlist(geneCNAlist, idcol = "sample")
CNs <- CNs %>% filter(gene == "FOLH1")
CNs$patient <- substr(CNs$sample, 1, 6)
sites <- CNs$sample
sites[grep("brain", sites)] <- "brain"
sites[grep("dura", sites)] <- "dura"
sites[grep("prostate", sites)] <- "prostate"
sites[grep("liver", sites)] <- "liver"
sites[grep("fat", sites)] <- "fat"
sites[grep("lymph|hilar", sites)] <- "LN"
sites[grep("rib|vertebra", sites)] <- "bone"
sites[grep("lung", sites)] <- "lung"
sites[grep("abdomen", sites)] <- "abdomen"
sites[grep("bladder", sites)] <- "bladder"
CNs$site <- sites
CNs$pathology <- ifelse(CNs$patient %in% c("CA0090", "CA0046"), "NE", ifelse(CNs$patient %in% c("CA0027", "CA0058"), "Mixed", "AD"))
CNs$patient <- gsub("00","", CNs$patient)
CNs <- CNs[-9, ]
CNs$sample[8] <- "CA0035_paraaortic_lymph_node_1"
CNs$sample[7] <- "CA0035_bladder_2" 


sites <- unique(CNs$sample)
sites[grep("brain", sites)] <- "brain"
sites[grep("dura", sites)] <- "dura"
sites[grep("prostate", sites)] <- "prostate"
sites[grep("liver", sites)] <- "liver"
sites[grep("fat", sites)] <- "fat"
sites[grep("lymph|hilar", sites)] <- "LN"
sites[grep("rib|vertebra", sites)] <- "bone"
sites[grep("lung", sites)] <- "lung"
sites[grep("abdomen", sites)] <- "abdomen"
sites[grep("bladder", sites)] <- "bladder"


meta <- readRDS("~/Fig2_signature/sc_integration/tumor_only_meta_expression.Rds")
meta3 <- meta %>% group_by(sample) %>% mutate(median_FOLH1 = median(FOLH1)) %>% select(sample,  median_FOLH1) %>% distinct() 

CNs <- merge(CNs, meta3, by = "sample")
CNs$median_FOLH1 <- CNs$median_FOLH1*3

pdf("PSMA_CN_check.pdf", width = 10, height = 4 )
p <- ggplot(CNs, aes(x=sample,fill = patient))+
  geom_bar(stat = "identity", aes (y = CN))+
  geom_point(aes( y = median_FOLH1))+
  scale_fill_manual(values = brewer.set1(9))+
  geom_hline(yintercept = 2) +
  scale_x_discrete(labels = sites, breaks = unique(CNs$sample))+
  scale_y_continuous("CN", breaks = 0:7,
                     sec.axis = sec_axis(~./3, "median experssion")) +
  facet_grid(.~pathology+patient, space = "free_x", scales = "free_x") + 
  theme_bw()+
  theme(axis.text.x = element_text(hjust =1, angle = 90, vjust = 0.5))
print(p)
dev.off()




# subclone difference -----------------------------------------------------
subclone_anno <- readRDS("~/Fig3_genomics_normal_cells/Infercnv/subclone_anno.Rds")
sig_region <- readRDS("~/Fig3_genomics_normal_cells/Infercnv/region_test/sig_region.Rds")
multi_g <- readRDS("~/genomic/exon_only_reference.Rds")
psma_reg <- multi_g[multi_g$gene.name == "FOLH1", ]

sig_region <- rbindlist(sig_region, idcol = "sample", fill = T)

sig_region1 <- sig_region %>% filter(chr == 11) %>% filter(start <= 49145092 & end >= 49208638)

sig_region1

paths <- system("realpath /trigos_team/CASCADE/Analysis/230113_seurat_intron/*/*/*non_tumour_removed.Rds", intern = T)
paths_split <- unlist(lapply(strsplit(paths, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
paths <- setNames(paths, paths_split)
exp <- list()
pdf("FOLH1_subclone_expression.pdf", width = 9, height = 4)
for (i in sig_region1$sample){
  srt <- readRDS(paths[i])
  srt$subclone <- subclone_anno[[i]][colnames(srt)]
  Idents(srt) <- srt$subclone
  subregion <- sig_region1 %>% filter(sample == i) %>% as.data.frame()
  subregion <- as.character(subregion[sort(paste0(unique(subclone_anno[[i]]), ".", unique(subclone_anno[[i]])))])
  h <- DimPlot(srt,label = T) + FeaturePlot(srt, "FOLH1") +VlnPlot(srt, "FOLH1")
  print(h+plot_annotation(title = i, subtitle = paste(paste(sort(unique(subclone_anno[[i]])), subregion, sep = ":"), collapse = " ")))
  exp[[i]] <- data.frame(subclone = srt$subclone, exp = srt[["RNA2"]]@data["FOLH1", ])
}
dev.off()
saveRDS(exp, "subclone_PSMA_exp.Rds")

library(ggpubr)
library(colorspace)
library(scales)
my_comparisons <- list( c("1", "2"), c("2", "3"), c("1", "3") )
ggboxplot(exp$CA0035_bladder_2, x = "subclone", y = "exp",
          color = "subclone", palette = hue_pal()(4),
          add = "jitter", 
          short.panel.labs = T, title = i) + 
  stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.x = 1, label.y = 6.5, method = "kruskal.test")+
  ggtitle("CA0035_bladder_2")

ggboxplot(exp$CA0035_paraaortic_lymph_node_1, x = "subclone", y = "exp",
          color = "subclone", palette = hue_pal()(4),
          add = "jitter", 
          short.panel.labs = T) + 
  stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.x = 1, label.y = 7, method = "kruskal.test")
