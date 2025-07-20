library(Seurat)
library(ggplot2)
library(patchwork)
library(stringr)
library(ggpubr)
library(dplyr)

source("~/cols.R")
source("~/functions.R")

normal3 <- readRDS("~/normal_cell_annotation/subtype/final_normalsrt.Rds")
normal3$cell_type <- ifelse(normal3$cell_anno%in% c("T cells","B cells","Plasma cells",  "Macrophages"), "immune cells", 
                            ifelse(normal3$cell_anno%in% c("Fibroblasts","Adipocytes","Pericytes", "Endothelial cells"), "stromal cells", "others"))
normal3$pathology <- ifelse(normal3$patient %in% c("CA0090", "CA0046"), "NE", ifelse(normal3$patient %in% c("CA0027", "CA0058"), "Mixed", "AD"))
##### A-------------
anno_marker <- list(
  Tcell= c("CD247", "CD4", "CD8A"), #"SKAP1"
  Bcell = c("CD19","BANK1", "CD83", "MS4A1"),
  Plasma= c("CD38","PDK1","FKBP11", "JCHAIN", "SDC1"),
  Macrophage = c("CD86","CD163","FMN1"),
  Neuron = c( "NRXN1", "CADM2", "CSMD1"),#"KCNIP4"
  Hepotocytes = c( "CP", "ALB", "SORBS2"),
  Epithelial = c("WFDC2","EPCAM", "EHF", "TMC5"), 
  Chondrocyte = c("IBSP", "SATB2", "INSC"),
  Endothelial = c( "VWF", "FLT1",  "BTNL9"), 
  LymEndo = c("PKHD1L1", "CD36", "PROX1"), 
  # Fibrobast = c("LAMA2","PRKG1" ,"CALD1","LHFPL6"),
  Fibroblast = c("LAMA2","DCN", "LUM", "CFD"),
  Adipocytes =  c("MIR99AHG","TRPC4", "KAZN"),
  # Adipose =  c("MIR99AHG","FATP1", "PAT2", "NRG4", "P2RX4"),
  
  Pericytes = c("ABCC9", "RGS5","EBF2")
  
)
labs <- levels(normal3$cell_anno)
labs <- str_wrap(labs, width = 10)
labs <- setNames(labs, levels(normal3$cell_anno))
g <- DimPlot(normal3, reduction = "css_umap", group.by = "cell_anno", cols = brewer.set1(14))
g <- coneraxes(g, labelx = "css_UMAP", labely = "css_UMAP")
g + scale_color_manual(labels = labs, values = normal_cols) + theme(plot.title = element_blank())
ggsave("pdf/cell_anno_dimplot.pdf", width = 6, height = 5)

FeaturePlot(normal3, "CD68", reduction = "css_umap")
##### B-------------
labs2 <- labs
labs2 <- sub("\n.*", "", labs)


g <- DotPlot(normal3, features = anno_marker, group.by = "cell_anno")+ theme_classic2(base_size = 18)+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
g + scale_y_discrete("Cell Identities", labels = labs2) + theme(legend.position = "top",
                                                               legend.margin = margin(0,0, -0.5,0, unit = "cm"), 
                                                               legend.text = element_text(size = 13),      # Decrease text size
                                                               legend.title = element_text(size = 13),   # Adjust title size (optional)
                                                               legend.key.size = unit(0.5, "cm"), 
                                                               axis.title.x = element_blank(),
                                                               axis.title.y = element_blank(), 
                                                               panel.spacing = unit(0.2, "cm"),
                                                               strip.text = element_text(size = 13))
ggsave("pdf/dotplot.pdf", width = 18, height = 5)

##### C-------------
# by patient samples
df <- normal3@meta.data %>% group_by(patient,site,sample, cell_anno) %>% summarise(n = n()) %>% mutate (freq = n/sum(n))
df$cell_type <- ifelse(df$cell_anno%in% c("T cells","B cells","Plasma cells",  "Macrophages"), "immune cells", 
                       ifelse(df$cell_anno%in% c("Fibroblasts","Adipocytes","Pericytes", "Endothelial cells"), "stromal cells", "others"))
df$pathology <- ifelse(df$patient %in% c("CA0090", "CA0046"), "NE", ifelse(df$patient %in% c("CA0027", "CA0058"), "Mixed", "AD"))
df$patient <- gsub("00", "", df$patient)
df$labs <- paste0(signif(df$freq*100, 2))
df$labs[df$freq < 0.03] <- NA
xlab <- setNames(df$site, df$sample)
xlab <- paste(xlab, sapply(strsplit(names(xlab), split = "_"), function(x) tail(x, 1)))
xlab <- setNames(xlab, df$sample)
df2 <- df %>% group_by(sample) %>% mutate(n = sum(n))

f <- ggplot(df, aes(x = sample, y = freq, fill = cell_anno)) + geom_bar(stat = "identity") + 
  geom_text(aes(label = labs),
            position = position_stack(vjust = 0.5),        # Position text in the middle of the bars
            color = "white",                               # Set text color
            size = 3)+
  geom_text(data=df2,aes(x=sample,y=1,label=n),hjust=0, angle = 90)+
  ggh4x::facet_nested(.~pathology+patient, scales = "free", space = "free") + 
  theme_bw(base_size = 18) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
                                   panel.spacing = unit(0.1, "cm"), 
                                   strip.background = element_rect(fill = "white"),
                                   legend.position = "none",
        plot.margin = margin(r = -5))+
  scale_fill_manual(values =normal_cols )+
  scale_x_discrete(labels= xlab)+
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0, 1.1))+
  labs(x = "Sample", y = "Cell Proportion")
f
ggsave("pdf/barplot_per_sample.pdf",f,  width = 18, height = 5)


tb<- compare_means(freq~patient, df, group.by = "cell_anno", p.adjust.method = "bonferroni") #all not_significant
write.table(tb, "cell_anno_freq_bwt_patient.csv", quote = F, row.names = F)
f2 <- ggplot(df, aes(x = patient, y = freq, fill = patient)) + geom_boxplot()+theme_minimal(base_size = 18)+
  facet_wrap(.~cell_anno, scales = "free")+
  labs(y = "Cell Proportion", x = "Patient")+
  theme(axis.text.x = element_text(angle = 90, vjust=1,hjust = 1))
ggsave("pdf/boxplot_interpatient.pdf", width = 18, height = 10)
f+f2+plot_annotation(tag_levels = 'A', tag_prefix = '(', tag_suffix = ')')+plot_layout(heights = c(1, 3), nrow = 2)
ggsave("pdf/S4.pdf", width = 18, height = 15)



tb<- compare_means(freq~pathology, df, group.by = "cell_anno", p.adjust.method = "bonferroni") # not significant
df2 <- df %>% group_by(pathology, cell_anno) %>% summarise(avg_freq = mean(freq))
df2$labs<- paste0(signif(df2$avg_freq*100, 2), "%")
df2$labs[df2$avg_freq < 0.03] <- NA

p1 <- ggplot(df2, aes(x = pathology, y = avg_freq, fill = cell_anno)) + geom_bar(stat = "identity") + 
  geom_text(aes(label = labs),
            position = position_stack(vjust = 0.5),        # Position text in the middle of the bars
            color = "black",                               # Set text color
            size = 4)+
  theme_bw(base_size = 18) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none", 
                                   panel.spacing = unit(0.1, "cm"), 
                                   strip.background = element_rect(fill = "white"))+
  scale_fill_manual(values =normal_cols )+
  labs(x = "Pathology", y = "Average Cell Proportion")

compdf <- compare_means(freq ~ pathology,group.by = "cell_anno",   data = df,method = "wilcox.test") # no significant results for both freq and n

ggplot(df %>% filter(cell_anno %in% c("T cells", "B cells", "Macrophages")), aes(x = pathology, y = n, fill = cell_anno)) + geom_boxplot()+
  scale_fill_manual(values =normal_cols )+theme_bw(base_size = 18) +
  ggh4x::facet_nested(.~cell_anno, scales = "free")+
  theme(legend.position = "none")+
  stat_compare_means(comparisons = my_comparisons, label =  "p.signif", label.y.npc =  0.5,
                     method = "wilcox.test",
                     hide.ns = T)

my_comparisons <- list( c("AD", "NE"), c("AD", "Mixed"), c("NE", "Mixed") )
ggboxplot(df %>% filter(cell_anno %in% c("T cells", "B cells", "Macrophages")), 
          x = "pathology", y = "freq",
         facet.by = "cell_anno")+ 
  stat_compare_means(comparisons = my_comparisons, label =  "p.signif", method = "wilcox.test")# Add pairwise comparisons p-value

ggplot(df, aes(x = pathology, y = freq, fill = cell_anno))+geom_bar(stat = "identity")+
  scale_fill_manual(values = normal_cols )

g <- ggboxplot(df, x = "pathology", y = "freq", 
          facet.by = "cell_type")+theme_classic(base_size = 18)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(x = "Pathology", y = "Cell Proportion")+
  stat_compare_means(comparisons = my_comparisons, label =  "p.signif", method = "wilcox.test")

ggboxplot(df, x = "pathology", y = "freq", 
          facet.by = "cell_anno")+
  stat_compare_means(comparisons = my_comparisons, label =  "p.signif", method = "wilcox.test")
p1+g+plot_layout(widths = c(2, 3))
ggsave("pdf/average_ferq_pathology.pdf", width = 9, height = 5)
# by site
df <- normal3@meta.data %>% group_by(site,cell_anno) %>% summarise(n = n()) %>% mutate (freq = n/sum(n))
df$labs <- paste0(signif(df$freq*100, 2))
df$labs[df$freq < 0.03] <- NA
df2 <- df %>% group_by(site) %>% mutate(n = sum(n))

f2 <- ggplot(df, aes(x = site, y = freq, fill = cell_anno)) + geom_bar(stat = "identity") + 
  geom_text(aes(label = labs),
            position = position_stack(vjust = 0.5),        # Position text in the middle of the bars
            color = "white",                               # Set text color
            size = 3)+
  theme_bw(base_size = 20) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none", 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(), 
        plot.margin = margin(l = -5))+
  scale_fill_manual(values =normal_cols)+
  geom_text(data=df2,aes(x=site,y=1,label=n),hjust=0, angle = 90)+
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0, 1.1))+
  labs(y = "Cell Proportion", x = "Sites")

f2
ggsave("pdf/barplots_per_site.pdf", f2,width = 6, height = 5)

f + f2 + plot_layout(widths = c(3, 1))
ggsave("pdf/F-barplot.pdf", width = 12, height = 5)

##### by cell anno color by patient and site --------

df <- normal3@meta.data %>% group_by(cell_anno,patient) %>% summarise(n = n()) %>% mutate (freq = n/sum(n))
df$pathology <- ifelse(df$patient %in% c("CA0090", "CA0046"), "NE", ifelse(df$patient %in% c("CA0027", "CA0058"), "Mixed", "AD"))
df$patient <- gsub("00", "", df$patient)
df$labs <- paste0(signif(df$freq*100, 2), "%")
df$labs[df$freq < 0.03] <- NA
f1 <- ggplot(df, aes(x = cell_anno, y = freq, fill = patient)) + geom_bar(stat = "identity") +
  geom_text(aes(label = labs),
            position = position_stack(vjust = 0.5),        # Position text in the middle of the bars
            color = "white",                               # Set text color
            size = 3)+
  theme_bw(base_size = 18) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  scale_fill_manual(values = patient_cols)
df <- normal3@meta.data %>% group_by(cell_anno,site) %>% summarise(n = n()) %>% mutate (freq = n/sum(n))
df$labs <- paste0(signif(df$freq*100, 2), "%")
df$labs[df$freq < 0.03] <- NA
f2 <- ggplot(df, aes(x = cell_anno, y = freq, fill = site)) + geom_bar(stat = "identity") +
  geom_text(aes(label = labs),
            position = position_stack(vjust = 0.5),        # Position text in the middle of the bars
            color = "black",                               # Set text color
            size = 3)+
  theme_bw(base_size = 18) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  scale_fill_manual(values = site_cols)

f1/f2 + plot_layout(axes = "collect", guides = "collect")
ggsave("pdf/bycellanno.pdf", height = 10, width = 10)

##### Dby pathology-------------
df <- normal3@meta.data %>% group_by(pathology,cell_anno ) %>% summarise(n = n()) %>% mutate (sum = sum(n), freq = n/sum(n))
df$labs <- paste0(signif(df$freq*100, 2), "%")
df$labs[df$freq < 0.03] <- NA
df2 <- df %>% group_by(pathology) %>% mutate(n = sum(n))

prop.test(x = c(df$n[df$pathology=="NE"&df$cell_anno == "T cells"],df$n[df$pathology=="AD"&df$cell_anno == "T cells"]), 
          n = c(df$sum[df$pathology=="NE"&df$cell_anno == "T cells"],df$sum[df$pathology=="AD"&df$cell_anno == "T cells"]), alternative = "greater")

# 1-sample proportions test with continuity correction
# 
# data:  1236 out of 4256, null probability 0.217839445
# X-squared = 131.14, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is greater than 0.2178394
# 95 percent confidence interval:
#   0.2789879 1.0000000
# sample estimates:
#   p 
# 0.2904135 

prop.test(x = c(df$n[df$pathology=="NE"&df$cell_anno == "Macrophages"],df$n[df$pathology=="AD"&df$cell_anno == "Macrophages"]), 
          n = c(df$sum[df$pathology=="NE"&df$cell_anno == "Macrophages"],df$sum[df$pathology=="AD"&df$cell_anno == "Macrophages"]), alternative = "less")

# 1-sample proportions test with continuity correction
# 
# data:  2433 out of 10090, null probability 0.1395676692
# X-squared = 865.82, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is greater than 0.1395677
# 95 percent confidence interval:
#   0.2341461 1.0000000
# sample estimates:
#   p 
# 0.2411298 

g <- ggplot(df, aes(x = pathology, y = freq, fill = cell_anno)) + geom_bar(stat = "identity") + 
  geom_text(aes(label = labs),
            position = position_stack(vjust = 0.5),        # Position text in the middle of the bars
            color = "white",                               # Set text color
            size = 4)+
  theme_bw(base_size = 20) + theme(legend.position = "none",plot.margin = margin())+
  scale_fill_manual(values =brewer.set1(14), name = "Annotation")+
  geom_text(data=df2,aes(x=pathology,y=1,label=n),hjust=0, angle = 90)+
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1.1))+
  labs(y = "Proportion", x = "Pathology Groups")
g

df <- normal3@meta.data %>% group_by(pathology,cell_type ) %>% summarise(n = n()) %>% mutate (sum = sum(n), freq = n/sum(n))
df$labs <- paste0(signif(df$freq*100, 2), "%")
df$labs[df$freq < 0.03] <- NA

df$cell_type <- factor(df$cell_type, levels = c("immune cells", "stromal cells", "others"))
df2 <- df %>% group_by(pathology) %>% mutate(n = sum(n))
p <- ggplot(df, aes(x = pathology, y = freq, fill = cell_type)) + geom_bar(stat = "identity") + 
  geom_text(aes(label = labs),
            position = position_stack(vjust = 0.5),        # Position text in the middle of the bars
            color = "white",                               # Set text color
            size = 4)+
  theme_bw(base_size = 20) + scale_fill_manual(values = c("tomato", "royalblue1", "lightgrey"), name = "Cell Type", labels = c("Immune Cells", "Stromal Cells","Others")) +
  theme(plot.margin = margin(), 
        legend.margin = margin(l = -0.5, unit = "cm"), 
        legend.text = element_text(size = 14),      # Decrease text size
        legend.title = element_text(size = 14),   # Adjust title size (optional)
        legend.key.size = unit(0.5, "cm"))+
  geom_text(data=df2,aes(x=pathology,y=1,label=n),hjust=0, angle = 90)+
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1.1))+
  labs(y = "Proportion", x = "Pathology Groups")

p
g+p+ plot_layout(axes  = "collect")
ggsave("pdf/barplots_pathology.pdf", height = 5, width = 6)
f2+g+p+ plot_layout(axes  = "collect", widths = c(3, 1, 1))
ggsave("pdf/barplots.pdf", height = 5, width = 12)
## p test 
head(df)

prop.test(3262,6516, p = 0.303, alternative = "greater") #Mix stromal vs AD stromal
# 1-sample proportions test with continuity correction
# 
# data:  3262 out of 6516, null probability 0.303
# X-squared = 1203.9, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is greater than 0.303
# 95 percent confidence interval:
#   0.4903506 1.0000000
# sample estimates:
#   p 
# 0.5006139 

prop.test(3262,6516, p = 0.309, alternative = "greater") #Mix stromal vs NE stromal
# 
# 1-sample proportions test with continuity correction
# 
# data:  3262 out of 6516, null probability 0.309
# X-squared = 1119.6, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is greater than 0.309
# 95 percent confidence interval:
#   0.4903506 1.0000000
# sample estimates:
#   p 
# 0.5006139 

prop.test(1319,6516, p = 0.587, alternative = "less") # Mix immune vs AD immune
# 1-sample proportions test with continuity correction
# 
# data:  1319 out of 6516, null probability 0.587
# X-squared = 3973.6, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is less than 0.587
# 95 percent confidence interval:
#   0.000000 0.210813
# sample estimates:
#   p 
# 0.2024248 

prop.test(1319,6516, p = 0.608, alternative = "less") # Mix immune vs N1-sample proportions test with continuity correction
# 
# data:  1319 out of 6516, null probability 0.608
# X-squared = 4495.4, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is less than 0.608
# 95 percent confidence interval:
#   0.000000 0.210813
# sample estimates:
#   p 
# 0.2024248 E immune


path <- "~/normal_cell_annotation/subtype/"

pdf("pdf/subtypes_overview.pdf", width = 10, height = 7)
for (t in c("Tcells", "Bcells", "Macro","Endo", "Epi", "Fibro")){
  subtype <- readRDS(paste0(path, t , ".Rds"))
  f1 <- barplot_bypatient(subtype)
  f2 <- barplot_bytissue_byrow(subtype)
  f3 <- barplot_bypathology(subtype)
  f4 <- coneraxes(DimPlot(subtype, reduction = "harmony_umap", label = T, group.by = "subtype", cols = ggsci::pal_observable()(8)) + 
                    ggtitle(t)+
                    theme(legend.position = "bottom")) +
    theme(legend.margin = margin(l = 50, t = -20))
  print(free(f4) +f3 +f1 +f2+guide_area()+
    plot_layout(ncol = 2, nrow = 3, guides = "collect", heights  = c(3, 3, 0.5)))
  print(t)
}
dev.off()

DefaultAssay(Fibro) <- "RNA"
FeaturePlot(Fibro, reduction = "harmony_umap", c("ATAC2", "IL6", "LY6G6C"))
