library(readxl)
library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
microscopy <- read_excel("~/CASCADEpaper/paper/normal_cells_202406/Microscopy_annotation_RY Feb2024.xlsx", sheet = 3)
microscopy <- as.data.frame(microscopy)
final_normalsrt <- readRDS("~/CASCADEpaper/paper/normal_cells_202406/subtype/final_normalsrt.Rds")
Tcells <- readRDS("~/CASCADEpaper/paper/normal_cells_202406/subtype/Tcells.Rds")
signature_meta_tumor_normal <- readRDS("~/CASCADEpaper/paper/Fig2_signature/signature_meta_tumor_normal.Rds")
temp <- signature_meta_tumor_normal %>% group_by(sample) %>% summarise(n = n())
tot_cll <- setNames(temp$n, temp$sample)

df2 <- Tcells@meta.data %>% group_by(sample, subtype) %>% summarise(n = n()) 

df <-  final_normalsrt@meta.data%>% group_by(sample, cell_anno) %>% summarise(n = n()) 
df <- df %>% filter(cell_anno %in% c("T cells", "B cells", "Macrophages"))
colnames(df2) <- colnames(df)
df <- rbind(df, df2)
df$patient <- gsub("00", "", substr(df$sample, 1, 6))
df$prop <- df$n/tot_cll[df$sample]


for (m in c("CD3", "CD8" ,"CD11c","CD20","CD68")){
  microscopy[, m] <- as.numeric(microscopy[,m])
  microscopy[, paste0("quan_",m)] <- ntile(microscopy[, m], 5)
  # microscopy[, paste0("quan_",m)] <- factor(microscopy[, paste0("quan_",m)])
}


matchl <- setNames(c("CD3", "CD8" ,"CD20","CD68"),c("T cells","CD8 memory/naive cells", "B cells", "Macrophages") )

pdf("~/CASCADEpaper/paper/normal_cells_202406/microscopy_cor_cellprop.pdf", width =9, height = 4)
for(i in c("T cells","CD8 memory/naive cells", "B cells", "Macrophages")){
  df2 <- df[df$cell_anno == i, ]
  m <- as.character(matchl[i])
    subdf <- merge(df2, microscopy[, c("srt obj sites", m, paste0("quan_",m))], by.x = "sample", by.y = "srt obj sites")
    colnames(subdf) <- c("sample","cell_anno","n","patient","prop", "marker", "quan")
    p <- DescTools::JonckheereTerpstraTest(subdf$prop, g = subdf$quan, alternative = "increasing", nperm = 5000)
    corR <-  cor.test(subdf$prop[!is.na(subdf$marker)], subdf$marker[!is.na(subdf$marker)], method = "spearman")
    f <- ggplot(subdf, aes(x = factor(quan), y = prop)) + 
      geom_boxplot(outlier.shape = NA) + geom_jitter(aes(colour = patient)) + 
      labs(x = paste("Bins of", m, "Percentage"), y = paste0(i, "\nProportion"), subtitle = paste("JT p value=",signif(p$p.value, 2)))+
      theme_bw(base_size = 18)
    
    g <- ggplot(subdf, aes(x = marker, y = prop ,colour = patient)) + geom_point() +  
      geom_smooth(method=lm, color="black", se=FALSE)+ 
      labs(x = paste( m, "Percentage"), y = paste0(i, "\nProportion"),subtitle = paste("rho=", signif(corR$estimate, 2), "p value=", signif(corR$p.value, 2)) )+
      theme_bw(base_size = 18)
    print(f+g+plot_layout(guides = "collect", axes = "collect"))
  
}
dev.off()

