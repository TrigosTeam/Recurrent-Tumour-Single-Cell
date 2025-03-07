library(StatCharrms)
library(readxl)
library(ggplot2)
microscopy <- read_excel("~/CASCADEpaper/paper/Fig6/microscopy estimation.xlsx")
phenotype_meta <- readRDS("~/integration/tumor_only_phenotype_meta.Rds")
b.cells <- normal2@meta.data[normal2$cell_id2 == "B cell",]
t.cells <- normal2@meta.data[normal2$cell_id2 == "T cell",]
tb.t <- as.data.frame(table(t.cells$sample)) %>% `colnames<-`(c("sample", "Tcell_n"))
tb.b <- as.data.frame(table(b.cells$sample)) %>% `colnames<-`(c("sample", "Bcell_n"))
tb2 <- as.data.frame(table(phenotype_meta$sample))%>% `colnames<-`(c("sample", "tumor_n"))
df <- merge(tb2, tb.t, by = "sample", all = T)
df <- merge(df, tb.b, by = "sample", all = T)
normal.tb <- as.data.frame(table(normal2$sample))%>% `colnames<-`(c("sample", "normal_n"))
df <- merge(df, normal.tb, by = "sample", all = T )

df$Tcell_n[is.na(df$Tcell_n)] <- 0
df$Bcell_n[is.na(df$Bcell_n)] <- 0
df <- df %>% mutate(sum = normal_n + tumor_n) %>% mutate(Tcell_prop = Tcell_n/sum, Bcell_prop = Bcell_n/sum)
df <- df[-28, ]# CA83 fat has no normal cells
microscopy <- microscopy %>% filter(FF_ID %in% df$sample)
microscopy <- microscopy %>% mutate(quan_CD3  = ntile(CD3, 5), quan_CD8  = ntile(CD8, 5), quan_CD20  = ntile(CD20, 5))

plotdf <- merge(df, microscopy, by.x = "sample", by.y ="FF_ID")
plotdf$patient <- substr(plotdf$sample, 1, 6)

pdf("microscopy_cor_cellprop.pdf", width = 6, height = 5)
p <- DescTools::JonckheereTerpstraTest(plotdf$Tcell_prop, g = plotdf$quan_CD3, alternative = "increasing", nperm = 5000)
ggplot(plotdf, aes(x = as.factor(quan_CD3), y = Tcell_prop)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(aes(colour = patient)) + scale_color_manual(labels  = gsub("00", "", unique(plotdf$patient)), values = brewer.set1(9))+
  labs(x = "bins of CD3 intensity", y = "T cell proportion", title = paste("correlation between CD3 stain and T cell\nproportion:", round(cor(plotdf$Tcell_prop, plotdf$CD3), 3)))+theme_bw()+
annotate("text", x = 3, y = 0.8*max(plotdf$Tcell_prop), label = paste0("JonckheereTerpstraTest p-value:", p$p.value))

ggplot(plotdf, aes(x = CD3, y = Tcell_prop,colour = patient)) + geom_point() +  geom_smooth(method=lm, color="black", se=FALSE)+ labs(x = "CD3 intensity", y = "T cell proportion", title = paste("correlation between CD3 stain and T cell proportion:", round(cor(plotdf$Tcell_prop, plotdf$CD3), 3)))+theme_bw()

p <- DescTools::JonckheereTerpstraTest(plotdf$Tcell_prop, g = plotdf$quan_CD8, alternative = "increasing", nperm = 5000)
ggplot(plotdf, aes(x = as.factor(quan_CD8), y = Tcell_prop)) + geom_boxplot(outlier.shape = NA) + geom_jitter(aes(colour = patient)) + labs(x = "bins of CD8 intensity", y = "T cell proportion", title = paste("correlation between CD8 stain and T cell\nproportion:", round(cor(plotdf$Tcell_prop, plotdf$CD8), 3)))+theme_bw()+scale_color_manual(labels  = gsub("00", "", unique(plotdf$patient)), values = brewer.set1(9))+
  annotate("text", x = 3, y = 0.8*max(plotdf$Tcell_prop), label = paste0("JonckheereTerpstraTest p-value:", p$p.value))

ggplot(plotdf, aes(x = CD8, y = Tcell_prop,colour = patient)) + geom_point() +  geom_smooth(method=lm, color="black", se=FALSE)+ labs(x = "CD8 intensity", y = "T cell proportion", title =paste("correlation between CD8 stain and T cell proportion:", round(cor(plotdf$Tcell_prop, plotdf$CD8), 3)))+theme_bw()


p <- DescTools::JonckheereTerpstraTest(plotdf$Bcell_prop, g = plotdf$quan_CD20, alternative = "increasing", nperm = 5000)
ggplot(plotdf, aes(x = as.factor(quan_CD20), y = Bcell_prop)) + geom_boxplot(outlier.shape = NA) + geom_jitter(aes(colour = patient)) + scale_color_manual(labels  = gsub("00", "", unique(plotdf$patient)), values = brewer.set1(9))+labs(x = "bins of CD20 intensity", y = "B cell proportion", title = paste("correlation between CD20 stain and B cell\nproportion:", round(cor(plotdf$Bcell_prop, plotdf$CD20), 3)))+theme_bw()+
  annotate("text", x = 3, y = 0.8*max(plotdf$Bcell_prop), label = paste0("JonckheereTerpstraTest p-value:", p$p.value))

ggplot(plotdf, aes(x = CD20, y = Bcell_prop,colour = patient)) + geom_point() +  geom_smooth(method=lm, color="black", se=FALSE)+ labs(x = "CD20 intensity", y = "B cell proportion", title =paste("correlation between CD20 stain and B cell proportion:", round(cor(plotdf$Bcell_prop, plotdf$CD20), 3)))+theme_bw()
dev.off()


library(DescTools)


