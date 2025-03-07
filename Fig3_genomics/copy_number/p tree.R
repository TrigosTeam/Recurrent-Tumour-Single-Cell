library(phangorn)
library(gridExtra)
library(dendextend)
library(ComplexHeatmap)
library("factoextra")
data <- readRDS("~/genomic/copy_number/geneCNAlist.Rds")
data <- rbindlist(data, idcol = "sample")

data <- data[data$sample != "CA0035_left_paraaortic_lymph_node_track_1B", ]
data$variant<- paste(data$chrom, data$g_start, data$g_end, data$CN, sep = ";")
data$patient <- substr(data$sample, 1, 6)


# per patient tree --------------------------------------------------------
ppdf <- split(data, data$patient)

pdf("per_patient_tree.pdf")
ppm <- lapply(names(ppdf), function(p){
  subdf <- ppdf[[p]]
  id <- unique(subdf$sample)
  m <- matrix(ncol = length(id), nrow = length(id), dimnames = list(id, id))
  ftab <- data.frame(table(subdf$variant,subdf$sample))
  for ( i in id){
    temp1 <- ftab$Freq[ftab$Var2==i]
    for (n in id[id!=i]){
      temp2 <- ftab$Freq[ftab$Var2==n]
      m[i, n] <- sum(temp1 != temp2)
    }
  }
  diag(m) <- 0
  
  treeUPGMA  <- as.dendrogram(hclust(as.dist(m), method = "average"))
  par(mar=c(5, 4, 4,20) + 0.1)
  plot(treeUPGMA,horiz = TRUE)
  return(m)
})
dev.off()

Reduce(rbind, sapply(ppm, colMeans))
dend_p = as.dendrogram(hclust(dist(rbind(colMeans(m1), colMeans(m2), colMeans(m3)))))
dend_m = merge_dendrogram(dend_p, list(dend1, dend2, dend3))

# all the patient ---------------------------------------------------------
id<- unique(data$sample)
m<- matrix(ncol = 43, nrow = 43)
colnames(m) <- id
rownames(m) <- id

## withMOIDFIER looks better 
ftab <- data.frame(table(data$variant,data$sample))

for ( i in id){
  temp1 <- ftab$Freq[ftab$Var2==i]
  for (n in id[id!=i]){
    temp2 <- ftab$Freq[ftab$Var2==n]
    m[i, n] <- sum(temp1 != temp2)
  }
}
diag(m) <- 0
m <- as.dist(m)

treeUPGMA  <- upgma(m)
dend <- as.dendrogram(treeUPGMA)

par(mar=c(5, 4, 4,20) + 0.1)
colors <- structure(brewer.dark2(12), names = unique(substr(unique(data$sample),1, 6)))
labs <- labels(dend)
labs <- gsub("_", " ", labs)
dend %>%  set("labels", labs) %>% set("labels_colors", as.character(colors[substr(labels(dend),1, 6)])) %>% plot(horiz = T)
colored_bars(colors = as.character(colors[substr(labels(dend),1, 6)]), horiz = T, y_shift = -2, rowLabels = "")
