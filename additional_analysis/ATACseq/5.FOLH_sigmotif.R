library(dplyr)
library(data.table)
library(ComplexHeatmap)

paths <- system("realpath ~/ATAC/*/FOLH1_sig.motifs.Rds", intern = T)
paths <- setNames(paths, sapply(strsplit(paths, split  = "/"), '[', 5))

paths1 <- system("realpath ~/ATAC/*/FOLH1_DE_peaks.Rds", intern = T)
paths1 <- setNames(paths1, sapply(strsplit(paths1, split  = "/"), '[', 5))
paths2 <- system("realpath ~/ATAC/*/*multi_srt.qs", intern = T)
paths2 <- setNames(paths2, sapply(strsplit(paths2, split  = "/"), '[', 5))


# "CA0027_dura_inner_skull_14" "CA0058_liver_38"   have no DE peaks


motif <- lapply(paths, readRDS)
motifdf <- lapply(motif, function(x) rbindlist(x, idcol = "group"))
motifdf <- rbindlist(motifdf, idcol = "sample", fill = T)
write.table(motifdf, "persample_FOLH1_sigmotif.csv", quote = F, row.names= F, sep = "\t")
head(motifdf)
motifdf <- motifdf[motifdf$p.adjust <0.05,]
countable <- as.data.frame.matrix(table(motifdf$group, motifdf$motif.name))

mat <- matrix("", ncol = n_distinct(motifdf$sample), nrow = n_distinct(motifdf$motif.name), dimnames = list(unique(motifdf$motif.name), unique(motifdf$sample)))
for (i in seq(nrow(motifdf))){
  mat[motifdf$motif.name[i], motifdf$sample[i]] <-  paste(mat[motifdf$motif.name[i], motifdf$sample[i]], motifdf$group[i], sep = ";")
}

mat <-  sub("^;", "", mat)


col = c(high = "salmon", low = "green4", neg = "royalblue")
source("~/CASCADEpaper/paper/cols.R")

alter_fun = list(
  background = function(x, y, w, h) grid.rect(x, y, w, h,  gp = gpar(fill = NA, col = NA)), 
  high = function(x, y, w, h) {
    grid.polygon(
      unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
      unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
      gp = gpar(fill = "salmon", col = "white"))
  },
  low = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
      unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
      gp = gpar(fill = "green4", col = "white"))
  }, 
  # neg = function(x, y, w, h) {
  #   grid.segments(x - w*0.4, y - h*0.4, x + w*0.4, y + h*0.4, gp = gpar(lwd = 2))
  #   grid.segments(x + w*0.4, y - h*0.4, x - w*0.4, y + h*0.4, gp = gpar(lwd = 2))
  # }
  neg = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.3, gp = gpar(fill = "royalblue", col = NA))
)

alter_fun = list(
  background = function(x, y, w, h) grid.rect(x, y, w, h,  gp = gpar(fill = NA, col = NA)), 
  # High: Left bar
  high = function(x, y, w, h) {
    # Shift x left by w/3
    grid.rect(x - w*1/3, y, w*1/3, h, 
              gp = gpar(fill = col["high"], col = NA))
  },
  
  # Low: Middle bar
  low = function(x, y, w, h) {
    # x stays at center
    grid.rect(x, y, w*1/3, h, 
              gp = gpar(fill = col["low"], col = NA))
  },
  
  # Neg: Right bar
  neg = function(x, y, w, h) {
    # Shift x right by w/3
    grid.rect(x + w*1/3, y, w*1/3, h, 
              gp = gpar(fill = col["neg"], col = NA))
  }
)

patients <- gsub("00", "", substr(colnames(mat), 1, 6))
meta <- readRDS("~/integration/202406/tumor_only_meta_expression.Rds")
FOLH1 <- meta %>% group_by(patient) %>% summarise(mean = mean(FOLH1))
FOLH1$patient <- gsub("00", "", FOLH1$patient)

ind <- rowSums(nchar(mat)>0)>2 & apply(nchar(mat)>0, 1, function(x) n_distinct(patients[x])>1)

ha <- HeatmapAnnotation(
  patient = patients, 
  FOLH1 = FOLH1$mean[match(patients, FOLH1$patient)],
  col = list(patient = patient_cols)
)

ht <- oncoPrint(mat[ind, ], alter_fun = alter_fun, col = col, 
          top_annotation = ha, 
          column_order = order(FOLH1$mean[match(patients, FOLH1$patient)]))
pdf("FOLH1_motif_full.pdf", width = 7, height = 20)
draw(ht, merge_legend = T)
dev.off()

#CA43 has some problems with lots of FOLH1 expression and conflict motif enrichment
#remove all coflicting motif 
ind2 <-  apply(mat ,1, function(x) !(any(grepl("high", x))&any(grepl("neg", x))))

ht2 <- oncoPrint(mat[ind&ind2, ], alter_fun = alter_fun, col = col, 
                top_annotation = ha, 
                column_order = order(FOLH1$mean[match(patients, FOLH1$patient)]))
pdf("FOLH1_motif.pdf", width = 7, height = 8)
draw(ht2, merge_legend = T)
dev.off()