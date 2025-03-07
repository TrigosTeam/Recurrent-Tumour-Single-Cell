#' ---
#' title: "oncoplots with limited genelist"
#' output: html_notebook
#' ---
#' 
#' 
## ----set up, include=FALSE--------------------------------------------------------------
geneCNA <- readRDS("~/CASCADEpaper/paper/Fig3_genomics/copy_number/geneCNAlist.Rds")
for (i in names(geneCNA)){
  geneCNA[[i]]$sample <- i 
}
library(dendextend)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(data.table)
library(readr)
library(pals)

sv <- readRDS("~/CASCADEpaper/paper/Fig3_genomics/structure_variants/sv_list.Rds")
SO_tablelist <- readRDS("~/CASCADEpaper/paper/Fig3_genomics/point_mutation/point_mutation_table_intersection.Rds")

#' 
#' ### genelist generated
## ----eval=FALSE, include=FALSE----------------------------------------------------------
DDR <- unique(c("BRCA1", "BRCA2", "ATM", "BRIP1", "BARD1", "CDK12", "CHEK1", "CHEK2", "FANCL", "PALB2", "PPP2R2A", "RAD51C", "RAD51D", "RAD54L", #HRR PROfound genes RAD51B
             "BRAC1", "BRAC2", "PALB2", "BRIP", "FANC", "RAD51B", "RAD51C", "RAD51D", "RAD54L", "CHEK1", "ATM", "CDK12", "PALB2", "RAD51", "FANCA", "HDAC2", "BRIP1", "MLH1", "MRE11A", "NBN"))
AR_regulation <- c("AR", "PTEN", "MYC", "NCOA2", "PTK2", "RB1", "APC", "FOXA1", "ZMIZ1", "CCND1", "SIRT1", "NCOR1", "ETV5", "FOXO1", "MDM2", "SPOP", "SOX2", "GSK3B", "PIAS3", "PIAS4", "RAD9A",  "CREB1",  "CREBBP", "CTNNB1","AKT1","AES", "FOLH1")
primary_matas <- c("TP53",  "ZMYM3","PTPRD", "ZFP36L2", "ADAM15",  "BRIP1", "APC", "KMT2C", "CCAR2", "NKX3-1", "C8orf58", "RYBP")
genelist <- c(AR_regulation, primary_matas, DDR)

#' "C11orf30" "FAM123B"  "FAM175A"  "FAM46C"   "GPR124"   "MLL"      "MRE11A"   "MYCL1"    "PAK7"     "PARK2"    "RFWD2"    "TCEB1"    "WISP3"
#' not included because lack CN annotation 
#' 
#' 
#' 
#' ## CN matrix
## ----CN matrix, include=FALSE-----------------------------------------------------------
CN.mat <- lapply(geneCNA, function(x){
  temp <- x[x$gene%in% genelist, c("gene", "CN")]
  if(any(temp$CN >=10)){
    temp$CN[temp$CN>=10] <- 10
  }
  colnames(temp) <- c("gene", unique(x$sample))
  return(temp)
})

CN.mat <- CN.mat %>% Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="gene"), .) 
CN.mat <- as.matrix(CN.mat)
rownames(CN.mat) <- CN.mat[,1]
CN.mat <- CN.mat[,2:ncol(CN.mat)]


rev(colorRampPalette(brewer.pal(9, "RdBu"))(9))
col.CN <- c("#5DA3CC", "#BCDAEA", "#F7F7F7", tail(rev(colorRampPalette(brewer.pal(9, "RdBu"))(17)), 8))
names(col.CN) <- c(as.character(0:10))
# col.CN = colorRamp2(c(0,2,25), c("#5DA3CC", "#F7F7F7", "#B2182B"))
class(CN.mat) <- "numeric"
gene_order <- order.dendrogram(as.dendrogram(hclust(dist(CN.mat))))
sample_order <- order.dendrogram(as.dendrogram(hclust(dist(t(CN.mat)))))
# CN.list <- list("0" = alter_graphic("rect", width = 0.9, height = 0.9, fill = col.CN["0"]),
#                 "1" = alter_graphic("rect", width = 0.9, height = 0.9, fill = col.CN["1"]),
#                 "2" = alter_graphic("rect", width = 0.9, height = 0.9, fill = col.CN["2"]),
#                 "3" = alter_graphic("rect", width = 0.9, height = 0.9, fill = col.CN["3"]),
#                 "4" = alter_graphic("rect", width = 0.9, height = 0.9, fill = col.CN["4"]),
#                 "5" = alter_graphic("rect", width = 0.9, height = 0.9, fill = col.CN["5"]),
#                 "6" = alter_graphic("rect", width = 0.9, height = 0.9, fill = col.CN["6"]))
CN.list <- lapply(col.CN, function(x) return(alter_graphic("rect", width = 0.9, height = 0.9, fill = x)))
# test_alter_fun(CN.list)

#' 
#' ## sv matrix 
## ----structure variant------------------------------------------------------------------
f_sv <- rbindlist(sv, idcol = "sample")
SV.mat <- CN.mat
class(SV.mat) <- "character"
for(i in rownames(CN.mat)){
  for( n in colnames(CN.mat)){
    cons <- f_sv %>% filter(sample == n) 
    cons <- cons[union(which(sapply(strsplit(cons$start_gene, ";"), function(x) i%in%x)), which(sapply(strsplit(cons$end_gene, ";"), function(x) i%in%x))),]
    if(nrow(cons)>0){
      SV.mat[i, n] <- paste(c(SV.mat[i,n], unique(cons$sv)), collapse  = ";")
    }
  }
}


SV_anno <- list(
  inversion = function(x, y, w, h) 
        grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = NA, lwd = 1)),
  translocation = function(x, y, w, h) 
        grid.segments(x + w*0.4, y , x - w*0.4, y , gp = gpar(lwd = 1)),
  tandem_duplication= function(x, y, w, h) 
        grid.segments(x - w*0.4, y - h*0.4, x + w*0.4, y + h*0.4, gp = gpar(lwd = 1)),
  deletion = function(x, y, w, h) 
       grid.segments(x + w*0.4, y - h*0.4, x - w*0.4, y + h*0.4, gp = gpar(lwd = 1)),
  insertion = function(x, y, w, h) 
        grid.segments(x , y - h*0.4, x , y + h*0.4, gp = gpar(lwd = 1))
)
test_alter_fun(SV_anno)



#' ## wholegenome doubling 
## ----message=FALSE, warning=FALSE, include=FALSE----------------------------------------
# genf_anno <- list("TURE"= alter_graphic("rect", width = 0.9, height = 0.9, fill = "black"),
#                   "FALSE" = alter_graphic("rect", width = 0.9, height = 0.9, fill = "grey"))
# paths <- system("realpath /trigos_team/CASCADE/Sequencing_updated_reference/CA*/*/purple/*purity.tsv", intern = T)
# paths_split <- data.table::transpose(strsplit(paths, "/"))
# ids <- paths_split[[5]]
# sites <- paths_split[[6]]
# purplereport <- lapply(paths, read_tsv)
# wgd <- unlist(lapply(purplereport , function(x) return(x$wholeGenomeDuplication)))
# wgd <- setNames(wgd, names(geneCNA))
# saveRDS(wgd, "wgd.Rds")
wgd <- readRDS("~/CASCADEpaper/paper/Fig3_genomics/wgd.Rds")
#' 
#' ## Viriant class matrix 
## ----point mutation, message=FALSE, warning=FALSE, include=FALSE------------------------
# tmb <- unlist(lapply(purplereport , function(x) return(x$tmbPerMb)))
# names(tmb) <- names(geneCNA)
SO_table <- rbindlist(SO_tablelist, idcol = "sample")


f_SO <- SO_table %>% filter(SYMBOL%in%genelist&IMPACT%in% c("HIGH", "MODERATE")) %>% filter(as.numeric(MAX_AF )<0.1 | is.na(as.numeric(MAX_AF)))
# f_SO$IMPACT <- factor(f_SO$IMPACT, levels = c("HIGH", "MODERATE", "LOW", "MODIFIER"))
SM.mat <- SV.mat
class(SM.mat) <- "character"
SOs <- c()
for(i in rownames(CN.mat)){
  for( n in colnames(CN.mat)){
    cons <- f_SO[f_SO$SYMBOL == i & f_SO$sample == n,]
      SM.mat[i, n] <- paste(SM.mat[i,n], names(which.max(table(cons$Consequence))), sep = ";")
      SOs <- c(SOs, names(which.max(table(cons$Consequence))))
  }
}

SO_col <- c("deeppink", "limegreen", "tan1", "deepskyblue", "slateblue1")#, "khaki1", "seagreen", "orangered2")
names(SO_col) <- unique(SOs)
SO_anno <- lapply(SO_col, function(x) return(alter_graphic("rect", width = 0.4, height = 0.4, fill = x)))

#' 
#' 
#' general matrix generation
## ----mat, fig.height=10, fig.width=6, warning=FALSE-------------------------------------
mat <- rbind(SM.mat)
#mat <- mat[gene_order, ]
row_col <- rep("black", nrow(mat))
row_col[which(rownames(mat)%in%DDR)] <- "red"
mat2 <- mat[genelist[genelist%in%rownames(mat)], ]
gene_order <- order.dendrogram(as.dendrogram(hclust(dist(CN.mat[genelist[genelist%in%rownames(mat)],]))))

alter_fun <- c(background = function(...) NULL,CN.list, SO_anno,SV_anno)
# gene_order <- c(gene_order, "TMPRSS2_ERG_fusion", "WGD")

patients <- substr(colnames(mat), 1, 6)

rowanno <- c(rep("AR regulation",sum(AR_regulation%in% rownames(mat))), rep("primary vs metastasis", sum(primary_matas%in%rownames(mat))),
             rep("DDR", sum(DDR%in%rownames(mat))))
colnames(mat2)


#' 
#' legend setting
## ---------------------------------------------------------------------------------------
lgd1 = Legend(at = c(as.character(0:9), ">=10"), title = "CNA", legend_gp = gpar(fill = col.CN))
lgd2 = Legend(at = c("Yes", "No"), title = "WGD", legend_gp = gpar(fill = c("Yes" = "black", "No" = "grey")))
lgd3 = Legend(labels = gsub("_", " ", names(SO_col)), title = "SNV consequence",graphics = SO_anno)
lgd4 = Legend(labels = c("inversion","translocation", "tandem duplication" ,"deletion" , "insertion"), title = "Structure variants", 
              graphics = SV_anno)
lgd5 = Legend(at = unique(patients), title = "Patients", legend_gp = gpar(fill = brewer.set3(length(unique(patients)))))
pd <- packLegend(lgd1, lgd3, lgd4)
draw(pd)

## ---------------------------------------------------------------------------------------
sites<- gsub("CA[0-9]+_","", colnames(mat2))

sites[grep("brain", sites)] <- "brain"
sites[grep("dura", sites)] <- "dura"
sites[grep("prostate", sites)] <- "prostate"
sites[grep("liver", sites)] <- "liver"
sites[grep("fat", sites)] <- "fat"
sites[grep("lymph|hilar|LN", sites)] <- "LN"
sites[grep("rib|vertebra|Vertebra", sites)] <- "bone"
sites[grep("lung", sites)] <- "lung"
sites[grep("abdomen", sites)] <- "abdomen"
sites[grep("bladder", sites)] <- "bladder"

#' 
#' right_annotation = NULL,
#'     left_annotation = rowAnnotation(legend = anno_empty(border = F), width = unit(5, "cm"), height = unit(10, "cm"))
## ----onco, fig.height=10, fig.width=10, message=FALSE, warning=FALSE--------------------
onc <- oncoPrint(mat2, alter_fun = alter_fun, name = "part1", 
          col = c(col.CN, SO_col), show_column_names = F, row_names_gp = gpar(fontsize = 7),show_pct = F, row_names_side = "left",
          row_order =  gene_order, column_names_rot = -90,column_order = colnames(mat2), 
          top_annotation = HeatmapAnnotation(WGD = wgd, col = list(WGD = c("TRUE" = "black", "FALSE" = "grey")),
                                             CN_composition = anno_oncoprint_barplot(as.character(0:10)),
                                             annotation_height =  unit(c(5,10), "mm") , annotation_name_gp = gpar(fontsize =7), show_annotation_name = T, gap = unit(1, "mm")), 
          bottom_annotation = HeatmapAnnotation(Patient = patients, Site = sites, col = list(Patient = structure(brewer.set1(length(unique(patients))), names = unique(patients)), Site = structure(brewer.set3(length(unique(sites))), names = unique(sites)))), 
          row_split = rowanno, border = T, row_title_gp = gpar(fontsize = 9, fontface = "bold")
          # top_annotation = HeatmapAnnotation( TMB_per_Mb =  anno_barplot(tmb), CN_composition = anno_oncoprint_barplot(as.character(0:6)),
          #   height = unit(2, "cm"),  annotation_name_gp = gpar(fontsize =8), show_annotation_name = T, gap = unit(2, "mm")),
    )



draw(onc, show_heatmap_legend = F,  ht_gap = unit(7, "mm"), annotation_legend_list = pd)

pdf("~/CASCADEpaper/paper/Fig3_genomics/oncoplot_allCASCADE.pdf", width = 13, height = 11)
draw(onc, show_heatmap_legend = F,  ht_gap = unit(7, "mm"), annotation_legend_list = pd, padding = unit(c(3, 1, 1, 1), "cm"))

dev.off()
#draw(pd, x = unit(7, "cm"), y = unit(20, "cm"),just = c("right"), test = F)
# dev.off()
ind <- grep("CA0063|CA0079", colnames(mat2))
mat3 <- mat2[ , ind]
onc2 <- oncoPrint(mat3, alter_fun = alter_fun, name = "variants", 
          col = c(col.CN, SO_col), show_column_names = T, row_names_gp = gpar(fontsize = 7),show_pct = F, row_names_side = "left",
          row_order =  gene_order, column_names_rot = -45,column_order = colnames(mat3), 
          top_annotation = HeatmapAnnotation(WGD = wgd[ind], col = list(WGD = c("TRUE" = "black", "FALSE" = "grey")),
                                             CN_composition = anno_oncoprint_barplot(as.character(0:10)),
                                             annotation_height =  unit(c(5,10), "mm") , annotation_name_gp = gpar(fontsize =7), show_annotation_name = T, gap = unit(1, "mm")), 
          bottom_annotation = HeatmapAnnotation(Patient = patients[ind], Site = sites[ind], col = list(Patient = structure(brewer.set1(length(unique(patients))), names = unique(patients)), Site = structure(brewer.set3(length(unique(sites))), names = unique(sites)))), 
          row_split = rowanno, border = T, row_title_gp = gpar(fontsize = 9, fontface = "bold")
          # top_annotation = HeatmapAnnotation( TMB_per_Mb =  anno_barplot(tmb), CN_composition = anno_oncoprint_barplot(as.character(0:6)),
          #   height = unit(2, "cm"),  annotation_name_gp = gpar(fontsize =8), show_annotation_name = T, gap = unit(2, "mm")),
)

pdf("~/genomic/oncoplot_CA6379.pdf", width = 7, height = 11)
draw(onc2, show_heatmap_legend = F,  ht_gap = unit(7, "mm"), annotation_legend_list = pd, padding = unit(c(3, 1, 1, 1), "cm"))

dev.off()
