library(phylogram)
library(dendextend)
library("beanplot")
library(mixtools)
library(pheatmap)
library(zoo)
library(squash)
library("biomaRt")
library(scran)
library(CONICSmat)
library(Seurat)
library(scales)
library(ggplot2)
# library(infercnv)
library(dplyr)

setwd("~/CASCADEpaper/paper/Fig3/Infercnv/CONIC_PURPLE_sigment")
# run the script of rewrite function 

# rewrite the function  ---------------------------------------------------
BIC.mix = function(out){
  bc=-2*out$loglik+log(length(out$x))*(length(unlist(out[2:4]))-1)
  return(bc)
}

plotChrEnichment2 <- function (expmat, chr, normFactor, gene_positions, n = 1, groups1 = NULL, 
                               groups2 = NULL, start = NULL, end = NULL, k = 2, vis = T, 
                               postProb = 0.95, repetitions = 5, my_para) 
{
  par(mfrow = c(2, 2))
  if (!is.null(groups1)) {
    cellcolor = rep("black", (length(groups1) + length(groups2)))
    cellcolor[groups2] = "red"
  }
  else {
    cellcolor = NULL
  }
  if (!is.null(start)) {
    chr_genes = gene_positions[which(gene_positions[, 3] == 
                                       chr & gene_positions[, 4] > start & gene_positions[, 5] < end), 2]
  }
  else {
    chr_genes = gene_positions[which(gene_positions[, 3] == chr), 2]
  }
  if(length(intersect(chr_genes, row.names(expmat)))>1){# new modification
    if (length(chr_genes) > my_para) {
      chr_exp = scale(colMeans(expmat[intersect(chr_genes, row.names(expmat)), ]) - normFactor)
      bestlog = (-Inf)
      bestmix = NULL
      loglik = NULL
      for (i in 1:repetitions) {
        print(paste("Fitting GMM for chr", chr, " ", start, 
                    ":", end, " iteration ", i, sep = ""))
        mixmdl = tryCatch(mixtools::normalmixEM(chr_exp, 
                                                k = k, maxit = 1000, maxrestarts = 10), error = function(e) {
                                                  print(paste("EM algorithm did not converge for region", 
                                                              chr, " ", start, " ", end))
                                                  mixmdl = NULL
                                                })
        if (!is.null(mixmdl)) {
          if (mixmdl$loglik > bestlog) {
            bestlog = mixmdl$loglik
            bestmix = mixmdl
          }
        }
      }
      if (is.null(bestmix)) {
        hist(chr_exp, breaks = 50, main = paste("Chr: ", 
                                                chr, ":", start, ":", end, "\n", "Unable to fit 2 component mixture model", 
                                                sep = ""), lwd2 = 3, xlab2 = "Expression z-score")
        plot(runif(length(chr_exp), 0, 100), chr_exp, pch = 16, 
             ylab = "Expression z-score", ylim = c(min(chr_exp), 
                                                   (max(chr_exp) + 2)), xlab = "Cells")
        hist(chr_exp, breaks = 50, main = paste("Chr: ", 
                                                chr, ":", start, ":", end, "\n", "Unable to fit 2 component mixture model", 
                                                sep = ""), lwd2 = 3, xlab2 = "Expression z-score")
        plot(runif(length(chr_exp), 0, 100), chr_exp, pch = 16, 
             ylab = "Expression z-score", ylim = c(min(chr_exp), 
                                                   (max(chr_exp) + 2)), xlab = "Cells")
      }
      else {
        out1 = list(x = chr_exp, mu = mean(chr_exp), sigma = sd(chr_exp), 
                    lambda = 1, loglik = sum(dnorm(chr_exp, mean(chr_exp), 
                                                   sd(chr_exp), log = TRUE)))
        bics = c(max(BIC.mix(out1), 1), max(BIC.mix(bestmix), 
                                            1))
        lrt = round(likelihoodRatioTest(out1$loglik, bestmix$loglik, 
                                        n), 6)
        bestmix$BIC = bics
        bestmix$lrt = lrt
        if (vis == T) {
          plot(bestmix, which = 2, breaks = 50, col1 = c("red", 
                                                         "green"), main2 = paste("Chr: ", chr, ":", 
                                                                                 start, ":", end, "\n", "Log likelihood ", round(bestmix$loglik, 
                                                                                                                                 1), sep = ""), lwd2 = 3, xlab2 = "Expression z-score")
        }
        if (length(cellcolor) > 1 & vis == T) {
          g1 = length(which(bestmix$posterior[groups1, 1] > postProb))/length(groups1) * 100
          g2 = length(which(bestmix$posterior[groups1, 2] > postProb))/length(groups1) * 100
          g3 = length(which(bestmix$posterior[groups1, 2] < postProb & bestmix$posterior[groups1,1] < postProb))/length(groups1) * 100
          g4 = length(which(bestmix$posterior[groups2, 1] > postProb))/length(groups2) * 100
          g5 = length(which(bestmix$posterior[groups2, 2] > postProb))/length(groups2) * 100
          g6 = length(which(bestmix$posterior[groups2, 2] < postProb & bestmix$posterior[groups2, 1] < postProb))/(length(groups2) + length(groups1)) * 100
          barplot(rbind(c(g1, g2, g3), c(g4, g5, g6)), 
                  ylim = c(0, 100), beside = T, ylab = "Percentage of cells", 
                  names = c("Cluster", "Cluster", "Ambigu"), 
                  legend = c("Non-malignant", "Malignant"), args.legend = list(title = "Pred. via transcript.", 
                                                                               x = "topright", cex = 0.65), xlab = "Predicted via transcriptomics")
          axis(1, at = c(0.5, 1, 2, 3, 3.3), line = 2, 
               tick = T, labels = rep("", 5), lwd = 3, lwd.ticks = 0, 
               col = "red")
          axis(1, at = c(3.5, 4, 5, 6, 6.5), line = 2, 
               tick = T, labels = rep("", 5), lwd = 3, lwd.ticks = 0, 
               col = "green")
          barplot(bics, names = c("1", "2"), ylab = "BIC", 
                  pch = 16, xlab = "Number of components", log = "y")
          plot(runif(length(chr_exp), 0, 100), chr_exp, 
               pch = 16, col = cellcolor, ylab = "Expression z-score", 
               ylim = c(min(chr_exp), (max(chr_exp) + 2)), 
               xlab = "Cells")
          legend("topright", col = c("black", "red"), c("Non-malignant", 
                                                        "Malignant"), bty = "o", box.col = "darkgreen", 
                 cex = 0.65, pch = 16, title = "Pred. via transcript.")
        }
        else {
          if (vis == T) {
            plot(runif(length(chr_exp), 0, 100), chr_exp, 
                 pch = 16, ylab = "Expression z-score", ylim = c(min(chr_exp), 
                                                                 (max(chr_exp) + 2)), xlab = "Cells")
            barplot(bics, names = c("1", "2"), ylab = "BIC", 
                    pch = 16, xlab = "Number of components", 
                    log = "y")
            hist(bestmix$posterior[, 1], main = "Posterior probablility distribution\n component 1", 
                 xlab = "Posterior probability", breaks = 20, 
                 xlim = c(0, 1))
          }
        }
        return(bestmix)
      }
    }
  }
}

plotAll2 <- function (mat, normFactor, regions, gene_pos, fname, normal = NULL, 
                      tumor = NULL, postProb = 0.8, repetitions = 4, my_para) 
{
  pdf(paste(fname, "_CNVs.pdf", sep = ""))
  loglik = c()
  bic = c()
  lrt = c()
  l = c()
  for (i in 1:nrow(regions)) {
    mixmdl = plotChrEnichment2(mat, regions[i, 1], normFactor, 
                               gene_pos, nrow(regions), normal, tumor, regions[i, 
                                                                               2], regions[i, 3], postProb = postProb, repetitions = repetitions, my_para = my_para)
    if (!is.null(mixmdl)) {
      loglik = c(loglik, mixmdl$loglik)
      bic = c(bic, mixmdl$BIC)
      lrt = c(lrt, mixmdl$lrt)
      names(loglik)[length(loglik)] = rownames(regions)[i]
      names(bic)[length(bic)] = paste(rownames(regions)[i], 
                                      "1_comp", sep = "_")
      names(bic)[length(bic) - 1] = paste(rownames(regions)[i], 
                                          "2_comp", sep = "_")
      if (mixmdl$mu[1] > mixmdl$mu[2]) {
        r = mixmdl$posterior[, 1]
      }
      else {
        r = mixmdl$posterior[, 2]
      }
      l = cbind(l, r)
      colnames(l)[ncol(l)] = rownames(regions)[i]
    }
  }
  par(mfrow = c(1, 1))
  barplot(sort(loglik), names = names(sort(loglik)), cex.axis = 0.8, 
          cex.names = 0.7, las = 2, ylab = "log-likelihood")
  dev.off()
  bicLRmat = matrix(ncol = 4, nrow = length(loglik))
  bicLRmat[, 1] = bic[seq(1, (length(bic) - 1), 2)]
  bicLRmat[, 2] = bic[seq(2, length(bic), 2)]
  bicLRmat[, 3] = bicLRmat[, 1] - bicLRmat[, 2]
  bicLRmat[, ncol(bicLRmat)] = lrt
  colnames(bicLRmat) = c("BIC 1 component", "BIC 2 components", 
                         "BIC difference", "LRT adj. p-val")
  rownames(bicLRmat) = names(loglik)
  rownames(l) = colnames(mat)
  write.table(bicLRmat, paste(fname, "BIC_LR.txt", sep = "_"), 
              sep = "\t")
  return(l)
}

# generate l  -------------------------------------------------------------
paths <- system("realpath /trigos_team/CASCADE/Analysis/230113_seurat_intron/*/*/*non_tumour_removed.Rds", intern = T)
paths_split <- unlist(lapply(strsplit(paths, split = "/"), function(x) paste(x[6], x[7], sep = "_")))
paths <- setNames(paths, paths_split)
cnlist <- readRDS("~/genomic/copy_number/cnlist.Rds")
multi_g <- readRDS("~/genomic/exon_only_reference.Rds")
geneCNAlist <- readRDS("~/genomic/copy_number/geneCNAlist.Rds")
multi_g$enID <- NA
multi_g <- multi_g[,c(5,1,2,3,4)]

names(cnlist)
cnlist <- cnlist[-9]
geneCNAlist <- geneCNAlist[-9]
names(cnlist) <- paths_split
names(geneCNAlist) <- paths_split


for ( i in paths_split){
  
  tumor <- readRDS(paths[i])
  DefaultAssay(tumor) <- "RNA"
  tumor = NormalizeData(tumor, normalization.method = "RC", scale.factor = 1e6)
  data <- GetAssayData(tumor, slot = "data")
  suva_expr <- log2((data/10) + 1)
  suva_expr [is.na(suva_expr)]=0
  # gene_pos=getGenePositions(rownames(suva_expr), ensembl_version = "https://www.ensembl.org")
  gene_pos = multi_g
  normFactor=calcNormFactors(suva_expr)
  suva_expr=filterMatrix(suva_expr,gene_pos$gene.name,minCells=5)
  regions <- cnlist[[i]]%>% dplyr::select(c(chromosome, start, end, bin_width, copyNumber))%>% as.data.frame()
  colnames(regions) <- c("Chrom", "Start", "End", "Length","CN")
  rownames(regions) <- paste(regions$Chrom, regions$Start, regions$End, sep = ";")
  # regions <- readRDS("~/CONICS/unmerged_regions.Rds")
  # regions<- regions[[paste(patient, sample, sep = "_")]]
  # regions <- regions[round(regions$CN)!=2, ]
  my_para = mean(table(geneCNAlist[[i]]$region))
  l <- try(plotAll2(suva_expr,normFactor,regions,gene_pos, paste0("~/CASCADEpaper/paper/Fig3/Infercnv/CONIC_PURPLE_sigment/", i), my_para = my_para), silent = F)
  while (inherits(l, "try-error")){
    my_para <- my_para + 10 
    l <- try(plotAll2(suva_expr,normFactor,regions,gene_pos, paste0("~/CASCADEpaper/paper/Fig3/Infercnv/CONIC_PURPLE_sigment/", i), my_para = my_para), silent = F)
  }
  write.table(my_para, file = paste0("~/CASCADEpaper/paper/Fig3/Infercnv/CONIC_PURPLE_sigment/",i, "_chr_gene.txt"), row.names = F, col.names = F)
  saveRDS(l, paste0("~/CASCADEpaper/paper/Fig3/Infercnv/CONIC_PURPLE_sigment/",i, "_count_mat.Rds"))
  print(i)
  print("DONE")
}


# hi = plotHistogram(l[, sig_regions],suva_expr,zscoreThreshold=2, clusters = 1, patients = subclone_anno$subclone)
# hi2 =  plotHistogram(CA0090_brain_2_count_mat,suva_expr,zscoreThreshold=2, clusters = 1, , patients = subclone_anno$subclone)
# lrbic=read.table(paste(i, "BIC_LR.txt", sep = "_"),sep="\t",header=T,row.names=1,check.names=F)
# candRegions=rownames(lrbic)[which(lrbic[,"BIC difference"]>100 & lrbic[,"LRT adj. p-val"]<0.01)]
# candRegions=rownames(lrbic)[which(lrbic[,"LRT adj. p-val"]<0.01)]
# hi = plotHistogram(l[, candRegions],suva_expr,zscoreThreshold=3, clusters = 2)
# class(hi)
# 
# dend <- read.dendrogram("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/unclustered/CA0027_dura_base_skull_13/infercnv.observations_dendrogram.txt")
# plot(dend)
# tree <- cutree(dend, k = 3)
# dend2 <- color_labels(dend, labels = names(hi[hi == 2]), col = 2) 
# plot(dend2)
# 
# r=generatePvalMat(suva_expr,regions[candRegions,],normFactor,normal,tumor,gene_pos,threshold=0.8)
# binr=ifelse(r>0.1,0,1)
# boxplot(r)

