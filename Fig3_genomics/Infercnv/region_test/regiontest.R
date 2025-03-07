library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(circlize)
library(scales)
library(pals)

paths <- system("realpath ~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/subclone/*/17_HMM_predHMMi6.hmm_mode-samples.pred_cnv_regions.dat", intern = T)
samples <- sapply(strsplit(paths, split = "/"), function(x){return(x[10])})
hmm <- lapply(paths, function(x){
  temp <- read.table(x,header = T, sep = "\t")
  temp <- temp[nchar(temp$cell_group_name)<10,]
  return(temp)
})

names(hmm) <- samples


range_create <- function(temp){
  ID <- unique(substr(temp$cell_group_name, 1, 3))
  x <- split(temp, temp$chr)
  mydf <- lapply(x, function(y){
    v = c(y[,"start"], y[,"end"])
    v = v[order(v, decreasing = F)]
    start <- c()
    end <- c()
    for (i in seq(length(v)-1)){
      start <- c(start, v[i])
      end <- c(end, v[i+1])
    }
    df <- data.frame(start = start, end = end)
    df <- df[(df$end- df$start)>0, ]
    return(df)
  })
  temp1 <- lapply(ID, function(c){
    return(rbindlist(mydf, idcol = "chr") %>% mutate(cell_group = c))
  })
  return(rbindlist(temp1))
}

diff_cn_region <- lapply(hmm, function(z){
  df <- range_create(z)
  for ( i in seq(nrow(df))){
    x <- df[i, ]
    temp <- z %>% filter(grepl(x$cell_group, z$cell_group_name)) %>% filter(chr == x$chr)
    cn <- temp$state[temp$start<= x$start&temp$end>=x$end]
    if( length(cn)>0){
      df$state[i] <- cn
      if(length(cn)>1){
        print(x)
        print(temp)
      }
    }else{
      df$state[i] <- 3
    }
  }
  df$region_state <- paste(df$chr,df$start, df$end, df$state, sep  = ";")
  rm_region <- names(which(table(df$region_state)==length(unique(df$cell_group))))
  df <- df[!df$region_state%in%rm_region,]
  df$region <- paste(df$chr, df$start, df$end, sep = ";")
  
  region_list <- split(df, df$region)
  region_list <- lapply(region_list, function(x){
    temp <- x[1, 1:3]
    for (i in seq(nrow(x))){
      temp[1, x$cell_group[i]]<- x$state[i]
    }
    return(temp)
  })
  final_df <- rbindlist(region_list, idcol = "region", fill = T) 
  final_df <- na.omit(final_df)
  final_df$length <- final_df$end-final_df$start
  return(final_df)
})

col_fun= colorRamp2(c(0, 3, 6), c("#00009D","white","#9D0002"))



regiontest_df <- list()
pdf("~/CASCADEpaper/paper/Fig3/Infercnv/region_test/regiontest.pdf", width = 8, height = 5)
for(i in samples){
  test <- diff_cn_region[[i]]
  test$chr <- as.numeric(test$chr)
  test$group <- "infercnv"
  test <- melt(test, id = c("region","chr", "start", "end", "length", "group"))
  col <- col_fun(seq(min(test$value), max(test$value)))
  
  lrbic=read.table(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/CONIC_PURPLE_sigment/",i, "_BIC_LR.txt"),sep="\t",header=T,row.names=1,check.names=F)
  sig_region =rownames(lrbic)[which(lrbic[,"LRT adj. p-val"]<0.01)]

  test2 <- as.data.frame(matrix(data = unlist(strsplit(sig_region, split = ";")), nrow = length(sig_region), ncol = 3, byrow = T))
  colnames(test2) <- c("chr", "start", "end")
  test2 <- as.data.frame(apply(test2, 2, as.numeric))
  test2$group <- "CONIC_region"
  test2$length <- test2$end - test2$start
  
  g <- rbindlist(list(test, test2), fill = T)
  g$variable <- as.character(g$variable)
  g$variable[which(is.na(g$variable))] <- "CONICS"
  p <- ggplot(g, aes(x=start, y = factor(variable), width = length)) + 
    geom_tile(aes(fill = factor(value)))+ 
    facet_grid(.~chr, scales = "free")+ggtitle(i)+
    scale_fill_manual(values = col, na.value = "yellow", labels = c(levels(factor(g$value)), "significant"))+
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(), 
          panel.spacing = unit(0.1, "lines"),
          panel.background = element_rect(fill = 'black', colour = 'black'),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())+labs(fill = "CN", x= "region", y = "subclone")
  print(p)
  print(i)
  regiontest_df[[i]] <- g
}
dev.off()

View(g)

# find overlaps between infercnv prediction and CONICS --------------------


final_region <- list()
for(i in samples){ 
  test <- diff_cn_region[[i]]
  test$chr <- as.numeric(test$chr)
  test$group <- "infercnv"

  lrbic=read.table(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/CONIC_PURPLE_sigment/",i, "_BIC_LR.txt"),sep="\t",header=T,row.names=1,check.names=F)
  sig_region =rownames(lrbic)[which(lrbic[,"LRT adj. p-val"]<0.01)]
  
  test2 <- as.data.frame(matrix(data = unlist(strsplit(sig_region, split = ";")), nrow = length(sig_region), ncol = 3, byrow = T))
  colnames(test2) <- c("chr", "start", "end")
  test2 <- as.data.frame(apply(test2, 2, as.numeric))
  test2$group <- "CONIC_region"
  test2$length <- test2$end - test2$start
  
  comb_test <- rbind(test[, c("chr", "start", "end")], test2[, c("chr", "start", "end")]) %>% arrange(chr)
  x <- split(comb_test, comb_test$chr)
  mydf <- lapply(x, function(y){
    v = c(y$start, y$end)
    v = v[order(v, decreasing = F)]
    start <- c()
    end <- c()
    for (i in seq(length(v)-1)){
      start <- c(start, v[i])
      end <- c(end, v[i+1])
    }
    df <- data.frame(start = start, end = end)
    df <- df[(df$end- df$start)>0, ]
    return(df)
  })
  temp <- rbindlist(mydf, idcol = "chr")
  temp <- apply(temp, 2, as.numeric)
  temp <- as.data.frame(temp)
  
  
  finaldf <- apply(temp, 1, function(temp_df){
    start <-  as.numeric(temp_df["start"])
    end <-as.numeric(temp_df["end"])
    subtest <- as.data.frame(test[test$chr == as.numeric(temp_df["chr"]),])
    subtest2 <- as.data.frame(test2[test2$chr == as.numeric(temp_df["chr"]),])
    ind1 <- (subtest$start <= start & subtest$end>= start)&(subtest$start <= end & subtest$end >= end)
    ind2 <- ( subtest2$start <= start & subtest2$end>= start)&(subtest2$start <= end & subtest2$end >= end)
    if(any(ind1)&any(ind2)){
      temp2 <- cbind(t(temp_df), subtest[ind1, !colnames(subtest)%in%c("region","chr", "start", "end", "length", "group")])
      
      return(temp2)
    }
  })
  finaldf <- rbindlist(finaldf)
  
  final_region[[i]] <- finaldf
  print(i)
}

subclone_anno <- readRDS("~/CASCADEpaper/paper/Fig3/Infercnv/subclone_anno.Rds")
multi_g <- readRDS("~/genomic/exon_only_reference.Rds")
multi_g <- multi_g[order(multi_g$chrom, multi_g$start, decreasing = F), ]




# comparing regions using wilcox test --------------------------------------
sig_regions <- list()
for( i in names(subclone_anno)){
  if( any(grepl(i , list.files("~/CASCADEpaper/paper/Fig3/Infercnv/CONIC_PURPLE_sigment")))){
    inf_obj <- readRDS(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/subclone/", i, "/preliminary.infercnv_obj" ))
    mat <- inf_obj@expr.data
    subclones <- subclone_anno[[i]]
    mat <- mat[,colnames(mat)%in%names(subclones)]
    candRegions= final_region[[i]]
    plot_list <- list()
    plot_fail_list <- list()
    sigregion <- apply(candRegions, 1, function(regions){
      ks_list <- list()
      chrom <- as.numeric(regions["chr"])
      start <- as.numeric(regions["start"])
      end <- as.numeric(regions["end"])
      CN <- as.numeric(regions[4: length(regions)])
      names(CN) <- substr(names(regions)[4:length(regions)], 1, 1)
      sbgenes <- unique(multi_g$gene.name[multi_g$chrom==chrom&multi_g$start>=start&multi_g$end<=end])
      if(length(which(rownames(mat)%in%sbgenes))>1){
        sbmat <- mat[rownames(mat)%in%sbgenes, ]
        sbdf <- reshape2::melt(sbmat)
        colnames(sbdf) <- c("gene", "ID", "score")
        sbdf$subclone <- subclones[sbdf$ID]
        # agg_df <- aggregate(score~gene+subclone, sbdf, mean)
        agg_df <- aggregate(score~ID+subclone, data = sbdf, mean)
        m_agg <- aggregate(score~subclone, agg_df, mean)
        m_agg$CN <-CN[as.character(m_agg$subclone)]
        r = paste(chrom, start, end, sep = ";")
        # g1<- ggplot(agg_df, aes(x = score, colour = as.factor(subclone))) + 
        #   geom_density(alpha = .3) + 
        #   ggtitle(r)+ 
        #   labs(colour = "subclone")+
        #   geom_vline(data = m_agg, aes(xintercept = score, colour = as.factor(subclone)), linetype = "dashed", size = 1) + 
        #   theme(legend.position = "none", plot.title = element_text(size = 9), axis.title = element_blank())
        # legend <- as_ggplot(get_legend(ggplot(agg_df, aes(x = score, colour = as.factor(subclone))) + 
        #                                  geom_density(alpha = .3) + 
        #                                  labs(colour = "subclone")+
        #                                  geom_vline(data = m_agg, aes(xintercept = score, colour = as.factor(subclone)), linetype = "dashed", size = 1)))
        ks_df <- split(agg_df, agg_df$subclone)
        comb <- combn(names(ks_df),2)
        p_value <- c()
        for (n in seq(ncol(comb))){ ## test whether the expression score of subclone cells agree with the CN prediction -- direction test
          clone1 <- comb[1, n]
          clone2 <- comb[2, n]
          if (CN[clone1] - CN[clone2] >0 ){
            k1<- wilcox.test(ks_df[[clone1]]$score, ks_df[[clone2]]$score, alternative = "greater")
            ks_list[[r]][[paste(comb[1,n], "vs", comb[2,n])]] <- k1
            p_value <- c(p_value, k1$p.value)}
          else if(CN[clone1] - CN[clone2] <0 ){
            k1<- wilcox.test(ks_df[[clone1]]$score, ks_df[[clone2]]$score, alternative = "less")
            ks_list[[r]][[paste(comb[1,n], "vs", comb[2,n])]] <- k1
            p_value <- c(p_value, k1$p.value)}
          else{
            k1<- wilcox.test(ks_df[[clone1]]$score, ks_df[[clone2]]$score, alternative = "two.sided")
            ks_list[[r]][[paste(comb[1,n], "vs", comb[2,n])]] <- k1
            p_value <- c(p_value, 1-k1$p.value) ## if the CNs are the same then the comparison result should be reversed
          }
        }
        
        adj_p_value <- p.adjust(p_value, method = "bonferroni")
        adj_p_value_great_less  <- adj_p_value[sapply(ks_list[[r]], function(x) x$alternative) %in% c("greater", "less")]
        if (all(adj_p_value_great_less<0.05)){
          return(cbind(data.frame(as.list(regions)), data.frame(comp = paste(names(ks_list[[r]])[sapply(ks_list[[r]], function(x) x$alternative) %in% c("greater", "less")], collapse = ";"))))
        
        }
      }
    })
    
    
    if (!is.null(sigregion)){
      sigregion <- sigregion[sapply(sigregion, length)>1]
      if (length(sigregion) > 1){
      sigregion <- rbindlist(sigregion)
      colnames(sigregion) <- c(colnames(candRegions), "comp")
        sig_regions[[i]]<- sigregion
        # saveRDS(plot_list, paste0("~/subclone/1_region_test/", i, "_plot_list.Rds"))
        # saveRDS(plot_fail_list, paste0("~/subclone/1_region_test/", i, "_plot_fail_list.Rds"))
        print(i)
        comb <- apply(combn(unique(sort(subclones)), 2), 2, function(x) paste(x, collapse = " vs "))
        combr <- unique(unlist(strsplit(sigregion$comp,split = ";"))) %in% comb
        names(combr) <- unique(unlist(strsplit(sigregion$comp,split = ";")))
        print(combr)
        
      }else{ print(paste(i, "no sig region"))}
      }else{print(paste(i, "no sig region"))}
    }
}


connect_neighbour <- function(ind,x){
  temp <- x[ind, ]
  temp2 <- x[ind+1, ]
  test1 <- paste(temp[, 4:grep("comp", colnames(temp))],collapse = ";")
  test2 <- paste(temp2[, 4:grep("comp", colnames(temp2))],collapse = ";")
  if (test1 == test2){
    new <- temp
    new$end <- temp2$end
    new$length <- new$end - new$start
    return(new)
  }else{
    return(rbind(temp, temp2))
  }
}

sig_region2 <- lapply(sig_regions, function(x){
  x$length <- x$end- x$start
  next_start <- x$start[2:nrow(x)]
  x$next_start <- c(next_start,0)
  ind <- which(x$end == x$next_start)
  if(length(ind) > 0){
  newdf <- mapply(connect_neighbour, as.list(ind), list(x), SIMPLIFY = F)
  newdf <- rbindlist(newdf) %>% distinct()
  x <- x[-unique(c(ind, ind+1)), ]
  x <- rbind(x, newdf)
  x <- x[order(x$chr, x$start),]}
  
  ind <- which(x$end == x$next_start-1)
  if(length(ind) > 0){
  newdf <- mapply(connect_neighbour, as.list(ind), list(x), SIMPLIFY = F)
  newdf <- rbindlist(newdf) %>% distinct()
  x <- x[-unique(c(ind, ind+1)), ]
  x <- rbind(x, newdf)
  x <- x[order(x$chr, x$start),]}
  return(x)
})
sigregiondf <- rbindlist(sig_region2, idcol = "sample")



saveRDS(sig_regions, "sig_region.Rds")
View(sig_regions)


