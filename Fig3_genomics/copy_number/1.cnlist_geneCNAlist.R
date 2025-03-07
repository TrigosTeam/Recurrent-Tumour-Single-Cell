library(readr)
library(purrr)
library(dplyr)
library(GenomicRanges)
library(data.table)
# mylist<- list( patient = c(rep("CA0090", 5), rep ("CA0083", 4), rep("CA0058", 3), rep("CA0043", 3), rep("CA0027", 3), rep("CA0034", 2)),
#                sample = c("liver_43", "abdomen_13", "brain_2","paraaortic_lymph_node_39", "porta_hepatis_lymph_node_52",
#                           "hilar_lymph_node_47", "liver_49", "paraaortic_lymph_node_12", "perinephric_fat_19",
#                           "hilar_50", "liver_29", "liver_38",
#                           "liver_7", "portal_lymph_node_4", "liver_12", 
#                           "dura_base_skull_13", "paraaortic_lymph_node_1", "prostate_9", 
#                           "liver_left_11", "liver_right_8"))
# 
# mylist<- list( patient = c(rep("CA0090", 5), rep("CA0058", 3)), 
#                               sample = c("liver_43", "abdomen_13", "brain_2","paraaortic_lymph_node_39", "porta_hepatis_lymph_node_52",
#                                          "hilar_50", "liver_29", "liver_38"))
# 
# cnlist<- list()
# for ( i in seq(length(mylist[[1]]))){
#   patient <-  mylist$patient[i]
#   sample <-  mylist$sample[i]
#   path <- grep(".purple.cnv.somatic.tsv",list.files(paste0("/trigos_team/CASCADE/Sequencing_updated_reference/", patient, "/", sample, "/purple/"), full.names = T), value = T)
#   cnlist[[paste(patient, sample, sep = "_")]]  <- read_tsv(path)
# }


gf = rtracklayer::import("/trigos_team/CASCADE/Analysis/reference_annotation/multiome/genes.gtf")


paths <- system("realpath /trigos_team/CASCADE/Sequencing_updated_reference/*/*/purple/*cnv.somatic.tsv", intern = T)
paths_split <- sapply(strsplit(paths, split = "/"), function(x) paste(x[5], x[6], sep = "_"))
paths <- setNames(paths, paths_split)
cnlist <- lapply(paths, read_tsv)

# cnlist <- lapply(cnlist, function(x) return(x[x$segmentStartSupport == "NONE"&x$segmentEndSupport == "NONE",])) 
cnlist <- lapply(cnlist, function(x) {x$bin_width <- x$end - x$start 
                 return(x)})
# cnlist <- lapply(cnlist, function(x) return(x[x$chromosome != "chrY", ]))
# # cnlist <- lapply(cnlist, function(x) {x[x$chromosome == "chrX", "chromosome"] <- "chr23"
# x$chromosome <- gsub("chr", "", x$chromosome)
# x$chromosome <- as.numeric(x$chromosome)
# return(x)})


saveRDS(cnlist, "cnlist.Rds")

multi_g <- gf[gf$type =="gene"]



# 2024_06 -----------------------------------------------------------------

geneCNA_list <- lapply(cnlist, function(x){
  temp <- as(x, "GRanges")
  hits <- findOverlaps(multi_g, temp, select = "all", minoverlap = 1)
  df <- lapply(seq(length(hits)), function(ind){
    n = from(hits)[ind]
    m = to(hits)[ind]
    g_start = start(multi_g)[n]
    g_end = end(multi_g)[n]
    f_start = start(temp)[m]
    f_end = end(temp)[m]
    len <- ((width(temp)[m]+width(multi_g)[n]) - (abs(g_start-f_start)+abs(g_end - f_end)))/2
    p_gene <- len/width(multi_g)[n]
    return(data.frame(gene=multi_g$gene_name[n], chrom = as.character(multi_g@seqnames[n]), g_start = g_start , g_end = g_end, 
    f_start = f_start, f_end = f_end, SCN = temp$copyNumber[m], CN =round(temp$copyNumber[m]), 
    p_gene = p_gene))
  })
  print("done")
  return(rbindlist(df))
})

saveRDS(geneCNA_list, "geneCNAlist.Rds")

# previous version --------------------------------------------------------
geneCNA_list <- list()
for (i in names(cnlist)){
  df <- cnlist[[i]]
  temp <- data.frame()
  for (n in seq(nrow(multi_g))){
    g_start <- multi_g[n, "start"]
    g_end <- multi_g[n, "end"]
    chrom <- multi_g[n, "chrom"]
    df_c <- subset(df, df$chromosome == chrom)
    end <- which(apply(df_c, 1, function(z) between(g_end, as.numeric(z["start"]), as.numeric(z["end"]))))
    start <- which(apply(df_c, 1, function(z) between(g_start, as.numeric(z["start"]), as.numeric(z["end"]))))
    
    if(length(end)>0 &length(start)>0){
      ind <- seq(start, end)
      if( length(ind) == 1){
        p_gene <- 1
      }
      if (length(ind)==2){
        g_length <- g_end - g_start
        g1 <- (df_c$end[head(ind, 1)]-g_start)/g_length
        g2 <- (g_end - df_c$start[tail(ind,1)])/g_length
        p_gene <- c(g1, g2)
      }
      if (length(ind)>2){
        g_length <- g_end - g_start
        g1 <- (df_c$end[head(ind, 1)]-g_start)/g_length
        g2 <- (g_end - df_c$start[tail(ind,1)])/g_length
        p_gene <- c(g1,df_c$bin_width[ind[2:(length(ind)-1)]]/g_length, g2)
      }
    }
    
    if (length(end)==0){
      ind <- seq(start, which(df_c$end == max(df_c$end)))
      g_length <- g_end - g_start
      if( length(ind) == 1){
        p_gene <- (df_c$end[tail(which(df_c$chromosome==chrom), 1)]-g_start)/g_length}
      if (length(ind)==2){
        g1 <- (df_c$end[head(ind, 1)]-g_start)/g_length
        g2 <- (df_c$end[tail(ind,1)] - df_c$start[tail(ind,1)])/g_length
        p_gene <- c(g1, g2)  
      }
      if (length(ind)>2){
        g1 <- (df_c$end[head(ind, 1)]-g_start)/g_length
        g2 <- (df_c$end[tail(ind,1)] - df_c$start[tail(ind,1)])/g_length
        p_gene <- c(g1,df_c$bin_width[ind[2:(length(ind)-1)]]/g_length,g2) 
      }
    }
    
    if (length(start)==0){
      ind <- seq(head(which(df_c$chromosome==chrom), 1), end)  
      g_length <- g_end - g_start
      if( length(ind) == 1){
        p_gene <- (g_end - df_c$start[head(which(df_c$chromosome==chrom), 1)])/g_length}
      if(length(ind) == 2){
        g1 <- (df_c$end[head(ind, 1)]-df_c$start[head(which(df_c$chromosome==chrom), 1)])/g_length
        g2 <- (g_end - df_c$start[tail(ind,1)])/g_length
        p_gene <- c(g1, g2)  
      }
      if(length(ind)>2){
        g1 <- (df_c$end[head(ind, 1)]-df_c$start[head(which(df_c$chromosome==chrom), 1)])/g_length
        g2 <- (g_end - df_c$start[tail(ind,1)])/g_length
        p_gene <- c(g1, df_c$bin_width[ind[2:(length(ind)-1)]]/g_length, g2)  
      }
    }
    
    temp <- rbind(temp, data.frame(gene=multi_g$gene.name[n], chrom = chrom, g_start = g_start, g_end = g_end, 
                                   f_start = df_c$start[ind], f_end = df_c$end[ind], SCN = df_c$copyNumber[ind], CN =round(df_c$copyNumber[ind]), p_gene = p_gene)) 
  }
  geneCNA_list[[i]]<- temp
  print(i)
}


# uncesssary step
# class.assign <- function(x) {
#   x$class[x$p_gene>=0.5]<- ">0.5"
#   x$class[x$p_gene<0.5&x$p_gene >0]<- "0~0.5"
#   x$class[x$p_gene< 0 & x$p_gene > -0.5]<- "-0.5~0"
#   x$class[x$p_gene< -0.5 & x$p_gene > -1]<- "-1~-0.5"
#   x$class[x$p_gene<=-1]<- "<-1"
#   return(x)
# }
# 
# geneCNA_list <- lapply(geneCNA_list, class.assign)
# gene_table <- lapply(geneCNA_list, function(x) return(as.data.frame.matrix(table(x$gene, x$class))))
# gene_fitlered <- lapply(gene_table, function(x) return(x[x$`>0.5`<1|x$`0~0.5`>0|x$`-0.5~0`>0,]))
saveRDS(geneCNA_list, "~/genomic/copy_number/geneCNAlist_preliminary.Rds")


gene_fitlered <- lapply(geneCNA_list, function(x) return(unique(x$gene[x$p_gene<1])))

#### This the version usued for all the downstream analysis
# genes disturbed by breakpoints just follow CN of segment the most gene body fall into same as the segmeng belonging

for ( i in names(geneCNA_list)){
  fil <- gene_fitlered[[i]]
  if(length(fil)>0){
    print(i)
    df <- geneCNA_list[[i]]
    for( g in fil){
      sbdf <- df[df$gene == g, ]
      max_p <- max(sbdf$p_gene)
      fdf <- sbdf[sbdf$p_gene == max_p, ]
      df<- df[df$gene != g, ]
      df <- rbind(df, fdf)
      print(g)
    }
    df$region <- paste(df$chrom, df$f_start, df$f_end, sep = ";")
    df$fragment_width <- df$f_end-df$f_start
    geneCNA_list[[i]]<-df
  }
  print(i)
}

# "CA0036_307"                                
# [11] "CA0036_310"  
names(geneCNA_list)[c(10, 11)] <- c("CA0036_right_liver_307", "CA0036_left_liver_310")
# "CA0063_26"  bladder_nodule                             
# [25] "CA0063_29"   bladder_base                               "CA0063_62"   R_lobe_liver_2                              
# [27] "CA0063_67"    R_lobe_liver_4                        "CA0063_69A"   Vertebra_A                             
# [29] "CA0063_97" soft_tissue_rib
names(geneCNA_list)[c(24:29)] <- paste("CA0063", c( "bladder_26","bladder_base_29","R_lobe_liver_62","R_lobe_liver_67","Vertebra_A_69A","soft_tissue_rib_97"),sep = "_")
# "CA0079_19"       R_paraaortic_LN_thoracic                           "CA0079_2"R_dural_tumour_1
names(geneCNA_list)[c(33:34)] <- paste("CA0079", c("R_paraaortic_LN_thoracic_19","R_dural_tumour_2"))

saveRDS(geneCNA_list, "~/genomic/copy_number/geneCNAlist.Rds")


#CN is recalculated bt percentage of gene segments and PURPLE CN, not used later
for ( i in names(geneCNA_list)){
  fil <- gene_fitlered[[i]]
  if(length(fil)>0){
    print(i)
    df <- geneCNA_list[[i]]
    for( g in rownames(fil)){
      sbdf <- df[df$gene == g, ]
      sbdf$temp_CN <- sbdf$SCN*sbdf$p_gene
      
      fdf <- cbind(sbdf[1, 1:4], data.frame(f_start = min(sbdf$f_start), f_end = max(sbdf$f_end), SCN = sum(sbdf$temp_CN), CN = round(sum(sbdf$temp_CN)), p_gene = sum(sbdf$p_gene), class = "mix"))
      df<- df[df$gene != g, ]
      df <- rbind(df, fdf)
      print(g)
    }
    df$region <- paste(df$chrom, df$f_start, df$f_end, sep = ";")
    df$fragment_width <- df$f_end-df$f_start
  geneCNA_list[[i]]<-df
  }
  print(i)
}
saveRDS(geneCNA_list, "~/PURPLE/genelist_updated.Rds")



# CN unchanged but region is mereged when assign the corresponding segments 
for ( i in names(geneCNA_list)){
  fil <- gene_fitlered[[i]]
  if(length(fil)>0){
    print(i)
    df <- geneCNA_list[[i]]
    for( g in rownames(fil)){
      sbdf <- df[df$gene == g, ]
      if (gene_table[[i]][g,'>0.5']==1){
        fdf <- sbdf[sbdf$class == ">0.5", ]
        fdf$class <- ">0.5 not complete"
        df<- df[df$gene != g, ]
        df <- rbind(df, fdf)
        print(g)
      }else{
        agCN <- aggregate(p_gene~CN, sbdf, sum)
        if( TRUE %in% agCN$p_gene>0.5){
          CN <- agCN$CN[agCN$p_gene>0.5]
          fdf <- cbind(sbdf[1, 1:4], data.frame(f_start = min(sbdf$f_start), f_end = max(sbdf$f_end), SCN = mean(sbdf$SCN[sbdf$CN ==CN]), CN = CN, p_gene = agCN$p_gene[agCN$CN==CN], class = "fragments"))
          df<- df[df$gene != g, ]
          df <- rbind(df, fdf)
          print(g)
        }else{
          df<- df[df$gene != g, ] 
          print(g)
          print("delete the gene")
                  }
      }
    }
    df$region <- paste(df$chrom, df$f_start, df$f_end, sep = ";")
    df$fragment_width <- df$f_end-df$f_start
    geneCNA_list[[i]]<-df
  }
  print(i)
}
saveRDS(geneCNA_list, "~/PURPLE/genelist_CN_updated.Rds")


# merged regions for CONICS run, not used later 
new_region <- list()
for ( i in names(geneCNA_list)){
  cnlist <- geneCNA_list[[i]]
  regions <- data.frame(Chrom = cnlist$chrom, Start = cnlist$f_start, End = cnlist$f_end, class = cnlist$class)
  regions <- distinct(regions)
  regions <- aggregate(End ~ Chrom + Start, regions, max)
  regions <- regions[order(regions$Chrom, regions$Start), ]
  sbregions <- data.frame()
  for (n in unique(regions$Chrom)){
    temp <- regions[regions$Chrom==n, ]
    temp2 <- temp[temp$Start%in% (temp$End+1)|temp$Start==1, ]
    if(nrow(temp2)==0){
      temp2 <- temp[temp$Start==min(temp$Start),]
    }
    if (length(which(duplicated(temp2$End)))>0){
      temp2 <- temp2[temp2$Start%in% (temp2$End+1)|temp2$Start==1, ]
    }
    sbregions <- rbind(sbregions, temp2)
  }
  sbregions$Length <- sbregions$End-sbregions$Start
  rownames(sbregions)<- paste(sbregions$Chrom, sbregions$Start, sbregions$End, sep = ";")
  new_region[[i]]<- sbregions
}
saveRDS(new_region, "~/CONICS/merged_regions.Rds")

for (i in names(geneCNA_list)){
  geneCNA_list[[i]]$region <- paste(geneCNA_list[[i]]$chrom, geneCNA_list[[i]]$f_start,geneCNA_list[[i]]$f_end, sep = ";")
}
regions <- list()
for (i in names(geneCNA_list)){
  temp <- data.frame(Chrom = geneCNA_list[[i]]$chrom, Start =geneCNA_list[[i]]$f_start, End = geneCNA_list[[i]]$f_end, Length = geneCNA_list[[i]]$f_end- geneCNA_list[[i]]$f_start)
  temp <- distinct(temp)
  rownames(temp)<- paste(temp$Chrom, temp$Start, temp$End, sep = ";")
  regions[[i]] <- temp 
}

saveRDS(regions, "~/CONICS/unmerged_regions.Rds")

