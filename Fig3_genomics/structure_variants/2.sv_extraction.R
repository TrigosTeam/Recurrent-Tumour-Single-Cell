library(rlist)
library(dplyr)
library(GenomicRanges)

gf = rtracklayer::import("/trigos_team/CASCADE/Analysis/reference_annotation/multiome/genes.gtf")
#previously is "exon"
gf_df <-as.data.frame(gf[gf$type!= "gene"])
c  = nchar(as.character(gf_df$seqnames))
gf_df <- gf_df[which(c==4|c==5),]
gf_df$seqnames <- as.character(gf_df$seqnames)
multi_g <- gf_df %>% filter(seqnames!= "chrM" &seqnames!= "chrY")
multi_g$seqnames <- gsub("chrX", "chr23", multi_g$seqnames)
multi_g$seqnames<- as.numeric(gsub("chr", "", multi_g$seqnames))
multi_g <- as(multi_g, "GRanges")


paths <- system("realpath /trigos_team/CASCADE/Sequencing_updated_reference/CA*/*/purple/circos/*link.circos", intern = T)
paths_split <- data.table::transpose(strsplit(paths, "/"))
ids <- paths_split[[5]]
sites <- paths_split[[6]]

colordf <- data.frame(color = c("color=black", "color=blue", "color=green", "color=red", "color=vdyellow"), sv = c("inversion", "translocation", "tandem_duplication", "deletion", "insertion"))

sv_list <- lapply(paths, function(x){
  df <- read.delim(x, header = F)
  colnames(df) <- c("start_chrom", "start_pos", "end_pos", "alt_chrom", "alt_start", "alt_end", "color")
  df <- merge(df, colordf, by = "color")
  df <- df %>% filter(start_chrom!="hsY") %>% filter(alt_chrom!="hsY") 
  df$start_chrom <- as.numeric(gsub("X", "23", gsub("hs", "", df$start_chrom)))
  df$alt_chrom <- as.numeric(gsub("X", "23", gsub("hs", "", df$alt_chrom)))
  return(list(ref = as(df[, c(2, 3, 4, 1, 8)] %>% `colnames<-`(c("seqnames", "start", "end", "color", "sv")), "GRanges"), 
              alt = as(df[, c(5, 6, 7, 1, 8)] %>% `colnames<-`(c("seqnames", "start", "end", "color", "sv")), "GRanges"), 
              df = df))
})

names(sv_list) <- paste(ids, sites, sep = "_")


#2024 version
sv_list2 <- lapply(sv_list, function(x){
  ind1 <- split(to(findOverlaps(x$ref, multi_g, minoverlap = 0, select = "all")), from(findOverlaps(x$ref, multi_g, minoverlap = 0, select = "all")))
  ind2 <- split(to(findOverlaps(x$alt, multi_g, minoverlap = 0, select = "all")), from(findOverlaps(x$alt, multi_g, minoverlap = 0, select = "all")))
  ref_gene <- sapply(ind1, function(y) paste(unique(multi_g$gene_name[y]), collapse = ";"))
  ref_type <- sapply(ind1, function(y) paste(unique(multi_g$type[y]), collapse = ";"))
  alt_gene <- sapply(ind2, function(y) paste(unique(multi_g$gene_name[y]), collapse = ";"))
  alt_type <- sapply(ind2, function(y) paste(unique(multi_g$type[y]), collapse = ";"))
  df <- cbind(x$df, data.frame(start_gene = NA, start_type = NA, end_gene = NA, end_type = NA))

  df$start_gene[as.numeric(names(ref_gene))] <- ref_gene
  df$start_type[as.numeric(names(ref_gene))] <- ref_type
  df$end_gene[as.numeric(names(ind2))] <- alt_gene
  df$end_type[as.numeric(names(ind2))] <- alt_type
  print("done")
  return(df)
})
saveRDS(sv_list2, "~/CASCADEpaper/paper/Fig3_genomics/structure_variants/sv_list.Rds")

#2023 version 
sv_list2 <- list()
for( i in names(sv_list)){
  x <- sv_list[[i]]
  df <- data.frame()
  for (n in seq(nrow(x))){
    y <- x[n, ]
    temp <- multi_g[multi_g$chrom==as.numeric(y["start_chrom"]), ]
    temp2 <- multi_g[multi_g$chrom==as.numeric(y["alt_chrom"]), ]
    start_gene <- paste(unique(temp$gene[which(apply(temp, 1, function(z) between(as.numeric(y["start_pos"]), as.numeric(z["start"]), as.numeric(z["end"]))))]), collapse = ";")
    end_gene <- paste(unique(temp2$gene[which(apply(temp2, 1, function(z) between(as.numeric(y["alt_start"]), as.numeric(z["start"]), as.numeric(z["end"]))))]), collapse = ";")
    if(length(start_gene)>0){
      start_type <- c()
      for(m in start_gene){
        temp <- gf %>% filter(gene_name == m)
        start_type <-c(start_type, paste(unique(temp$type[which(apply(temp, 1, function(z) between(as.numeric(y["start_pos"]), as.numeric(z["start"]), as.numeric(z["end"]))))]), collapse = ";")) 
      }
      start_type <- paste(start_type, collapse = ";")
    }else { 
      start_type <- NA
      start_gene <- NA}
    
    if(length(end_gene)>0){
      end_type <- c()
      for (e in end_gene){
        temp2 <- gf %>% filter(gene_name ==e)
        end_type <- c(end_type, paste(unique(temp$type[which(apply(temp, 1, function(z) between(as.numeric(y["start_pos"]), as.numeric(z["start"]), as.numeric(as.numeric(z["end"])))))]), collapse = ";")) 
      }
      end_type <- paste(end_type, collapse = ";")  
    }else {
    end_type <- NA
    end_gene <- NA}
    # n <- list(start_gene = start_gene, start_type= start_type, end_gene= end_gene, end_type=end_type)
    # ml <- max(sapply(n, function(x) return(length(x))))
    # lapply(n, function(x) {
    #   length(x) <- ml
    #   return(x)
    #   })
    df <- rbind(df, cbind(y, start_gene = start_gene, start_type= start_type, end_gene= end_gene, end_type=end_type))
  }
  sv_list2[[i]] <- df
  print(i)
}
saveRDS(sv_list2, "~/genomic/structure_variants/sv_list.Rds")
 
  