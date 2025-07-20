library(ggvenn)
source("~/cols.R")
library(VennDiagram)
library(UpSetR)
library(ComplexUpset)
library(ggplot2)
library(patchwork)
 data <- readRDS("~/Fig3_genomics_normal_cells/point_mutation/point_mutation_table_intersection.Rds")
 variants <- lapply(data, function(x) {
   temp <- x %>% filter(IMPACT != "MODIFIER") %>% filter(as.numeric(MAX_AF )<0.1 | is.na(as.numeric(MAX_AF)))
   return(unique(paste(temp$CHROM, temp$POS, temp$REF, temp$ALT, sep = ";")))
 })
 
impacts <- lapply(data, function(x) {
   temp <- x %>% filter(IMPACT != "MODIFIER") %>% filter(as.numeric(MAX_AF )<0.1 | is.na(as.numeric(MAX_AF)))
   temp$variant <- paste(temp$CHROM, temp$POS, temp$REF, temp$ALT, sep = ";")
   impact <- sapply(split(temp$IMPACT, temp$variant), function(y)  head(intersect(c("HIGH","MODERATE","LOW","MODIFIER"), y), 1))
   return(impact)
 })
 sites <- names(variants)
 sites[grep("brain", sites)] <- "brain"
 sites[grep("dura", sites)] <- "dura"
 sites[grep("prostate", sites)] <- "prostate"
 sites[grep("liver", sites)] <- "liver"
 sites[grep("fat", sites)] <- "fat"
 sites[grep("lymph|hilar", sites)] <- "LN"
 sites[grep("rib|vertebra", sites)] <- "bone"
 sites[grep("lung", sites)] <- "lung"
 sites[grep("abdomen", sites)] <- "abdomen"
 sites[grep("bladder", sites)] <- "bladder"
 pa <- substr(names(data), 1, 6)
 names(variants) <- paste(pa, sites, sapply(strsplit(names(variants), split = "_"), function(x) tail(x, 1)))
names(impacts) <- NULL
 patent_snp <- split(variants, pa)
 patient_impact <- split(impacts, pa)
 
 plotlist <- list()
 pdf("upset_snv.pdf", width = 6, height = 4)
 for (p in names(patent_snp)){
   plotl <- patent_snp[[p]]
   site <- split(sites, pa)[[p]]
   impact <- unlist(patient_impact[[p]]) 
   test <- fromList(plotl)
   rownames(test) <- unique(unlist(plotl))
   test$impact <- factor(impact[rownames(test)], levels = c("HIGH","MODERATE","LOW","MODIFIER"))
   
   f <- ComplexUpset::upset(test, 
                            intersect =  names(plotl),
                            base_annotations=list(
                              'Intersection of SNVs'=intersection_size(
                                counts=T,
                                mapping=aes(fill=impact)
                              )
                              +scale_fill_manual(values=c(
                                'HIGH'='#E41A1C', 'LOW'='#377EB8',
                                'MODIFIER'='#4DAF4A', 'MODERATE'='#FF7F00'
                              )) + ggtitle(paste(p, "SNVs intersection"))),
                            set_sizes=(
                              upset_set_size(
                                geom=geom_bar(
                                  aes(fill=impact),
                                  width=0.8
                                ),
                                position='right'
                              )
                            )+ geom_text(aes(label=..count..), hjust=1, stat='count')
                            # you can also add annotations on top of bars:
                            + theme(axis.text.x=element_text(angle=90), 
                                    legend.position = "none", 
                                    plot.margin = margin(r = 20))
                            +scale_fill_manual(values=c(
                              'HIGH'='#E41A1C', 'LOW'='#377EB8',
                              'MODIFIER'='#4DAF4A', 'MODERATE'='#FF7F00'
                            )),
                            
                            # moves legends over the set sizes
                            guides='over', 
                            queries= lapply(seq_along(plotl), function(ind){
                              upset_query(
                                set=names(plotl)[ind],
                                color=site_cols[site][ind],
                                only_components=c('intersections_matrix')
                              )
                            })
   )
   print(f)
   plotlist[[p]] <- f
 }
dev.off()
 


