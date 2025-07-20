clean_module <- readRDS("~/Fig4-5_archetype/clean_module.Rds")
final_signature <- clean_module
names(final_signature) <- paste0("Module", 1:6)

paths <- list.files("~/Fig4-5_archetype/4_srt_meta", full.names = T)

annos <- lapply(paths, function(i){
  meta <- readRDS(i)
  df <- lapply(seq(from = 45, to = 70, by =1), function(p){
    group<- apply(meta, 2, function(x){
      bins <- cut(x, breaks = 100)
      group <- bins %in% levels(bins)[p:100]
    })
    ind_high <- apply(group, 1, function(x) sum(x))
    return(as.data.frame(table(ind_high)) %>% mutate(cutoff = p))
  })
  df <- rbindlist(df)
  df$cutoff <- factor(df$cutoff, levels = df %>% filter(ind_high ==1) %>% arrange(Freq) %>% pull(cutoff))
  df %>% filter(ind_high == "1") -> subdf
  
  cutoff <- as.numeric(as.character(subdf$cutoff[which.max(subdf$Freq)]))
  
  
  group<-  apply(meta, 2, function(x){
    bins <- cut(x, breaks = 100)
    group <- bins %in% levels(bins)[cutoff:100]
  })
  ind_high <- apply(group, 1, function(x) if(sum(x) > 1) data.frame( group = paste(names(clean_module)[x], collapse = "&"), freq = sum(x)) else if(sum(x)==1) data.frame(group = names(clean_module)[x], freq = 1) else data.frame(group = "background", freq = 0))
  ind_high <- rbindlist(ind_high)
  meta$oldgroup <- ind_high$group
  ind <- grepl("AR&", ind_high$group)&ind_high$freq == 2
  ind_high$group[ind] <- gsub("AR&", "",ind_high$group[ind])
  ind_high$freq[ind] <- ind_high$freq[ind]-1
  
  meta$group <- ind_high
  meta$cutoff <- cutoff
  return(meta)
})


names(annos) <- gsub(".Rds", "", list.files("~/Fig4-5_archetype/4_srt_meta"))

saveRDS(annos,"~/Fig4-5_archetype/per_sample_anno.Rds" )

