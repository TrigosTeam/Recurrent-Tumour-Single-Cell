data <- readRDS("~/CASCADEpaper/paper/Fig3_genomics/point_mutation/point_mutation_table_intersection.Rds")
data <-rbindlist(data, idcol = "sample")

data %>% filter(SYMBOL == "FOLH1") -> psma
psma <- psma[, 1:25]
psma$variant <- paste(psma$CHROM, psma$POS, psma$REF, psma$ALT)

psma <- split(psma, psma$sample)
psma <- lapply(psma, function(x) return(split(x, x$variant)))

psma <- lapply(psma, lapply, function(cons){
  if(length(unique(cons$IMPACT)) == 1){
  temp <- unlist(strsplit(cons$Consequence, split = "&"))
  temp <- names(which.max(table(temp)))
  return(data.frame(IMPACT = unique(cons$IMPACT[cons$Consequence==temp]), consequence = temp, variant = unique(cons$variant)))}
 else{
    ind <-  which(c("HIGH", "MODERATE", "MODIFIER") %in% cons$IMPACT)
    temp <- unlist(strsplit(cons$Consequence[cons$IMPACT == c("HIGH", "MODERATE", "MODIFIER")[min(ind)]], split = "&"))
    temp <- names(which.max(table(temp)))
    return(data.frame(IMPACT = unique(cons$IMPACT[cons$Consequence==temp]), consequence = temp, variant = unique(cons$variant)[cons$Consequence==temp]))
  }
  
})

psma <- rbindlist(lapply(psma, rbindlist), idcol = "sample")
View(psma)


sv <- readRDS("~/CASCADEpaper/paper/Fig3/structure_variants/sv_list.Rds")

psma <- lapply(sv, function(x){
  return(x %>% filter(grepl("FOLH1", fixed=T, start_gene)|grepl("FOLH1", fixed = T, end_gene)))
})
