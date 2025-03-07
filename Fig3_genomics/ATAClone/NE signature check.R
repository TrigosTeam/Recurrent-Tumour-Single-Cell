load("/trigos_team/CASCADE/Analysis/signatures_analysis/objects/all_gene_sets.rda")
NE_gens <- genelist$NE_signature
NE_gens <-all_gene_sets$NE
NE_gens <- NE_gens[1:9]
NE_gens <- unique(unlist(NE_gens))

DE_region <- readRDS("~/CASCADEpaper/paper/Fig3_genomics/ATAClone/DE_region.Rds")

lapply(DE_region, function(x) sum(NE_gens %in% unlist(x$genes)))
NE_gens %in% unlist(DE_region$CA0058_hilar_50$genes)
DE_region$CA0058_hilar_50$region[DE_region$CA0058_hilar_50$region$gene %in% names(which(sapply(DE_region$CA0058_hilar_50$genes, function(x) any(x %in% NE_gens)))), ] %>%
  filter(p_val_adj < 0.05) %>% View()


      