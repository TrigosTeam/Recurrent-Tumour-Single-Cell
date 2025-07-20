library(data.table)
paths <- list.files("~/Fig2_signature/sc_integration/DEG/inter_lesion", pattern = "CA", full.names = T)
inter_lesionDEG <- lapply(paths, readRDS)
names(inter_lesionDEG) <- gsub(".Rds", "",list.files("~/Fig2_signature/sc_integration/DEG/inter_lesion", pattern = "CA"))
inter_lesionDEG <- rbindlist(inter_lesionDEG, idcol = "patient")
inter_lesionDEG <- inter_lesionDEG %>% filter(p_val_adj <0.05) %>% mutate(mergetag  = paste(cluster, gene))

###### inter-organ DEG --------------------
inter_organ_DEG <- readRDS("~/Fig2_signature/sc_integration/DEG/inter_organ_DEG.Rds")
inter_organ_DEG <- inter_organ_DEG %>% filter(p_val_adj <0.05)
write.table(inter_organ_DEG,"~/Supp_tables/inter_organ_DEG.csv", quote = F, row.names = F, col.names = T, sep = "\t")

#####inter-lesion within patients CNA influence check
inter_interlesion_comp <- readRDS("~/Fig2_signature/sc_integration/DEG/inter_interlesion_comp.Rds")
df <- lapply(inter_interlesion_comp, lapply, function(x) x %>% rownames_to_column("gene"))
df <- lapply(df, rbindlist)
df <- rbindlist(df)
df$mergetag <- paste(df$patient, df$gene)
head(inter_lesionDEG)

df2 <- merge(inter_lesionDEG, df, by = "mergetag", all = T)
df2 <- df2[, c(2:9,12:14,16:19)]
colnames(df2)[1] <- "patient"
write.table(df2,"~/Supp_tables/inter_lesion_CNA_comp.csv", quote = F, row.names = F, col.names = T, sep = "\t")

##### DEG DCNR between subclone
DEGs <- readRDS("~/Fig3_genomics_normal_cells/ATAClone/DEGs.Rds")
DE_region <- readRDS("~/Fig3_genomics_normal_cells/ATAClone/DE_region.Rds")

df <- lapply(DEGs[1:22], function(x) if(length(x)>1) rbindlist(x, fill = T) else x[[1]])
df[23] <- DEGs[23]
df <- rbindlist(df, idcol = "sample", fill = T)
df <- df[df$p_val_adj<0.05, ]
View(df)
write.table(df,"~/Supp_tables/inter_subclone_DEG.csv", quote = F, row.names = F, col.names = T, sep = "\t")


df <- lapply(DE_region,function(x){
  temp <- x[["region"]]
  colnames(temp)[1] <- "region"
  genes <- sapply(x[["genes"]], paste, collapse = ",")
  temp$gene <- genes[temp$region]
  temp <- temp[temp$p_val_adj<0.05, ]
  return(temp)
})
df <- rbindlist(df, idcol = "sample")
head(df)
write.table(df,"~/Supp_tables/inter_subclone_DCNR.csv", quote = F, row.names = F, col.names = T, sep = "\t")

##### FOLH1 variants-------
FOLH1_SNP <- readRDS("~/Fig6_PSMA/FOLH1_SNP.Rds")
write.table(FOLH1_SNP,"~/Supp_tables/FOLH1_SNP.csv", quote = F, row.names = F, col.names = T, sep = "\t")

