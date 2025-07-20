library(data.table)
sv <- readRDS("~/Fig3_genomics_normal_cells/structure_variants/sv_list.Rds")
fgene <- read.delim("~/Fig3_genomics_normal_cells/structure_variants/fusiongenes.txt", header = F)
fgene <- as.vector(fgene$V1)
gene_sv <- lapply(sv, function(x) {
  df <- data.frame(start_gene = x$start_gene, end_gene=x$end_gene, type = x$sv)
  df <- df%>% filter(start_gene != "") %>% filter(end_gene != "") %>% filter(start_gene!=end_gene)
  
  ind1 <- which(sapply(strsplit(df$start_gene, split = ";"), function(y) any(y%in%fgene)))
  ind2 <- which(sapply(strsplit(df$end_gene, split = ";"), function(y) any(y%in%fgene)))
  ind <- union(ind1,ind2)
  df <- df[ind, ]
  return(df)
})

genef_table <- rbindlist(gene_sv, idcol = "sample")
saveRDS(genef_table, "~/Fig3_genomics_normal_cells/structure_variants/genef_table.Rds")

del_g <- c(";NRG1-IT1", "SH2D2A;")
for ( i in del_g){
  genef_table$end_gene <- gsub(i, '', genef_table$end_gene)
  genef_table$start_gene <- gsub(i, '', genef_table$start_gene)
}
unique(paste(genef_table$start_gene, genef_table$end_gene, sep = "-"))

genef_table$pair <- paste(genef_table$start_gene, genef_table$end_gene, sep = "-")

