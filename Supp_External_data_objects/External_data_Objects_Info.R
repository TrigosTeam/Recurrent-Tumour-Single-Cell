# Savings External Dataset Objects

#Molecular profiling stratifies diverse phenotypes of treatment-refractory metastatic castration-resistant prostate cancer
#GSE126078
local_path <- "~/ExternalDatasets/Public_prostate_single_cell/Molecular_profiling_mCRPC/"
temp <- read.table(gzfile(paste("raw", "GSM3590965_C42B_siCont_1_processed_data.txt.gz", sep="/")),sep="\t")

saveRDS(mol_prof, file = "~/External_data_objects/mol_prof.Rds")
saveRDS(mol_prof_specimen, file = "~/External_data_objects/mol_prof_specimen.Rds")

# To add more genes to the mol_prof_specimen dataset (Patient Sample only) use the following code
gene_add = function(gene, mol_prof_xeno1, mol_prof1) {
  #Check if gene already added
  if (gene %in% colnames(mol_prof_xeno1)) {
    print("Gene already in data")
    return(mol_prof_xeno1)
  }
  library(dplyr)
  mol_prof_temp = filter(mol_prof, Gene_Names %in% c(gene))
  mol_prof_temp <- unlist(mol_prof_temp)
  mol_prof_xeno1$gene <- mol_prof_temp[match(rownames(mol_prof_xeno1), names(mol_prof_temp))]
  mol_prof_xeno1$gene <- as.numeric(as.character(mol_prof_xeno1$gene))
  names(mol_prof_xeno1)[length(names(mol_prof_xeno1))] = gene
  
  return(mol_prof_xeno1)
}
mol_prof_specimen = gene_add("EZH2", mol_prof_specimen, mol_prof)
for (variable in pre_NE_genes) {
  mol_prof_specimen = gene_add(variable, mol_prof_specimen, mol_prof)
}


# Roxanne's data set, as suggested by Luc Furic
library(ggplot2)
library(reshape2)
exp <- read.delim("~/ExternalDatasets/Mouse_models/GSE92721_RNAseq_normalized.txt")
sample_map <- read.delim("~/ExternalDatasets/Mouse_models/Sample_map_GSE92721.txt", sep="\t")
# N: control
# NP: PTEN-null
# NPp53: PTEN and p53-null (best resembles full mCRPC)
saveRDS(exp, file = "~/External_data_objects/exp.Rds")
saveRDS(sample_map, file = "~/External_data_objects/sample_map.Rds")
# To plot a gene use the following
rox_get_gene_plots = function(gene_name, exp_matrix, sample_mapping) {
  exp_gene <- exp_matrix[rownames(exp_matrix) == gene_name,]
  exp_gene <- melt(exp_gene)
  colnames(exp_gene) <- c("Sample", paste0(gene_name, "_level"))
  
  
  exp_gene$Condition <- sample_mapping$description[match(exp_gene$Sample, sample_mapping$ID1)]
  exp_gene$Mouse <- sample_mapping$Mouse[match(exp_gene$Sample, sample_mapping$ID1)]
  
  exp_gene$Condition <- factor(exp_gene$Condition,
                               levels=c("prostate tissue",
                                        "intact prostate tumor tissue treated with vehicle",
                                        "castrated prostate tumor tissue treated with vehicle",
                                        "castrated prostate tumor tissue treated with vehicle (outlier)",
                                        "castrated prostate tumor tissue treated with abiraterone",
                                        "castrated prostate tumor tissue treated with abiraterone (exp)"))
  return(exp_gene)
}
# This code uses the above function to get a gene for plotting
exp_HEPACAM2 = rox_get_gene_plots("Hepacam2", exp, sample_map)
ggplot(exp_HEPACAM2, aes(x=Mouse, y=Hepacam2_level))+
  geom_boxplot(aes(fill=Condition)) +
  theme(legend.position="none")

# Prostate Cancer Genome Atlas
# Prostate cancer atlas (bulk RNAseq merged from different sources)
# Gene names here are in 
atlas <- read.delim("~/ExternalDatasets/Public_prostate_single_cell/Dynamic_prostate_cancer_transcriptome/pca_atlas_dataset/pcatlas_dataset_vst_normalized.txt")
atlas_meta <- read.delim("~/ExternalDatasets/Public_prostate_single_cell/Dynamic_prostate_cancer_transcriptome/pca_atlas_dataset/pcatlas_dataset_annotations.txt")

saveRDS(atlas_meta, file = "~/External_data_objects/atlas_meta.Rds")
saveRDS(atlas_meta_NORMAL,file = "~/External_data_objects/atlas_meta_normal.Rds")
saveRDS(atlas_meta_PRIMARY,file = "~/External_data_objects/atlas_meta_primary.Rds")
saveRDS(atlas_meta_CRPC,file = "~/External_data_objects/atlas_meta_crpc.Rds")
saveRDS(atlas_meta_NEPC,file = "~/External_data_objects/atlas_meta_nepc.Rds")

saveRDS(prostate_atlas_genes_added, file = "~/External_data_objects/prostate_atlas_genes_added.Rds")
# To add more genes use the following
# Install the package if you have not installed by running this command: 
BiocManager::install("EnsDb.Hsapiens.v79")

library(EnsDb.Hsapiens.v86)

# 1. Convert from ensembl.gene to gene.symbol
ensembl.genes <- c("ENSG00000150676", "ENSG00000099308", "ENSG00000142676", "ENSG00000180776", "ENSG00000108848", "ENSG00000277370", "ENSG00000103811", "ENSG00000101473")

geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))

# 2. Convert from gene.symbol to ensembl.gene
geneSymbols <-  c(pre_NE_genes_filtered, cluster4_dura_marker_genes_filtered, cluster7_dura_from_4_genes_filtered, cluster7_dura_marker_genes_filtered, cluster89_from_7_genes_filtered, NE_Dura_Top_Genes_filtered)

geneIDs2 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= geneSymbols, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))

geneIDs2 = geneIDs2[-c(5, 53, 69, 71), ]
library(foreach)

atlas_gene_add = function(i, j) {
  if (j %in% colnames(atlas_meta)) {
    print("Gene already in data")
    return(atlas_meta)
  }
  PCDH9 = unlist(atlas[rownames(atlas) == i,])
  atlas_meta$PCD = PCDH9
  colnames(atlas_meta)[ncol(atlas_meta)] = j
  return(atlas_meta)
}

foreach(i = geneIDs2$GENEID, j = geneIDs2$SYMBOL) %do% {
  atlas_meta = atlas_gene_add(i, j)
}

# Primary Prostate Cancer Dataset - Erho
##Feng dataset #For investigating whether genes mark differences in metastasis or gleason grade
expr <- read.delim("~/ExternalDatasets/GSE46691/GSE46691_expression_aggregated.txt")
metadata <- read.csv("~/ExternalDatasets/GSE46691/GSE46691_series_matrix_for_R.csv")
saveRDS(metadata, file = "~/External_data_objects/metadata.Rds") #The dataframe for plotting
#To add more genes to the dataframe use
for (gene in NE_Dura_Top_Genes_filtered) { #this is a list of gene names
  g_exp = expr[rownames(expr) == gene,]
  names(g_exp) = gsub("\\.CEL", "", names(g_exp))
  metadata$gene = unlist(g_exp[match(metadata$Sample, names(g_exp))])
  colnames(metadata)[ncol(metadata)] = gene
}

#Living Tumour Laboratory PDX's
#Polycomb-mediated silencing in neuroendocrine prostate cancer
#GSE41193
##They have a website with more information:
#http://www.livingtumorlab.com/PDX_Prostate.html

##This information was derived from a PDF in the GSE entry
###Patient 927 -> 331 + 331R
##Patient 946 and 972 (NE metastasis of same patient) -> 352+370
##Patient 1005 -> 412
##Patient 1015 -> 418
#Living Tumour Lab Dataset
living_dataset <- read.delim("~/ExternalDatasets/Public_prostate_single_cell/Living_tumour_laboratory/Living_Tumour_Laboratory_samples.txt")
further_info <- read.delim("~/ExternalDatasets/Public_prostate_single_cell/Living_tumour_laboratory/Living_Tumour_Laboratory_data_sheets.txt")
living_GSE_samples <- read.delim("~/ExternalDatasets/Public_prostate_single_cell/Living_tumour_laboratory/Living_Tumour_Laboratory_patient_sample_key.txt")

living_GSE_samples$File <- NA
for(sample in living_GSE_samples$Sample_ID){
  temp <- grep(sample, list.files("~/ExternalDatasets/Public_prostate_single_cell/Living_tumour_laboratory/"))
  temp_name <- list.files("~/ExternalDatasets/Public_prostate_single_cell/Living_tumour_laboratory/")[temp]
  living_GSE_samples$File[living_GSE_samples$Sample_ID == sample] <- paste("~/ExternalDatasets/Public_prostate_single_cell/Living_tumour_laboratory/", temp_name, sep="") 
}

saveRDS(living_data, file = "~/External_data_objects/living_data.Rds")
saveRDS(living_tumour_dataset, file = "~/External_data_objects/living_tumour_dataset.Rds")
#Create a subset of the dataset to plot different genes
# If you want to add more genes to the dataset use this code
add_to_living = function(gene_name, dataset, living_d) {
  #Check if gene already added
  if (gene_name %in% colnames(dataset)) {
    print("Gene already in data")
    return(dataset)
  }
  gene_to_add = living_d[rownames(living_d) == gene_name,]
  if(is.matrix(gene_to_add)) {
    gene_to_add <- apply(gene_to_add, 2, median)
  }
  dataset$gene_name <- gene_to_add
  names(dataset)[length(names(dataset))] = gene_name
  
  return(dataset)
}

living_tumour_dataset = add_to_living("AR", living_tumour_dataset, living_data)


# Dynamic prostate cancer transcriptome analysis delineates the trajectory to disease progression 
# Paper has the Prostate Cancer Atlas in it as well
# PDX's Pre and Post Castration # This is single cell data from PDX's see paper for more info
path <- "~/ExternalDatasets/Public_prostate_single_cell/Dynamic_prostate_cancer_transcriptome/E-MTAB-9903/"
# The two objects for analysis are PDX_merge and PDX_mergeNE. PDX_merge is all samples and PDX_mergeNE are only the samples that expressed HDAC9 or ASCL1
saveRDS(PDX_merge, file = "~/External_data_objects/PDX_merge.Rds")
saveRDS(PDX_mergeNE, file = "~/External_data_objects/PDX_mergeNE.Rds")
