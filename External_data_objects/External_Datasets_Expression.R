" 
 Title: Script to Explore Bulk and single cell RNA-seq Data from external data sets.
 
 Purpose: To explore and visualise the differential expression of identified genes from custom analysis 
          in external RNA-seq data sets.
 Author: James Comben
 Project: Honours Research Project Investigating Neuroendocrine Transformation in Late stage Prostate Cancer
 Date script created: 27/03/2023
 Date script finished: -
 
 # Package Dependencies #
 - tidyverse
 - ggplot2
 - dplyr
 
 Data Location: mCRPC Molecular profiling - /trigos_team/CASCADE/ExternalDatasets/Public_prostate_single_cell/Molecular_profiling_mCRPC/
 
 "
library(tidyverse)
library(dplyr)
library(ggplot2)

#Molecular profiling stratifies diverse phenotypes of treatment-refractory metastatic castration-resistant prostate cancer
#GSE126078
local_path <- "/trigos_team/CASCADE/ExternalDatasets/Public_prostate_single_cell/Molecular_profiling_mCRPC/"

temp <- read.table(gzfile(paste(local_path, "GSM3590965_C42B_siCont_1_processed_data.txt.gz", sep="/")),sep="\t")
temp <- temp[-1,]
mol_prof <- temp$V2
mol_prof <- data.frame(Gene = mol_prof)

for(temp_file in list.files(local_path)[2:length(list.files(local_path))]){
  temp <- read.table(gzfile(paste(local_path, temp_file, sep="/")),sep="\t")
  colnames(temp) <- temp[1,]
  temp <- temp[-1,]
  temp$Sample <- temp_file
  mol_prof <- cbind(mol_prof, temp$FPKM)
  local_name <- temp_file
  local_name <- gsub("_processed_data.txt.gz", "", local_name)
  colnames(mol_prof)[length(colnames(mol_prof))] <- local_name
}

colnames(mol_prof) <- substr(colnames(mol_prof), 12, nchar(colnames(mol_prof)))
#Metadata

ARnegNEpos <- c("03-192C3_LN",
                "05-144E6_LIVER",
                "05-144V3_ADRENAL",
                "13-084H5_LIVER",
                "15-010H1_LIVER",
                "16-071H10_LIVER",
                "16-071I2_LUNG",#
                "16-080L2_LIVER",
                "17-017H8_LIVER",
                "17-017N1_LUNG")

which(ARnegNEpos %in% colnames(mol_prof))



ARnegNEneg <- c("11-028G3_ADRENAL",
                "11-028H1_LIVER",
                "11-028L3_LUNG",
                "12-021FF6_BONE",
                "12-021H6_LIVER",
                "12-021I1_LUNG",
                "13-099H7_LIVER",
                "13-117H7_LIVER")


ARlowNEneg <- c("03-163S5_LIVER",
                "13-099J2_LUNG",
                "14-039E3_LUNG",
                "14-039K4_LIVER",
                "16-052H1_LIVER",
                "16-052K1_LN",
                "17-043I1_LN",
                "17-043J1_PERIRECTAL_OMENTAL",
                "17-043K1_PELVIC_PERITONEAL")

ARposNEpos <- c("06-131TT1_BONE",
                "07-042K2_LN",
                "07-042L3_ADRENAL",
                "10-056J1_LN",
                "10-056K1_LIVER",
                "12-011I8_LN",
                "13-042M3_LN",
                "14-096J1_LEFT_POSTERIOR",
                "15-096L2_PERIAORTIC",
                "15-096M1_PERIAORTIC_DIAPHRAGM",
                "17-033J1_PANCREAS")

which(ARposNEpos %in% colnames(mol_prof))
sort(colnames(mol_prof)[grep("LIVER", colnames(mol_prof))])


mol_prof_HDAC9 <- mol_prof[mol_prof[,1] == "HDAC9",]
mol_prof_HDAC9 <- mol_prof_HDAC9[1,]
mol_prof_HDAC9 <- unlist(mol_prof_HDAC9)
mol_prof_HDAC9 <- mol_prof_HDAC9[-1]
mol_prof_HDAC9 <- data.frame(HDAC9=mol_prof_HDAC9)

mol_prof_HDAC9$Classification <- "AR+NE-"
mol_prof_HDAC9$Classification[rownames(mol_prof_HDAC9) %in% ARnegNEpos] <- "AR-NE+"
mol_prof_HDAC9$Classification[rownames(mol_prof_HDAC9) %in% ARnegNEneg] <- "AR-NE-"

mol_prof_HDAC9$Classification[rownames(mol_prof_HDAC9) %in% ARlowNEneg] <- "ARlowNE-"
mol_prof_HDAC9$Classification[rownames(mol_prof_HDAC9) %in% ARposNEpos] <- "AR+NE+"

mol_prof_HDAC9$HDAC9 <- as.numeric(as.character(mol_prof_HDAC9$HDAC9))
mol_prof_HDAC9_LuCap <- mol_prof_HDAC9[grep("LuCaP", rownames(mol_prof_HDAC9)),]

mol_prof_HDAC9_xeno <- mol_prof_HDAC9[-grep("LuCaP", rownames(mol_prof_HDAC9)),]

# Use mol_prof_specimen as it contains only the classified patient samples
mol_prof_specimen = mol_prof_HDAC9_xeno[-grep("C42B", rownames(mol_prof_HDAC9_xeno)),]
mol_prof_specimen = mol_prof_specimen[-grep("Pac", rownames(mol_prof_specimen)),]
mol_prof_specimen = mol_prof_specimen[-grep("PC3", rownames(mol_prof_specimen)),]

ggplot(mol_prof_HDAC9_xeno, aes(x=Classification, y=HDAC9))+
  geom_boxplot()


mol_prof_AR <- mol_prof[mol_prof[,1] == "AR",]
mol_prof_AR <- mol_prof_AR[1,]
mol_prof_AR <- unlist(mol_prof_AR)
mol_prof_AR <- mol_prof_AR[-1]
mol_prof_HDAC9_xeno$AR <- mol_prof_AR[match(rownames(mol_prof_HDAC9_xeno), names(mol_prof_AR))]
mol_prof_HDAC9_xeno$AR <- as.numeric(as.character(mol_prof_HDAC9_xeno$AR))

ggplot(mol_prof_HDAC9_xeno, aes(x=log2(AR), y=log2(HDAC9)))+
  geom_point()

cor.test(mol_prof_HDAC9_xeno$AR, mol_prof_HDAC9_xeno$HDAC9, method="sp")

mol_prof_SYP <- mol_prof[mol_prof[,1] == "SYP",]
mol_prof_SYP <- mol_prof_SYP[1,]
mol_prof_SYP <- unlist(mol_prof_SYP)
mol_prof_SYP <- mol_prof_SYP[-1]
mol_prof_HDAC9_xeno$SYP <- mol_prof_SYP[match(rownames(mol_prof_HDAC9_xeno), names(mol_prof_SYP))]
mol_prof_HDAC9_xeno$SYP <- as.numeric(as.character(mol_prof_HDAC9_xeno$SYP))

ggplot(mol_prof_HDAC9_xeno, aes(x=log2(SYP), y=log2(HDAC9)))+
  geom_point()

cor.test(mol_prof_HDAC9_xeno$SYP, mol_prof_HDAC9_xeno$HDAC9, method="sp")

#ASCL1 & HEPACAM2 - Added Manually no function
mol_prof_xeno = mol_prof_HDAC9_xeno

mol_prof_ASCL1 <- mol_prof[mol_prof[,1] == "ASCL1",]
mol_prof_ASCL1 <- mol_prof_ASCL1[1,]
mol_prof_ASCL1 <- unlist(mol_prof_ASCL1)
mol_prof_ASCL1 <- mol_prof_ASCL1[-1]
mol_prof_xeno$ASCL1 <- mol_prof_ASCL1[match(rownames(mol_prof_xeno), names(mol_prof_ASCL1))]
mol_prof_xeno$ASCL1 <- as.numeric(as.character(mol_prof_xeno$ASCL1))

ggplot(mol_prof_xeno, aes(x=Classification, y=ASCL1))+
  geom_boxplot()

mol_prof_HEPACAM2 <- mol_prof[mol_prof[,1] == "HEPACAM2",]
mol_prof_HEPACAM2 <- mol_prof_HEPACAM2[4,]
mol_prof_HEPACAM2 <- unlist(mol_prof_HEPACAM2)
mol_prof_HEPACAM2 <- mol_prof_HEPACAM2[-1]
mol_prof_xeno$HEPACAM2 <- mol_prof_HEPACAM2[match(rownames(mol_prof_xeno), names(mol_prof_HEPACAM2))]
mol_prof_xeno$HEPACAM2 <- as.numeric(as.character(mol_prof_xeno$HEPACAM2))

ggplot(mol_prof_xeno, aes(x=Classification, y=HEPACAM2))+
  geom_boxplot()

# NEED TO ADDRESS THIS ASAP!!!! Some of the samples are siRNA REST Knock Down
########################################################################################
#  THE SAMPLE SET IN XENO IS NOT JUST XENOGRAFTS - XENO IS LUCAP SYMBOLS, ONLY THE     #
#  PATIENT SAMPLES HAVE BEEN CLASSIFIED FOR THE GRAPH, THIS IS DESTORTING RESULTS!!!!  #
########################################################################################
# Update: Has been addressed instead use the mol_prof_specimen data object

#Some more cleaning up of the dataset
n_occur <- data.frame(table(mol_prof[,1]))
n_occur[n_occur$Freq > 1,]
vocabulary[vocabulary$id %in% n_occur$Var1[n_occur$Freq > 1],]
mol_prof <- mol_prof[-c(10521,10502,15695),]

colnames(mol_prof)[1] ="Gene_Names"
mol_prof_gene_names = mol_prof[,1]
row.names(mol_prof) = mol_prof_gene_names
# Ended up using dplyr since the other method doesnt work well for the a function
library(dplyr)

mol_prof_temp = filter(mol_prof, Gene_Names %in% c(gene))
mol_prof_temp <- as.matrix(mol_prof[mol_prof[,1] == gene,])

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

mol_prof_xeno = gene_add("NTN4", mol_prof_xeno, mol_prof)

pre_NE_genes = c("HEPACAM2", "NTN4", "BLCAP", "RUNX1", "NRXN3", "GALNTL6")
cluster4_dura_marker_genes = c("HEPACAM2", "NTN4", "RUNX1", "NRXN3", "BLCAP", "AC083837.1")

cluster7_dura_marker_genes = c("ST7", "COL25A1","MID1","BLCAP","NTN4","FER1L6","PTPRE","EPHA5",
                               "NRXN3", "FOXP2", "GALNTL6", "PBX1") # This is 7 from 8 and 9
cluster7_dura_from_4_genes = c("PCDH9","KCNJ6","AC093916.1","TNIK","ST7","PPFIA2","KIAA0825","PTPRN2","GAS2L3","ODF2L","NPAS3","PCED1B","CIT","HEPACAM2",
                               "THSD7A","SOX4","PODXL","MTSS1","AMIGO2","SRGAP1","CDH6","GPR137C","SLC18A1","LINC02466","EPHA5","PTPRE","PBX1",
                               "TGFBR3","SDK2","ERBB4","GRIK2","SCGN","SORBS1","SOX5","ST7-AS2","MID1")
cluster89_from_7_genes = c("HDAC9", "CNTNAP2","SYT1","NELL2","ST18","ASCL1","CACNA2D1","SMOC2","SDK1","RIMBP2","RALYL","MIAT","THSD7A","ADCY2","TRPM3","KCNMA1","PEX5L","DIRAS2","LINC00511","SLCO3A1","TRPM8",
                           "CAMK1D","DMD","CACNA1A","CHGA","AMACR","MGAT4C","ATP8A2","NRXN1","MIR924HG","RIMS2","RASSF6","TNS3","NKAIN2","SDK2","BCL2","LHFPL6","DPYSL3","AC105031.2","SPOCK1",
                           "JAZF1","NLGN1","STXBP5L","PAM","DPP6","BRINP2","NFATC2","ERC2","NFASC","NPNT", "ESRRG","JAKMIP2","SRRM4")
NE_Dura_Top_Genes = c("PCDH9","THSD7A","CACNA2D1","RALYL","CACNA1A","KCNJ6","ASCL1","PPFIA2","NELL2","CNTNAP2","PTPRN2","SDK2","RIMBP2","ST18","KCNMA1","HEPACAM2","AC099520.1","HDAC9","SLCO3A1","TNIK","SCGN","MIAT","TRPM3","RIMS2","KCNMB2","ADCY2")
mol_prof_specimen = gene_add("EZH2", mol_prof_specimen, mol_prof)
for (variable in pre_NE_genes) {
  mol_prof_specimen = gene_add(variable, mol_prof_specimen, mol_prof)
}
for (variable in cluster4_dura_marker_genes) {
  mol_prof_specimen = gene_add(variable, mol_prof_specimen, mol_prof)
}
for (variable in cluster7_dura_marker_genes) {
  mol_prof_specimen = gene_add(variable, mol_prof_specimen, mol_prof)
}
for (variable in cluster7_dura_from_4_genes) {
  mol_prof_specimen= gene_add(variable, mol_prof_specimen, mol_prof)
}
for (variable in cluster89_from_7_genes) {
  mol_prof_specimen = gene_add(variable, mol_prof_specimen, mol_prof)
}
for (variable in NE_Dura_Top_Genes) {
  mol_prof_specimen = gene_add(variable, mol_prof_specimen, mol_prof)
}

# Plot the Genes now that they have been added
for (variable in pre_NE_genes) {
  print(ggplot(mol_prof_specimen, aes(x=Classification, y=eval(as.symbol(variable)))) + ylab(variable) +
          geom_boxplot())
  Sys.sleep(2)
}

for (variable in cluster4_dura_marker_genes) {
  print(ggplot(mol_prof_specimen, aes(x=Classification, y=eval(as.symbol(variable)))) + ylab(variable) +
          geom_boxplot())
  Sys.sleep(2)
}

for (variable in cluster7_dura_marker_genes) {
  print(ggplot(mol_prof_specimen, aes(x=Classification, y=eval(as.symbol(variable)))) + ylab(variable) +
          geom_boxplot())
  Sys.sleep(2)
}

for (variable in cluster7_dura_from_4_genes) {
  print(ggplot(mol_prof_specimen, aes(x=Classification, y=eval(as.symbol(variable)))) + ylab(variable) +
          geom_boxplot())
  Sys.sleep(2)
}

for (variable in cluster89_from_7_genes) {
  print(ggplot(mol_prof_specimen, aes(x=Classification, y=eval(as.symbol(variable)))) + ylab(variable) +
          geom_boxplot())
  Sys.sleep(2)
}

for (variable in NE_Dura_Top_Genes) {
  print(ggplot(mol_prof_specimen, aes(x=Classification, y=eval(as.symbol(variable)))) + ylab(variable) +
          geom_boxplot())
  Sys.sleep(2)
}

ggplot(mol_prof_specimen, aes(x=Classification, y=SPOCK1)) + ylab("SPOCK1") +
  geom_boxplot()

ggplot(mol_prof_specimen, aes(x=Classification, y=SRRM4)) + ylab("SRRM4") +
  geom_boxplot()

ggplot(mol_prof_specimen, aes(x=Classification, y=PPFIA2)) + ylab("PPFIA2") +
  geom_boxplot()

ggplot(mol_prof_specimen, aes(x=Classification, y=PPFIA2)) + ylab("KCNMA1") +
  geom_boxplot()

ggplot(mol_prof_specimen, aes(x=Classification, y=HEPACAM2)) + ylab("HEPACAM2") +
  geom_boxplot() + geom_signif(test = wilcox.test, comparisons = list(c("AR+NE+", "AR-NE-"), c("AR+NE+", "ARlowNE-"),c("AR+NE+", "AR-NE+"), c("AR+NE+", "AR+NE-"), c("AR-NE+", "ARlowNE-" ), c("AR-NE+", "AR-NE-"), c("AR-NE+", "AR+NE-" )), step_increase = 0.1, map_signif_level=TRUE) +
ggplot(mol_prof_specimen, aes(x=Classification, y=NTN4)) + ylab("NTN4") +
  geom_boxplot() +
ggplot(mol_prof_specimen, aes(x=Classification, y=BLCAP)) + ylab("BLCAP") +
  geom_boxplot() + geom_signif(test = wilcox.test, comparisons = list(c("AR+NE+", "AR-NE+"), c("AR-NE+", "ARlowNE-" ), c("AR-NE+", "AR+NE-" )), step_increase = 0.1, map_signif_level=TRUE) +
ggplot(mol_prof_specimen, aes(x=Classification, y=RUNX1)) + ylab("RUNX1") +
  geom_boxplot() + geom_signif(test = wilcox.test, comparisons = list(c("AR+NE+", "AR-NE-"), c("AR+NE+", "ARlowNE-"), c("AR-NE+", "ARlowNE-" ), c("AR-NE+", "AR-NE-")), step_increase = 0.1, map_signif_level=TRUE) +
ggplot(mol_prof_specimen, aes(x=Classification, y=NRXN3)) + ylab("NRXN3") +
  geom_boxplot() + geom_signif(test = wilcox.test, comparisons = list(c("AR+NE+", "AR-NE+"), c("AR-NE+", "ARlowNE-" ), c("AR-NE+", "AR+NE-" )), step_increase = 0.1, map_signif_level=TRUE) +
ggplot(mol_prof_specimen, aes(x=Classification, y=GALNTL6)) + ylab("GALNTL6") +
  geom_boxplot() + geom_signif(test = wilcox.test, comparisons = list(c("AR+NE+", "AR-NE-"), c("AR+NE+", "ARlowNE-"),c("AR+NE+", "AR-NE+"), c("AR+NE+", "AR+NE-"), c("AR-NE+", "AR+NE-" )), step_increase = 0.1, map_signif_level=TRUE)

# Wilcoxin ranked sum test
# Use AR+ as the base
ggplot(mol_prof_specimen, aes(x=Classification, y=HEPACAM2)) + ylab("HEPACAM2") +
  geom_boxplot()
class_filt = subset(mol_prof_specimen, mol_prof_specimen$Classification %in% c("AR-NE+", "AR+NE-") ) # Filter for 2 classes
class_filt$Classification = factor(class_filt$Classification)
class_filt$Classification <- relevel(class_filt$Classification,"AR-NE+") # Set which category will be the base
wilcox.test(class_filt$GALNTL6 ~ class_filt$Classification, conf.int = T) # Two sided Test between the specified groups


wilcox.test(HEPACAM2 ~ Classification, data = class_filt, alternative = "greater", conf.int = T) # One-sided Test
# Multiple comparisons correction for each set of genes marking a trajectory state
pvalue_vector = c(1.71E-05, 9.08E-07, 1.25E-01, 3.78E-01, 0.7567, 4.22E-03, 4.41E-01, 7.06E-01,
                  2.01E-01, 4.08E-06, 9.82E-03, 5.29E-04)
p.adjust(pvalue_vector, method = "BH")

#GSVA Gene Signature Plots
install.packages("BiocManager")
BiocManager::install("GSVA")
library(GSVA)
pre_NE_genes_filtered = c("HEPACAM2", "NTN4", "BLCAP", "RUNX1", "NRXN3", "GALNTL6")
cluster4_dura_marker_genes_filtered = c("HEPACAM2", "NTN4", "RUNX1", "NRXN3", "BLCAP")

cluster7_dura_marker_genes_filtered = c("ST7", "COL25A1","MID1","BLCAP","NTN4","FER1L6","PTPRE","EPHA5",
                               "NRXN3", "FOXP2", "GALNTL6", "PBX1")
cluster7_dura_from_4_genes_filtered = c("PCDH9","KCNJ6","TNIK","ST7","PPFIA2","KIAA0825","PTPRN2","GAS2L3","ODF2L","NPAS3","PCED1B","CIT","HEPACAM2",
                               "THSD7A","SOX4","PODXL","MTSS1","AMIGO2","SRGAP1","CDH6","GPR137C","SLC18A1","EPHA5","PTPRE","PBX1",
                               "TGFBR3","SDK2","ERBB4","GRIK2","SCGN","SORBS1","SOX5","ST7-AS2","MID1")
cluster89_from_7_genes_filtered = c("HDAC9", "CNTNAP2","SYT1","NELL2","ST18","ASCL1","CACNA2D1","SMOC2","SDK1","RIMBP2","RALYL","MIAT","THSD7A","ADCY2","TRPM3","KCNMA1","PEX5L","DIRAS2","LINC00511","SLCO3A1","TRPM8",
                           "CAMK1D","DMD","CACNA1A","CHGA","AMACR","MGAT4C","ATP8A2","NRXN1","MIR924HG","RIMS2","RASSF6","TNS3","NKAIN2","SDK2","BCL2","DPYSL3","SPOCK1",
                           "JAZF1","NLGN1","STXBP5L","PAM","DPP6","BRINP2","NFATC2","ERC2","NFASC","NPNT", "ESRRG","JAKMIP2","SRRM4")
NE_Dura_Top_Genes_filtered = c("PCDH9","THSD7A","CACNA2D1","RALYL","CACNA1A","KCNJ6","ASCL1","PPFIA2","NELL2","CNTNAP2","PTPRN2","SDK2","RIMBP2","ST18","KCNMA1","HEPACAM2","HDAC9","SLCO3A1","TNIK","SCGN","MIAT","TRPM3","RIMS2","KCNMB2","ADCY2")

# Get new data frame with correct format rows are genes and columns are samples
gene_sig_mol_prof_specimen = mol_prof_specimen[,-c(2, 10, 21, 40, 82, 84, 99)]
gene_sig_mol_prof_specimen = t(gene_sig_mol_prof_specimen)

gene_sig_mol_prof_specimen = as.matrix(gene_sig_mol_prof_specimen)
gene_sigs = list(pre_NE_genes_filtered, cluster7_dura_from_4_genes_filtered, cluster89_from_7_genes_filtered, NE_Dura_Top_Genes_filtered)

gene_sig_mol_prof_specimen[1:5, 1:5]
gsva.es <- gsva(gene_sig_mol_prof_specimen, gene_sigs, verbose=FALSE, method = "ssgsea") # Look into what NA deprecated means

gsva.es_df = as.data.frame(t(gsva.es))
colnames(gsva.es_df) = c("pre_NE_sig", "cluster7_from_4_sig", "cluster89_from_7_sig", "Top_NE_Dura_sig")

gsva.es_df$Classification = "AR+NE-"
gsva.es_df$Classification[rownames(gsva.es_df) %in% ARnegNEpos] = "AR-NE+"
gsva.es_df$Classification[rownames(gsva.es_df) %in% ARnegNEneg] = "AR-NE-"
gsva.es_df$Classification[rownames(gsva.es_df) %in% ARlowNEneg] = "ARlowNE-"
gsva.es_df$Classification[rownames(gsva.es_df) %in% ARposNEpos] = "AR+NE+"

ggplot(gsva.es_df, aes(x=Classification, y=pre_NE_sig)) + ylab("Pre-NE Signature") +
        geom_boxplot() + geom_signif(test = wilcox.test, comparisons = list(c("AR+NE+", "AR-NE+"), c("AR-NE+", "ARlowNE-" ), c("AR-NE+", "AR-NE-"), c("AR-NE+", "AR+NE-" )), step_increase = 0.1, map_signif_level=TRUE) +

ggplot(gsva.es_df, aes(x=Classification, y=cluster7_from_4_sig)) + ylab("Partial-NE Signature") +
  geom_boxplot() + geom_signif(test = wilcox.test, comparisons = list( c("AR+NE+", "ARlowNE-"),c("AR+NE+", "AR-NE+"), c("AR-NE+", "ARlowNE-" ), c("AR-NE+", "AR-NE-"), c("AR-NE+", "AR+NE-" )), step_increase = 0.1, map_signif_level=TRUE) +

ggplot(gsva.es_df, aes(x=Classification, y=cluster89_from_7_sig)) + ylab("Full-NE Signature") +
  geom_boxplot() + geom_signif(test = wilcox.test, comparisons = list(c("AR+NE+", "AR-NE-"), c("AR+NE+", "ARlowNE-"),c("AR+NE+", "AR-NE+"), c("AR+NE+", "AR+NE-"), c("AR-NE+", "ARlowNE-" ), c("AR-NE+", "AR-NE-"), c("AR-NE+", "AR+NE-" )), step_increase = 0.1, map_signif_level=TRUE)

ggplot(gsva.es_df, aes(x=Classification, y=Top_NE_Dura_sig)) + ylab("Top_NE_Dura_sig") +
  geom_boxplot()

# It seems as though in terms of signatures, a pre-NE signature does not appear to be highly expressed in NE+ samples
# relative to other samples which suggests that these genes may contribute but to trajectory switching but not be required to 
# maintain NEPC lineage 

#PCDH9 Does it split the NEPC Samples?
NE_only_mol_prof_specimen = filter(mol_prof_specimen, Classification == "AR-NE+")
NE_mol_prof_specimen = filter(mol_prof_specimen, Classification == c("AR-NE+", "AR+NE+"))
NotAR_mol_prof_specimen = filter(mol_prof_specimen, Classification %in% c("AR-NE+", "AR+NE+", "ARlowNE-"))

ggplot(NE_only_mol_prof_specimen, aes(x=ASCL1, y=PCDH9)) + ylab("PCDH9") +
  geom_point() + geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE)

ggplot(NE_only_mol_prof_specimen, aes(x=ASCL1, y=PCDH9)) + ylab("PCDH9") +
  geom_point() + geom_smooth(method = "lm", formula = y~log(x), se = FALSE)

cor.test(NE_only_mol_prof_specimen$ASCL1, NE_only_mol_prof_specimen$PCDH9, method = "spearman")

# Correlation between ASCL1 and EZH2 - Investigating myself after it was shown in the paper
ggplot(NE_only_mol_prof_specimen, aes(x=ASCL1, y=EZH2)) + ylab("EZH2") +
  geom_point() + geom_smooth(method = "lm", formula = y~x, se = FALSE)

ggplot(NE_mol_prof_specimen, aes(x=ASCL1, y=EZH2)) + ylab("EZH2") +
  geom_point() + geom_smooth(method = "lm", formula = y~x, se = FALSE)

ggplot(NotAR_mol_prof_specimen, aes(x=ASCL1, y=EZH2)) + ylab("EZH2") +
  geom_point() + geom_smooth(method = "lm", formula = y~x, se = FALSE)

cor(NE_only_mol_prof_specimen$ASCL1, NE_only_mol_prof_specimen$EZH2)
cor(mol_prof_specimen$ASCL1, mol_prof_specimen$EZH2)
cor(NE_mol_prof_specimen$ASCL1, NE_mol_prof_specimen$EZH2)
cor(NotAR_mol_prof_specimen$ASCL1, NotAR_mol_prof_specimen$EZH2, method = "sp")

# Z score transformation - Testing out code, data is not normally distributed
library(tibble)
column <- NotAR_mol_prof_specimen$Classification
NotAR_mol_prof_specimen = add_column(NotAR_mol_prof_specimen, column, .before = 1)
NotAR_mol_prof_specimen = NotAR_mol_prof_specimen[-c(3)]
colnames(NotAR_mol_prof_specimen)[1] = "Classification"
NotAR_mol_prof_specimen = NotAR_mol_prof_specimen[-c(1)]
znorm <- function (data) 
  {
    d = dim(data)
    c = class(data)
    cnames = colnames(data)
    rnames = rownames(data)
    classes = data[, d[2]]
    zdata = scale(data)
    attributes(zdata) = NULL
    zdata = matrix(zdata, dim(data)[1], dim(data)[2])
    zdata = cbind(zdata, classes)
    zdata = as.data.frame(zdata)
    colnames(zdata) = cnames
    rownames(zdata) = rnames
    return(zdata)
  }

NotAR_mol_prof_specimen = sapply(NotAR_mol_prof_specimen[, 1:100],as.numeric)
NotAR_mol_prof_specimen = znorm(NotAR_mol_prof_specimen)
rownames(NotAR_mol_prof_specimen) = rownames(NotAR_mol_prof_specimen1)

ggplot(NotAR_mol_prof_specimen, aes(x=ASCL1, y=EZH2)) + ylab("EZH2") +
  geom_point() + geom_smooth(method = "lm", formula = y~x, se = FALSE)

# Roxanne's data set, as suggested by Luc Furic
library(ggplot2)
library(reshape2)
exp <- read.delim("/trigos_team/CASCADE/ExternalDatasets/Mouse_models/GSE92721_RNAseq_normalized.txt")
rownames(exp) <- exp$GeneName
exp$GeneName <- NULL
boxplot(exp)

sample_map <- read.delim("/trigos_team/CASCADE/ExternalDatasets/Mouse_models/Sample_map_GSE92721.txt", sep="\t")
sample_map$Mouse <- sapply(strsplit(sample_map$ID, "_"), function(x){ return(x[2]) })
sample_map$ID1 <- sapply(strsplit(sample_map$ID, "_"), function(x){ return(x[1]) })


pca <- prcomp(exp)
pca <- as.data.frame(pca$rotation[,c(1:2)])
pca$Sample <- rownames(pca)
pca$Condition <- sample_map$description[match(pca$Sample, sample_map$ID1)]
pca$Mouse <- sample_map$Mouse[match(pca$Sample, sample_map$ID1)]

ggplot(pca, aes(x=PC1, y=PC2))+
  geom_point(aes(shape=Mouse, colour=Condition), size=3)

exp_Hdac9 <- exp[rownames(exp) == "Hdac9",]
exp_Hdac9 <- melt(exp_Hdac9)
colnames(exp_Hdac9) <- c("Sample", "Hdac9_level")


exp_Hdac9$Condition <- sample_map$description[match(exp_Hdac9$Sample, sample_map$ID1)]
exp_Hdac9$Mouse <- sample_map$Mouse[match(exp_Hdac9$Sample, sample_map$ID1)]

exp_Hdac9$Condition <- factor(exp_Hdac9$Condition,
                              levels=c("prostate tissue",
                                       "intact prostate tumor tissue treated with vehicle",
                                       "castrated prostate tumor tissue treated with vehicle",
                                       "castrated prostate tumor tissue treated with vehicle (outlier)",
                                       "castrated prostate tumor tissue treated with abiraterone",
                                       "castrated prostate tumor tissue treated with abiraterone (exp)"))
ggplot(exp_Hdac9, aes(x=Mouse, y=Hdac9_level))+
  geom_boxplot(aes(fill=Condition))

ggplot(exp_Hdac9, aes(x=Mouse, y=Hdac9_level))+
  geom_boxplot()

##Correlation with SOX11 #associated with NE
exp_Sox11 <- exp[rownames(exp) == "Sox11",]
exp_Sox11 <- melt(exp_Sox11)
colnames(exp_Sox11) <- c("Sample", "Sox11_level")

exp_Hdac9$Sox11 <- exp_Sox11$Sox11_level

ggplot(exp_Hdac9, aes(x=Hdac9_level, y = Sox11))+
  geom_point(aes(shape=Mouse, colour=Condition), size=3)

cor.test(exp_Hdac9$Hdac9_level, exp_Hdac9$Sox11, method="sp")


exp_Ar <- exp[rownames(exp) == "Ar",]
exp_Ar <- melt(exp_Ar)
colnames(exp_Ar) <- c("Sample", "Ar_level")

exp_Hdac9$Ar <- exp_Ar$Ar_level

ggplot(exp_Hdac9, aes(x=Hdac9_level, y = Ar))+
  geom_point(aes(shape=Mouse, colour=Condition), size=3)

cor.test(exp_Hdac9$Hdac9_level, exp_Hdac9$Ar, method="sp")

# Running Genes in dataset extracted from Analysis of CA0027
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

# Plot after the function returns the data frame for specified gene
#HDAC9

# N: control
# NP: PTEN-null
# NPp53: PTEN and p53-null (best resembles full mCRPC)
# Don't know what the outlier samples are/mean

geneSymbols = unique(geneSymbols)


gene_list = geneSymbols
exp_matrix = exp 
rownames(exp_matrix) = toupper(rownames(exp_matrix))
sample_mapping = sample_map

exp_gene <- exp_matrix[gene_list,]
exp_gene = na.omit(exp_gene)
exp_gene$gene = rownames(exp_gene)
exp_gene <- melt(exp_gene)
colnames(exp_gene) <- c("Gene", "Sample", "Expression")


exp_gene$Condition <- sample_mapping$description[match(exp_gene$Sample, sample_mapping$ID1)]
exp_gene$Mouse <- sample_mapping$Mouse[match(exp_gene$Sample, sample_mapping$ID1)]
exp_gene$Condition <- factor(exp_gene$Condition,
                             levels=c("prostate tissue",
                                      "intact prostate tumor tissue treated with vehicle",
                                      "castrated prostate tumor tissue treated with vehicle",
                                      "castrated prostate tumor tissue treated with vehicle (outlier)",
                                      "castrated prostate tumor tissue treated with abiraterone",
                                      "castrated prostate tumor tissue treated with abiraterone (exp)"))
# Gene Signature Generation
rox_all_genes = exp_gene
library(reshape2)
rox_gene_sigs = dcast(formula = Sample~Gene, data = rox_all_genes, value.var = "Expression")

rownames(rox_gene_sigs) = rox_gene_sigs$Sample
rox_gene_sigs = rox_gene_sigs[,-c(1)]
# Need to get the sample ID's to be the row names and make sure its unique first
rox_gene_sigs = t(rox_gene_sigs)

rox_gene_sigs = as.matrix(rox_gene_sigs)
gene_sigs = list(pre_NE_genes_filtered, cluster7_dura_from_4_genes_filtered, cluster89_from_7_genes_filtered, NE_Dura_Top_Genes_filtered)

rox_gene_sigs[1:5, 1:5]
gsva.es_rox <- gsva(rox_gene_sigs, gene_sigs, verbose=FALSE, method = "ssgsea") # Look into what NA deprecated means

gsva.es_rox_df = as.data.frame(t(gsva.es_rox))
colnames(gsva.es_rox_df) = c("pre_NE_sig", "cluster7_from_4_sig", "cluster89_from_7_sig", "Top_NE_Dura_sig")

gsva.es_rox_df$Mouse = "AR+NE-"
rox_N = dcast(formula = Sample + Mouse ~Gene, data = rox_all_genes, value.var = "Expression")
rox_N_mouse = rox_N[rox_N$Mouse == "N",]
rox_NP_mouse = rox_N[rox_N$Mouse == "NP",]
rox_NPp53_mouse = rox_N[rox_N$Mouse == "NPp53",]

gsva.es_rox_df$Mouse[rownames(gsva.es_rox_df) %in% rox_N_mouse$Sample] = "N"
gsva.es_rox_df$Mouse[rownames(gsva.es_rox_df) %in% rox_NP_mouse$Sample] = "NP"
gsva.es_rox_df$Mouse[rownames(gsva.es_rox_df) %in% rox_NPp53_mouse$Sample] = "NPp53"

#Add in tissue types
rox_N_tissue = dcast(formula = Sample + Mouse + Condition ~Gene, data = rox_all_genes, value.var = "Expression")

rox_C1 = rox_N_tissue[rox_N_tissue$Condition == "prostate tissue",]
rox_C2 = rox_N_tissue[rox_N_tissue$Condition == "intact prostate tumor tissue treated with vehicle",]
rox_C3 = rox_N_tissue[rox_N_tissue$Condition == "castrated prostate tumor tissue treated with vehicle",]
rox_C4 = rox_N_tissue[rox_N_tissue$Condition == "castrated prostate tumor tissue treated with vehicle (outlier)",]
rox_C5 = rox_N_tissue[rox_N_tissue$Condition == "castrated prostate tumor tissue treated with abiraterone",]
rox_C6 = rox_N_tissue[rox_N_tissue$Condition == "castrated prostate tumor tissue treated with abiraterone (exp)",]

gsva.es_rox_df$Condition = "prostate tissue"
gsva.es_rox_df$Condition[rownames(gsva.es_rox_df) %in% rox_C1$Sample] = "prostate tissue"
gsva.es_rox_df$Condition[rownames(gsva.es_rox_df) %in% rox_C2$Sample] = "intact prostate tumor tissue treated with vehicle"
gsva.es_rox_df$Condition[rownames(gsva.es_rox_df) %in% rox_C3$Sample] = "castrated prostate tumor tissue treated with vehicle"
gsva.es_rox_df$Condition[rownames(gsva.es_rox_df) %in% rox_C4$Sample] = "castrated prostate tumor tissue treated with vehicle (outlier)"
gsva.es_rox_df$Condition[rownames(gsva.es_rox_df) %in% rox_C5$Sample] = "castrated prostate tumor tissue treated with abiraterone"
gsva.es_rox_df$Condition[rownames(gsva.es_rox_df) %in% rox_C6$Sample] = "castrated prostate tumor tissue treated with abiraterone (exp)"

gsva.es_rox_df$Condition = factor(gsva.es_rox_df$Condition,
       levels=c("prostate tissue",
                "intact prostate tumor tissue treated with vehicle",
                "castrated prostate tumor tissue treated with vehicle",
                "castrated prostate tumor tissue treated with vehicle (outlier)",
                "castrated prostate tumor tissue treated with abiraterone",
                "castrated prostate tumor tissue treated with abiraterone (exp)"))

ggplot(gsva.es_rox_df, aes(x=Mouse, y=pre_NE_sig)) + ylab("pre_NE_sig") +
  geom_boxplot(aes(fill=Condition))

ggplot(gsva.es_rox_df, aes(x=Mouse, y=cluster7_from_4_sig)) + ylab("cluster7_from_4_sig") +
  geom_boxplot(aes(fill=Condition))

ggplot(gsva.es_rox_df, aes(x=Mouse, y=cluster89_from_7_sig)) + ylab("cluster89_from_7_sig") +
  geom_boxplot(aes(fill=Condition))

ggplot(gsva.es_rox_df, aes(x=Mouse, y=Top_NE_Dura_sig)) + ylab("Top_NE_Dura_sig") +
  geom_boxplot(aes(fill=Condition))





exp_HDAC9 = rox_get_gene_plots("Hdac9", exp, sample_map)
ggplot(exp_HDAC9, aes(x=Mouse, y=Hdac9_level))+
  geom_boxplot(aes(fill=Condition))

ggplot(exp_HDAC9, aes(x=Mouse, y=Hdac9_level))+
  geom_boxplot()

#HEPACAM2
exp_HEPACAM2 = rox_get_gene_plots("Hepacam2", exp, sample_map)
exp_NTN4 = rox_get_gene_plots("Ntn4", exp, sample_map)
exp_BLCAP = rox_get_gene_plots("Blcap", exp, sample_map)
exp_RUNX1 = rox_get_gene_plots("Runx1", exp, sample_map)
exp_NRXN3 = rox_get_gene_plots("Nrxn3", exp, sample_map)
exp_GALNTL6 = rox_get_gene_plots("Galntl6", exp, sample_map)
ggplot(exp_HEPACAM2, aes(x=Mouse, y=Hepacam2_level))+
  geom_boxplot(aes(fill=Condition)) +
  theme(legend.position="none") +
ggplot(exp_NTN4, aes(x=Mouse, y=Ntn4_level))+
  geom_boxplot(aes(fill=Condition)) +
  theme(legend.position="none") +
ggplot(exp_BLCAP, aes(x=Mouse, y=Blcap_level))+
  geom_boxplot(aes(fill=Condition)) +
  theme(legend.position="none") +
ggplot(exp_RUNX1, aes(x=Mouse, y=Runx1_level))+
  geom_boxplot(aes(fill=Condition)) + theme(legend.position="none") +
ggplot(exp_NRXN3, aes(x=Mouse, y=Nrxn3_level))+
  geom_boxplot(aes(fill=Condition)) + theme(legend.position="none") +
ggplot(exp_GALNTL6, aes(x=Mouse, y=Galntl6_level))+
  geom_boxplot(aes(fill=Condition)) + theme(legend.position="none")

ggplot(exp_HEPACAM2, aes(x=Mouse, y=Hepacam2_level))+
  geom_boxplot()

#SPOCK1
exp_SPOCK1 = rox_get_gene_plots("Spock1", exp, sample_map)
ggplot(exp_SPOCK1, aes(x=Mouse, y=Spock1_level))+
  geom_boxplot(aes(fill=Condition))

ggplot(exp_SPOCK1, aes(x=Mouse, y=Spock1_level))+
  geom_boxplot()


# Prostate Cancer Genome Atlas
# Prostate cancer atlas (bulk RNAseq merged from different sources)

atlas <- read.delim("/trigos_team/CASCADE/ExternalDatasets/Public_prostate_single_cell/Dynamic_prostate_cancer_transcriptome/pca_atlas_dataset/pcatlas_dataset_vst_normalized.txt")

rownames(atlas) <- atlas$X
atlas$X <- NULL

atlas_meta <- read.delim("/trigos_team/CASCADE/ExternalDatasets/Public_prostate_single_cell/Dynamic_prostate_cancer_transcriptome/pca_atlas_dataset/pcatlas_dataset_annotations.txt")

HDAC9 <- unlist(atlas[rownames(atlas) == "ENSG00000048052",])
AR <- unlist(atlas[rownames(atlas) == "ENSG00000169083",])
SYP <- unlist(atlas[rownames(atlas) == "ENSG00000102003",])
ASCL1 <- unlist(atlas[rownames(atlas) == "ENSG00000139352",])
INSM1 = unlist(atlas[rownames(atlas) == "ENSG00000173404",])

atlas_meta$HDAC9 <- HDAC9
atlas_meta$AR <- AR
atlas_meta$SYP <- SYP
atlas_meta$ASCL1 <- ASCL1
atlas_meta$INSM1 = INSM1

atlas_meta$Sample.Type <- factor(atlas_meta$Sample.Type,
                                 levels=c("NORMAL", "PRIMARY", "CRPC", "NEPC"))

ggplot(atlas_meta, aes(x=Sample.Type, y=HDAC9))+
  geom_boxplot()

ggplot(atlas_meta, aes(x=Sample.Type, y=AR))+
  geom_boxplot()

ggplot(atlas_meta, aes(x=Sample.Type, y=SYP))+
  geom_boxplot()

ggplot(atlas_meta, aes(x=Sample.Type, y=ASCL1))+
  geom_boxplot()


cor.test(atlas_meta$HDAC9, atlas_meta$AR)
cor.test(atlas_meta$HDAC9, atlas_meta$SYP)
cor.test(atlas_meta$HDAC9, atlas_meta$ASCL1)


atlas_meta_NORMAL <- atlas_meta[atlas_meta$Sample.Type == "NORMAL",]
atlas_meta_PRIMARY <- atlas_meta[atlas_meta$Sample.Type == "PRIMARY",]
atlas_meta_CRPC <- atlas_meta[atlas_meta$Sample.Type == "CRPC",]
atlas_meta_NEPC <- atlas_meta[atlas_meta$Sample.Type == "NEPC",]

cor.test(atlas_meta_NORMAL$HDAC9, atlas_meta_NORMAL$AR)
cor.test(atlas_meta_NORMAL$HDAC9, atlas_meta_NORMAL$SYP)
cor.test(atlas_meta_NORMAL$HDAC9, atlas_meta_NORMAL$ASCL1)

cor.test(atlas_meta_PRIMARY$HDAC9, atlas_meta_PRIMARY$AR)
cor.test(atlas_meta_PRIMARY$HDAC9, atlas_meta_PRIMARY$SYP)
cor.test(atlas_meta_PRIMARY$HDAC9, atlas_meta_PRIMARY$ASCL1)

cor.test(atlas_meta_CRPC$HDAC9, atlas_meta_CRPC$AR)
cor.test(atlas_meta_CRPC$HDAC9, atlas_meta_CRPC$SYP)
cor.test(atlas_meta_CRPC$HDAC9, atlas_meta_CRPC$ASCL1)

cor.test(atlas_meta_NEPC$HDAC9, atlas_meta_NEPC$AR)
cor.test(atlas_meta_NEPC$HDAC9, atlas_meta_NEPC$SYP)
cor.test(atlas_meta_NEPC$HDAC9, atlas_meta_NEPC$ASCL1)

# Prostate Cancer Atlas Plots
HDAC9 <- unlist(atlas[rownames(atlas) == "ENSG00000048052",])
AR <- unlist(atlas[rownames(atlas) == "ENSG00000169083",])
ASCL1 <- unlist(atlas[rownames(atlas) == "ENSG00000139352",])

atlas_meta$HDAC9 <- HDAC9
atlas_meta$AR <- AR
atlas_meta$ASCL1 <- ASCL1

ggplot(atlas_meta, aes(x=Sample.Type, y=HDAC9))+
  geom_boxplot()

HEPACAM2 = unlist(atlas[rownames(atlas) == "ENSG00000188175",])
atlas_meta$HEPACAM2 = HEPACAM2
ggplot(atlas_meta, aes(x=Sample.Type, y=HEPACAM2))+
  geom_boxplot()

SPOCK1 = unlist(atlas[rownames(atlas) == "ENSG00000152377",])
atlas_meta$SPOCK1 = SPOCK1
ggplot(atlas_meta, aes(x=Sample.Type, y=SPOCK1))+
  geom_boxplot()

KCNMA1 = unlist(atlas[rownames(atlas) == "ENSG00000156113",])
atlas_meta$KCNMA1 = KCNMA1
ggplot(atlas_meta, aes(x=Sample.Type, y=KCNMA1))+
  geom_boxplot()

SRRM4 = unlist(atlas[rownames(atlas) == "ENSG00000139767",])
atlas_meta$SRRM4 = SRRM4
ggplot(atlas_meta, aes(x=Sample.Type, y=SRRM4))+
  geom_boxplot()

PPFIA2 = unlist(atlas[rownames(atlas) == "ENSG00000139220",])
atlas_meta$PPFIA2 = PPFIA2
ggplot(atlas_meta, aes(x=Sample.Type, y=PPFIA2))+
  geom_boxplot()

PCDH9 = unlist(atlas[rownames(atlas) == "ENSG00000184226",])
atlas_meta$PCDH9 = PCDH9
ggplot(atlas_meta, aes(x=Sample.Type, y=PCDH9))+
  geom_boxplot()

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

prostate_atlas_genes_added <- atlas_meta[,colnames(atlas_meta) %in% c("Sample.ID", "Sample.Type", geneSymbols)]
saveRDS(prostate_atlas_genes_added, file = "/home/jcomben/NE_Transformation/prostate_atlas_genes_added.Rds")

ggplot(atlas_meta, aes(x=Sample.Type, y=PCDH9))+
  geom_boxplot()
atlas_meta$Sample.Type = factor(atlas_meta$Sample.Type, levels = c("NORMAL", "PRIMARY", "CRPC", "NEPC") )

library(ggsignif)
# annotations = c("*") to specify your own
ggplot(atlas_meta, aes(x=Sample.Type, y=HEPACAM2))+ xlab("Classification") +
  geom_boxplot() + geom_signif(test = wilcox.test, comparisons = list(c("NEPC", "CRPC"), c("NEPC", "PRIMARY"),c("NEPC", "NORMAL") ), step_increase = 0.1, map_signif_level=TRUE) +
  ggplot(atlas_meta, aes(x=Sample.Type, y=NTN4))+ xlab("Classification") +
  geom_boxplot() + geom_signif(test = wilcox.test, comparisons = list(c("NEPC", "PRIMARY"),c("NEPC", "NORMAL") ), step_increase = 0.1, map_signif_level=TRUE) +
  ggplot(atlas_meta, aes(x=Sample.Type, y=BLCAP))+ xlab("Classification") +
  geom_boxplot() + geom_signif(test = wilcox.test, comparisons = list(c("NEPC", "PRIMARY")), step_increase = 0.1, map_signif_level=TRUE) +
  ggplot(atlas_meta, aes(x=Sample.Type, y=RUNX1))+ xlab("Classification") +
  geom_boxplot() + geom_signif(test = wilcox.test, comparisons = list(c("NEPC", "PRIMARY"),c("NEPC", "NORMAL") ), step_increase = 0.1, map_signif_level=TRUE) +
  ggplot(atlas_meta, aes(x=Sample.Type, y=NRXN3))+ xlab("Classification") +
  geom_boxplot() + geom_signif(test = wilcox.test, comparisons = list(c("NEPC", "CRPC"), c("NEPC", "PRIMARY")), step_increase = 0.1, map_signif_level=TRUE) +
  ggplot(atlas_meta, aes(x=Sample.Type, y=GALNTL6))+ xlab("Classification") +
  geom_boxplot() + geom_signif(test = wilcox.test, comparisons = list(c("NEPC", "PRIMARY")), step_increase = 0.1, map_signif_level=TRUE)





plot(atlas_meta_NEPC$ASCL1, atlas_meta_NEPC$PCDH9)

ggplot(atlas_meta_NEPC, aes(x=ASCL1, y=PCDH9)) + ylab("PCDH9") +
  geom_point() + geom_smooth(method = "lm", formula = y~log(x), se = FALSE)

ggplot(atlas_meta_NEPC, aes(x=ASCL1, y=PCDH9)) + ylab("PCDH9") +
  geom_point() + geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE)

ggplot(atlas_meta_NEPC, aes(x=ASCL1, y=PCDH9)) + ylab("PCDH9") +
  geom_point()

# Gene Signatures in Prostate Cancer Atlas
# Get new data frame with correct format rows are genes and columns are samples
gene_sig_atlas = prostate_atlas_genes_added[,-c(2)]
gene_sig_atlas = gene_sig_atlas[!duplicated(gene_sig_atlas$Sample.ID), ]
rownames(gene_sig_atlas) = gene_sig_atlas$Sample.ID
gene_sig_atlas = gene_sig_atlas[,-c(1)]
# Need to get the sample ID's to be the row names and make sure its unique first
gene_sig_atlas = t(gene_sig_atlas)

#Genes should be rows
gene_sig_atlas = as.matrix(gene_sig_atlas)
gene_sigs = list(pre_NE_genes_filtered, cluster7_dura_from_4_genes_filtered, cluster89_from_7_genes_filtered, NE_Dura_Top_Genes_filtered)

gene_sig_atlas[1:5, 1:5]
gsva.es_pa <- gsva(gene_sig_atlas, gene_sigs, verbose=FALSE, method = "ssgsea") # Look into what NA deprecated means

gsva.es_pa_df = as.data.frame(t(gsva.es_pa))
colnames(gsva.es_pa_df) = c("pre_NE_sig", "cluster7_from_4_sig", "cluster89_from_7_sig", "Top_NE_Dura_sig")

gsva.es_pa_df$Classification = "AR+NE-"
gsva.es_pa_df$Classification[rownames(gsva.es_pa_df) %in% atlas_meta_NORMAL$Sample.ID] = "NORMAL"
gsva.es_pa_df$Classification[rownames(gsva.es_pa_df) %in% atlas_meta_PRIMARY$Sample.ID] = "PRIMARY"
gsva.es_pa_df$Classification[rownames(gsva.es_pa_df) %in% atlas_meta_CRPC$Sample.ID] = "CRPC"
gsva.es_pa_df$Classification[rownames(gsva.es_pa_df) %in% atlas_meta_NEPC$Sample.ID] = "NEPC"

gsva.es_pa_df$Classification = factor(gsva.es_pa_df$Classification, levels = c("NORMAL", "PRIMARY", "CRPC", "NEPC") )

ggplot(gsva.es_pa_df, aes(x=Classification, y=pre_NE_sig)) + ylab("Pre-NE Signature") +
  geom_boxplot() + geom_signif(test = wilcox.test, comparisons = list(c("NEPC", "CRPC"), c("NEPC", "PRIMARY"),c("NEPC", "NORMAL") ), step_increase = 0.1, map_signif_level=TRUE) +

ggplot(gsva.es_pa_df, aes(x=Classification, y=cluster7_from_4_sig)) + ylab("Partial-NE Signature") +
  geom_boxplot() + geom_signif(test = wilcox.test, comparisons = list(c("NEPC", "CRPC"), c("NEPC", "PRIMARY"),c("NEPC", "NORMAL") ), step_increase = 0.1, map_signif_level=TRUE) +

ggplot(gsva.es_pa_df, aes(x=Classification, y=cluster89_from_7_sig)) + ylab("Full-NE Signature") +
  geom_boxplot() + geom_signif(test = wilcox.test, comparisons = list(c("NEPC", "CRPC"), c("NEPC", "PRIMARY"),c("NEPC", "NORMAL") ), step_increase = 0.1, map_signif_level=TRUE)

ggplot(gsva.es_pa_df, aes(x=Classification, y=Top_NE_Dura_sig)) + ylab("Top_NE_Dura_sig") +
  geom_boxplot()

# Wilcoxin ranked sum test
# Use CRPC as the base
ggplot(class_filt_pa, aes(x=Sample.Type, y=HEPACAM2))+
  geom_boxplot()
class_filt_pa = subset(prostate_atlas_genes_added, prostate_atlas_genes_added$Sample.Type %in% c("NORMAL", "NEPC") ) # Filter for 2 classes
class_filt_pa$Sample.Type = factor(class_filt_pa$Sample.Type)
class_filt_pa$Sample.Type <- relevel(class_filt_pa$Sample.Type,"NEPC") # Set which category will be the base
wilcox.test(class_filt_pa$HEPACAM2 ~ class_filt_pa$Sample.Type, conf.int = T) # Two sided Test between the specified groups
wilcox.test(class_filt_pa$NTN4 ~ class_filt_pa$Sample.Type, conf.int = T) # Two sided Test between the specified groups
wilcox.test(class_filt_pa$BLCAP ~ class_filt_pa$Sample.Type, conf.int = T) # Two sided Test between the specified groups
wilcox.test(class_filt_pa$RUNX1 ~ class_filt_pa$Sample.Type, conf.int = T) # Two sided Test between the specified groups
wilcox.test(class_filt_pa$NRXN3 ~ class_filt_pa$Sample.Type, conf.int = T) # Two sided Test between the specified groups
wilcox.test(class_filt_pa$GALNTL6 ~ class_filt_pa$Sample.Type, conf.int = T) # Two sided Test between the specified groups



wilcox.test(HEPACAM2 ~ Sample.Type, data = class_filt_pa, alternative = "greater", conf.int = T) # One-sided Test
# Multiple comparisons correction for each set of genes marking a trajectory state
pvalue_vector = c(2.04E-08, 0.177, 0.1283, 0.06755, 7.36E-05, 0.5823) # NEPC vs CRPC
pvalue_vector = c(1.41E-08, 2.52E-09, 0.02816, 0.03705, 2.56E-05, 0.02864) # NEPC vs PRIMARY
pvalue_vector = c(9.94E-10, 1.06E-11, 0.9638, 0.0001634, 0.7604, 0.4810191) # NEPC vs NORMAL
p.adjust(pvalue_vector, method = "BH")

cor.test(atlas_meta_NEPC$ASCL1, atlas_meta_NEPC$PCDH9, method = "spearman")

# Primary Prostate Cancer Dataset - Erho
##Feng dataset #For investigating whether genes mark differences in metastasis or gleason grade
expr <- read.delim("/trigos_team/CASCADE/ExternalDatasets/GSE46691/GSE46691_expression_aggregated.txt")
rownames(expr) <- expr$ID_REF
expr$ID_REF <- NULL

metadata <- read.csv("/trigos_team/CASCADE/ExternalDatasets/GSE46691/GSE46691_series_matrix_for_R.csv")

metadata <- metadata[,colnames(metadata) %in% c("X.Sample_description", "X.Sample_characteristics_ch1",
                                                "X.Sample_characteristics_ch1.1")]
metadata <- metadata[,c(3,1,2)]
colnames(metadata) <- c("Sample", "Gleason", "Metastasis")

metadata$Metastasis <- gsub("metastatic event: ", "", metadata$Metastasis)
metadata$Gleason <- gsub("gleason score: ", "", metadata$Gleason)

HDAC9_exp <- expr[rownames(expr) == "HDAC9",]
names(HDAC9_exp) <- gsub("\\.CEL", "", names(HDAC9_exp))
metadata$HDAC9 <- unlist(HDAC9_exp[match(metadata$Sample, names(HDAC9_exp))])

AR_exp <- expr[rownames(expr) == "AR",]
names(AR_exp) <- gsub("\\.CEL", "", names(AR_exp))
metadata$AR <- unlist(AR_exp[match(metadata$Sample, names(AR_exp))])

ASCL1_exp <- expr[rownames(expr) == "ASCL1",]
names(ASCL1_exp) <- gsub("\\.CEL", "", names(ASCL1_exp))
metadata$ASCL1 <- unlist(ASCL1_exp[match(metadata$Sample, names(ASCL1_exp))])

SYP_exp <- expr[rownames(expr) == "SYP",]
names(SYP_exp) <- gsub("\\.CEL", "", names(SYP_exp))
metadata$SYP <- unlist(SYP_exp[match(metadata$Sample, names(SYP_exp))])

INSM1_exp <- expr[rownames(expr) == "INSM1",]
names(INSM1_exp) <- gsub("\\.CEL", "", names(INSM1_exp))
metadata$INSM1 <- unlist(INSM1_exp[match(metadata$Sample, names(INSM1_exp))])

HEPACAM2_exp = expr[rownames(expr) == "HEPACAM2",]
names(HEPACAM2_exp) = gsub("\\.CEL", "", names(HEPACAM2_exp))
metadata$HEPACAM2 = unlist(HEPACAM2_exp[match(metadata$Sample, names(HEPACAM2_exp))])

PCDH9_exp = expr[rownames(expr) == "PCDH9",]
names(PCDH9_exp) = gsub("\\.CEL", "", names(PCDH9_exp))
metadata$PCDH9 = unlist(PCDH9_exp[match(metadata$Sample, names(PCDH9_exp))])

SCGN_exp = expr[rownames(expr) == "SCGN",]
names(SCGN_exp) = gsub("\\.CEL", "", names(SCGN_exp))
metadata$SCGN = unlist(SCGN_exp[match(metadata$Sample, names(SCGN_exp))])

for (gene in pre_NE_genes_filtered) {
  g_exp = expr[rownames(expr) == gene,]
  names(g_exp) = gsub("\\.CEL", "", names(g_exp))
  metadata$gene = unlist(g_exp[match(metadata$Sample, names(g_exp))])
  colnames(metadata)[ncol(metadata)] = gene
}

for (gene in cluster4_dura_marker_genes_filtered) {
  g_exp = expr[rownames(expr) == gene,]
  names(g_exp) = gsub("\\.CEL", "", names(g_exp))
  metadata$gene = unlist(g_exp[match(metadata$Sample, names(g_exp))])
  colnames(metadata)[ncol(metadata)] = gene
}

for (gene in cluster7_dura_from_4_genes_filtered) {
  g_exp = expr[rownames(expr) == gene,]
  names(g_exp) = gsub("\\.CEL", "", names(g_exp))
  metadata$gene = unlist(g_exp[match(metadata$Sample, names(g_exp))])
  colnames(metadata)[ncol(metadata)] = gene
}

for (gene in cluster7_dura_marker_genes_filtered) {
  g_exp = expr[rownames(expr) == gene,]
  names(g_exp) = gsub("\\.CEL", "", names(g_exp))
  metadata$gene = unlist(g_exp[match(metadata$Sample, names(g_exp))])
  colnames(metadata)[ncol(metadata)] = gene
}

for (gene in cluster89_from_7_genes_filtered) {
  g_exp = expr[rownames(expr) == gene,]
  names(g_exp) = gsub("\\.CEL", "", names(g_exp))
  metadata$gene = unlist(g_exp[match(metadata$Sample, names(g_exp))])
  colnames(metadata)[ncol(metadata)] = gene
}

for (gene in NE_Dura_Top_Genes_filtered) {
  g_exp = expr[rownames(expr) == gene,]
  names(g_exp) = gsub("\\.CEL", "", names(g_exp))
  metadata$gene = unlist(g_exp[match(metadata$Sample, names(g_exp))])
  colnames(metadata)[ncol(metadata)] = gene
}

library(ggplot2)

ggplot(metadata, aes(x=Metastasis, y=AR))+
  geom_boxplot()

ggplot(metadata, aes(x=Metastasis, y=HEPACAM2))+
  geom_boxplot() +
ggplot(metadata, aes(x=Metastasis, y=NTN4))+
  geom_boxplot() +
ggplot(metadata, aes(x=Metastasis, y=BLCAP))+
geom_boxplot() +
ggplot(metadata, aes(x=Metastasis, y=RUNX1))+
  geom_boxplot() +
ggplot(metadata, aes(x=Metastasis, y=NRXN3))+
geom_boxplot() +
  ggplot(metadata, aes(x=Metastasis, y=GALNTL6))+
  geom_boxplot()
ggplot(metadata, aes(x=Metastasis, y=PCDH9))+
  geom_boxplot()

ggplot(metadata, aes(x=Metastasis, y=SCGN))+
  geom_boxplot()

ggplot(metadata, aes(x=Metastasis, y=HDAC9))+
  geom_boxplot()

##Samples more likely to have metastasis within 5 years have higher HDAC9
wilcox.test(metadata$HDAC9[metadata$Metastasis == 1],
            metadata$HDAC9[metadata$Metastasis == 0], alternative="greater")
wilcox.test(metadata$GALNTL6[metadata$Metastasis == 1],
            metadata$GALNTL6[metadata$Metastasis == 0], alternative="greater")

ggplot(metadata, aes(x=Gleason, y=HDAC9))+
  geom_boxplot()

ggplot(metadata, aes(x=HDAC9, y=AR))+
  geom_point()

cor(metadata$HDAC9, metadata$AR, method="sp")
cor(metadata$HDAC9, metadata$SYP, method="sp")
cor(metadata$HDAC9, metadata$ASCL1, method="sp")

cor(metadata$ASCL1, metadata$HEPACAM2, method = "sp")
cor(metadata$ASCL1, metadata$PCDH9, method = "sp")


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

##The downloaded dataset does not only include gene expression data.
#Will include only GE samples were.
#I got the sample sheet from the GSE entry
living_dataset <- read.delim("/trigos_team/CASCADE/ExternalDatasets/Public_prostate_single_cell/Living_tumour_laboratory/Living_Tumour_Laboratory_samples.txt")
further_info <- read.delim("/trigos_team/CASCADE/ExternalDatasets/Public_prostate_single_cell/Living_tumour_laboratory/Living_Tumour_Laboratory_data_sheets.txt")
living_GSE_samples <- read.delim("/trigos_team/CASCADE/ExternalDatasets/Public_prostate_single_cell/Living_tumour_laboratory/Living_Tumour_Laboratory_patient_sample_key.txt")

living_GSE_samples$File <- NA
for(sample in living_GSE_samples$Sample_ID){
  temp <- grep(sample, list.files("/trigos_team/CASCADE/ExternalDatasets/Public_prostate_single_cell/Living_tumour_laboratory/"))
  temp_name <- list.files("/trigos_team/CASCADE/ExternalDatasets/Public_prostate_single_cell/Living_tumour_laboratory/")[temp]
  living_GSE_samples$File[living_GSE_samples$Sample_ID == sample] <- paste("/trigos_team/CASCADE/ExternalDatasets/Public_prostate_single_cell/Living_tumour_laboratory/", temp_name, sep="") 
}

library(limma)
x <- read.maimages(living_GSE_samples$File,
                   source="agilent", green.only=TRUE, other.columns="gIsWellAboveBG")

agilent <- read.delim("/trigos_team/CASCADE/ExternalDatasets/Public_prostate_single_cell/Living_tumour_laboratory/GPL14550-9757.txt",
                      skip=17)

x$genes$Gene_name <- agilent$GENE_SYMBOL[match(x$genes$ProbeName, agilent$SPOT_ID)]

y <- backgroundCorrect(x, method="normexp")
y <- normalizeBetweenArrays(y, method="quantile")

Control <- y$genes$ControlType==1L
NoSymbol <- y$genes$Gene_name == ""
IsExpr <- rowSums(y$other$gIsWellAboveBG > 0) >= 5 ##highly expressed above background in at least 5 samples
yfilt <- y[!Control & !NoSymbol & IsExpr, ]
dim(yfilt)

living_data <- yfilt$E

living_GSE_samples$Description <- gsub("Patient Derived Xenograft ", "", living_GSE_samples$Description)
living_GSE_samples$Description <- gsub("Clinical ", "", living_GSE_samples$Description)
living_GSE_samples$File <- gsub(".gz", "", living_GSE_samples$File)

colnames(living_data) <- living_GSE_samples$Description[match(colnames(living_data), living_GSE_samples$File)]

rownames(living_data) <- yfilt$genes$Gene_name


living_metadata <- further_info
living_metadata$Metastasis <- living_dataset$Metastasis[match(living_metadata$Sample, living_dataset$Case.No.)]

living_metadata$Hormone_sensitivity <- gsub("Androgen -independent", "Androgen-independent", living_metadata$Hormone_sensitivity)

living_metadata$Derived_from <- gsub("LTL-313B NA", "LTL-313B", living_metadata$Derived_from)

living_metadata[living_metadata$Sample == "LTL-313D","Pathology"] <- "High grade adenocarcinoma"
living_metadata[living_metadata$Sample == "LTL-313D","Derived_from"] <- "Primary"
living_metadata[living_metadata$Sample == "LTL-313D","Doubling_time"] <- "11"

living_metadata <- rbind(living_metadata,
                         unlist(c("LTL-311-G8", living_metadata[living_metadata$Sample == "LTL-311",2:6])),
                         unlist(c("LTL-313B-G5", living_metadata[living_metadata$Sample == "LTL-313B",2:6])),
                         unlist(c("LTL-313B-G8", living_metadata[living_metadata$Sample == "LTL-313B",2:6])),
                         unlist(c("LTL-313H-G16", living_metadata[living_metadata$Sample == "LTL-313H",2:6])),
                         unlist(c("LTL-313H-G7", living_metadata[living_metadata$Sample == "LTL-313H",2:6])),
                         unlist(c("LTL-352-G13", living_metadata[living_metadata$Sample == "LTL-352",2:6])),
                         unlist(c("LTL-412-G5", living_metadata[living_metadata$Sample == "LTL-412",2:6])),
                         unlist(c("LTL-418-G6", living_metadata[living_metadata$Sample == "LTL-418",2:6])))


living_metadata <- rbind(living_metadata,
                         unlist(c("LTL-311-G8", living_metadata[living_metadata$Sample == "LTL-311",2:6])),
                         unlist(c("LTL-310", living_metadata[living_metadata$Sample == "LTL-310F",2:6])),
                         unlist(c("LTL-310R", "Androgen-independent", "CRPC", "LTL-310", NA, NA)),
                         unlist(c("LTL-311R", "Androgen-independent", NA, "LTL-311", NA, NA)),
                         unlist(c("LTL-313HR", "Androgen-independent", "CRPC", "LTL-313H", NA, NA)),
                         unlist(c("LTL-418R", "Androgen-independent", NA, "LTL-418", NA, NA)),
                         unlist(c("LTL-467R", "Androgen-independent", NA, "LTL-467", NA, NA)),
                         unlist(c("LTL-484R", "Androgen-independent", NA, "LTL-484", NA, NA)))


living_data_HDAC9 <- living_data[rownames(living_data) == "HDAC9",]
living_data_HDAC9 <- apply(living_data_HDAC9, 2, median)

living_data_AR <- living_data[rownames(living_data) == "AR",]
living_data_AR <- apply(living_data_AR, 2, median)

living_data_SYP <- living_data[rownames(living_data) == "SYP",]
living_data_SYP <- apply(living_data_SYP, 2, median)

living_data_ASCL1 <- living_data[rownames(living_data) == "ASCL1",]
living_data_ASCL1 <- apply(living_data_ASCL1, 2, median)

living_data_EZH2 <- living_data[rownames(living_data) == "EZH2",]
living_data_EZH2 <- apply(living_data_EZH2, 2, median)

living_data_INSM1 <- living_data[rownames(living_data) == "INSM1",]
living_data_INSM1 <- apply(living_data_INSM1, 2, median)


library(reshape2)
living_data_HDAC9 <- melt(living_data_HDAC9)
colnames(living_data_HDAC9) <- "HDAC9"
living_data_HDAC9$AR <- living_data_AR
living_data_HDAC9$SYP <- living_data_SYP
living_data_HDAC9$ASCL1 <- living_data_ASCL1
living_data_HDAC9$EZH2 <- living_data_EZH2
living_data_HDAC9$INSM1 <- living_data_INSM1
living_data_HDAC9$Sample <- rownames(living_data_HDAC9)

living_data_HDAC9 <- melt(living_data_HDAC9)
colnames(living_data_HDAC9) <- c("Sample", "Gene", "Expr")

living_data_HDAC9$Hormone_sensitivty <- living_metadata$Hormone_sensitivity[match(living_data_HDAC9$Sample, living_metadata$Sample)]
living_data_HDAC9$Pathology <- living_metadata$Pathology[match(living_data_HDAC9$Sample, living_metadata$Sample)]
living_data_HDAC9$Doubling_time <- living_metadata$Doubling_time[match(living_data_HDAC9$Sample, living_metadata$Sample)]
living_data_HDAC9$Metastasis <- living_metadata$Metastasis[match(living_data_HDAC9$Sample, living_metadata$Sample)]

living_data_HDAC9$Pathology[living_data_HDAC9$Pathology == "Castration-resistant prostate adenocarcinoma"] <- "CRPC"
living_data_HDAC9$Pathology[living_data_HDAC9$Pathology == "Poorly differentiated neuroendo- crine carcinoma of the prostate"] <- "PD NE carcinoma"
living_data_HDAC9$Pathology[living_data_HDAC9$Pathology == "Metastatic small cell carcinoma of the prostate"] <- "Metastatic small cell carcinoma"


living_data_HDAC9$Patient <- substr(living_data_HDAC9$Sample, 1, 8)
living_data_HDAC9$Patient <- gsub("-$", "", living_data_HDAC9$Patient)
living_data_HDAC9$Patient <- gsub("R$", "", living_data_HDAC9$Patient)
living_data_HDAC9$Patient[living_data_HDAC9$Sample == "Patient-927"] <- "LTL-331"
living_data_HDAC9$Patient[living_data_HDAC9$Sample == "Patient-1005"] <- "LTL-412"
living_data_HDAC9$Patient[living_data_HDAC9$Sample == "Patient-1015"] <- "LTL-418"

living_data_HDAC9$Sample <- factor(living_data_HDAC9$Sample,
                                   levels=c("Patient-927","LTL-331","LTL-331R",
                                            "Patient-1005","LTL-412","LTL-412-G5",
                                            "Patient-1015","LTL-418","LTL-418-G6","LTL-418R",
                                            "LTL-310","LTL-310R",
                                            "LTL-311","LTL-311-G8","LTL-311R",
                                            "LTL-313A",
                                            "LTL-313B","LTL-313B-G5","LTL-313B-G8","LTL-313BR",
                                            "LTL-313C",
                                            "LTL-313D",    
                                            "LTL-313H","LTL-313H-G7","LTL-313H-G16","LTL-313HR",
                                            "LTL-352-G13",
                                            "LTL-370",
                                            "LTL-467","LTL-467R",
                                            "LTL-484","LTL-484R"))

living_data_HDAC9$Metastasis[living_data_HDAC9$Metastasis == "*"] <- NA
living_data_HDAC9$Metastasis[is.na(living_data_HDAC9$Metastasis)] <- "Unknown"

library(ggplot2)
ggplot(living_data_HDAC9, aes(x=Sample, y=Expr))+
  geom_point(aes(colour=Pathology, shape=Metastasis))+
  geom_line(aes(group = Patient))+
  geom_point(aes(colour=Pathology, shape=Metastasis, size=2))+
  facet_grid(Gene~Patient, scales="free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggplot(living_data_HDAC9, aes(x=Sample, y=Expr))+
  geom_point(aes(colour=Hormone_sensitivty, shape=Metastasis))+
  geom_line(aes(group = Patient))+
  geom_point(aes(colour=Hormone_sensitivty, shape=Metastasis, size=2))+
  facet_grid(Gene~Patient, scales="free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Add new genes to the Living Tumour Dataset
library(reshape2)
living_tumour_dataset <- living_data[rownames(living_data) == "HEPACAM2",]
living_tumour_dataset = melt(living_tumour_dataset)
colnames(living_tumour_dataset) = "HEPACAM2"
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

for (variable in pre_NE_genes) {
  living_tumour_dataset = add_to_living(variable, living_tumour_dataset, living_data)
}
for (variable in cluster4_dura_marker_genes) {
  living_tumour_dataset = add_to_living(variable, living_tumour_dataset, living_data)
}
for (variable in cluster7_dura_marker_genes) {
  living_tumour_dataset = add_to_living(variable, living_tumour_dataset, living_data)
}
for (variable in cluster7_dura_from_4_genes) {
  living_tumour_dataset = add_to_living(variable, living_tumour_dataset, living_data)
}
for (variable in cluster89_from_7_genes) {
  living_tumour_dataset = add_to_living(variable, living_tumour_dataset, living_data)
}
for (variable in NE_Dura_Top_Genes) {
  living_tumour_dataset = add_to_living(variable, living_tumour_dataset, living_data)
}

add_meta_living = function(living_data_HDAC9) {
  living_data_HDAC9$Sample <- rownames(living_data_HDAC9)
  
  living_data_HDAC9 <- melt(living_data_HDAC9)
  colnames(living_data_HDAC9) <- c("Sample", "Gene", "Expr")
  
  living_data_HDAC9$Hormone_sensitivty <- living_metadata$Hormone_sensitivity[match(living_data_HDAC9$Sample, living_metadata$Sample)]
  living_data_HDAC9$Pathology <- living_metadata$Pathology[match(living_data_HDAC9$Sample, living_metadata$Sample)]
  living_data_HDAC9$Doubling_time <- living_metadata$Doubling_time[match(living_data_HDAC9$Sample, living_metadata$Sample)]
  living_data_HDAC9$Metastasis <- living_metadata$Metastasis[match(living_data_HDAC9$Sample, living_metadata$Sample)]
  
  living_data_HDAC9$Pathology[living_data_HDAC9$Pathology == "Castration-resistant prostate adenocarcinoma"] <- "CRPC"
  living_data_HDAC9$Pathology[living_data_HDAC9$Pathology == "Poorly differentiated neuroendo- crine carcinoma of the prostate"] <- "PD NE carcinoma"
  living_data_HDAC9$Pathology[living_data_HDAC9$Pathology == "Metastatic small cell carcinoma of the prostate"] <- "Metastatic small cell carcinoma"
  
  
  living_data_HDAC9$Patient <- substr(living_data_HDAC9$Sample, 1, 8)
  living_data_HDAC9$Patient <- gsub("-$", "", living_data_HDAC9$Patient)
  living_data_HDAC9$Patient <- gsub("R$", "", living_data_HDAC9$Patient)
  living_data_HDAC9$Patient[living_data_HDAC9$Sample == "Patient-927"] <- "LTL-331"
  living_data_HDAC9$Patient[living_data_HDAC9$Sample == "Patient-1005"] <- "LTL-412"
  living_data_HDAC9$Patient[living_data_HDAC9$Sample == "Patient-1015"] <- "LTL-418"
  
  living_data_HDAC9$Sample <- factor(living_data_HDAC9$Sample,
                                     levels=c("Patient-927","LTL-331","LTL-331R",
                                              "Patient-1005","LTL-412","LTL-412-G5",
                                              "Patient-1015","LTL-418","LTL-418-G6","LTL-418R",
                                              "LTL-310","LTL-310R",
                                              "LTL-311","LTL-311-G8","LTL-311R",
                                              "LTL-313A",
                                              "LTL-313B","LTL-313B-G5","LTL-313B-G8","LTL-313BR",
                                              "LTL-313C",
                                              "LTL-313D",    
                                              "LTL-313H","LTL-313H-G7","LTL-313H-G16","LTL-313HR",
                                              "LTL-352-G13",
                                              "LTL-370",
                                              "LTL-467","LTL-467R",
                                              "LTL-484","LTL-484R"))
  
  living_data_HDAC9$Metastasis[living_data_HDAC9$Metastasis == "*"] <- NA
  living_data_HDAC9$Metastasis[is.na(living_data_HDAC9$Metastasis)] <- "Unknown"
  return(living_data_HDAC9)
}

living_tumour_dataset = add_meta_living(living_tumour_dataset)

# Need to subset by genes first before plotting
living_subset = filter(living_tumour_dataset, Gene %in% c("HEPACAM2", "NTN4", "BLCAP", "NRXN3", "RUNX1", "GALNTL6"))
library(ggplot2)
ggplot(living_subset, aes(x=Sample, y=Expr))+ ylab("Gene Expression") +
  geom_point(aes(colour=Pathology, shape=Metastasis))+
  geom_line(aes(group = Patient))+
  geom_point(aes(colour=Pathology, shape=Metastasis, size=2))+
  facet_grid(Gene~Patient, scales="free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggplot(living_subset, aes(x=Sample, y=Expr))+
  geom_point(aes(colour=Hormone_sensitivty, shape=Metastasis))+
  geom_line(aes(group = Patient))+
  geom_point(aes(colour=Hormone_sensitivty, shape=Metastasis, size=2))+
  facet_grid(Gene~Patient, scales="free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Focus on the mouse the best mouse model that developed NEPC - Generate Gene Signatures 
NE_mouse_data = dplyr::filter(living_tumour_dataset, Patient %in% c("LTL-331"))
NE_mouse_subset = dplyr::filter(NE_mouse_data, Gene %in% c("HEPACAM2", "NTN4", "BLCAP", "NRXN3", "RUNX1", "GALNTL6"))
NE_mouse_subset = dplyr::filter(NE_mouse_data, Gene %in% c("PCDH9","KCNJ6","TNIK","ST7","PPFIA2","KIAA0825"))
NE_mouse_subset = dplyr::filter(NE_mouse_data, Gene %in% c("PCDH9","THSD7A","CACNA2D1","RALYL","CACNA1A","KCNJ6"))
NE_mouse_subset = dplyr::filter(NE_mouse_data, Gene %in% c("HDAC9", "CNTNAP2", "SYT1", "NELL2", "ST18", "ASCL1"))
ggplot(NE_mouse_subset, aes(x=Sample, y=Expr))+ ylab("Gene Expression") +
  geom_point(aes(colour=Pathology, shape=Metastasis))+
  geom_line(aes(group = Patient))+
  geom_point(aes(colour=Pathology, shape=Metastasis), size=3)+
  facet_grid(Gene~Patient, scales="free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggplot(NE_mouse_subset, aes(x=Sample, y=Expr))+
  geom_point(aes(colour=Hormone_sensitivty, shape=Metastasis))+
  geom_line(aes(group = Patient))+
  geom_point(aes(colour=Hormone_sensitivty, shape=Metastasis, size=2))+
  facet_grid(Gene~Patient, scales="free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Gene Sig for the Mouse Model that developed NEPC
# Gene Signature Generation
library(reshape2)
NE_mouse_data_sigs = dcast(formula = Sample~Gene, data = NE_mouse_data, value.var = "Expr")

rownames(NE_mouse_data_sigs) = NE_mouse_data_sigs$Sample
NE_mouse_data_sigs = NE_mouse_data_sigs[,-c(1)]
# Need to get the sample ID's to be the row names and make sure its unique first
NE_mouse_data_sigs = t(NE_mouse_data_sigs)

NE_mouse_data_sigs = as.matrix(NE_mouse_data_sigs)
gene_sigs = list(pre_NE_genes_filtered, cluster7_dura_from_4_genes_filtered, cluster89_from_7_genes_filtered, NE_Dura_Top_Genes_filtered)

NE_mouse_data_sigs[1:5, 1:5]
library(GSVA)
row_names_df_to_remove<-c("AC083837.1","FER1L6","AC093916.1","PCED1B","LINC02466", "ST7-AS2", "PEX5L", "MIR924HG", "LHFPL6", "AC105031.2", "BRINP2", "AC099520.1", "SOX5")
NE_mouse_data_sigs = NE_mouse_data_sigs[!(row.names(NE_mouse_data_sigs) %in% row_names_df_to_remove),]
gsva.es_NE_mouse_data <- gsva(NE_mouse_data_sigs, gene_sigs, verbose=FALSE, method = "ssgsea") # Look into what NA deprecated means

gsva.es_NE_mouse_data_df = as.data.frame(t(gsva.es_NE_mouse_data))
colnames(gsva.es_NE_mouse_data_df) = c("pre_NE_sig", "cluster7_from_4_sig", "cluster89_from_7_sig", "Top_NE_Dura_sig")
gsva.es_NE_mouse_data_df$Sample = rownames(gsva.es_NE_mouse_data_df)
gsva.es_NE_mouse_data_df$Patient = "LTL-331"
gsva.es_NE_mouse_data_df$Sample = factor(gsva.es_NE_mouse_data_df$Sample, c("Patient-927", "LTL-331", "LTL-331R"))
ggplot(gsva.es_NE_mouse_data_df, aes(x=Sample, y=pre_NE_sig)) + ylab("Pre-NE Signature") +
  geom_point(size=4)+
  geom_line(aes(group = Patient))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
ggplot(gsva.es_NE_mouse_data_df, aes(x=Sample, y=cluster7_from_4_sig)) + ylab("Partial-NE Signature") +
  geom_point(size=4)+
  geom_line(aes(group = Patient))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
ggplot(gsva.es_NE_mouse_data_df, aes(x=Sample, y=cluster89_from_7_sig)) + ylab("Full-NE Signature") +
  geom_point(size=4)+
  geom_line(aes(group = Patient))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(gsva.es_NE_mouse_data_df, aes(x=Sample, y=Top_NE_Dura_sig)) + ylab("Top_NE_Dura_sig") +
  geom_point()+
  geom_line(aes(group = Patient))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# Dynamic prostate cancer transcriptome analysis delineates the trajectory to disease progression #Paper has the Prostate Cancer Atlas in it as well
# PDX's Pre and Post Castration 
library(Seurat)
library(Matrix)
path <- "/trigos_team/CASCADE/ExternalDatasets/Public_prostate_single_cell/Dynamic_prostate_cancer_transcriptome/E-MTAB-9903/"

Read10X_h52 = function (filename, use.names = TRUE, unique.features = TRUE) 
{
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Please install hdf5r to read HDF5 files")
  }
  if (!file.exists(filename)) {
    stop("File not found")
  }
  infile <- hdf5r::H5File$new(filename = filename, mode = "r")
  genomes <- names(x = infile)
  output <- list()
  if (hdf5r::existsGroup(infile, "matrix")) {
    if (use.names) {
      feature_slot <- "features/name"
    }
    else {
      feature_slot <- "features/id"
    }
  }
  else {
    if (use.names) {
      feature_slot <- "gene_names"
    }
    else {
      feature_slot <- "genes"
    }
  }
  for (genome in genomes) {
    counts <- infile[[paste0(genome, "/data")]]
    indices <- infile[[paste0(genome, "/indices")]]
    indptr <- infile[[paste0(genome, "/indptr")]]
    shp <- infile[[paste0(genome, "/shape")]]
    features <- infile[[paste0(genome, "/", feature_slot)]][]
    barcodes <- infile[[paste0(genome, "/barcodes")]]
    sparse.mat <- sparseMatrix(i = indices[] + 1, p = indptr[], 
                               x = as.numeric(x = counts[]), dims = shp[], repr = "T")
    if (unique.features) {
      features <- make.unique(names = features)
    }
    rownames(x = sparse.mat) <- features
    colnames(x = sparse.mat) <- barcodes[]
    sparse.mat <- as(object = sparse.mat, Class = "CsparseMatrix")
    if (infile$exists(name = paste0(genome, "/features"))) {
      types <- infile[[paste0(genome, "/features/feature_type")]][]
      types.unique <- unique(x = types)
      if (length(x = types.unique) > 1) {
        message("Genome ", genome, " has multiple modalities, returning a list of matrices for this genome")
        sparse.mat <- sapply(X = types.unique, FUN = function(x) {
          return(sparse.mat[which(x = types == x), ])
        }, simplify = FALSE, USE.NAMES = TRUE)
      }
    }
    output[[genome]] <- sparse.mat
  }
  infile$close_all()
  if (length(x = output) == 1) {
    return(output[[genome]])
  }
  else {
    return(output)
  }
} # Edited since the deprecation of the matrix function

file1 <- Read10X_h52(paste(path, "LNCaP_castr_residual_filtered_feature_bc_matrix.h5", sep="/"))
file1_Seurat <- CreateSeuratObject(counts = file1, project = "LNCaP_castr_residual")  
file2 <- Read10X_h52(paste(path, "LNCaP_castr_residual_GSK_filtered_feature_bc_matrix.h5", sep="/"))
file2_Seurat <- CreateSeuratObject(counts = file2, project = "LNCaP_castr_residual_GSK")  
file3 <- Read10X_h52(paste(path, "pdx_LNCaP_castr_regrow_filtered_feature_bc_matrix.h5", sep="/"))
file3_Seurat <- CreateSeuratObject(counts = file3, project = "pdx_LNCaP_castr_regrow")  
file4 <- Read10X_h52(paste(path, "pdx_LNCaP_castr_regrow_GSK_filtered_feature_bc_matrix.h5", sep="/"))
file4_Seurat <- CreateSeuratObject(counts = file4, project = "pdx_LNCaP_castr_regrow_GSK_filtered")  
file5 <- Read10X_h52(paste(path, "pdx_LNCaP_filtered_feature_bc_matrix.h5", sep="/"))
file5_Seurat <- CreateSeuratObject(counts = file5, project = "pdx_LNCaP")  
file6 <- Read10X_h52(paste(path, "pdx_PDX145_2_filtered_feature_bc_matrix.h5", sep="/"))
file6_Seurat <- CreateSeuratObject(counts = file6, project = "pdx_PDX145_2")  
file7 <- Read10X_h52(paste(path, "pdx_PDX147_filtered_feature_bc_matrix.h5", sep="/"))
file7_Seurat <- CreateSeuratObject(counts = file7, project = "pdx_PDX147")  
file8 <- Read10X_h52(paste(path, "pdx_PDX147_regrow_filtered_feature_bc_matrix.h5", sep="/"))
file8_Seurat <- CreateSeuratObject(counts = file8, project = "pdx_PDX147_regrow")  
file9 <- Read10X_h52(paste(path, "pdx_PDX35_filtered_feature_bc_matrix.h5", sep="/"))
file9_Seurat <- CreateSeuratObject(counts = file9, project = "pdx_PDX35")  
file10 <- Read10X_h52(paste(path, "pdx_PDX78_filtered_feature_bc_matrix.h5", sep="/"))
file10_Seurat <- CreateSeuratObject(counts = file10, project = "pdx_PDX78")  
file11 <- Read10X_h52(paste(path, "pdx_PNPCa_filtered_feature_bc_matrix.h5", sep="/"))
file11_Seurat <- CreateSeuratObject(counts = file11, project = "pdx_PNPCa")  

PDX_merge <- merge(file1_Seurat, y=c(file2_Seurat, file3_Seurat, file4_Seurat, file5_Seurat,
                                     file6_Seurat, file7_Seurat, file8_Seurat, file9_Seurat,
                                     file10_Seurat, file11_Seurat),
                   add.cell.ids=c("LNCaP_castr_residual","LNCaP_castr_residual_GSK","pdx_LNCaP_castr_regrow",
                                  "pdx_LNCaP_castr_regrow_GSK_filtered", 
                                  "pdx_LNCaP","pdx_PDX145_2",  "pdx_PDX147",
                                  "pdx_PDX147_regrow","pdx_PDX35","pdx_PDX78","pdx_PNPCa"))


PDX_merge$PDX <- unname(PDX_merge$orig.ident)

PDX_merge <- NormalizeData(PDX_merge)
PDX_merge <- ScaleData(PDX_merge)
PDX_merge <- FindVariableFeatures(PDX_merge)
PDX_merge <- RunPCA(PDX_merge)
PDX_merge <- FindNeighbors(PDX_merge)
PDX_merge <- FindClusters(PDX_merge)
PDX_merge <- RunUMAP(PDX_merge, dims=1:20)

DimPlot(PDX_merge, group="PDX", split.by="PDX")

DimPlot(PDX_merge, group.by = "PDX") # Added this line

FeaturePlot(PDX_merge, c("GRCh38-AR", "GRCh38-SYP", "GRCh38-ASCL1", "GRCh38-HDAC9"))

#Only samples with HDAC9

PDX_merge2 <- merge(file6_Seurat, y=c(file7_Seurat, file8_Seurat, file9_Seurat),
                    add.cell.ids=c("pdx_PDX145_2",  "pdx_PDX147",
                                   "pdx_PDX147_regrow","pdx_PDX35"))


PDX_merge2$PDX <- unname(PDX_merge2$orig.ident)

PDX_merge2 <- NormalizeData(PDX_merge2)
PDX_merge2 <- ScaleData(PDX_merge2)
PDX_merge2 <- FindVariableFeatures(PDX_merge2)
PDX_merge2 <- RunPCA(PDX_merge2)
PDX_merge2 <- FindNeighbors(PDX_merge2)
PDX_merge2 <- FindClusters(PDX_merge2)
PDX_merge2 <- RunUMAP(PDX_merge2, dims=1:20)

DimPlot(PDX_merge2, group="PDX")

FeaturePlot(PDX_merge2, c("GRCh38-AR", "GRCh38-SYP", "GRCh38-ASCL1", "GRCh38-HDAC9"))

# How does the expression of Pre-NE Genes look in these PDX's
PDX_mergeNE = PDX_merge2
FeaturePlot(PDX_mergeNE, c("GRCh38-HEPACAM2", "GRCh38-ASCL1", "GRCh38-PCDH9", "GRCh38-AR")) # Need GRCh38- in the name

FeaturePlot(PDX_mergeNE, features = c("GRCh38-HEPACAM2", "GRCh38-NTN4", "GRCh38-BLCAP", "GRCh38-RUNX1", "GRCh38-GALNTL6", "GRCh38-NRXN3")) # The 6 Pre-NE candidates from CA0027 Dura

FeaturePlot(PDX_mergeNE, features = c("GRCh38-PCDH9", "GRCh38-KCNJ6", "GRCh38-TNIK", "GRCh38-ST7", "GRCh38-PPFIA2", "GRCh38-KIAA0825")) # Top 6 Cluster 7 from 4 Dura

FeaturePlot(PDX_mergeNE, features = c("GRCh38-HDAC9", "GRCh38-CNTNAP2", "GRCh38-SYT1", "GRCh38-NELL2", "GRCh38-ST18", "GRCh38-ASCL1")) # Top 6 Cluster 8 & 9 from 7 Dura

FeaturePlot(PDX_mergeNE, features = c("GRCh38-PCDH9", "GRCh38-THSD7A", "GRCh38-CACNA2D1", "GRCh38-RALYL", "GRCh38-CACNA1A", "GRCh38-KCNJ6")) # Top 6 Dura NE above 2.5 LFC

FeaturePlot(PDX_mergeNE, c("GRCh38-AR","GRCh38-ASCL1"))

DimPlot(PDX_mergeNE, group.by = "seurat_clusters", label = T)

FeaturePlot(PDX_mergeNE, c("GRCh38-MKI67")) # Cycling Cells

PDX_35_all = FindMarkers(PDX_mergeNE, ident.1 = c(4,6,10,17), min.pct = 0.25)

PDX_35_adeno_from_NE = FindMarkers(PDX_mergeNE, ident.1 = c(4,10,17), ident.2 = c(6), min.pct = 0.25)

PDX_35_NE_from_adeno = FindMarkers(PDX_mergeNE, ident.1 = c(6), ident.2 = c(4,10,17), min.pct = 0.25)


Find_best_genes = function(gene_data) {
  genes_list = c()
  gene_data = gene_data[order(gene_data$avg_log2FC, decreasing = TRUE),]
  for(i in 1 : nrow(gene_data)) {
    if(gene_data[i, 2] > 1.5 & (gene_data[i, 3] - gene_data[i, 4]) > 0.2) {
      genes_list = c(genes_list, rownames(gene_data[i,]))
    }
  }
  if (is.null(genes_list) == TRUE) {
    genes_list = c(genes_list, "NA")
  }
  return(genes_list)
}

PDX_35_adeno_from_NE_best = Find_best_genes(PDX_35_adeno_from_NE)

PDX_35_NE_from_adeno_best = Find_best_genes(PDX_35_NE_from_adeno)

Cluster10_from_4_PDX35 = FindMarkers(PDX_mergeNE, ident.1 = c(10), ident.2 = c(4), min.pct = 0.25)

Cluster10_from_4_PDX35_best = Find_best_genes(Cluster10_from_4_PDX35)

PDX_147_from_adeno = FindMarkers(PDX_mergeNE, ident.1 = c(1,3,5,7,13), ident.2 = c(4,10,17), min.pct = 0.25)

PDX_147_from_adeno_best = Find_best_genes(PDX_147_from_adeno)

#Is Cluster 0 Pre-NE then? or just ASCL1 negative PCDH9 positive NEPC
cluster0_from_adeno = FindMarkers(PDX_mergeNE, ident.1 = c(0), ident.2 = c(4,10,17), min.pct = 0.25)
cluster0_from_adeno_best = Find_best_genes(cluster0_from_adeno)

# Generate Gene Signatures - Probably less useful here
pre_NE_genes_filtered_PDX = c("GRCh38-HEPACAM2", "GRCh38-NTN4", "GRCh38-BLCAP", "GRCh38-RUNX1", "GRCh38-NRXN3", "GRCh38-GALNTL6")
cluster7_dura_from_4_genes_filtered_PDX = c("GRCh38-PCDH9","GRCh38-KCNJ6","GRCh38-TNIK","GRCh38-ST7","GRCh38-PPFIA2","GRCh38-KIAA0825","GRCh38-PTPRN2","GRCh38-GAS2L3","GRCh38-ODF2L","GRCh38-NPAS3","GRCh38-PCED1B","GRCh38-CIT","GRCh38-HEPACAM2",
                                        "GRCh38-THSD7A","GRCh38-SOX4","GRCh38-PODXL","GRCh38-MTSS1","GRCh38-AMIGO2","GRCh38-SRGAP1","GRCh38-CDH6","GRCh38-GPR137C","GRCh38-SLC18A1","GRCh38-EPHA5","GRCh38-PTPRE","GRCh38-PBX1",
                                        "GRCh38-TGFBR3","GRCh38-SDK2","GRCh38-ERBB4","GRCh38-GRIK2","GRCh38-SCGN","GRCh38-SORBS1","GRCh38-SOX5","GRCh38-ST7-AS2","GRCh38-MID1")
cluster89_from_7_genes_filtered_PDX = c("GRCh38-HDAC9", "GRCh38-CNTNAP2","GRCh38-SYT1","GRCh38-NELL2","GRCh38-ST18","GRCh38-ASCL1","GRCh38-CACNA2D1","GRCh38-SMOC2","GRCh38-SDK1","GRCh38-RIMBP2","GRCh38-RALYL","GRCh38-MIAT","GRCh38-THSD7A","GRCh38-ADCY2","GRCh38-TRPM3","GRCh38-KCNMA1","GRCh38-PEX5L","GRCh38-DIRAS2","GRCh38-LINC00511","GRCh38-SLCO3A1","GRCh38-TRPM8",
                                    "GRCh38-CAMK1D","GRCh38-DMD","GRCh38-CACNA1A","GRCh38-CHGA","GRCh38-AMACR","GRCh38-MGAT4C","GRCh38-ATP8A2","GRCh38-NRXN1","GRCh38-MIR924HG","GRCh38-RIMS2","GRCh38-RASSF6","GRCh38-TNS3","GRCh38-NKAIN2","GRCh38-SDK2","GRCh38-BCL2","GRCh38-DPYSL3","GRCh38-SPOCK1",
                                    "GRCh38-JAZF1","GRCh38-NLGN1","GRCh38-STXBP5L","GRCh38-PAM","GRCh38-DPP6","GRCh38-BRINP2","GRCh38-NFATC2","GRCh38-ERC2","GRCh38-NFASC","GRCh38-NPNT", "GRCh38-ESRRG","GRCh38-JAKMIP2","GRCh38-SRRM4")


PDX_mergeNE = AddModuleScore(
  object = PDX_mergeNE,
  features = list(pre_NE_genes_filtered_PDX),
  name = 'pre_diff_NE_sig_PDX')

PDX_mergeNE = AddModuleScore(
  object = PDX_mergeNE,
  features = list(cluster7_dura_from_4_genes_filtered_PDX),
  name = 'partial_diff_NE_sig_PDX')

PDX_mergeNE = AddModuleScore(
  object = PDX_mergeNE,
  features = list(cluster89_from_7_genes_filtered_PDX),
  name = 'full_diff_NE_sig_PDX')

FeaturePlot(PDX_mergeNE, features = "pre_diff_NE_sig_PDX1")
FeaturePlot(PDX_mergeNE, features = "partial_diff_NE_sig_PDX1")
FeaturePlot(PDX_mergeNE, features = "full_diff_NE_sig_PDX1")


