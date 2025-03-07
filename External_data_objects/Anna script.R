##HDAC9 investigation
#Can we find HDAC9 in other datasets (patients and also in vivo models)?

library(Seurat)
#Single-cell analysis supports a luminal-neuroendocrine transdifferentiation in human prostate cancer
##https://www.nature.com/articles/s42003-020-01476-1

##Patients 2, 5 and 6 have NE

##Patient 2
raw_counts_o<-read.delim(file="/trigos_team/CASCADE/ExternalDatasets/Public_prostate_single_cell/Luminal_neuro_transdifferentiation/GSE137829/GSM4089152_P2_gene_cell_exprs_table.txt")
duplicates <- c("CCDC7", "CYB561D2", "LINC01481", "MATR3", "PGM5-AS1", "PIK3R3", "RAET1E-AS1", "RGS5", "SPATA13", "TMEM256-PLSCR3")
raw_counts2 <- raw_counts_o[!(raw_counts_o$Symbol %in% duplicates),]
rownames(raw_counts2) <- raw_counts2[,2]
raw_counts <- raw_counts2[, -c(1,2)]
patient2 <- CreateSeuratObject(counts = raw_counts, min.cells = 3, project = "GSE137829")
patient2[["percent.mt"]] <- PercentageFeatureSet(patient2, pattern = "^MT-")
patient2[["percent.rp"]] <- PercentageFeatureSet(patient2, pattern = "^RP")
VlnPlot(patient2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",
                               "percent.rp"), ncol = 4)
patient2_filtered <- subset(patient2, subset = nFeature_RNA > 200  & percent.mt < 10)

patient2_filtered <- NormalizeData(patient2_filtered)
patient2_filtered <- ScaleData(patient2_filtered)
patient2_filtered <- FindVariableFeatures(patient2_filtered)
patient2_filtered <- RunPCA(patient2_filtered)
patient2_filtered <- FindNeighbors(patient2_filtered)
patient2_filtered <- FindClusters(patient2_filtered)
patient2_filtered <- RunUMAP(patient2_filtered, dims=1:20)

DimPlot(patient2_filtered)

FeaturePlot(patient2_filtered, c("AR", "SYP", "CHGA", "ASCL1", "HDAC9"))



##Patient 5
raw_counts_o<-read.delim(file="/trigos_team/CASCADE/ExternalDatasets/Public_prostate_single_cell/Luminal_neuro_transdifferentiation/GSE137829/GSM4711414_P5_gene_cell_exprs_table.txt")
duplicates <- c("ABCF2", "ATXN7", "CCDC39", "COG8", "CYB561D2", "EMG1", "HSPA14", "LINC01238", "MATR3", "POLR2J3", "RGS5", "TBCE", "TMSB15B")
raw_counts2 <- raw_counts_o[!(raw_counts_o$Symbol %in% duplicates),]
rownames(raw_counts2) <- raw_counts2[,2]
raw_counts <- raw_counts2[, -c(1,2)]

patient5 <- CreateSeuratObject(counts = raw_counts, min.cells = 3, project = "GSE137829")
patient5[["percent.mt"]] <- PercentageFeatureSet(patient5, pattern = "^MT-")
patient5[["percent.rp"]] <- PercentageFeatureSet(patient5, pattern = "^RP")

VlnPlot(patient5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",
                               "percent.rp"), ncol = 4)
patient5_filtered <- subset(patient5, subset = nFeature_RNA > 200  & percent.mt < 5)

patient5_filtered <- NormalizeData(patient5_filtered)
patient5_filtered <- ScaleData(patient5_filtered)
patient5_filtered <- FindVariableFeatures(patient5_filtered)
patient5_filtered <- RunPCA(patient5_filtered)
patient5_filtered <- FindNeighbors(patient5_filtered)
patient5_filtered <- FindClusters(patient5_filtered)
patient5_filtered <- RunUMAP(patient5_filtered, dims=1:20)

DimPlot(patient5_filtered)

FeaturePlot(patient5_filtered, c("AR", "SYP", "CHGA", "ASCL1", "HDAC9"))



##Patient 6
raw_counts_o<-read.delim(file="/trigos_team/CASCADE/ExternalDatasets/Public_prostate_single_cell/Luminal_neuro_transdifferentiation/GSE137829/GSM4711415_P6_gene_cell_exprs_table.txt")
duplicates <- c("ABCF2", "ATXN7", "COG8", "CYB561D2", "DIABLO", "EMG1", "HSPA14", "LINC01238", "MATR3", "POLR2J3", "RGS5", "SOD2", "TBCE", "TMSB15B")
raw_counts2 <- raw_counts_o[which(!raw_counts_o$Symbol %in% duplicates),]
rownames(raw_counts2) <- raw_counts2[,2]
raw_counts <- raw_counts2[, -c(1,2)]
patient6 <- CreateSeuratObject(counts = raw_counts, min.cells = 3, project = "GSE137829")
patient6[["percent.mt"]] <- PercentageFeatureSet(patient6, pattern = "^MT-")
patient6[["percent.rp"]] <- PercentageFeatureSet(patient6, pattern = "^RP")
VlnPlot(patient6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",
                               "percent.rp"), ncol = 4)
patient6_filtered <- subset(patient6, subset = nFeature_RNA > 200  & percent.mt < 10)

patient6_filtered <- NormalizeData(patient6_filtered)
patient6_filtered <- ScaleData(patient6_filtered)
patient6_filtered <- FindVariableFeatures(patient6_filtered)
patient6_filtered <- RunPCA(patient6_filtered)
patient6_filtered <- FindNeighbors(patient6_filtered)
patient6_filtered <- FindClusters(patient6_filtered)
patient6_filtered <- RunUMAP(patient6_filtered, dims=1:20)

DimPlot(patient6_filtered)

FeaturePlot(patient6_filtered, c("AR", "SYP", "CHGA", "ASCL1", "HDAC9"))



# metadata <- read.delim("/trigos_team/CASCADE/ExternalDatasets/Public_prostate_single_cell/Luminal_neuro_transdifferentiation/Metadata.txt")
#
# metadata$Patient <- metadata$orig.ident
# metadata$Patient <- gsub(" #", "", metadata$Patient)
# metadata$Cell_ID <- paste(metadata$Patient, metadata$X, sep="_")
# #data1_filtered$Patient <- metadata$Patient[match(colnames(data1_filtered), metadata$X)]
#
#
# data1 <- merge(patient2_filtered, y=c(patient5_filtered, patient6_filtered),
#                add.cell.ids = c("patient2", "patient5", "patient6"),
#                project = "NE_patients")
#
# data1$Patient <- substr(names(data1$orig.ident), 1, 8)
# data1 <- NormalizeData(data1)
# data1 <- ScaleData(data1)
# data1 <- FindVariableFeatures(data1)
# data1 <- RunPCA(data1)
# data1 <- FindNeighbors(data1)
# data1 <- FindClusters(data1)
# data1 <- RunUMAP(data1, dims=1:20)
#
# DimPlot(data1, group="Patient", split.by = "Patient")
#
# FeaturePlot(data1, c("AR", "SYP", "CHGA"), split.by = "Patient")




#Transcriptional mediators of treatment resistance in lethal prostate cancer
#https://www.nature.com/articles/s41591-021-01244-6

metadata <- read.delim("/trigos_team/CASCADE/ExternalDatasets/Public_prostate_single_cell/Treatment_resistance_in_lethal_prostate_cancer/scp_clustering.tsv")
metadata <- metadata[-1,]

tissue_of_origin <- read.delim("/trigos_team/CASCADE/ExternalDatasets/Public_prostate_single_cell/Treatment_resistance_in_lethal_prostate_cancer/scp_metadata.tsv")

metadata$Tissue_of_origin <- tissue_of_origin$organ__ontology_label[match(metadata$NAME,
                                                                          tissue_of_origin$NAME)]
metadata$NAME <- paste("Cell", metadata$NAME, sep="_")
colnames(metadata)[2:3] <- c("Cell_type", "Cluster_name")

treatments <- read.delim("/trigos_team/CASCADE/ExternalDatasets/Public_prostate_single_cell/Treatment_resistance_in_lethal_prostate_cancer/Table_S1.txt")

treatments$Abi_enza <- paste(treatments$prior.enzalutamide, treatments$prior.abiraterone)
treatments$Abi_enza[treatments$Abi_enza == "True True"] <- "Abi_enza"
treatments$Abi_enza[treatments$Abi_enza == "True False"] <- "Abi_enza"
treatments$Abi_enza[treatments$Abi_enza == "False True"] <- "Abi_enza"
treatments$Abi_enza[treatments$Abi_enza == "False False"] <- "Naive"

#CHGA marks small cell carcinoma cells
#09171135 has a small cell carcinoma histology -> abi enza naive, but given platinum

treatments$Histology <- "Adenocarcinoma"
treatments$Histology[treatments$biopsy == "09171135"] <- "Small_cell"

participant <- read.delim("/trigos_team/CASCADE/ExternalDatasets/Public_prostate_single_cell/Treatment_resistance_in_lethal_prostate_cancer/Table_S4.txt")
participant$cell_id <- paste("Cell", participant$cell_id, sep="_")

scp_tpm.tsv. #gene counts, but TPM because continuous
cells <- read.delim("/trigos_team/CASCADE/ExternalDatasets/Public_prostate_single_cell/Treatment_resistance_in_lethal_prostate_cancer/scp_tpm.tsv")

rownames(cells) <- cells[,1]
cells <- cells[,-1]
colnames(cells) <- paste("Cell", colnames(cells),sep="_")
colnames(cells) <- gsub("X", "", colnames(cells))

SMART_lethal <- CreateSeuratObject(counts = cells, min.cells = 3, min.genes = 200,
                                   project = "SMART_lethal")

SMART_lethal[["percent.mt"]] <- PercentageFeatureSet(SMART_lethal, pattern = "^MT-")

VlnPlot(SMART_lethal, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

SMART_lethal <- NormalizeData(SMART_lethal)
SMART_lethal <- ScaleData(SMART_lethal)
SMART_lethal <- FindVariableFeatures(SMART_lethal, selection.method = "vst", nfeatures = 2000)
SMART_lethal <- RunPCA(SMART_lethal)
SMART_lethal <- FindNeighbors(SMART_lethal, dims = 1:10)
SMART_lethal <- FindClusters(SMART_lethal, resolution = 0.8)
SMART_lethal <- RunUMAP(SMART_lethal, dims = 1:10)

colnames(metadata)[2:3] <- c("Cell_type", "Cluster_name")

SMART_lethal@meta.data$Cell_type <- metadata$Cell_type[match(rownames(SMART_lethal@meta.data), metadata$NAME)]
SMART_lethal@meta.data$Cluster_name <- metadata$Cluster_name[match(rownames(SMART_lethal@meta.data), metadata$NAME)]

SMART_lethal@meta.data$Patient <- participant$participant[match(rownames(SMART_lethal@meta.data), participant$cell_id)]
SMART_lethal@meta.data$Biopsy <- participant$biopsy[match(rownames(SMART_lethal@meta.data), participant$cell_id)]

SMART_lethal@meta.data$Abi_enza <- treatments$Abi_enza[match(SMART_lethal@meta.data$Biopsy, treatments$biopsy)]
SMART_lethal@meta.data$Histology <- treatments$Histology[match(SMART_lethal@meta.data$Biopsy, treatments$biopsy)]


DimPlot(SMART_lethal, reduction = "umap", group.by="Cluster_name")

DimPlot(SMART_lethal, reduction = "umap", group.by="Cell_type")

DimPlot(SMART_lethal, reduction = "umap", group.by="Patient")

DimPlot(SMART_lethal, reduction = "umap", group.by="Biopsy")

DimPlot(SMART_lethal, reduction = "umap", group.by="Abi_enza")

DimPlot(SMART_lethal, reduction = "umap", group.by="Abi_enza", split.by="Abi_enza")
FeaturePlot(SMART_lethal, c("AR", "SYP", "CHGA", "HDAC9"))

DimPlot(SMART_lethal, reduction = "umap", group.by="Histology", split.by="Histology")
FeaturePlot(SMART_lethal, c("AR", "SYP", "CHGA", "HDAC9"), split.by="Histology")


##Feng dataset
expr <- read.delim("/home/atrigos/GSE46691/GSE46691_expression_aggregated.txt")
rownames(expr) <- expr$ID_REF
expr$ID_REF <- NULL

metadata <- read.csv("/home/atrigos/GSE46691/GSE46691_series_matrix_for_R.csv")

metadata <- metadata[,colnames(metadata) %in% c("X.Sample_description", "X.Sample_characteristics_ch1",
                                                "X.Sample_characteristics_ch1.1")]
metadata <- metadata[,c(3,1,2)]
colnames(metadata) <- c("Sample", "Gleason", "Metastasis")

metadata$Metastasis <- gsub("metastatic event: ", "", metadata$Metastasis)
metadata$Gleason <- gsub("gleason score: ", "", metadata$Gleason)

HDAC9_exp <- expr[rownames(expr) == "HDAC9",]
names(HDAC9_exp) <- gsub("\\.CEL", "", names(HDAC9_exp))

metadata$HDAC9 <- HDAC9_exp[match(metadata$Sample, names(HDAC9_exp))]

AR_exp <- expr[rownames(expr) == "AR",]
names(AR_exp) <- gsub("\\.CEL", "", names(AR_exp))
metadata$AR <- unlist(AR_exp[match(metadata$Sample, names(AR_exp))])

ASCL1_exp <- expr[rownames(expr) == "ASCL1",]
names(ASCL1_exp) <- gsub("\\.CEL", "", names(ASCL1_exp))
metadata$ASCL1 <- unlist(ASCL1_exp[match(metadata$Sample, names(ASCL1_exp))])

SYP_exp <- expr[rownames(expr) == "SYP",]
names(SYP_exp) <- gsub("\\.CEL", "", names(SYP_exp))
metadata$SYP <- unlist(SYP_exp[match(metadata$Sample, names(SYP_exp))])


library(ggplot2)

ggplot(metadata, aes(x=Metastasis, y=HDAC9))+
  geom_boxplot()

##Samples more likely to have metastasis within 5 years have higher HDAC9
wilcox.test(metadata$HDAC9[metadata$Metastasis == 1],
            metadata$HDAC9[metadata$Metastasis == 0], alternative="greater")

ggplot(metadata, aes(x=Gleason, y=HDAC9))+
  geom_boxplot()

ggplot(metadata, aes(x=HDAC9, y=AR))+
  geom_point()

cor(metadata$HDAC9, metadata$AR, method="sp")
cor(metadata$HDAC9, metadata$SYP, method="sp")
cor(metadata$HDAC9, metadata$ASCL1, method="sp")



#Fraser dataset -. 73 sampels
norm_data <- read.delim("/home/atrigos/Prostate_TILs_new_data3/Fraser_paper/GSE84042_Normalized_data_with_annotation.txt")
norm_data <- norm_data[!is.na(norm_data$Symbol_UCSC),]
norm_data$Entrez_Gene_ID <- NULL
norm_data$Name_UCSC <- NULL
norm_data$Chr_UCSC <- NULL
norm_data$Start_UCSC <- NULL
norm_data$End_UCSC <- NULL
rownames(norm_data) <- norm_data$Symbol_UCSC
norm_data$Symbol_UCSC <- NULL

clinical_data <- read.delim("/home/atrigos/Prostate_TILs_new_data3/Fraser_paper/Fraser_clinical_data.txt")

clinical_data <- clinical_data[,c("SampleID", "Gleason.Score", "T.Category", "PreTreatment.PSA", "Age.at.Treatment", "BCR", "Time.to.BCR..years.")]

clinical_data$SampleID <- gsub("-F1", "", clinical_data$SampleID)
clinical_data <- clinical_data[clinical_data$SampleID %in% colnames(norm_data),]

clinical_data$Grade_group <- clinical_data$Gleason.Score
clinical_data$Grade_group[clinical_data$Grade_group == "3+3"] <- "GG1"
clinical_data$Grade_group[clinical_data$Grade_group == "3+4"] <- "GG2"
clinical_data$Grade_group[clinical_data$Grade_group == "4+3"] <- "GG3"

HDAC9_exp <- unlist(norm_data[rownames(norm_data) == "HDAC9",])
clinical_data$HDAC9 <- HDAC9_exp[match(clinical_data$SampleID, names(HDAC9_exp))]

ggplot(clinical_data, aes(x=Grade_group, y=HDAC9))+
  geom_boxplot()

library(survival)
library(survminer)

clinical_data$Censored_BCR <- ifelse(clinical_data$BCR == "Yes", 1, 0)

median_HDAC9 <- median(clinical_data$HDAC9)

clinical_data$HDAC9_categorical <- ifelse(clinical_data$HDAC9 > median_HDAC9, "High_HDAC9", "Low_HDAC9")

surv <- Surv(time = clinical_data$Time.to.BCR..years.,
             event = clinical_data$Censored_BCR)
fit1 <- survfit(as.formula(surv ~ HDAC9_categorical), data = clinical_data)

ggsurvplot(fit1, data = clinical_data, pval=TRUE, palette = c("red3", "deepskyblue3"))


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
##ASCL1


##Dynamic prostate cancer transcriptome analysis delineates the trajectory to disease progression
#Single-cell data

path <- "/trigos_team/CASCADE/ExternalDatasets/Public_prostate_single_cell/Dynamic_prostate_cancer_transcriptome/E-MTAB-9903/"

file1 <- Read10X_h5(paste(path, "LNCaP_castr_residual_filtered_feature_bc_matrix.h5", sep="/"))
file1_Seurat <- CreateSeuratObject(counts = file1, project = "LNCaP_castr_residual") 
file2 <- Read10X_h5(paste(path, "LNCaP_castr_residual_GSK_filtered_feature_bc_matrix.h5", sep="/"))
file2_Seurat <- CreateSeuratObject(counts = file2, project = "LNCaP_castr_residual_GSK") 
file3 <- Read10X_h5(paste(path, "pdx_LNCaP_castr_regrow_filtered_feature_bc_matrix.h5", sep="/"))
file3_Seurat <- CreateSeuratObject(counts = file3, project = "pdx_LNCaP_castr_regrow") 
file4 <- Read10X_h5(paste(path, "pdx_LNCaP_castr_regrow_GSK_filtered_feature_bc_matrix.h5", sep="/"))
file4_Seurat <- CreateSeuratObject(counts = file4, project = "pdx_LNCaP_castr_regrow_GSK_filtered") 
file5 <- Read10X_h5(paste(path, "pdx_LNCaP_filtered_feature_bc_matrix.h5", sep="/"))
file5_Seurat <- CreateSeuratObject(counts = file5, project = "pdx_LNCaP") 
file6 <- Read10X_h5(paste(path, "pdx_PDX145_2_filtered_feature_bc_matrix.h5", sep="/"))
file6_Seurat <- CreateSeuratObject(counts = file6, project = "pdx_PDX145_2") 
file7 <- Read10X_h5(paste(path, "pdx_PDX147_filtered_feature_bc_matrix.h5", sep="/"))
file7_Seurat <- CreateSeuratObject(counts = file7, project = "pdx_PDX147") 
file8 <- Read10X_h5(paste(path, "pdx_PDX147_regrow_filtered_feature_bc_matrix.h5", sep="/"))
file8_Seurat <- CreateSeuratObject(counts = file8, project = "pdx_PDX147_regrow") 
file9 <- Read10X_h5(paste(path, "pdx_PDX35_filtered_feature_bc_matrix.h5", sep="/"))
file9_Seurat <- CreateSeuratObject(counts = file9, project = "pdx_PDX35") 
file10 <- Read10X_h5(paste(path, "pdx_PDX78_filtered_feature_bc_matrix.h5", sep="/"))
file10_Seurat <- CreateSeuratObject(counts = file10, project = "pdx_PDX78") 
file11 <- Read10X_h5(paste(path, "pdx_PNPCa_filtered_feature_bc_matrix.h5", sep="/"))
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


##Prostate cancer atlas (bulk RNAseq merged from different sources)

atlas <- read.delim("/trigos_team/CASCADE/ExternalDatasets/Public_prostate_single_cell/Dynamic_prostate_cancer_transcriptome/pca_atlas_dataset/pcatlas_dataset_vst_normalized.txt")

rownames(atlas) <- atlas$X
atlas$X <- NULL




atlas_meta <- read.delim("/trigos_team/CASCADE/ExternalDatasets/Public_prostate_single_cell/Dynamic_prostate_cancer_transcriptome/pca_atlas_dataset/pcatlas_dataset_annotations.txt")
no_primary <- atlas_meta[atlas_meta$Sample.Type %in% c("CRPC", "NEPC"),]

atlas_no_primary <- atlas[,colnames(atlas) %in% no_primary$Entry]
atlas_no_primary <- cbind(atlas[,1], atlas_no_primary)
colnames(atlas_no_primary)[1] <- "ENSEMBL_ID"

saveRDS(atlas_no_primary, file="/trigos_team/CASCADE/ExternalDatasets/Public_prostate_single_cell/Dynamic_prostate_cancer_transcriptome/pca_atlas_dataset/pcatlas_no_primary.Rds")

HDAC9 <- unlist(atlas[rownames(atlas) == "ENSG00000048052",])
AR <- unlist(atlas[rownames(atlas) == "ENSG00000169083",])
SYP <- unlist(atlas[rownames(atlas) == "ENSG00000102003",])
ASCL1 <- unlist(atlas[rownames(atlas) == "ENSG00000139352",])

atlas_meta$HDAC9 <- HDAC9
atlas_meta$AR <- AR
atlas_meta$SYP <- SYP
atlas_meta$ASCL1 <- ASCL1

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



library(pheatmap)


###Roxane's dataset, as suggested by Luc
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