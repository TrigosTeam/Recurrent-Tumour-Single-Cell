# pathtable <- read.csv("~/CASCADEpaper/paper/Fig1/pathology review summary.csv", sep = ",")
# pathtable$patient <- gsub(" $", "", pathtable$patient)
# pathology_review <- readRDS("~/CASCADEpaper/paper/Fig1/pathology_review.Rds")
# pathology_review$patient <- gsub(" $", "", pathology_review$patient)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(pals)
library(ggsci)
library(ggplot2)

pathology_review <- data.frame(
   stringsAsFactors = FALSE,
        check.names = FALSE,
                          srt.obj.sites = c("CA0027_paraaortic_lymph_node_1","CA0027_prostate_9",
                                            "CA0027_dura_base_skull_13",
                                            "CA0027_dural_inner_skull_14",
                                            "CA0034_paraaortic_lymph_node_2",
                                            "CA0034_liver_right_8","CA0034_liver_left_11",
                                            "CA0035_paraaortic_lymph_node_1",
                                            "CA0035_bladder_2",
                                            "CA0043_portal_lymph_node_4","CA0043_liver_7",
                                            "CA0043_liver_12","CA0046_hilar_lymph_node_5",
                                            "CA0046_lung_7","CA0046_liver_12",
                                            "CA0046_paraaortic_lymph_node_15",
                                            "CA0046_prostate_17",
                                            "CA0046_pelvic_lymph_node_19","CA0058_hilar_50",
                                            "CA0058_liver_29","CA0058_liver_38",
                                            "CA0076_left_rib_17",
                                            "CA0076_vertebra_25","CA0076_liver_right_40",
                                            "CA0083_paraaortic_lymph_node_12",
                                            "CA0083_perinephric_fat_19",
                                            "CA0083_hilar_lymph_node_47","CA0083_liver_49",
                                            "CA0083_lung_55","CA0090_brain_2",
                                            "CA0090_abdomen_13",
                                            "CA0090_paraaortic_lymph_node_39",
                                            "CA0090_liver_43",
                                            "CA0090_porta_hepatis_lymph_node_52"),
                                patient = c("CA0027","CA0027",
                                            "CA0027","CA0027","CA0034","CA0034",
                                            "CA0034","CA0035","CA0035",
                                            "CA0043","CA0043","CA0043","CA0046",
                                            "CA0046","CA0046","CA0046","CA0046",
                                            "CA0046","CA0058","CA0058",
                                            "CA0058","CA0076","CA0076","CA0076",
                                            "CA0083","CA0083","CA0083",
                                            "CA0083","CA0083","CA0090","CA0090",
                                            "CA0090","CA0090","CA0090"),
                                 sample = c(1L,9L,13L,14L,2L,8L,
                                            11L,1L,2L,4L,7L,12L,5L,7L,
                                            12L,15L,17L,19L,49L,28L,37L,
                                            18L,26L,41L,11L,18L,46L,48L,54L,
                                            3L,14L,40L,44L,53L),
                                   site = c("paraaortic lymph nodes (left upper)",
                                            "prostate - apex left side","dural base skull",
                                            "dural inner skull",
                                            "Right infra-renal paraaortic node","Liver: right lobe A",
                                            "Liver: left lobe A",
                                            "Left para-aortic lymph node track",
                                            "Bladder base 1","lymph node - portal",
                                            "Liver: right lobe nodule 2",
                                            "Liver: left lobe nodule 3",
                                            "Hilar lymph node - left","Lung: left lower lobe nodule",
                                            "Liver: right lobe nodule (2)",
                                            "Left para-aortic node",
                                            "Prostate - posterior","Pelvic lymph node (1)",
                                            "Left Lung - Hilar Lymph Node",
                                            "Liver - Left Lobe Site 2",
                                            "Liver - Right Lobe Site 1",
                                            "left 6th rib soft tumour tissue",
                                            "Anterior surface of L5 vertebra soft tumour",
                                            "Right lobe liver, discrete tumour nodule",
                                            "right para-aortic LN 1",
                                            "right side peri-nephric tissue (?LN)",
                                            "right hilar LN 4",
                                            "left lobe liver Tumour 1","left upper lobe lung",
                                            "Brain Right temporal fossa",
                                            "Left anterior abdomen, subcutaneous tumour",
                                            "Para aortic LN ( right side) 1",
                                            "liver, right lobe tumour 1",
                                            "LN, porta hepatus"),
                             `PSMA.%.0` = c(0L,0L,10L,0L,40L,50L,
                                            30L,5L,5L,98L,10L,0L,15L,
                                            40L,20L,50L,40L,30L,5L,100L,
                                            100L,90L,90L,10L,50L,50L,60L,90L,
                                            5L,100L,100L,100L,100L,100L),
                             `PSMA.%1+` = c(20L,20L,30L,0L,40L,
                                            10L,50L,5L,5L,0L,60L,0L,10L,
                                            40L,40L,10L,40L,50L,30L,0L,0L,
                                            5L,5L,20L,30L,20L,20L,5L,10L,
                                            0L,0L,0L,0L,0L),
                             `PSMA.%2+` = c(50L,20L,30L,50L,20L,
                                            40L,15L,40L,60L,2L,30L,20L,
                                            15L,20L,40L,40L,20L,15L,40L,0L,
                                            0L,5L,5L,70L,15L,20L,10L,5L,
                                            50L,0L,0L,0L,0L,0L),
                             `PSMA.%3+` = c(30L,60L,30L,50L,0L,
                                            0L,5L,50L,30L,0L,0L,80L,60L,
                                            0L,0L,0L,0L,5L,25L,0L,0L,0L,
                                            0L,0L,5L,10L,10L,0L,35L,0L,
                                            0L,0L,0L,0L),
                             PSMA.score = c(0.7,0.8,0.6,
                                            0.834343434,0.266666667,0.3,0.316666667,
                                            0.783434343,0.716666667,0.013434343,
                                            0.4,0.934343434,0.734343434,
                                            0.266666667,0.4,0.3,0.266666667,
                                            0.316666667,0.616666667,0,0,0.05,
                                            0.05,0.534343434,0.25,0.3,
                                            0.234343434,0.05,0.716666667,0,0,0,0,
                                            0),
                      PSMA.Localization = c("cytoplasmic & membranous","cytoplasmic & membranous",
                                            "cytoplasmic & membranous",
                                            "cytoplasmic & membranous","cytoplasmic",
                                            "cytoplasmic & membranous",
                                            "cytoplasmic & membranous",
                                            "cytoplasmic & membranous","cytoplasmic & membranous",
                                            "cytoplasmic, tiny focus of PSMA-reactive tumour, separate to main tissue mass",
                                            "cytoplasmic, ? Fixation artefact","cytoplasmic & membranous",
                                            "cytoplasmic & membranous",
                                            "cytoplasmic","cytoplasmic & membranous",
                                            "cytoplasmic & membranous",
                                            "cytoplasmic & membranous",
                                            "cytoplasmic & membranous","cytoplasmic & membranous",
                                            NA,NA,"cytoplasmic & luminal",
                                            "cytoplasmic",
                                            "cytoplasmic & luminal","cytoplasmic (perinuclear dot)",
                                            "cytoplasmic (perinuclear dot)",
                                            "cytoplasmic","cytoplasmic",
                                            "cytoplasmic & membranous",NA,NA,NA,NA,
                                            NA),
                          AR_percentage = c(10L,5L,10L,10L,75L,
                                            50L,30L,90L,90L,5L,30L,65L,0L,
                                            0L,0L,0L,0L,0L,50L,0L,0L,5L,
                                            0L,5L,30L,5L,1L,0L,90L,0L,
                                            0L,0L,0L,0L),
                           AR_intensity = c(1L,1L,1L,1L,3L,3L,
                                            2L,2L,2L,1L,1L,2L,0L,0L,0L,
                                            0L,0L,0L,2L,0L,0L,1L,0L,1L,
                                            2L,2L,1L,0L,3L,0L,0L,0L,0L,
                                            0L),
                               AR.score = c(0.034343434,0.016666667,
                                            0.034343434,0.034343434,0.75,0.5,
                                            0.2,0.6,0.6,0.016666667,0.1,
                                            0.434343434,0,0,0,0,0,0,
                                            0.343434343,0,0,0.016666667,0,
                                            0.016666667,0.2,0.034343434,0.003434343,
                                            0,0.9,0,0,0,0,0),
                              `INSM1.%` = c(1L,2L,2L,1L,0L,1L,
                                            1L,0L,0L,0L,0L,0L,80L,70L,50L,
                                            80L,30L,50L,1L,5L,10L,0L,0L,
                                            0L,0L,0L,0L,0L,0L,50L,30L,
                                            40L,15L,30L),
                        INSM1.intensity = c(1L,3L,2L,2L,0L,2L,
                                            2L,0L,0L,0L,0L,0L,3L,3L,3L,
                                            3L,1L,2L,1L,1L,2L,0L,0L,0L,
                                            0L,0L,0L,0L,0L,3L,3L,3L,3L,
                                            3L),
                            INSM1.score = c(0.003434343,0.02,
                                            0.013434343,0.006666667,0,0.006666667,
                                            0.006666667,0,0,0,0,0,0.8,0.7,
                                            0.5,0.8,0.1,0.343434343,
                                            0.003434343,0.016666667,0.066666667,0,
                                            0,0,0,0,0,0,0,0.5,0.3,0.4,
                                            0.15,0.3),
                                 `SYP%` = c(1L,0L,0L,0L,10L,1L,
                                            1L,2L,1L,0L,0L,0L,95L,95L,
                                            95L,95L,95L,95L,1L,70L,10L,NA,
                                            NA,0L,0L,0L,0L,0L,1L,80L,95L,
                                            95L,95L,95L),
                                   CNGA = c("1","1","2","0","1",
                                            "1","1","0","0","1","0","0",
                                            "1","1","1","1","0","1",
                                            "2","10","60",NA,NA,"1","0",
                                            "0","0","0","1","50","40",
                                            "35","60","20")
                    )
pathology_review$CNGA <- as.numeric(pathology_review$CNGA)
sites <- pathology_review$srt.obj.sites
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
pathology_review$tissue <- sites

pathology_review$patient

mat1 <- matrix(pathology_review$patient, ncol = 34, byrow = T)
mat1 <- gsub("00", "", mat1)
#colnames(mat1) <- pathology_review$srt.obj.sites
rownames(mat1) <- "patient"
colors1 <- structure(brewer.pal(9, "Set3"), names = unique(mat1[1,]))
h1 <- Heatmap(mat1, name = "patient", col = colors1, cluster_columns = F,
              row_names_side = "left", rect_gp = gpar(col = "white", lwd = 2))

mat2 <- matrix(pathology_review$tissue, ncol = 34, byrow = T)
#colnames(mat1) <- pathology_review$srt.obj.sites
rownames(mat2) <- "sites"
colors2 <- structure(pal_d3()(10), names = unique(pathology_review$tissue))
h2 <- Heatmap(mat2, name = "sites", col = colors2, cluster_columns = F,
              row_names_side = "left", rect_gp = gpar(col = "white", lwd = 2))

mat3 <- matrix(pathology_review$PSMA.score, ncol = 34, byrow = T)
#colnames(mat1) <- pathology_review$srt.obj.sites
rownames(mat3) <- "PSMA"
colors3 <- colorRamp2(c(0,1), c("snow2", "darkgreen"))
h3 <- Heatmap(mat3, name = "PSMA\nscore", col = colors3, cluster_columns = F,
              row_names_side = "left", rect_gp = gpar(col = "white", lwd = 2))

mat4 <- matrix(pathology_review$AR.score, ncol = 34, byrow = T)
#colnames(mat1) <- pathology_review$srt.obj.sites
rownames(mat4) <- "AR"
colors4 <- colorRamp2(c(0,1), c("snow2", "navy"))
h4 <- Heatmap(mat4, name = "AR\nscore", col = colors4, cluster_columns = F,
              row_names_side = "left", rect_gp = gpar(col = "white", lwd = 2))

mat5 <- matrix(pathology_review$INSM1.score, ncol = 34, byrow = T)
#colnames(mat1) <- pathology_review$srt.obj.sites
rownames(mat5) <- "INSM1"
colors5 <- colorRamp2(c(0,1), c("snow2", "orangered3"))
h5 <- Heatmap(mat5, name = "INSM1\nscore", col = colors5, cluster_columns = F,
              row_names_side = "left", rect_gp = gpar(col = "white", lwd = 2))

mat6 <- matrix(pathology_review$`SYP%`, ncol = 34, byrow = T)
#colnames(mat1) <- pathology_review$srt.obj.sites
rownames(mat6) <- "SYP%"
colors6 <- colorRamp2(c(0,100), c("snow2", "darkgreen"))
h6 <- Heatmap(mat6, name = "SYP%", col = colors6, cluster_columns = F,
              row_names_side = "left", rect_gp = gpar(col = "white", lwd = 2))


mat7 <- matrix(pathology_review$CNGA, ncol = 34, byrow = T)
#colnames(mat1) <- pathology_review$srt.obj.sites
rownames(mat7) <- "CNGA%"
colors7 <- colorRamp2(c(0,100), c("snow2", "deeppink"))
h7 <- Heatmap(mat7, name = "CNGA%", col = colors7, cluster_columns = F,
              row_names_side = "left", rect_gp = gpar(col = "white", lwd = 2))
p
 h1 %v% h2 %v% h3 %v% h4 %v% h5 %v% h6 %v% h7 
ggsave("patients overview.pdf", width = 19 , height = 5)
ggsave("patients overview2.pdf", width = 15 , height = 5)
