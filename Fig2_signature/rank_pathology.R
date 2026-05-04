library(ggplot2)
library(patchwork)
library(ggpubr)
source("~/cols.R")
source("~/CASCADEpaper/paper/cols.R")
phenotype_meta <- readRDS("~/CASCADEpaper/paper/PSMA/phenotype_meta.Rds")
phenotype_meta$pathology <- ifelse(phenotype_meta$patien %in% c("CA27", "CA58"), "Mixed", ifelse(phenotype_meta$patient %in% c("CA46", "CA90"), "NE", "AD"))

data  <- phenotype_meta %>% group_by(pathology, sample) %>% summarise(AR_mean = mean(AR_exp), ASCL1_mean = mean(ASCL1_exp))
ar_score <- c(0.033333333, 0.016666667, 0.033333333, 0.033333333, 0.75, 0.5, 0.2, 0.6, 0.6, 0.016666667, 0.1, 0.433333333, 0, 0, 0, 0, 0, 0, 0.333333333, 0, 0, 0.016666667, 0, 0.016666667, 0.2, 0.033333333, 0.003333333, 0, 0.9, 0, 0, 0, 0, 0)
sample <- c("CA0027_paraaortic_lymph_node_1", "CA0027_prostate_9", "CA0027_dura_base_skull_13", "CA0027_dura_inner_skull_14", "CA0034_paraaortic_lymph_node_2", "CA0034_liver_right_8", "CA0034_liver_left_11", "CA0035_paraaortic_lymph_node_1", "CA0035_bladder_2", "CA0043_portal_lymph_node_4", "CA0043_liver_7", "CA0043_liver_12", "CA0046_hilar_lymph_node_5", "CA0046_lung_7", "CA0046_liver_12", "CA0046_paraaortic_lymph_node_15", "CA0046_prostate_17", "CA0046_pelvic_lymph_node_19", "CA0058_hilar_50", "CA0058_liver_29", "CA0058_liver_38", "CA0076_left_rib_17", "CA0076_vertebra_25", "CA0076_liver_right_40", "CA0083_paraaortic_lymph_node_12", "CA0083_perinephric_fat_19", "CA0083_hilar_lymph_node_47", "CA0083_liver_49", "CA0083_lung_55", "CA0090_brain_2", "CA0090_abdomen_13", "CA0090_paraaortic_lymph_node_39", "CA0090_liver_43", "CA0090_porta_hepatis_lymph_node_52")
ar_score <- setNames(ar_score, sample)
ar_score <- sort(ar_score)
data$ar_score <- ar_score[data$sample]
data$ar_path <- ifelse(data$ar_score, "positive", "negative")
data$ar_path <- paste(data$ar_path, "\n", "n =", table(data$ar_path)[data$ar_path])

insm1_score <- c(0.003333333, 0.02, 0.013333333, 0.006666667, 0, 0.006666667, 0.006666667, 0, 0, 0, 0, 0, 0.8, 0.7, 0.5, 0.8, 0.1, 0.333333333, 0.003333333, 0.016666667, 0.066666667, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0.3, 0.4, 0.15, 0.3)
insm1_score <- sort(setNames(insm1_score, sample))
data$insm1_score <- insm1_score[data$sample]
data$insm1_path <- ifelse(data$insm1_score > 0, "positive", "negative")
data$insm1_path <- paste(data$insm1_path, "\n", "n =", table(data$insm1_path)[data$insm1_path])
data$pathology <- paste(data$pathology, "\n", "n =", table(data$pathology)[data$pathology])

compare_means(AR_mean~ar_path, data = data, p.adjust.method = "bonferroni")
g <- ggplot(data, aes(x = ar_path, y = AR_mean)) + 
  geom_boxplot()+
  geom_jitter(position = position_jitter(width = 0.1, height = 0.1), alpha = 0.3, color = "blue")+
  stat_compare_means( method = "wilcox", label = "p.format", label.x.npc = 0.5)+
  theme_minimal(base_size = 18)+
  labs(x = "AR IHC score group", y = "AR mean expression")
data %>% group_by(ar_path) %>% summarise(AR = mean(AR_mean))

mycomp <- list(c("AD \n n = 16",  "NE \n n = 11" ), c( "AD \n n = 16" ,"Mixed \n n = 7"), c("Mixed \n n = 7",  "NE \n n = 11" ))
p <- ggplot(data, aes(x = pathology, y = AR_mean)) + 
  geom_boxplot()+
  geom_jitter(position = position_jitter(width = 0.1, height = 0.1), alpha = 0.3, color = "blue")+
  stat_compare_means( method = "wilcox", comparisons = mycomp, label = "p.format")+
  theme_minimal(base_size = 18)+
  labs(x = "Pathology group", y = "AR mean expression")
data %>% group_by(pathology) %>% summarise(AR = mean(AR_mean))

compare_means(AR_mean~ar_path, data = data, p.adjust.method = "bonferroni")



compare_means(ASCL1_exp~insm1_path, data = phenotype_meta, p.adjust.method = "bonferroni")
g2 <- ggplot(data, aes(x = insm1_path, y = ASCL1_mean)) + 
  geom_boxplot()+
  geom_jitter(position = position_jitter(width = 0.1, height = 0.1), alpha = 0.3, color = "blue")+
  stat_compare_means( method = "wilcox", label = "p.format", label.x.npc = 0.5)+
  theme_minimal(base_size = 18)+
  labs(x = "INSM1 IHC score group", y= "ASCL1 mean expression")
data %>% group_by(insm1_path) %>% summarise(AR = mean(ASCL1_mean))
p2 <- ggplot(data, aes(x = pathology, y = ASCL1_mean)) + 
  geom_boxplot()+
  geom_jitter(position = position_jitter(width = 0.1, height = 0.1), alpha = 0.3, color = "blue")+
  stat_compare_means( method = "wilcox", comparisons = mycomp, label = "p.format")+
  theme_minimal(base_size = 18)+
  labs(x = "Pathology group", y= "ASCL1 expression")
data %>% group_by(pathology) %>% summarise(AR = mean(ASCL1_mean))

saveRDS(data, "~/CASCADEpaper/sourcedata/figs1.Rds")

compare_means(ASCL1_exp~pathology, data = phenotype_meta, p.adjust.method = "bonferroni")


g+g2 +p +p2 +plot_layout(nrow = 1)+plot_annotation(tag_levels = 'A', tag_prefix = '(', tag_suffix = ')')& 
  theme(plot.tag = element_text(size = 18))
ggsave("pdf/S1pothology_vs_expression.pdf", width = 18, height = 6)

