library(ggplot2)
library(patchwork)
library(ggpubr)
source("~/CASCADEpaper/paper/cols.R")
setwd("~/CASCADEpaper/paper/Fig2_signature")
phenotype_meta <- readRDS("~/CASCADEpaper/paper/PSMA/phenotype_meta.Rds")
phenotype_meta$pathology <- ifelse(phenotype_meta$patien %in% c("CA27", "CA58"), "Mixed", ifelse(phenotype_meta$patient %in% c("CA46", "CA90"), "NE", "AD"))

ar_score <- c(0.033333333, 0.016666667, 0.033333333, 0.033333333, 0.75, 0.5, 0.2, 0.6, 0.6, 0.016666667, 0.1, 0.433333333, 0, 0, 0, 0, 0, 0, 0.333333333, 0, 0, 0.016666667, 0, 0.016666667, 0.2, 0.033333333, 0.003333333, 0, 0.9, 0, 0, 0, 0, 0)
sample <- c("CA0027_paraaortic_lymph_node_1", "CA0027_prostate_9", "CA0027_dura_base_skull_13", "CA0027_dura_inner_skull_14", "CA0034_paraaortic_lymph_node_2", "CA0034_liver_right_8", "CA0034_liver_left_11", "CA0035_paraaortic_lymph_node_1", "CA0035_bladder_2", "CA0043_portal_lymph_node_4", "CA0043_liver_7", "CA0043_liver_12", "CA0046_hilar_lymph_node_5", "CA0046_lung_7", "CA0046_liver_12", "CA0046_paraaortic_lymph_node_15", "CA0046_prostate_17", "CA0046_pelvic_lymph_node_19", "CA0058_hilar_50", "CA0058_liver_29", "CA0058_liver_38", "CA0076_left_rib_17", "CA0076_vertebra_25", "CA0076_liver_right_40", "CA0083_paraaortic_lymph_node_12", "CA0083_perinephric_fat_19", "CA0083_hilar_lymph_node_47", "CA0083_liver_49", "CA0083_lung_55", "CA0090_brain_2", "CA0090_abdomen_13", "CA0090_paraaortic_lymph_node_39", "CA0090_liver_43", "CA0090_porta_hepatis_lymph_node_52")
ar_score <- setNames(ar_score, sample)
ar_score <- sort(ar_score)

phenotype_meta$ar_path <- ifelse(ar_score[phenotype_meta$sample] > 0, "positive", "negative")

compare_means(AR_exp~ar_path, data = phenotype_meta, p.adjust.method = "bonferroni")
g <- ggplot(phenotype_meta, aes(x = ar_path, y = AR_exp)) + 
  geom_jitter(position = position_jitter(width = 0.1, height = 0.1), alpha = 0.3, color = "grey")+
  geom_violin()+
  stat_compare_means( method = "wilcox", label = "p.signif", label.x.npc = 0.5)+
  theme_minimal(base_size = 18)+
  labs(x = "AR H score", y = "AR expression")

phenotype_meta %>% group_by(ar_path) %>% summarise(mean = mean(AR_exp))

mycomp <- list(c("AD", "NE"), c("AD", "Mixed"), c("Mixed", "NE"))
p <- ggplot(phenotype_meta, aes(x = pathology, y = AR_exp)) + 
  geom_jitter(position = position_jitter(width = 0.1, height = 0.1), alpha = 0.3, color = "grey")+
  geom_violin(width = 1.2)+
  stat_compare_means( method = "wilcox", comparisons = mycomp, label = "p.signif")+
  theme_minimal(base_size = 18)+
  labs(x = "Pathology group", y = "AR expression")
phenotype_meta %>% group_by(pathology) %>% summarise(mean = mean(AR_exp))

compare_means(AR_exp~pathology, data = phenotype_meta, p.adjust.method = "bonferroni")
insm1_score <- c(0.003333333, 0.02, 0.013333333, 0.006666667, 0, 0.006666667, 0.006666667, 0, 0, 0, 0, 0, 0.8, 0.7, 0.5, 0.8, 0.1, 0.333333333, 0.003333333, 0.016666667, 0.066666667, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0.3, 0.4, 0.15, 0.3)
insm1_score <- sort(setNames(insm1_score, sample))
phenotype_meta$insm1_path <- ifelse(insm1_score[phenotype_meta$sample] > 0, "positive", "negative")

compare_means(ASCL1_exp~insm1_path, data = phenotype_meta, p.adjust.method = "bonferroni")
g2 <- ggplot(phenotype_meta, aes(x = insm1_path, y = ASCL1_exp)) + 
  geom_jitter(position = position_jitter(width = 0.1, height = 0.1), alpha = 0.3, color = "grey")+
  geom_violin()+
  stat_compare_means( method = "wilcox", label = "p.signif", label.x.npc = 0.5)+
  theme_minimal(base_size = 18)+
  labs(x = "INSM1 H score", y= "ASCL1 expression")
phenotype_meta %>% group_by(insm1_path) %>% summarise(mean = mean(ASCL1_exp))
p2 <- ggplot(phenotype_meta, aes(x = pathology, y = ASCL1_exp)) + 
  geom_jitter(position = position_jitter(width = 0.1, height = 0.1), alpha = 0.3, color = "grey")+
  geom_violin(width = 1.2)+
  stat_compare_means( method = "wilcox", comparisons = mycomp, label = "p.signif")+
  theme_minimal(base_size = 18)+
  labs(x = "Pathology group", y= "ASCL1 expression")
compare_means(ASCL1_exp~pathology, data = phenotype_meta, p.adjust.method = "bonferroni")
phenotype_meta %>% group_by(pathology) %>% summarise(mean = mean(ASCL1_exp))

g+g2 +p +p2 +plot_layout(nrow = 1)+plot_annotation(tag_levels = 'A', tag_prefix = '(', tag_suffix = ')')& 
  theme(plot.tag = element_text(size = 18))
ggsave("pdf/pothology_vs_expression.pdf", width = 18, height = 6)
