
library(lme4)
library(lmerTest)      # p-values
library(emmeans)
library(dplyr)
library(purrr)
library(broom.mixed)
library(ggplot2)

# public signature GLM test - assess organ, pathology effect =============
load("~/CASCADE/Analysis/signatures_analysis/objects/all_gene_sets.rda")
genelist <- lapply(all_gene_sets, function(x) unique(unlist(x))) # union of different signatures

tang_2022 <- all_gene_sets$tang_2022
names(tang_2022) <- paste0("tang_2022_", names(tang_2022))

temp <- lapply(all_gene_sets, function(x){
  names(which(table(unlist(x)) >= length(x)/2))
})

temp <- temp[!names(temp)%in% c("CRPC", "cancer", "hillock", "metastasis", "club", "invasion", "tang_2022")]
genelist <- c(temp, tang_2022)
names(genelist) <- paste0(names(genelist), "_signature")
meta <- readRDS("~/CASCADEpaper/paper/Fig2_signature/signature_meta_tumor_only.Rds")


signature_score <- meta[, names(genelist)[1:20]]
apply(signature_score, 2, var)

signature_score$sample <- meta$sample

sigature_list <- split(signature_score, signature_score$sample)
signature_mat <- lapply(sigature_list, function(x){
  return(as.data.frame(colMeans(x[, 1:20])) %>% `colnames<-`(unique(x[,21])))
})

signature_mat <- as.matrix(Reduce(cbind, signature_mat))
rownames(signature_mat) <- gsub("_signature", "", rownames(signature_mat))

calc_zscore <- function(x){
  (x-mean(x))/sd(x)
}

sites <- colnames(signature_mat)
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


patients <- gsub("00", "", substr(colnames(signature_mat), 1, 6))
classes <- ifelse(patients %in% c("CA27", "CA58"), "Mixed", ifelse(patients %in% c("CA46", "CA90"), "NE", "AD"))

test <- t(apply(signature_mat, 1, calc_zscore))
df <- as.data.frame(t(signature_mat)) %>% mutate(patient = factor(patients), organ = factor(sites), class = factor(classes))
# Vector of signature column names
signatures <- rownames(signature_mat)  

results <- purrr::map_df(signatures, function(sig) {
  formula <- as.formula(paste0(sig, " ~ patient"))
  fit <- lm(formula, data = df)
  formula2 <- as.formula(paste0(sig, " ~ patient + organ"))
  fit2 <- lm(formula2, data = df)
  
  lrt_df <- as.data.frame(anova(fit, fit2))
  lrt_df$model     <- c("~patient", "~patient+organ")
  lrt_df$signature <- sig
  return(lrt_df)
})

results$p_BH <- p.adjust(results$`Pr(>F)`)
write.table(results, "GLME_result_Fig2.siangures.csv", quote = F, row.names = F, sep  =",")


#results added as supplementary table

# normal cell number enrichment test ====================

### Test with distribution and model ========
library(glmmTMB)
final_normalsrt <- readRDS("~/CASCADEpaper/paper/normal_cells_202406/subtype/final_normalsrt.Rds")
final_normalsrt$cell_type <- ifelse(final_normalsrt$cell_anno%in% c("T cells","B cells","Plasma cells",  "Macrophages"), "immune cells", 
                            ifelse(final_normalsrt$cell_anno%in% c("Fibroblasts","Adipocytes","Pericytes", "Endothelial cells"), "stromal cells", "others"))
final_normalsrt$pathology <- ifelse(final_normalsrt$patient %in% c("CA0090", "CA0046"), "NE", ifelse(final_normalsrt$patient %in% c("CA0027", "CA0058"), "Mixed", "AD"))
meta <- final_normalsrt@meta.data
meta$patient <- gsub("00", "", meta$patient)

df <- meta %>% group_by(patient, pathology, site, cell_type, sample, cell_anno) %>% 
  summarise(n = n()) %>% group_by( sample) %>% mutate(total_normal = sum(n))
df_expanded <- df %>% 
  dplyr::select(sample, cell_anno, n) %>%
  ungroup() %>%
  tidyr::complete(
    sample,
    cell_anno,
    fill = list(n = 0)
  )

sample_meta <- df  %>%
  dplyr::select(sample, patient, pathology, site, total_normal) %>%
  distinct()
df_expanded <- df_expanded %>%
  left_join(sample_meta, by = "sample")



celltypes <- unique(meta$cell_anno)

test_model <- lapply(celltypes, function(anno){
  cat(anno, "\n")
  fit1 <- glmmTMB(
    cbind(n, total_normal - n) ~ site + pathology + (1 | patient),
    data = df %>% filter(cell_anno == anno),
    family = betabinomial()
  )
  fit2 <- glmmTMB(
    cbind(n, total_normal - n) ~ pathology + patient +  (1 | site),
    data = df %>% filter(cell_anno == anno),
    family = betabinomial()
  )
  fit3 <- glmmTMB(
    cbind(n, total_normal - n) ~ pathology + (1 | patient/site),
    data = df %>% filter(cell_anno == anno),
    family = betabinomial()
  )
  fit4 <- glmmTMB(
    cbind(n, total_normal - n) ~ pathology,
    data = df %>% filter(cell_anno == anno),
    family = betabinomial()
  )
  fit5 <- glmmTMB(
    cbind(n, total_normal - n) ~ pathology+ (1|patient),
    data = df %>% filter(cell_anno == anno),
    family = betabinomial()
  )
  fit6 <- glmmTMB(
    cbind(n, total_normal - n) ~ pathology+site,
    data = df %>% filter(cell_anno == anno),
    family = betabinomial()
  )
  fit7<- glmmTMB(
    cbind(n, total_normal - n) ~ pathology+site +patient,
    data = df %>% filter(cell_anno == anno),
    family = betabinomial()
  )
  
  fit8 <- glmmTMB(
    cbind(n, total_normal - n) ~ site,
    data = df %>% filter(cell_anno == anno),
    family = betabinomial()
  )

  fit <- glmmTMB(
    cbind(n, total_normal - n) ~  pathology + (1 | site),
    data = df %>% filter(cell_anno == anno),
    family = betabinomial()
  )
  
  fit9 <- glmmTMB(
    cbind(n, total_normal - n) ~ pathology + patient,
    data = df %>% filter(cell_anno == anno),
    family = betabinomial()
  )
  
  fixed_only          = fit4 
  fixed_site          = fit6
  fixed_patient       = fit9
  fixed_both          = fit7
  site_only = fit8
  re_site             = fit
  re_patient          = fit5
  re_patient_site     = fit3
  re_site_fix_patient = fit2
  re_patient_fix_site = fit1
  print(anova(site_only,fixed_only,fixed_site,fixed_patient,fixed_both,re_site,re_patient,re_patient_site,re_site_fix_patient,re_patient_fix_site))
  return(anova(site_only,fixed_only,fixed_site,fixed_patient,fixed_both,re_site,re_patient,re_patient_site,re_site_fix_patient,re_patient_fix_site))
})

### centered log-ratio transform proportion and use distance===========
set.seed(123)
celltypes <- unique(df$cell_anno)

# Metadata
sample_meta <- df  %>%
  dplyr::select(sample, patient, pathology, site, total_normal) %>%
  distinct()


counts <- reshape2::dcast(data = df %>% dplyr::select(sample, cell_anno, n), sample~ cell_anno, fill = 0 )
rownames(counts) <- counts$sample
counts <- counts[, -1]
counts <- counts[sample_meta$sample, ]
# Proportions (compositional data)
props <- counts / rowSums(counts)

#========== PERMANOVA ==========#

# Aitchison distance (recommended for compositional data)
# Replace zeros first with a small pseudocount
counts_nz <- counts + 0.5
clr_mat   <- as.matrix(clr(counts_nz))   # centered log-ratio transform

aitchison_dist <- dist(clr_mat, method = "euclidean")  # Euclidean on CLR = Aitchison

perm_patient_strata <- adonis2(
  aitchison_dist ~ pathology + site ,
  data         = sample_meta,
  strata       = sample_meta$patient,
  permutations = 999,
  by           = "margin"
)

print(perm_patient_strata)
# Permutation test for adonis under reduced model
# Marginal effects of terms
# Blocks:  strata 
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = aitchison_dist ~ pathology + site, data = sample_meta, permutations = 999, by = "margin", strata = sample_meta$patient)
# Df SumOfSqs      R2      F Pr(>F)    
# pathology  2    57.44 0.06891 1.6798  0.044 *  
#   site       8   360.15 0.43207 2.6331  0.001 ***
#   Residual  22   376.14 0.45125                  
# Total     32   833.54 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#========== BETADISPER (Homogeneity of Dispersion) ==========#

bd_disease <- betadisper(aitchison_dist, sample_meta$pathology)
print(permutest(bd_disease, permutations = 999))
# 
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)
# Groups     2  0.902 0.45118 0.1468    999  0.856
# Residuals 30 92.216 3.07386  

bd_organ <- betadisper(aitchison_dist, sample_meta$site)
print(permutest(bd_organ, permutations = 999))
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)  
# Groups     8 43.879  5.4849 2.6415    999  0.018 *
#   Residuals 24 49.834  2.0764                       
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
plot(bd_organ, main = "Dispersion by metastatic site")


cat("\n========== Linear Mixed Models (per cell type) ==========\n")

clr_df <- as.data.frame(clr_mat)
clr_df <- cbind(clr_df, sample_meta)

lmm_results <- lapply(celltypes, function(ct) {
  formula_str <- paste0("`", ct, "` ~  site + (1 | patient)")
  fit <- lmer(as.formula(formula_str), data = clr_df, REML = FALSE)
  
  # Type III ANOVA-style F-tests
  an  <- anova(fit)
  data.frame(
    cell_type     = ct,
    term          = rownames(an),
    F_value       = an$`F value`,
    p_value       = an$`Pr(>F)`
  )
})

lmm_df <- do.call(rbind, lmm_results)

# FDR correction separately for each term
lmm_df <- lmm_df %>%
  group_by(term) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup()

print(lmm_df %>% filter(p_adj_BH < 0.05) %>% arrange(term, p_adj_BH))


sig_celltypes <- lmm_df %>%
  filter(term == "site", p_adj_BH < 0.05) %>%
  pull(cell_type)

# Run emmeans contrasts for each significant cell type
posthoc_results <- lapply(celltypes, function(ct) {
  
  formula_str <- paste0("`", ct, "` ~ site + (1 | patient)")
  fit <- lmer(as.formula(formula_str), data = clr_df, REML = T)
  
  # Estimated marginal means per site
  emm  <- emmeans(fit, ~ site)
  
  # Each site vs grand mean (enrichment/depletion)
  cont <- contrast(emm, method = "eff") %>%
    as.data.frame() %>%
    mutate(
      cell_type = ct,
      direction = ifelse(estimate > 0, "Enriched", "Depleted")
    )
  cont
})

posthoc_df <- do.call(rbind, posthoc_results)

# FDR correction ACROSS all cell types and sites together
posthoc_df$p_adj_BH <- p.adjust(posthoc_df$p.value, method = "BH")

# Significant enrichments
sig_enrichments <- posthoc_df %>%
  filter( p.value < 0.05) %>%
  arrange(cell_type,  p.value)

print(sig_enrichments)
write.table(posthoc_df, "GLM_emm_result_celltype_site.csv", quote = F, row.names = F, sep = ",")

props_long <- as.data.frame(props) %>%
  mutate(sample_id = sample_meta$sample_id,
         pathology = sample_meta$pathology,
         site = sample_meta$site) %>%
  pivot_longer(cols = all_of(celltypes), names_to = "cell_anno", values_to = "proportion")

p_bar_disease <- props_long %>%
  group_by(pathology, cell_anno) %>%
  summarise(mean_prop = mean(proportion), .groups = "drop") %>%
  mutate(labs =signif(mean_prop, 3)) %>% 
  mutate(labs = ifelse(labs>0.05, labs, NA)) %>%
  ggplot(aes(pathology, mean_prop, fill = cell_anno)) +
  geom_col(position = "stack") +
  geom_text(aes(label = labs),
            position = position_stack(vjust = 0.5),        # Position text in the middle of the bars
            color = "white",                               # Set text color
            size = 4)+
  labs(x = "Disease Group", y = "Mean Cell Proportion", fill = "Cell Types") +
  scale_fill_manual(values = normal_cols)+
  theme_bw(base_size = 18) +
  theme(legend.position = "right")
p_bar_disease
ggsave("celltype_bypathology.pdf", width = 8, height = 4)


## Compare Neuron cells between patients -----------
ct = "Neurons"
formula_str <- paste0("`", ct, "` ~ patient + (1 | site)")
fit <- lmer(as.formula(formula_str), data = clr_df, REML = TRUE)

# Estimated marginal means per patient 
emm  <- emmeans(fit, ~ patient)

# Each patient vs grand mean (enrichment/depletion)
cont <- contrast(emm, method = "eff") %>%
  as.data.frame() %>%
  mutate(
    cell_type = ct,
    direction = ifelse(estimate > 0, "Enriched", "Depleted")
  )
cont


ct = "Macrophages"
formula_str <- paste0("`", ct, "` ~ patient + (1 | site)")
fit <- lmer(as.formula(formula_str), data = clr_df, REML = TRUE)

# Estimated marginal means per patient 
emm  <- emmeans(fit, ~ patient)

# Each patient vs grand mean (enrichment/depletion)
cont <- contrast(emm, method = "eff") %>%
  as.data.frame() %>%
  mutate(
    cell_type = ct,
    direction = ifelse(estimate > 0, "Enriched", "Depleted")
  )
cont %>% filter(p.value <0.05)

# Compare archetype signature expression variance between subclones
subclones <- readRDS("~/CASCADEpaper/paper/Fig3_genomics/ATAClone/subclones.Rds")
paths <- system(paste0("realpath ~/CASCADEpaper/paper/Fig5_archetype/4_srt_meta/*.Rds"), intern = T)
names(paths) <- gsub(".Rds", "", sapply(strsplit(paths, split = "/"), '[', 8))


signclonelist <- lapply(names(subclones), function(i){
meta <- readRDS(paths[i])
colnames(meta) <- c("AR","Inflammation", "NE1","NE2", "Cycling","Glycolysis")
subclone <- subclones[[i]]
meta$sample <- i
meta$subclone <- as.character(subclone[rownames(meta)])
meta <- meta[!is.na(meta$subclone), ]
meta$sample_subclone <- paste(meta$sample, meta$subclone)
return(meta)
})

df <- rbindlist(signclonelist)
df$patient <- substr(df$sample, 1, 6)
df$patient_sample <- paste(df$patient, df$sample)
df$patient_sample_subclone <- paste(df$patient, df$sample_subclone)
sites <-df$sample
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
df$site <- sites

reslist <- lapply(c("AR","Inflammation", "NE1","NE2", "Cycling","Glycolysis"), function(sig){
  formula_str1 <- paste0("`", sig, "` ~ (1 | patient) + (1 | patient_sample)")
  # Null model: patient + sample explain variance, but not subclone
  m0 <- lmer(as.formula(formula_str1),
  data = df, REML = FALSE)
  formula_str2 <- paste0("`", sig, "` ~ (1 | patient) + (1 | patient_sample) + (1 | patient_sample_subclone)")
  m1 <- lmer(as.formula(formula_str2),
  data = df, REML = FALSE)
  formula_str3 <- paste0("`", sig, "` ~ (1 | patient) + (1 | patient_sample) + (1 | patient_sample_subclone) + (1 | site)")
  m2 <- lmer(as.formula(formula_str3),
             data = df, REML = FALSE)
  lrt <- anova(m0, m1, m2)
  lrt_df <- as.data.frame(lrt)
  lrt_df$model     <- rownames(lrt_df)
  lrt_df$signature <- sig
  lrt_df$model <- c("NULL", "subclone", "subclone+organ")
  colnames(lrt_df) <- c("npar", "AIC", "BIC", "logLik", "deviance", 
                        "Chisq", "Chi_df", "p_value", "model", "signature")
  lrt_df$is.singular <- c(isSingular(m0), isSingular(m1), isSingular(m2))
  
  # Variance components
  vc <- as.data.frame(VarCorr(m1))
  vc$pct <- vc$vcov / sum(vc$vcov) * 100
  
  lrt_df$var_pct[lrt_df$model =="subclone"] <- vc$pct[vc$grp == "patient_sample_subclone"]
  return(lrt_df)
})
anovadf<- rbindlist(reslist)
anovadf$p_BH <- p.adjust(anovadf$p_value)
write.table(anovadf, "anova_module_subclone.csv", quote = F, row.names = F, sep = ",")

