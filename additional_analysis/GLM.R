
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
counts <- reshape2::dcast(data = df %>% dplyr::select(sample, cell_anno, n), sample~ cell_anno, fill = 0 )
rownames(counts) <- counts$sample
counts <- counts[, -1]
counts <- counts[sample_meta$sample, ]
# Proportions (compositional data)
props <- counts / rowSums(counts)
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

library(speckle)
library(limma)
library(ggplot2)

final_normalsrt <- readRDS("~/CASCADEpaper/paper/normal_cells_202406/subtype/final_normalsrt.Rds")
final_normalsrt$cell_type <- ifelse(final_normalsrt$cell_anno%in% c("T cells","B cells","Plasma cells",  "Macrophages"), "immune cells", 
                                    ifelse(final_normalsrt$cell_anno%in% c("Fibroblasts","Adipocytes","Pericytes", "Endothelial cells"), "stromal cells", "others"))
final_normalsrt$pathology <- ifelse(final_normalsrt$patient %in% c("CA0090", "CA0046"), "NE", ifelse(final_normalsrt$patient %in% c("CA0027", "CA0058"), "Mixed", "AD"))

props <- getTransformedProps(clusters = final_normalsrt$cell_anno,sample = final_normalsrt$sample, transform="logit")
propeller(x = final_normalsrt, clusters = final_normalsrt$cell_anno,sample = final_normalsrt$sample, group = final_normalsrt$site,  transform="logit")
my_table <- props$TransformedProps
propmat <- matrix(my_table, 
                  nrow = nrow(my_table), 
                  ncol = ncol(my_table),
                  dimnames = dimnames(my_table))


patient <- substr(colnames(props$Proportions), 1, 6)
sites <- colnames(props$Proportions)
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
pathology <-  ifelse(patient %in% c("CA0090", "CA0046"), "NE", ifelse(patient %in% c("CA0027", "CA0058"), "Mixed", "AD"))


design <- model.matrix(~patient)
fit0 <- lmFit(props$TransformedProps, design=design)
fit0 <- eBayes(fit0)


design <- model.matrix(~ sites+ patient)
fit1<- lmFit(props$TransformedProps, design=design)
fit1 <- eBayes(fit1)


des.tech <- model.matrix(~sites)
dupcor <- duplicateCorrelation(props$TransformedProps, design=des.tech,
                               block=patient)
dupcor
fit1d <- lmFit(props$TransformedProps, design=des.tech, block=patient, 
               correlation=dupcor$consensus)
fit1d <- eBayes(fit1d)


design <- model.matrix(~patient + sites + pathology )
fit2<- lmFit(props$TransformedProps, design=design)
fit2 <- eBayes(fit2)

# Coefficients not estimable: pathologyMixed pathologyNE 
# Warning: Partial NA coefficients for 14 probe(s)

des.tech2 <- model.matrix(~ sites + pathology)

dupcor <- duplicateCorrelation(props$TransformedProps, design=des.tech2,
                               block=patient)
dupcor
fit2d<- lmFit(props$TransformedProps, design=des.tech2, block=patient, 
              correlation=dupcor$consensus)
fit2d<- eBayes(fit2d)

aic_model <- function(fit, n) {
  k <- ncol(fit$design)          # number of parameters
  rss <- fit$sigma^2 * fit$df.residual
  aic <- n * log(rss / n) + 2 * k
  return(aic)
}

compare_model <- sapply(list(fit0, fit1, fit1d, fit2, fit2d), function(fit){
  data.frame(median_sigma = median(fit$sigma), 
             mean_AIC = mean(aic_model(fit, 33)))
})
colnames(compare_model) <- c("~ patient", "~ patient + sites", "~ sites + block(patient)", "~ patient + sites + pathology", "~ patient + sites + block(pathology)")
compare_model

# best fitting model ~ 0 + patient + sites
write.table(topTable(fit1, coef= 2:9, number = Inf), "limma_celltype_site.csv", sep = "," ) 



## Compare Neuron cells between patients -----------
design <- model.matrix(~ 0+ patient)
cont.matrix <- makeContrasts(
  CA90_vs_CA27 = patientCA0090 - patientCA0034,
  CA90_vs_CA34 = patientCA0090 - patientCA0034,
  CA90_vs_CA35 = patientCA0090 - patientCA0035,
  CA90_vs_CA43 = patientCA0090 - patientCA0043,
  CA90_vs_CA46 = patientCA0090 - patientCA0046,
  CA90_vs_CA58 = patientCA0090 - patientCA0058,
  CA90_vs_CA76 = patientCA0090 - patientCA0076,
  CA90_vs_CA83 = patientCA0090 - patientCA0083,
  levels = design
)
fitp <- lmFit(props$TransformedProps, design=design)
fitp <- contrasts.fit(fitp, cont.matrix)
topTable(eBayes(fitp), number = Inf)
# CA90_vs_CA27 CA90_vs_CA34 CA90_vs_CA35 CA90_vs_CA43 CA90_vs_CA46 CA90_vs_CA58 CA90_vs_CA76 CA90_vs_CA83   AveExpr         F     P.Value  adj.P.Val
# Endothelial cells             -2.5199156   -2.5199156  -1.73908572   -2.6270029   -2.9575086   -3.3209973  -2.00663270  -0.63396677 -2.201612 3.5111344 0.005495571 0.04898542
# Neurons                        3.8987880    3.8987880   3.21590340    4.0314680    3.0307559    3.4643485   4.90614646   3.48011916 -5.938957 3.3719327 0.006997917 0.04898542
# Pericytes                     -1.2738429   -1.2738429  -1.39604541   -1.6497529   -2.2873180   -0.2361942  -0.25524445   0.26768844 -3.399101 2.4955119 0.033267857 0.15525000
# Fibroblasts                   -1.7544596   -1.7544596  -1.61174734   -1.0486089   -0.9350013    0.6009924  -0.18498027  -0.31619621 -1.637859 1.7244610 0.133426276 0.38333635
# T cells                        0.5557080    0.5557080  -0.63264510    0.0359266    1.2142092    1.4797940   0.09761292  -0.97024369 -2.253133 1.7099374 0.136905839 0.38333635
# Lymphatic endothelial cells   -0.4024341   -0.4024341   0.95605414    0.1700341   -0.4240143    1.1055200   1.74331705   1.06129677 -4.985943 1.3837625 0.241258689 0.56293694
# Adipocytes                     0.6044301    0.6044301  -0.50302759    1.7753021    0.6043327    1.4151116   1.30152085   1.02053036 -4.659972 1.0338552 0.424640966 0.77593417
# Macrophages                   -1.0830686   -1.0830686   0.04459526   -0.9538206   -0.9832777   -1.4334289  -0.66701718  -1.26354296 -1.618441 0.9669673 0.469278309 0.77593417
# Chondrocytes                  -0.1896550   -0.1896550  -0.69016457    0.3103482   -0.5096072   -0.5408318  -2.05914297  -0.11270080 -5.814471 0.8729933 0.536762653 0.77593417
# Plasma cells                   1.5574820    1.5574820  -0.08408156    0.0744516    0.4052557    1.1230425   0.22497801  -0.18014676 -5.127399 0.7979835 0.594003624 0.77593417
# EndoMT                        -1.2595707   -1.2595707  -1.66813152   -1.0347518   -1.0056734   -1.4234303  -0.56502062  -0.04034088 -4.907418 0.7779829 0.609662559 0.77593417
# B cells                        1.4131165    1.4131165   0.71019027    1.2623198    0.7446202    2.0150578  -0.26372095   0.41581555 -4.435993 0.7081994 0.665170791 0.77603259
# Hepotocytes                   -1.0465846   -1.0465846   0.08327947   -1.3312275   -0.1018680   -1.4702409  -0.72198448  -1.07968347 -5.482341 0.4032363 0.894139039 0.92294200
# Epithelial cells               0.0611833    0.0611833   1.60696835    0.8842529    0.9814196    0.3800160   0.09944445   0.37464922 -4.842821 0.3536568 0.922941997 0.92294200
# Compare archetype signature expression variance between subclones =====================
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

