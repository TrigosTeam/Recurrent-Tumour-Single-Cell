library(pals)
library(scales)
library(RColorBrewer)
meta <- readRDS("~/CASCADEpaper/paper/Fig2_signature/signature_meta_tumor_only.Rds")
# sites <- colnames(signature_mat)
# sites[grep("brain", sites)] <- "brain"
# sites[grep("dura", sites)] <- "dura"
# sites[grep("prostate", sites)] <- "prostate"
# sites[grep("liver", sites)] <- "liver"
# sites[grep("fat", sites)] <- "fat"
# sites[grep("lymph|hilar", sites)] <- "LN"
# sites[grep("rib|vertebra", sites)] <- "bone"
# sites[grep("lung", sites)] <- "lung"
# sites[grep("abdomen", sites)] <- "abdomen"
# sites[grep("bladder", sites)] <- "bladder"
site_cols <- setNames(brewer.set3(10), unique(meta$site))
patient_cols <- setNames(hue_pal()(9), gsub("00", "", unique(meta$patient)))
patient_cols2 <- setNames(hue_pal()(9), unique(meta$patient))

final_normal_anno <- readRDS("~/CASCADEpaper/paper/normal_cells_202406/subtype/final_normal_anno.Rds")
normal_cols <- setNames(brewer.set1(14), levels(final_normal_anno))
module_cols <- setNames( c("dodgerblue3", brewer.dark2(6)[2:6]),c("Module1", "Module2", "Module3", "Module4", "Module5", "Module6"))

module_cols2 <- setNames( c("dodgerblue3", brewer.dark2(6)[2:6]),c("AR","Inflammation", "NE1","NE2", "Cycling","Glycolysis"  ))
library(colorspace)
pal_names <- ls("package:pals")
discrete_colors <- lapply(pal_names, function(p) {
  try(get(p, envir = asNamespace("pals"))(8), silent = T)  # Get 8 colors from each palette
})

randcol <- function(n){
  colors <- unlist(discrete_colors, use.names = F)
  colors <- unique(colors[nchar(colors) == 7])
  set.seed(123)
  sample(colors, size = n)
}

