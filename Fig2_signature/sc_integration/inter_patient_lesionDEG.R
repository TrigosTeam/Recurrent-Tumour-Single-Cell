library(Seurat)
library(MAST)
library(tidyverse)

srt <- readRDS("~/Fig2_signature/sc_integration/tumor_only_ruvsrt.Rds")
Idents(srt) <- srt$patient

# groups <- combn(unique(srt$patient), 2, simplify = F)
# markerlist <- lapply(groups, function(x){
#   marker <- FindMarkers(srt, ident.1 = x[1], ident.2 = x[2], test.use = "MAST")
#   print(x[1])
#   print(x[2])
#   return(marker)
#   })
# names(markerlist) <- sapply(groups, function(x) paste(x, collapse = "-"))


### do 1 vs the rest
markers <- FindAllMarkers(srt, test.use = "MAST")
saveRDS(markers, "~/Fig2_signature/sc_integration/DEG/inter_patient_DEG.Rds")
rm(list = ls())
gc()

paths <- system("realpath ~/Fig2_signature/sc_integration/perpatient/tumor_only/*_ruvsrt.Rds", intern = T)
paths_split <- substr(unlist(lapply(strsplit(paths, split = "/"), function(x) x[8])),1 ,6)
paths <- setNames(paths, paths_split)

### intra-patient , inter-lesion DEG
for(i in paths_split){
  srt <- readRDS(paths[i])
  Idents(srt) <- srt$sample
  # groups <- combn(unique(srt$sample), 2, simplify = F)
  # markerlist <- lapply(groups, function(x){
  #   marker <- FindMarkers(srt, ident.1 = x[1], ident.2 = x[2], test.use = "MAST")
  #   print(x[1])
  #   print(x[2])
  #   return(marker)
  # })
  # names(markerlist) <- sapply(groups, function(x) paste(x, collapse = "-"))
  markers <- FindAllMarkers(srt, test.use = "MAST")
  saveRDS(markers, paste0("~/Fig2_signature/sc_integration/DEG/inter_lesion/", i, ".Rds"))
  print(i)
}

meta <- readRDS("~/Fig4-5_archetype/5_integrated_final_module_meta.Rds")
srt <- AddMetaData(srt, meta)
Idents(srt) <- srt$site
markers <- FindAllMarkers(srt, test.use = "MAST")
saveRDS(markers, "~/Fig2_signature/sc_integration/DEG/inter_organ_DEG.Rds")
