pkgs <- c(
  "broom.mixed",
  "BSgenome.Hsapiens.UCSC.hg38",
  "CellChat",
  "chromVAR",
  "circlize",
  "ComplexHeatmap",
  "cowplot",
  "data.table",
  "DelayedArray",
  "dorothea",
  "dplyr",
  "emmeans",
  "EnsDb.Hsapiens.v86",
  "future",
  "GenomicRanges",
  "ggplot2",
  "ggpubr",
  "ggrastr",
  "ggridges",
  "glmmTMB",
  "harmony",
  "infercnv",
  "JASPAR2020",
  "JASPAR2024",
  "limma",
  "lme4",
  "lmerTest",
  "Matrix",
  "motifmatchr",
  "pals",
  "patchwork",
  "purrr",
  "qs",
  "ruvIIInb",
  "scran",
  "Seurat",
  "Signac",
  "speckle",
  "TFBSTools",
  "tidyverse"
)

lines <- sapply(pkgs, function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    sprintf("%s v%s", pkg, packageVersion(pkg))
  } else {
    sprintf("%s NOT INSTALLED", pkg)
  }
})

out <- "additional_analysis/package_versions.txt"
writeLines(lines, con = out)
cat("Written to", out, "\n")
