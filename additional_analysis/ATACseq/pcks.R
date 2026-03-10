## CRAN packages ----
cran_pkgs <- c(
  "Seurat", "Matrix", "ggplot2", "cowplot", "tidyverse",
  "dplyr", "patchwork", "future"
)

install.packages(cran_pkgs)

## Bioconductor packages ----
bioc_pkgs <- c(
  "chromVAR",
  "motifmatchr",
  "scran",
  "Signac",
  "EnsDb.Hsapiens.v86",
  "BSgenome.Hsapiens.UCSC.hg38",
  "JASPAR2020",
  "TFBSTools",
  "dorothea",
  "ComplexHeatmap"
)

BiocManager::install(bioc_pkgs, ask = FALSE, update = TRUE)

## Extra step: dorothea (sometimes requires this explicit install)
# BiocManager::install("saezlab/dorothea")  # only if CRAN/Bioc version fails
module load conda/24.7.1-2
module load python/3.9.18-gcc-13.2.0
# 1. Create a new environment (optional but recommended)
conda create -n macs3_env python=3.9

# 2. Activate the environment
conda activate macs3_env

# 3. Install MACS3 from the bioconda channel
conda install -c bioconda -c conda-forge macs3