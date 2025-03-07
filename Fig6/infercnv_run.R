library(Seurat)
library(infercnv)
library(dplyr)

infercnvobject = CreateInfercnvObject(raw_counts_matrix= readRDS("~/CASCADEpaper/paper/Fig2/normal_cells/normal_inferncv_clean_countsv1.Rds"),
                                      annotations_file=readRDS("~/CASCADEpaper/paper/Fig2/normal_cells/normal_inferncv_clean_annov1.Rds"),
                                      delim="\t",
                                      gene_order_file="~/subclone/Infercnv/multi_g.txt",
                                      min_max_counts_per_cell = c(100, +Inf),
                                      ref_group_names= c("ref"))

infercnvobject1 = infercnv::run(infercnvobject,
                                cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                out_dir= "~/CASCADEpaper/paper/Fig2/normal_cells/infercnv/clustered",
                                cluster_by_groups=T,
                                denoise=T,
                                HMM=F, num_threads = 6)

infercnvobject1 = infercnv::run(infercnvobject,
                                cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                out_dir= "~/CASCADEpaper/paper/Fig2/normal_cells/infercnv/unclustered",
                                cluster_by_groups=F,
                                denoise=T,
                                HMM=F, num_threads = 6)