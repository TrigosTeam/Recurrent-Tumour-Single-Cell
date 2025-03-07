library(infercnv)
library(dplyr)
my_path <- list.files("~/CASCADEpaper/paper/Fig3/Infercnv/count_matrix/RNA2/", pattern = ".Rds") %>% gsub(".Rds", "", .)

for (i in my_path[2]){
  infercnvobject = CreateInfercnvObject(raw_counts_matrix=readRDS(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/count_matrix/RNA2/", i, ".Rds")),
                                        annotations_file=readRDS(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/anno/RNA2/", i, ".Rds")),
                                        delim="\t",
                                        gene_order_file="~/subclone/Infercnv/multi_g.txt",
                                        ref_group_names=readRDS("~/CASCADEpaper/paper/Fig3/Infercnv/ref_group.Rds"))
  
  dir.create(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/RNA2/unclustered/",i))
  infercnvobject1 = infercnv::run(infercnvobject,
                                  cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                  out_dir= paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/RNA2/unclustered/",i),
                                  cluster_by_groups=F,
                                  denoise=F,
                                  HMM=F, 
                                  num_threads = 8)
  rm(list = c("infercnvobject", "infercnvobject1"))
  gc()
  # 
  # dir.create(paste0("~/subclone/Infercnv/unclustered_run/",i))
  # infercnvobject2 = infercnv::run(infercnvobject,
  #                                 cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
  #                                 out_dir= paste0("~/subclone/Infercnv/unclustered_run/",i),
  #                                 cluster_by_groups=F,
  #                                 denoise=FALSE,
  #                                 HMM=F)
  # rm(list = c("infercnvobject","infercnvobject2" ))
  # gc()
  
  # if (file.exists(paste0("~/subclone/Infercnv/comb_annotation/","subclone_",i, "_annotation.Rds"))){
  # infercnvobject3 = CreateInfercnvObject(raw_counts_matrix=readRDS(paste0("~/subclone/Infercnv/comb_matrix/",i, "_matrix.Rds")),
  #                                       annotations_file=readRDS(paste0("~/subclone/Infercnv/comb_annotation/","subclone_",i, "_annotation.Rds")),
  #                                       delim="\t",
  #                                       gene_order_file="~/subclone/Infercnv/multi_g.txt",
  #                                       ref_group_names=readRDS(paste0("~/subclone/Infercnv/comb_annotation/",substr(i, 1, 6), "_n_anno.Rds")))
  # 
  # dir.create(paste0("~/subclone/Infercnv/subclone_HMM_run/",i))
  # infercnvobject4 = infercnv::run(infercnvobject3,
  #                                 cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
  #                                 out_dir= paste0("~/subclone/Infercnv/subclone_HMM_run/",i),
  #                                 cluster_by_groups=T,
  #                                 denoise=FALSE,
  #                                 HMM=T)
  # }
}