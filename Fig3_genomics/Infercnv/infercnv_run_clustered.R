library(infercnv)
library(dplyr)
my_path <- list.files("~/CASCADEpaper/paper/Fig3/Infercnv/count_matrix/RNA2/", pattern = ".Rds") %>% gsub(".Rds", "", .)

for (i in my_path[2]){
  infercnvobject = CreateInfercnvObject(raw_counts_matrix=readRDS(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/count_matrix/RNA2/", i, ".Rds")),
                                        annotations_file=readRDS(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/anno/RNA2/", i, ".Rds")),
                                        delim="\t",
                                        gene_order_file="~/subclone/Infercnv/multi_g.txt",
                                        ref_group_names=readRDS("~/CASCADEpaper/paper/Fig3/Infercnv/ref_group.Rds"))
  
  dir.create(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/RNA2/clustered/",i))
  infercnvobject1 = infercnv::run(infercnvobject,
                                  cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                  out_dir= paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/RNA2/clustered/",i),
                                  cluster_by_groups=T,
                                  denoise=F,
                                  HMM=F, 
                                  num_threads = 8)
  rm(list = c("infercnvobject", "infercnvobject1"))
  gc()
}