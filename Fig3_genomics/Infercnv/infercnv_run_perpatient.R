library(infercnv)
library(dplyr)
my_path <- list.files("~/CASCADEpaper/paper/Fig3/Infercnv/count_matrix", pattern = ".Rds") %>% gsub(".Rds", "", .)
my_path <-my_path[nchar(my_path)==6]
subclone_anno <- readRDS("~/CASCADEpaper/paper/Fig3/Infercnv/subclone_anno.Rds")

for (i in my_path[4:9]){
  anno <- readRDS(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/anno/", i, ".Rds"))

  infercnvobject = CreateInfercnvObject(raw_counts_matrix=readRDS(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/count_matrix/", i, ".Rds")),
                                        annotations_file=anno,
                                        delim="\t",
                                        gene_order_file="~/subclone/Infercnv/multi_g.txt",
                                        ref_group_names=readRDS("~/CASCADEpaper/paper/Fig3/Infercnv/ref_group.Rds"))
  
  # lapply(sample, function(x) dir.create(paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/subclone/",x)))
  infercnvobject1 = infercnv::run(infercnvobject,
                                  cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                  out_dir= paste0("~/CASCADEpaper/paper/Fig3/Infercnv/infercnv_run/perpatient_subclone/",i),
                                  cluster_by_groups=T,
                                  denoise=F,
                                  HMM=T, 
                                  num_threads = 8)
  rm(list = c("infercnvobject", "infercnvobject1"))
  gc()
}