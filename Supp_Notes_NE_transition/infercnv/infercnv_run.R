library(infercnv)
library(dplyr)
ref_group <- c("Adipocytes","Fibroblasts","Macrophages","EndoMT", "Endothelial cells","Pericytes" ,"T cells" ,"Epithelial cells",           
  "Lymphatic endothelial cells","Hepotocytes","B cells","Plasma cells","Chondrocytes","Neurons" )

infercnvobject = CreateInfercnvObject(raw_counts_matrix=readRDS("~/Supp_Notes_NE_transition/infercnv/p13/CA90_p13_counts.Rds"),
                                      annotations_file=readRDS("~/Supp_Notes_NE_transition/infercnv/p13/CA90_p13_anno.Rds"),
                                      delim="\t",
                                      gene_order_file="~/subclone/Infercnv/multi_g.txt",
                                      ref_group_names=ref_group)


infercnvobject1 = infercnv::run(infercnvobject,
                                cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                out_dir= "~/Supp_Notes_NE_transition/infercnv/p13/",
                                cluster_by_groups=T,
                                denoise=T,
                                HMM=F, 
                                num_threads = 8)

infercnvobject = CreateInfercnvObject(raw_counts_matrix=readRDS("~/Supp_Notes_NE_transition/infercnv/p18/CA90_p18_counts.Rds"),
                                      annotations_file=readRDS("~/Supp_Notes_NE_transition/infercnv/p18/CA90_p18_anno.Rds"),
                                      delim="\t",
                                      gene_order_file="~/subclone/Infercnv/multi_g.txt",
                                      ref_group_names=ref_group)


infercnvobject1 = infercnv::run(infercnvobject,
                                cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                out_dir= "~/Supp_Notes_NE_transition/infercnv/p18/",
                                cluster_by_groups=T,
                                denoise=T,
                                HMM=F, 
                                num_threads = 8)