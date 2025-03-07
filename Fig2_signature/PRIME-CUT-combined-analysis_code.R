#####PRIME-CUT combined analysis
source('/Users/aleksandar/Documents/Documents/MED SCHOOL/summer 2019/OncoTreat_ForAleks/Functions_WeightedOncoTreat.R')

#Function to draw heatmap of genes in seurat object dat, with cluster labels clust
#seurat object must have labels "tissue" and "patient" 
geneHeatmap_primecut=function(dat,clust,genes,genes_by_cluster=T,n_top_genes_per_cluster=5,viper=F,color_palette=NA,scaled=F,fontsize=8){
  identities <- levels(clust)
  if(is.na(color_palette)){my_color_palette <- hue_pal()(length(identities))}
  else{my_color_palette=color_palette}
  features=genes
  i=sample(1:ncol(dat),min(10000,ncol(dat)),replace = F)
  if(viper==F){x=dat[["integrated"]]@data[features,i]}
  if(viper==T){x=dat[["RNA"]]@data[features,i]}
  df <- data.frame(clust[i],dat$treatment[i])
  rownames(df)=colnames(x)
  colnames(df)=c("cluster","treatment")
  anno_colors <- list(cluster = my_color_palette,treatment=c("cornflowerblue","darkgoldenrod1","darkgoldenrod4","coral3"))
  names(anno_colors$cluster) <- levels(df$cluster)
  names(anno_colors$treatment)<-c("Pre-Treatment","On-ADT","On-Combination","Recurrence")
  o=order(df$cluster,df$treatment)
  x=x[,o]
  df=df[o,]
  quantile_breaks <- function(xs, n = 10) {
    breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
    breaks[!duplicated(breaks)]
  }
  if(scaled==F){t=as.matrix(apply(x,1,function(x){(x-mean(x))/sd(x)}))}
  if(scaled==T){t=as.matrix(x)}
  mat_breaks <- c(quantile_breaks(t[which(t<0)], n = 10),0,quantile_breaks(t[which(t>0)], n = 10))
  mat_breaks=mat_breaks[2:(length(mat_breaks)-1)] #restrict range of data to quantiles 5%-95%, extreme values excluded
  if(genes_by_cluster){
    anno_colors$group=anno_colors$cluster
    anno_row=data.frame(group=unlist(lapply(unique(df$cluster),function(x){rep(x,n_top_genes_per_cluster)})))
    gene_names=rownames(x)
    rownames(x)=1:nrow(x)
    if(!scaled){pheatmap(x, cluster_rows=FALSE,show_rownames=T,cluster_cols=FALSE, annotation_row = anno_row,annotation_col=df,breaks=mat_breaks,color = colorRampPalette(colors = c('blue', 'white', 'red'))(length(mat_breaks)),fontsize_row = fontsize,show_colnames = F,annotation_colors = anno_colors,scale="row",gaps_row=(2:length(unique(clust))-1)*n_top_genes_per_cluster,annotation_names_row = F,labels_row=gene_names,row_annotation_legend=F)}
    if(scaled){pheatmap(x, cluster_rows=FALSE,show_rownames=T,cluster_cols=FALSE, annotation_row = anno_row,annotation_col=df,breaks=mat_breaks,color = colorRampPalette(colors = c('blue', 'white', 'red'))(length(mat_breaks)),fontsize_row = fontsize,show_colnames = F,annotation_colors = anno_colors,gaps_row=(2:length(unique(clust))-1)*n_top_genes_per_cluster,annotation_names_row = F,labels_row=gene_names,row_annotation_legend=F)}
  }
  else{
    if(!scaled){pheatmap(x, cluster_rows=FALSE,show_rownames=T,cluster_cols=FALSE, annotation_col=df,breaks=mat_breaks,color = colorRampPalette(colors = c('blue', 'white', 'red'))(length(mat_breaks)),fontsize_row = fontsize,show_colnames = F,annotation_colors = anno_colors,scale="row")}
    if(scaled){pheatmap(x, cluster_rows=FALSE,show_rownames=T,cluster_cols=FALSE, annotation_col=df,breaks=mat_breaks,color = colorRampPalette(colors = c('blue', 'white', 'red'))(length(mat_breaks)),fontsize_row = fontsize,show_colnames = F,annotation_colors = anno_colors)}
  }
}


####Patient1 pre-tx vs on-tx (bone met. androgen deprivation + anti-pd1)
P1_BL<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/190603_CHARLES_NIVI_9_HUMAN_10X/CM2/outs/filtered_feature_bc_matrix")
P1_BL<-CreateSeuratObject(counts = P1_BL,min.features = 500,min.cells=50)
P1_BL$patient="Patient1"
P1_BL$treatment="Pre-Treatment"
P1_BL$tissue="Bone"
P1_W10<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/190825_VICTOR_VICTOR_1_HUMAN_10X/S1863-PRIMECUT_P001_W10_bone-bx_GEX/outs/filtered_feature_bc_matrix")
P1_W10<-CreateSeuratObject(counts = P1_W10,min.features = 500,min.cells=50)
P1_W10$patient="Patient1"
P1_W10$treatment="On-Combination"
P1_W10$tissue="Bone"
P1<-merge(x=P1_BL,y=P1_W10,project="PRIME-CUT")

####Patient3 pre-tx vs on-tx (bone met. androgen deprivation only)
P3_BL<-Read10X("/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/191118_VICTOR_YA_16_HUMAN_10X/S1863_PRIMECUT_P003_BL_bone-bx_GEX/filtered_feature_bc_matrix")
P3_BL<-CreateSeuratObject(counts = P3_BL,min.features = 500,min.cells=50)
P3_BL$patient="Patient3"
P3_BL$treatment="Pre-Treatment"
P3_BL$tissue="Bone"
P3_W4<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/191118_VICTOR_YA_16_HUMAN_10X/S1863_PRIMECUT_P003_W4_bone-bx_GEX/filtered_feature_bc_matrix")
P3_W4<-CreateSeuratObject(counts = P3_W4,min.features = 500,min.cells=50)
P3_W4$patient="Patient3"
P3_W4$treatment="On-ADT"
P3_W4$tissue="Bone"
P3_R<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/victor_nygc_samples/SJ15-S1863-PRIMECUT-P003-End-of-treatment-Bone-bx-GEX/filtered_feature_bc_matrix")
P3_R<-CreateSeuratObject(counts = P3_R,min.features = 500,min.cells=50)
P3_R$patient="Patient3"
P3_R$treatment="Recurrence"
P3_R$tissue="Bone"
P3<-merge(x=P3_BL,y=c(P3_W4,P3_R),project="PRIME-CUT")

####Patient8 pre-tx vs on-tx (bone met. androgen deprivation only)
P8_BL<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/191118_VICTOR_YA_16_HUMAN_10X/PRIMECUT_P008_BL_lumbar-spine/filtered_feature_bc_matrix")
P8_BL<-CreateSeuratObject(counts = P8_BL,min.features = 500,min.cells=50)
P8_BL$patient="Patient8"
P8_BL$treatment="Pre-Treatment"
P8_BL$tissue="Bone"
P8_W4<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/191118_VICTOR_YA_16_HUMAN_10X/PRIMECUT_P008_W4_lumbar-spine/filtered_feature_bc_matrix")
P8_W4<-CreateSeuratObject(counts = P8_W4,min.features = 500,min.cells=50)
P8_W4$patient="Patient8"
P8_W4$treatment="On-ADT"
P8_W4$tissue="Bone"
P8<-merge(x=P8_BL,y=P8_W4,project="PRIME-CUT")

####Patient5 pre-tx vs on-tx (lymph node met. androgen deprivation only)
P5_BL<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/191118_VICTOR_YA_16_HUMAN_10X/PRIMECUT_P005_BL_neck-LN-bx/filtered_feature_bc_matrix")
P5_BL<-CreateSeuratObject(counts = P5_BL,min.features = 500,min.cells=50)
P5_BL$patient="Patient5"
P5_BL$treatment="Pre-Treatment"
P5_BL$tissue="LN"
P5_W4<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/191118_VICTOR_YA_16_HUMAN_10X/PRIMECUT_P005_W4_LN-bx/filtered_feature_bc_matrix")
P5_W4<-CreateSeuratObject(counts = P5_W4,min.features = 500,min.cells=50)
P5_W4$patient="Patient5"
P5_W4$treatment="On-ADT"
P5_W4$tissue="LN"
P5<-merge(x=P5_BL,y=P5_W4,project="PRIME-CUT")

#Patient 7 pretx vs on-tx (lung metastasis, adt+anti-pd1)
#/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/191118_VICTOR_YA_16_HUMAN_10X/PRIMECUT_P007_BL_lung-bx/filtered_feature_bc_matrix
#/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/191118_VICTOR_YA_16_HUMAN_10X/PRIMECUT_P007_W10_lung-bx/filtered_feature_bc_matrix
P7_BL<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/191118_VICTOR_YA_16_HUMAN_10X/PRIMECUT_P007_BL_lung-bx/filtered_feature_bc_matrix")
P7_BL<-CreateSeuratObject(counts = P7_BL,min.features = 500,min.cells=50)
P7_BL$patient="Patient7"
P7_BL$treatment="Pre-Treatment"
P7_BL$tissue="Lung"
P7_W10<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/191118_VICTOR_YA_16_HUMAN_10X/PRIMECUT_P007_W10_lung-bx/filtered_feature_bc_matrix")
P7_W10<-CreateSeuratObject(counts = P7_W10,min.features =500,min.cells=50)
P7_W10$patient="Patient7"
P7_W10$treatment="On-Combination"
P7_W10$tissue="Lung"
P7<-merge(x=P7_BL,y=P7_W10,project="PRIME-CUT")

#Patient 6 liver met-- on-treatment only: W9 adt+pd1
#/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/191118_VICTOR_YA_16_HUMAN_10X/PRIMECUT_Pt6_W9_liver/filtered_feature_bc_matrix
P6_W10<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/191118_VICTOR_YA_16_HUMAN_10X/PRIMECUT_Pt6_W9_liver/filtered_feature_bc_matrix")
P6_W10<-CreateSeuratObject(counts = P6_W10,min.features = 500,min.cells=50)
P6_W10$patient="Patient6"
P6_W10$treatment="On-Combination"
P6_W10$tissue="Liver"
P6<-P6_W10

#Patient 10 lymph node met ADT+PD1
P10_BL<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/victor_nygc_samples/S1863-PRIMECUT-P010-BL-bone-bx-GEX/outs/filtered_feature_bc_matrix/")
P10_BL<-CreateSeuratObject(counts = P10_BL,min.features = 500,min.cells=50)
P10_BL$patient="Patient10"
P10_BL$treatment="Pre-Treatment"
P10_BL$tissue="LN"
P10_W10<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/victor_nygc_samples/S1863-PRIMECUT-P010-W10-bone-bx-GEX/outs/filtered_feature_bc_matrix/")
P10_W10<-CreateSeuratObject(counts = P10_W10,min.features = 500,min.cells=50)
P10_W10$patient="Patient10"
P10_W10$treatment="On-Combination"
P10_W10$tissue="LN"
P10<-merge(x=P10_BL,y=P10_W10,project="PRIME-CUT")

#Patient 11 lymph node met ADT only
P11_BL<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/victor_nygc_samples/SJ14-S1863-PRIMECUT-P011-BL-LN-bx-GEX/filtered_feature_bc_matrix")
P11_BL<-CreateSeuratObject(counts = P11_BL,min.features = 500,min.cells=50)
P11_BL$patient="Patient11"
P11_BL$treatment="Pre-Treatment"
P11_BL$tissue="LN"
P11_W4<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/victor_nygc_samples/SJ16-S1863-PRIMECUT-P011-W4-LN-bx-GEX/filtered_feature_bc_matrix")
P11_W4<-CreateSeuratObject(counts = P11_W4,min.features = 500,min.cells=50)
P11_W4$patient="Patient11"
P11_W4$treatment="On-ADT"
P11_W4$tissue="LN"
P11<-merge(x=P11_BL,y=P11_W4,project="PRIME-CUT")

###NEW SAMPLES
#/Users/aleksandar/genomecenter.columbia.edu/primecut_S1863-012-W4 Bone met /201118_ANDREA_ADINA_1_HUMAN_10X-CUAC2278SC-cellranger-count-default/CUAC2278SC_cellranger_count_outs/filtered_feature_bc_matrix
#/Users/aleksandar/genomecenter.columbia.edu/primecut-S1863-013-B Retroperitoneal biopsy./201214_ANDREA_ADINA_1_HUMAN_10X-CUAC2299SC-cellranger-count-default/CUAC2299SC_cellranger_count_outs/filtered_feature_bc_matrix
#/Users/aleksandar/genomecenter.columbia.edu/primecut_S1863-014 Liver met biopsy pretreatmen/CUAC2401SC_cellranger_count_outs/filtered_feature_bc_matrix

P12_W4<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/primecut_S1863-012-W4 Bone met /201118_ANDREA_ADINA_1_HUMAN_10X-CUAC2278SC-cellranger-count-default/CUAC2278SC_cellranger_count_outs/filtered_feature_bc_matrix/")
P12_W4<-CreateSeuratObject(counts = P12_W4,min.features = 500,min.cells=50)
P12_W4$patient="Patient12"
P12_W4$treatment="On-ADT"
P12_W4$tissue="Bone"
P12=P12_W4

P13_BL<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/primecut-S1863-013-B Retroperitoneal biopsy./201214_ANDREA_ADINA_1_HUMAN_10X-CUAC2299SC-cellranger-count-default/CUAC2299SC_cellranger_count_outs/filtered_feature_bc_matrix/")
P13_BL<-CreateSeuratObject(counts = P13_BL,min.features = 500,min.cells=50)
P13_BL$patient="Patient13"
P13_BL$treatment="Pre-Treatment"
P13_BL$tissue="Liver"
P13=P13_BL

P14_BL<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/primecut_S1863-014 Liver met biopsy pretreatmen/CUAC2401SC_cellranger_count_outs/filtered_feature_bc_matrix/")
P14_BL<-CreateSeuratObject(counts = P14_BL,min.features = 500,min.cells=50)
P14_BL$patient="Patient14"
P14_BL$treatment="Pre-Treatment"
P14_BL$tissue="Liver"
P14_W4<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/primecut_014_livermet_w4/CUAC2422SC_cellranger_count_outs/filtered_feature_bc_matrix")
P14_W4<-CreateSeuratObject(counts = P14_W4,min.features = 500,min.cells=50)
P14_W4$patient="Patient14"
P14_W4$treatment="On-ADT"
P14_W4$tissue="Liver"
P14=merge(x=P14_BL,y=P14_W4,project="PRIME-CUT")

P16_BL<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/primecut_016_bonemet_bl/CUAC2423SC_cellranger_count_outs/filtered_feature_bc_matrix")
P16_BL<-CreateSeuratObject(counts = P16_BL,min.features = 500,min.cells=50)
P16_BL$patient="Patient16"
P16_BL$treatment="Pre-Treatment"
P16_BL$tissue="Bone"
P16=P16_BL
####


###NEW DATA
#P17 Baseline, site: ABDOMEN/RETRO.
#/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/220510_ANDREA_VICTOR_1_HUMAN_10X/CUAC2517SC/analysis/220510_ANDREA_VICTOR_1_HUMAN_10X-CUAC2517SC-cellranger-count-default/CUAC2517SC_cellranger_count_outs/filtered_feature_bc_matrix

P17_BL<-Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/220510_ANDREA_VICTOR_1_HUMAN_10X/CUAC2517SC/analysis/220510_ANDREA_VICTOR_1_HUMAN_10X-CUAC2517SC-cellranger-count-default/CUAC2517SC_cellranger_count_outs/filtered_feature_bc_matrix")
P17_BL<-CreateSeuratObject(counts = P17_BL,min.features = 500,min.cells=50)
P17_BL$patient="Patient17"
P17_BL$treatment="Pre-Treatment"
P17_BL$tissue="Abdomen"
P17=P17_BL

#primecut_list=list(P1_BL,P1_W10,P3_BL,P3_W4,P3_R,P5_BL,P5_W4,P6_W10,P7_BL,P7_W10,P8_BL,P8_W4,P10_BL,P11_BL,P11_W4)
primecut_list=list(P1,P3,P5,P6,P7,P8,P10,P11,P12,P13,P14,P16)
for(i in 1:length(primecut_list)){
  p=primecut_list[[i]]
  p <- PercentageFeatureSet(p, pattern = "^MT-", col.name = "percent.mt")
  p <- subset(p, subset = percent.mt < 10 & nCount_RNA > 1000 & nCount_RNA < 15000)
  p <- SCTransform(p,vars.to.regress = c("nCount_RNA","percent.mt"),return.only.var.genes = F,verbose = T,conserve.memory = T)
  p.singler=CreateSinglerObject(p[["SCT"]]@counts, annot = NULL, 
                                              project.name = "primecut", min.genes = 0, 
                                              technology = "10X", species = "Human", citation = "", 
                                              do.signatures = F, clusters = NULL, numCores = numCores,
                                              fine.tune=F,temp.dir = "/Users/aleksandar/Downloads",
                                              variable.genes = "de",reduce.file.size = T,do.main.types = T)
  p$hpca_labels=p.singler$singler[[1]][[1]][[2]]
  p$hpca_main_labels=p.singler$singler[[1]][[4]][[2]]
  p$blueprint_labels=p.singler$singler[[2]][[1]][[2]]
  p$blueprint_main_labels=p.singler$singler[[2]][[4]][[2]]
  p$hpca_pvals=p.singler$singler[[1]][[1]][[3]]
  p$hpca_main_pvals=p.singler$singler[[1]][[4]][[3]]
  p$blueprint_pvals=p.singler$singler[[2]][[1]][[3]]
  p$blueprint_main_pvals=p.singler$singler[[2]][[4]][[3]]
  primecut_list[[i]]=p
}

features <- SelectIntegrationFeatures(object.list = primecut_list, nfeatures = 2000)
primecut_list <- PrepSCTIntegration(object.list = primecut_list, anchor.features = features, verbose = T)
anchors <- FindIntegrationAnchors(object.list = primecut_list, normalization.method = "SCT", anchor.features = features, verbose = T,reference = 2)
rm(primecut_list,features)
primecut.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = T)
rm(anchors)
primecut.integrated$treatment=factor(primecut.integrated$treatment, levels = c("Pre-Treatment","On-ADT","On-Combination","Recurrence"))

primecut.integrated <- RunPCA(primecut.integrated, features = VariableFeatures(object = primecut.integrated))
primecut.integrated <- RunUMAP(primecut.integrated, dims = 1:50, verbose = FALSE,umap.method="umap-learn",metric="correlation")
primecut.integrated <- FindNeighbors(primecut.integrated, dims = 1:50, verbose = FALSE)
primecut.integrated <- FindClusters(primecut.integrated, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) 
clust=primecut.integrated@meta.data[,which(grepl("integrated_snn_res.",colnames(primecut.integrated@meta.data)))]
mat=as.data.frame(t(primecut.integrated$pca@cell.embeddings))
out=sil_subsample(mat,clust)
means=out[[1]]
sd=out[[2]]
x=seq(0.01,1,by=0.01)
errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
lines(x,means)
best=tail(x[which(means==max(means))],n=1)
legend("topright",paste("Best",best,sep = " = "))
primecut.integrated$seurat_clusters=primecut.integrated@meta.data[,which(colnames(primecut.integrated@meta.data)==paste("integrated_snn_res.",best,sep=""))]
Idents(primecut.integrated) <- "seurat_clusters"
# iterative louvain sub-clustering
primecut.integrated$original_seurat_clusters=primecut.integrated$seurat_clusters
cluster_output=as.character(primecut.integrated$seurat_clusters)
for(i in unique(as.character(primecut.integrated$seurat_clusters))){
  print(i)
  s <- subset(primecut.integrated, subset = seurat_clusters == i)
  if(ncol(s)/ncol(primecut.integrated)>0.01){
    s <- RunPCA(s, features = VariableFeatures(object = s))
    s <- RunUMAP(s, dims = 1:50, verbose = FALSE,umap.method="umap-learn",metric="correlation")
    s <- FindNeighbors(s, dims = 1:50, verbose = FALSE)
    s <- FindClusters(s, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) 
    clust_s=s@meta.data[,which(grepl("integrated_snn_res.",colnames(s@meta.data)))]
    clust_s=clust_s[,10:ncol(clust_s)]
    mat_s=as.data.frame(t(s$pca@cell.embeddings))
    out=sil_subsample(mat_s,clust_s)
    means=out[[1]]
    sd=out[[2]]
    x=seq(0.1,1,by=0.01)
    best=tail(x[which(means==max(means))],n=1)
    s$seurat_clusters=s@meta.data[,which(colnames(s@meta.data)==paste("integrated_snn_res.",best,sep=""))]
    if(max(means)>0.1 & min(table(s$seurat_clusters))/sum(table(s$seurat_clusters))>0.05){cluster_output[which(cluster_output==i)]=paste(i,s$seurat_clusters)}
  }
}
primecut.integrated$seurat_clusters=as.factor(as.numeric(as.factor(cluster_output)))
Idents(primecut.integrated) <- "seurat_clusters"
#end iterative louvain sub-clustering
primecut.integrated$seurat_clusters=mapvalues(primecut.integrated$seurat_clusters, from = 1:length(unique(primecut.integrated$seurat_clusters)), to = c("CD4 T-cell 1","CD4 T-cell 2","CD4 T-cell 3","CD8 T-cell 1", "B-cell 1", "Monocyte.Macrophage","Epithelial","Treg 1","Treg 2","CD8 T-cell 2","Fibroblast","Endothelial","Plasma Cell 1","Plasma Cell 2","B-cell 2"))
Idents(primecut.integrated) <- "seurat_clusters"
plot(DimPlot(primecut.integrated,reduction="umap",group.by="seurat_clusters",split.by="treatment") + NoLegend())
plot(DimPlot(primecut.integrated,reduction="umap",group.by="seurat_clusters",split.by="tissue") + NoLegend())
plot(DimPlot(primecut.integrated,reduction="umap",group.by="seurat_clusters",split.by="patient",ncol=3))
plot(DimPlot(primecut.integrated, reduction = "umap",label = TRUE,label.size=7,repel=T) + NoLegend())
markers <- FindAllMarkers(primecut.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "MAST")
top10 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
geneHeatmap_primecut(primecut.integrated,primecut.integrated$seurat_clusters,top10$gene)
l=primecut.integrated$blueprint_labels
l[which(primecut.integrated$blueprint_pvals>0.1)]=NA
l[which(l %in% names(which(table(l)<150)))]=NA
primecut.integrated$l=l
Idents(primecut.integrated) <- "l"
plot(DimPlot(primecut.integrated, reduction = "umap",label=TRUE,repel=T,label.size=5)+NoLegend())
saveRDS(primecut.integrated, file = "~/Documents/Documents/MED SCHOOL/summer 2019/primecut.integrated.rds")
VlnPlot(primecut.integrated, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3,pt.size=0,group.by="patient")
plot(DimPlot(primecut.integrated[,which(primecut.integrated$patient=="Patient1")],reduction="umap",group.by="seurat_clusters",split.by="treatment") + NoLegend())
plot(DimPlot(primecut.integrated[,which(primecut.integrated$patient=="Patient3")],reduction="umap",group.by="seurat_clusters",split.by="treatment") + NoLegend())
plot(DimPlot(primecut.integrated[,which(primecut.integrated$patient=="Patient5")],reduction="umap",group.by="seurat_clusters",split.by="treatment") + NoLegend())
plot(DimPlot(primecut.integrated[,which(primecut.integrated$patient=="Patient6")],reduction="umap",group.by="seurat_clusters",split.by="treatment") + NoLegend())
plot(DimPlot(primecut.integrated[,which(primecut.integrated$patient=="Patient7")],reduction="umap",group.by="seurat_clusters",split.by="treatment") + NoLegend())
plot(DimPlot(primecut.integrated[,which(primecut.integrated$patient=="Patient8")],reduction="umap",group.by="seurat_clusters",split.by="treatment") + NoLegend())
plot(DimPlot(primecut.integrated[,which(primecut.integrated$patient=="Patient10")],reduction="umap",group.by="seurat_clusters",split.by="treatment") + NoLegend())
plot(DimPlot(primecut.integrated[,which(primecut.integrated$patient=="Patient11")],reduction="umap",group.by="seurat_clusters",split.by="treatment") + NoLegend())
##############
table(primecut.integrated$patient,primecut.integrated$seurat_clusters,primecut.integrated$treatment)
pre=table(primecut.integrated$patient,primecut.integrated$seurat_clusters,primecut.integrated$treatment)[,c(1,3,5,8,9,11),1]+1
adt=table(primecut.integrated$patient,primecut.integrated$seurat_clusters,primecut.integrated$treatment)[,c(1,3,5,8,9,11),2]+1
comb=table(primecut.integrated$patient,primecut.integrated$seurat_clusters,primecut.integrated$treatment)[,c(1,3,5,8,9,11),3]+1
recurr=table(primecut.integrated$patient,primecut.integrated$seurat_clusters,primecut.integrated$treatment)[,c(1,3,5,8,9,11),4]+1
pre=apply(pre,1,function(x){x/sum(x)})
adt=apply(adt,1,function(x){x/sum(x)})
comb=apply(comb,1,function(x){x/sum(x)})
recurr=apply(recurr,1,function(x){x/sum(x)})
adt.pre=adt[,c("Patient11","Patient3","Patient5","Patient8")]/pre[,c("Patient11","Patient3","Patient5","Patient8")]
comb.pre=comb[,c("Patient1","Patient10","Patient7")]/pre[,c("Patient1","Patient10","Patient7")]
recurr.pre=recurr[,c("Patient3")]/pre[,c("Patient3")]
combined=cbind(adt.pre,comb.pre,recurr.pre)
colnames(combined)=c("Patient11","Patient3","Patient5","Patient8","Patient1","Patient10","Patient7","Patient3_Recurrence")
y=melt(combined)
colnames(y)=c("cluster","patient","logFC")
y$patient=factor(y$patient, levels = c("Patient11","Patient3","Patient5","Patient8","Patient1","Patient10","Patient7","Patient3_Recurrence"))
y$cluster=as.factor(y$cluster)
y$logFC=log2(y$logFC)
p=ggplot(y, aes(x=cluster, y=logFC,fill=patient)) +
  #geom_boxplot() +
  geom_col(position="dodge")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Post/Pre Cluster Frequency")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+scale_fill_manual(values=c("cyan3","dodgerblue3","darkturquoise","deepskyblue2","darkgoldenrod2","gold1","lightgoldenrod1","firebrick"))


##make metaCells
primecut.integrated.meta=MakeCMfA(dat.mat=as.matrix(primecut.integrated[["SCT"]]@counts),clustering=primecut.integrated$seurat_clusters,out.dir="/Users/aleksandar/Documents/Documents/MED SCHOOL/summer 2019/prime-cut-metacells/",out.name="primecut_integrated",sizeThresh = 50)
meta=primecut.integrated.meta[[1]]
for(i in 2:length(primecut.integrated.meta)){
  meta=merge(meta,primecut.integrated.meta[[i]],by=0,all=T)
  rownames(meta)=meta[,1]
  meta=meta[,2:ncol(meta)]
}
meta[is.na(meta)] <- 0
primecut.integrated.meta<-CreateSeuratObject(counts = meta)
primecut.integrated.meta$patient=primecut.integrated$patient[colnames(primecut.integrated.meta)]
primecut.integrated.meta$treatment=primecut.integrated$treatment[colnames(primecut.integrated.meta)]
primecut.integrated.meta$tissue=primecut.integrated$tissue[colnames(primecut.integrated.meta)]
primecut.integrated.meta$blueprint_labels=primecut.integrated$blueprint_labels[colnames(primecut.integrated.meta)]
primecut.integrated.meta$blueprint_pvals=primecut.integrated$blueprint_pvals[colnames(primecut.integrated.meta)]
primecut.integrated.meta$l=primecut.integrated$l[colnames(primecut.integrated.meta)]
primecut_list_meta=SplitObject(primecut.integrated.meta,split.by="patient")
for(i in 1:length(primecut_list_meta)){
  p=primecut_list_meta[[i]]
  p <- PercentageFeatureSet(p, pattern = "^MT-", col.name = "percent.mt")
  p <- PercentageFeatureSet(p, pattern = "^RP", col.name = "percent.ribosomal")
  #p <- subset(p, subset = percent.mt < 25 & nCount_RNA > 1000 & nCount_RNA < 25000)
  p <- SCTransform(p,vars.to.regress = c("nCount_RNA","percent.mt"),return.only.var.genes = F,verbose = T,conserve.memory = T)
  primecut_list_meta[[i]]=p
}
features <- SelectIntegrationFeatures(object.list = primecut_list_meta, nfeatures = 2000)
primecut_list_meta <- PrepSCTIntegration(object.list = primecut_list_meta, anchor.features = features, verbose = T)
anchors <- FindIntegrationAnchors(object.list = primecut_list_meta, normalization.method = "SCT", anchor.features = features, verbose = T,reference = 2)
rm(primecut_list_meta,features)
primecut.integrated.meta <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = T)
rm(anchors)
saveRDS(primecut.integrated.meta,file = "~/Documents/Documents/MED SCHOOL/summer 2019/primecut.integrated.meta.rds")



##run VIPER
##load aracne networks
primecut.integrated.meta=readRDS("~/Documents/Documents/MED SCHOOL/summer 2019/primecut.integrated.meta.rds")
filenames <- list.files("/Users/aleksandar/Documents/Documents/MED SCHOOL/summer 2019/prime-cut-metacells/sc_nets/", pattern="*.rds", full.names=TRUE)
nets=lapply(filenames,readRDS)
#filenames <- list.files("/Users/aleksandar/Documents/Documents/MED SCHOOL/summer 2019/neoredp_aracnenets", pattern="*.rds", full.names=TRUE)
#nets2=lapply(filenames,readRDS)
#nets=c(nets,nets2)
dat=primecut.integrated.meta@assays$integrated@scale.data
#dat=as.data.frame(primecut.integrated.meta@assays$SCT@counts)
#g=rownames(primecut.integrated.meta@assays$SCT@counts)
rm(primecut.integrated.meta)
#dat=log10(apply(dat,2,function(x){(x+1)/(sum(x)+length(x))*1000000}))
#means=apply(dat,1,mean)
#sds=apply(dat,1,sd)
#dat=(dat-means)/sds
#dat=t(apply(dat,1,function(x){(x-mean(x))/sd(x)}))
#rownames(dat)=g
#dat=dat[complete.cases(dat),]
vp_meta_1<-viper(dat[,1:4000], nets, method = 'none')
dat=dat[,4001:ncol(dat)]
vp_meta_2<-viper(dat[,1:4000], nets, method = 'none')
dat=dat[,4001:ncol(dat)]
vp_meta_3<-viper(dat[,1:4000], nets, method = 'none')
dat=dat[,4001:ncol(dat)]
vp_meta_4<-viper(dat[,1:4000], nets, method = 'none')
dat=dat[,4001:ncol(dat)]
vp_meta_5<-viper(dat[,1:4000], nets, method = 'none')
dat=dat[,4001:ncol(dat)]
vp_meta_6<-viper(dat[,1:4000], nets, method = 'none')
dat=dat[,4001:ncol(dat)]
vp_meta_7<-viper(dat[,1:4000], nets, method = 'none')
dat=dat[,4001:ncol(dat)]
vp_meta_8<-viper(dat[,1:4000], nets, method = 'none')
dat=dat[,4001:ncol(dat)]
vp_meta_9<-viper(dat[,1:4000], nets, method = 'none')
dat=dat[,4001:ncol(dat)]
vp_meta_10<-viper(dat[,1:4000], nets, method = 'none')
dat=dat[,4001:ncol(dat)]
vp_meta_11<-viper(dat, nets, method = 'none')
vp_meta=cbind(vp_meta_1,vp_meta_2,vp_meta_3,vp_meta_4,vp_meta_5,vp_meta_6,vp_meta_7,vp_meta_8,vp_meta_9,vp_meta_10,vp_meta_11)
rm(vp_meta_1,vp_meta_2,vp_meta_3,vp_meta_4,vp_meta_5,vp_meta_6,vp_meta_7,vp_meta_8,vp_meta_9,vp_meta_10)
saveRDS(vp_meta, "~/Documents/Documents/MED SCHOOL/summer 2019/primecut.integrated.meta.vp.rds")

##VIPER clustering
cbcMRs <- CBCMRs(vp_meta)
primecut.integrated.meta.vp <- CreateSeuratObject(counts = vp_meta[cbcMRs,])
primecut.integrated.meta.vp@assays$RNA@scale.data=as.matrix(primecut.integrated.meta.vp@assays$RNA@data)
primecut.integrated.meta.vp <- RunPCA(primecut.integrated.meta.vp,features=rownames(primecut.integrated.meta.vp))
primecut.integrated.meta.vp <- RunUMAP(primecut.integrated.meta.vp, dims = 1:50, verbose = FALSE,umap.method="umap-learn",metric="correlation")
primecut.integrated.meta.vp <- FindNeighbors(primecut.integrated.meta.vp, dims = 1:50, verbose = FALSE)
primecut.integrated.meta.vp <- FindClusters(primecut.integrated.meta.vp, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) 
clust=primecut.integrated.meta.vp@meta.data[,which(grepl("RNA_snn_res.",colnames(primecut.integrated.meta.vp@meta.data)))]
mat=primecut.integrated.meta.vp@assays$RNA@scale.data
out=sil_subsample(mat,clust)
means=out[[1]]
sd=out[[2]]
x=seq(0.01,1,by=0.01)
errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
lines(x,means)
best=tail(x[which(means[5:length(means)]==max(means[5:length(means)]))+4],n=1)
legend("topright",paste("Best",best,sep = " = "))
primecut.integrated.meta.vp$seurat_clusters=primecut.integrated.meta.vp@meta.data[,which(colnames(primecut.integrated.meta.vp@meta.data)==paste("RNA_snn_res.",best,sep=""))]
Idents(primecut.integrated.meta.vp) <- "seurat_clusters"
primecut.integrated.meta.vp$patient=primecut.integrated$patient[colnames(primecut.integrated.meta.vp)]
primecut.integrated.meta.vp$treatment=primecut.integrated$treatment[colnames(primecut.integrated.meta.vp)]
primecut.integrated.meta.vp$tissue=primecut.integrated$tissue[colnames(primecut.integrated.meta.vp)]
primecut.integrated.meta.vp$blueprint_labels=primecut.integrated$blueprint_labels[colnames(primecut.integrated.meta.vp)]
primecut.integrated.meta.vp$blueprint_pvals=primecut.integrated$blueprint_pvals[colnames(primecut.integrated.meta.vp)]
primecut.integrated.meta.vp$l=primecut.integrated$l[colnames(primecut.integrated.meta.vp)]
primecut.integrated.meta.vp$treatment=ordered(primecut.integrated.meta.vp$treatment,levels=c("Pre-Treatment","On-ADT","On-Combination","Recurrence"))
# iterative louvain sub-clustering
primecut.integrated.meta.vp$original_seurat_clusters=primecut.integrated.meta.vp$seurat_clusters
cluster_output=as.character(primecut.integrated.meta.vp$seurat_clusters)
for(i in unique(as.character(primecut.integrated.meta.vp$seurat_clusters))){
  print(i)
  s <- subset(primecut.integrated.meta.vp, subset = seurat_clusters == i)
  if(ncol(s)/ncol(primecut.integrated.meta.vp)>0.05){
    s <- RunPCA(s, features = rownames(s))
    s <- RunUMAP(s, dims = 1:50, verbose = FALSE,umap.method="umap-learn",metric="correlation")
    s <- FindNeighbors(s, dims = 1:50, verbose = FALSE)
    s <- FindClusters(s, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) 
    clust=s@meta.data[,which(grepl("RNA_snn_res.",colnames(s@meta.data)))]
    mat=s@assays$RNA@scale.data
    out=sil_subsample(mat,clust)
    means=out[[1]]
    sd=out[[2]]
    x=seq(0.01,1,by=0.01)
    best=tail(x[5:length(x)][which(means[5:length(means)]==max(means[5:length(means)]))],n=1)
    #best=tail(x[which(means==max(means))],n=1)
    s$seurat_clusters=s@meta.data[,which(colnames(s@meta.data)==paste("RNA_snn_res.",best,sep=""))]
    if(max(means)>0.15 & min(table(s$seurat_clusters))/sum(table(s$seurat_clusters))>0.01){cluster_output[which(cluster_output==i)]=paste(i,s$seurat_clusters)}
  }
}
primecut.integrated.meta.vp$seurat_clusters=as.factor(as.numeric(as.factor(cluster_output)))
Idents(primecut.integrated.meta.vp) <- "seurat_clusters"
#end iterative louvain sub-clustering
#primecut.integrated.meta.vp$seurat_clusters=mapvalues(primecut.integrated.meta.vp$seurat_clusters, from = 1:length(unique(primecut.integrated.meta.vp$seurat_clusters)), to = c("Treg","CD4 T-cell","B-cells 1","B-cells 2","Neutrophils","MEP","CD8 T-cell 1","CD8 T-cell 2","Epithelial 1","Epithelial 2","Monocyte 1","Monocyte 2","Erythrocytes","Fibroblast","Macrophages","Plasma cells","Memory B-cells"))
primecut.integrated.meta.vp$seurat_clusters=mapvalues(primecut.integrated.meta.vp$seurat_clusters, from = 1:length(unique(primecut.integrated.meta.vp$seurat_clusters)), to = c("CD4 T-cell 1","CD4 T-cell 2","CD8 T-cell 1","CD8 T-cell 2","CD8 T-cell 3","B-cells 1","B-cells 2","B-cells 3","Monocytes 1","Monocytes 2","Macrophages","Monocytes 3","Treg 1","Neutrophils","Epithelial 1","Epithelial 2","Erythrocytes 1","Treg 2","B-cells 4","Epithelial 3","Treg 3","Endothelial","Plasma cells","B-cells 5"))
Idents(primecut.integrated.meta.vp) <- "seurat_clusters"
plot(DimPlot(primecut.integrated.meta.vp,reduction="umap",group.by="seurat_clusters",split.by="treatment") + NoLegend())
plot(DimPlot(primecut.integrated.meta.vp,reduction="umap",group.by="seurat_clusters",split.by="tissue") + NoLegend())
plot(DimPlot(primecut.integrated.meta.vp,reduction="umap",group.by="seurat_clusters",split.by="patient",ncol=3))
plot(DimPlot(primecut.integrated.meta.vp, reduction = "umap",label = TRUE,label.size=7,repel=T) + NoLegend())
markers.vp <- FindAllMarkers(primecut.integrated.meta.vp, only.pos = TRUE, min.pct = 0, logfc.threshold = 0,test.use = "t")
top10 <- markers.vp %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
geneHeatmap_primecut(primecut.integrated.meta.vp,primecut.integrated.meta.vp$seurat_clusters,top10$gene,viper = T,scaled = T,n_top_genes_per_cluster = 3)
l=primecut.integrated.meta.vp$blueprint_labels
l[which(primecut.integrated.meta.vp$blueprint_pvals>0.1)]=NA
l[which(l %in% names(which(table(l)<150)))]=NA
primecut.integrated.meta.vp$l=l
Idents(primecut.integrated.meta.vp) <- "l"
plot(DimPlot(primecut.integrated.meta.vp, reduction = "umap",label=TRUE,repel=T,label.size=5)+NoLegend())
saveRDS(primecut.integrated.meta.vp, file = "~/Documents/Documents/MED SCHOOL/summer 2019/primecut.integrated.meta.vp.rds")
plot(DimPlot(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$patient=="Patient1")],reduction="umap",group.by="seurat_clusters",split.by="treatment") + NoLegend())
plot(DimPlot(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$patient=="Patient3")],reduction="umap",group.by="seurat_clusters",split.by="treatment") + NoLegend())
plot(DimPlot(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$patient=="Patient5")],reduction="umap",group.by="seurat_clusters",split.by="treatment") + NoLegend())
plot(DimPlot(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$patient=="Patient6")],reduction="umap",group.by="seurat_clusters",split.by="treatment") + NoLegend())
plot(DimPlot(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$patient=="Patient7")],reduction="umap",group.by="seurat_clusters",split.by="treatment") + NoLegend())
plot(DimPlot(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$patient=="Patient8")],reduction="umap",group.by="seurat_clusters",split.by="treatment") + NoLegend())
plot(DimPlot(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$patient=="Patient10")],reduction="umap",group.by="seurat_clusters",split.by="treatment") + NoLegend())
plot(DimPlot(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$patient=="Patient11")],reduction="umap",group.by="seurat_clusters",split.by="treatment") + NoLegend())
plot(DimPlot(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$patient=="Patient12")],reduction="umap",group.by="seurat_clusters",split.by="treatment") + NoLegend())
plot(DimPlot(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$patient=="Patient13")],reduction="umap",group.by="seurat_clusters",split.by="treatment") + NoLegend())
plot(DimPlot(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$patient=="Patient14")],reduction="umap",group.by="seurat_clusters",split.by="treatment") + NoLegend())
######
x=primecut.integrated
x=x[,colnames(primecut.integrated.meta.vp)]
x$seurat_clusters=primecut.integrated.meta.vp$seurat_clusters
Idents(x)="seurat_clusters"
markers <- FindAllMarkers(x, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "MAST")
top10 <- markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
geneHeatmap_primecut(x,tumorcells$seurat_clusters,top10$gene,viper = F,scaled = F,n_top_genes_per_cluster = 3)
primecut.integrated.meta.vp$KLK3_geneexp=log10(x@assays$SCT@counts["KLK3",]+1)
VlnPlot(primecut.integrated.meta.vp,"KLK3_geneexp",pt.size=0)+NoLegend()

##InferCNV
write.table(x@assays$RNA@counts,file="~/Downloads/raw_counts_matrix.txt",sep="\t",quote=F)
write.table(x$seurat_clusters,file="~/Downloads/annotations.txt",sep="\t",quote=F,col.names = F)
infercnv_prostate = CreateInfercnvObject(raw_counts_matrix= '~/Downloads/raw_counts_matrix.txt',
                                         annotations_file='~/Downloads/annotations.txt',
                                         gene_order_file='/Users/aleksandar/Documents/Documents/MED SCHOOL/summer 2019/chromosome_locations_noDupGenename_armloc.txt',
                                         delim = '\t',
                                         ref_group_names=c("CD4 T-cell 1","CD4 T-cell 2","CD8 T-cell 1","CD8 T-cell 2","CD8 T-cell 3","B-cells 1","B-cells 2","B-cells 3","Monocytes 1","Monocytes 2","Macrophages","Monocytes 3","Treg 1","Neutrophils","Treg 2","B-cells 4","Treg 3","Plasma cells","B-cells 5", "Erythrocytes 1"))
infercnv_prostate = infercnv::run(infercnv_prostate,
                                  cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                  out_dir='/Users/aleksandar/Documents/Documents/MED SCHOOL/summer 2019/primecut_integrated_infercnv', 
                                  cluster_by_groups=TRUE, 
                                  denoise=TRUE,
                                  noise_logistic = TRUE, 
                                  scale_data = FALSE,
                                  HMM=F, 
                                  up_to_step = 50, 
                                  k_obs_groups = 1, 
                                  HMM_type = 'i6',
                                  analysis_mode = 'subclusters', 
                                  no_plot = F, 
                                  plot_steps = T,
                                  tumor_subcluster_partition_method = 'qnorm', 
                                  num_threads = 11, 
                                  debug = F)


####PLOTTING COMMANDS
primecut.integrated.meta.vp=readRDS("~/Documents/Documents/MED SCHOOL/summer 2019/primecut.integrated.meta.vp.rds")
primecut.integrated=readRDS("~/Documents/Documents/MED SCHOOL/summer 2019/primecut.integrated.rds")
tumorcells.vp=readRDS("~/Documents/Documents/MED SCHOOL/summer 2019/primecut.integrated.tumorcells.vp.rds")

rgb2col = function(rgbmat){
  ProcessColumn = function(col){
    rgb(rgbmat[1, col], 
        rgbmat[2, col], 
        rgbmat[3, col], 
        maxColorValue = 255)
  }
  sapply(1:ncol(rgbmat), ProcessColumn)
}
# Colors to darken
ColorsHex = c("#D55E00","#009E73")
# Convert to rgb
ColorsRGB = col2rgb(ColorsHex)
# Darken colors by lowering values of RGB
ColorsRGBDark = round(ColorsRGB*.7)
# Convert back to hex
ColorsHexDark = rgb2col(ColorsRGBDark)

vipercolors=c("#e5f5e0","#c7e9c0","#a1d99b","#74c476","#41ab5d","#c6dbef","#9ecae1","#6baed6","#fa9fb5","#f768a1","#ae017e","#dd3497","#238b45","#7a0177","#fd8d3c","#f03b20","#49006a","#006d2c","#4292c6","#bd0026","#00441b","#fecc5c","#08519c","#2171b5")
genecolors=c("#e5f5e0","#c7e9c0","#c7e9c0","#a1d99b","#c6dbef","#fa9fb5","#fd8d3c","#238b45","#006d2c","#74c476","#49006a","#fecc5c","#08519c","#08519c","#9ecae1")
vipercolors=rgb2col(round(col2rgb(vipercolors)*0.85))
genecolors=rgb2col(round(col2rgb(genecolors)*0.85))

plot(DimPlot(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")], reduction = "umap",label = TRUE,label.size=7,repel=T,cols=vipercolors) + NoLegend())
plot(DimPlot(primecut.integrated[,which(primecut.integrated$treatment=="Pre-Treatment")], reduction = "umap",label = TRUE,label.size=7,repel=T,cols=genecolors) + NoLegend())

plot(DimPlot(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")], reduction = "umap",label = TRUE,label.size=5,repel=T,split.by="tissue",ncol=2,cols=vipercolors) + NoLegend())
plot(DimPlot(primecut.integrated[,which(primecut.integrated$treatment=="Pre-Treatment")], reduction = "umap",label = TRUE,label.size=5,repel=T,split.by="tissue",ncol=2,cols=genecolors) + NoLegend())

primecut.integrated.meta.vp$cell_lineages=as.character(primecut.integrated.meta.vp$seurat_clusters)
primecut.integrated.meta.vp$cell_lineages[which(primecut.integrated.meta.vp$seurat_clusters %in% c("CD4 T-cell 1","CD4 T-cell 2"))]="CD4nonTreg"
primecut.integrated.meta.vp$cell_lineages[which(primecut.integrated.meta.vp$seurat_clusters %in% c("CD8 T-cell 1","CD8 T-cell 2","CD8 T-cell 3"))]="CD8"
primecut.integrated.meta.vp$cell_lineages[which(primecut.integrated.meta.vp$seurat_clusters %in% c("Treg 1", "Treg 2", "Treg 3"))]="Treg"
primecut.integrated.meta.vp$cell_lineages[which(primecut.integrated.meta.vp$seurat_clusters %in% c("B-cells 1","B-cells 2","B-cells 3","B-cells 4","B-cells 5","Plasma cells"))]="B-cell"
primecut.integrated.meta.vp$cell_lineages[which(primecut.integrated.meta.vp$seurat_clusters %in% c("Monocytes 1","Monocytes 2","Monocytes 3","Macrophages","Neutrophils"))]="Myeloid"
primecut.integrated.meta.vp$cell_lineages[which(primecut.integrated.meta.vp$seurat_clusters %in% c("Epithelial 1","Epithelial 2","Epithelial 3"))]="Tumor"
primecut.integrated.meta.vp$cell_lineages[which(primecut.integrated.meta.vp$seurat_clusters %in% c("Erythrocytes 1"))]="Erythrocyte"
table(primecut.integrated.meta.vp$cell_lineages)
primecut.integrated.meta.vp$cell_lineages_toplevel="Immune"
primecut.integrated.meta.vp$cell_lineages_toplevel[which(primecut.integrated.meta.vp$cell_lineages %in% c("Tumor","Endothelial"))]="Non-Immune"
lineage_colors=vipercolors[c(6,1,3,22,17,11,13,15)]
plot(DimPlot(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")], group.by="cell_lineages",reduction = "umap",label = TRUE,label.size=7,repel=T,cols = lineage_colors) + NoLegend())
plot(DimPlot(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")], group.by="cell_lineages", reduction = "umap",label = TRUE,label.size=5,repel=T,split.by="tissue",ncol=2,cols=lineage_colors) + NoLegend())

markers.vp <- FindAllMarkers(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")], only.pos = TRUE, min.pct = 0, logfc.threshold = 0,test.use = "t")
top10 <- markers.vp %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
geneHeatmap_primecut(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")],primecut.integrated.meta.vp$seurat_clusters[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")],top10$gene,viper = T,scaled = T,n_top_genes_per_cluster = 5,color_palette = vipercolors)

markers <- FindAllMarkers(primecut.integrated[,which(primecut.integrated$treatment=="Pre-Treatment")], only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "wilcox")
top10 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
geneHeatmap_primecut(primecut.integrated[,which(primecut.integrated$treatment=="Pre-Treatment")],primecut.integrated$seurat_clusters[which(primecut.integrated$treatment=="Pre-Treatment")],top10$gene,viper = F,scaled = F,n_top_genes_per_cluster = 5,color_palette = genecolors)

primecut.integrated.meta.vp$tissue_treatment<-paste(primecut.integrated.meta.vp$tissue, primecut.integrated.meta.vp$treatment, sep = "_")
primecut.integrated$tissue_treatment<-paste(primecut.integrated$tissue, primecut.integrated$treatment, sep = "_")
primecut.integrated.meta.vp$tissue_treatment=factor(primecut.integrated.meta.vp$tissue_treatment,levels=c("Bone_Pre-Treatment","Bone_On-ADT","Bone_On-Combination","Liver_Pre-Treatment","Liver_On-ADT","Liver_On-Combination","LN_Pre-Treatment","LN_On-ADT","LN_On-Combination","Lung_Pre-Treatment","Lung_On-ADT","Lung_On-Combination","Bone_Recurrence"))
#primecut.integrated.meta.vp$tissue_treatment=mapvalues(primecut.integrated.meta.vp$tissue_treatment, from = levels(primecut.integrated.meta.vp$tissue_treatment), to = c("Bone_Pre-Treatment","Bone_On-ADT","Bone_On-Combination","Viscera_Pre-Treatment","Viscera_On-ADT","Viscera_On-Combination","LN_Pre-Treatment","LN_On-ADT","LN_On-Combination","Viscera_Pre-Treatment","Viscera_On-ADT","Viscera_On-Combination","Bone_Recurrence"))
#primecut.integrated.meta.vp$tissue_treatment=factor(primecut.integrated.meta.vp$tissue_treatment,levels=c("Bone_Pre-Treatment","Bone_On-ADT","Bone_On-Combination","Viscera_Pre-Treatment","Viscera_On-ADT","Viscera_On-Combination","LN_Pre-Treatment","LN_On-ADT","LN_On-Combination","Bone_Recurrence"))
names(vipercolors)=levels(primecut.integrated.meta.vp)
plot(DimPlot(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$tissue_treatment!="Bone_Recurrence" & primecut.integrated.meta.vp$tissue=="Bone")], reduction = "umap",label = TRUE,label.size=4,repel=T,split.by="tissue_treatment",ncol=3,cols=vipercolors[levels(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$tissue_treatment!="Bone_Recurrence" & primecut.integrated.meta.vp$tissue=="Bone")])]) + NoLegend())
plot(DimPlot(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$tissue=="Bone")], reduction = "umap",label = TRUE,label.size=4,repel=T,split.by="tissue_treatment",ncol=4,cols=vipercolors[levels(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$tissue=="Bone")])]) + NoLegend())
plot(DimPlot(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$tissue=="LN")], reduction = "umap",label = TRUE,label.size=4,repel=T,split.by="tissue_treatment",ncol=3,cols=vipercolors[levels(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$tissue=="LN")])]) + NoLegend())
plot(DimPlot(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$tissue=="Liver")], reduction = "umap",label = TRUE,label.size=4,repel=T,split.by="tissue_treatment",ncol=3,cols=vipercolors[levels(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$tissue=="Liver")])]) + NoLegend())
plot(DimPlot(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$tissue=="Lung")], reduction = "umap",label = TRUE,label.size=4,repel=T,split.by="tissue_treatment",ncol=3,cols=vipercolors[levels(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$tissue=="Lung")])]) + NoLegend())
#plot(DimPlot(primecut.integrated, reduction = "umap",label = TRUE,label.size=5,repel=T,split.by="tissue_treatment",ncol=4,cols=genecolors) + NoLegend())
plot(DimPlot(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$tissue_treatment!="Bone_Recurrence")], reduction = "umap",label = TRUE,label.size=4,repel=T,split.by="treatment",ncol=3,cols=vipercolors) + NoLegend())

markers.vp.all <- FindAllMarkers(primecut.integrated.meta.vp, only.pos = TRUE, min.pct = 0, logfc.threshold = 0,test.use = "t")
top10 <- markers.vp.all %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
geneHeatmap_primecut(primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$tissue_treatment!="Bone_Recurrence")],primecut.integrated.meta.vp$seurat_clusters[which(primecut.integrated.meta.vp$tissue_treatment!="Bone_Recurrence")],top10$gene,viper = T,scaled = T,n_top_genes_per_cluster = 5,color_palette = vipercolors,fontsize=5)

#Early responder	#91bfdb
#Stable disease	#ffffbf
#Late progressor	#fc8d59

Bone=as.data.frame.matrix(table(primecut.integrated.meta.vp$patient[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")],primecut.integrated.meta.vp$seurat_clusters[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")],primecut.integrated.meta.vp$tissue[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")])[c(1,6,9),,1])
Liver=as.data.frame.matrix(table(primecut.integrated.meta.vp$patient[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")],primecut.integrated.meta.vp$seurat_clusters[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")],primecut.integrated.meta.vp$tissue[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")])[c(4,5),,2])
LN=as.data.frame.matrix(table(primecut.integrated.meta.vp$patient[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")],primecut.integrated.meta.vp$seurat_clusters[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")],primecut.integrated.meta.vp$tissue[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")])[c(2,3,7),,3])
Lung=as.data.frame.matrix(table(primecut.integrated.meta.vp$patient[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")],primecut.integrated.meta.vp$seurat_clusters[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")],primecut.integrated.meta.vp$tissue[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")])[c(8,8),,4])
Bone=apply(Bone,1,function(x){x/sum(x)})
Liver=apply(Liver,1,function(x){x/sum(x)})
LN=apply(LN,1,function(x){x/sum(x)})
Lung=apply(Lung,1,function(x){x/sum(x)})
combined=as.data.frame(t(cbind(Bone,Liver,LN,Lung)))
combined$tissue=c(rep("Bone",3),rep("Liver",2),rep("LN",3),rep("Lung",2))
combined=combined[c(1:6,8:10),]
y <- melt(combined)
colnames(y)=c("Tissue","Cluster","Frequency")
p=ggplot(y, aes(x=Cluster, y=Frequency,fill=Tissue)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Cluster Frequency By Tissue Site")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))

Bone=as.data.frame.matrix(table(primecut.integrated.meta.vp$patient[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")],primecut.integrated.meta.vp$cell_lineages[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")],primecut.integrated.meta.vp$tissue[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")])[c(1,6,9),,1])
Liver=as.data.frame.matrix(table(primecut.integrated.meta.vp$patient[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")],primecut.integrated.meta.vp$cell_lineages[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")],primecut.integrated.meta.vp$tissue[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")])[c(4,5),,2])
LN=as.data.frame.matrix(table(primecut.integrated.meta.vp$patient[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")],primecut.integrated.meta.vp$cell_lineages[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")],primecut.integrated.meta.vp$tissue[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")])[c(2,3,7),,3])
Lung=as.data.frame.matrix(table(primecut.integrated.meta.vp$patient[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")],primecut.integrated.meta.vp$cell_lineages[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")],primecut.integrated.meta.vp$tissue[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")])[c(8,8),,4])
Bone=apply(Bone,1,function(x){x/sum(x)})
Liver=apply(Liver,1,function(x){x/sum(x)})
LN=apply(LN,1,function(x){x/sum(x)})
Lung=apply(Lung,1,function(x){x/sum(x)})
combined=as.data.frame(t(cbind(Bone,Liver,LN,Lung)))
combined$tissue=c(rep("Bone",3),rep("Liver",2),rep("LN",3),rep("Lung",2))
combined$patient=paste(combined$tissue,rownames(combined),sep="_")
combined=combined[1:9,]
combined$patient=c("Bone_1","Bone_2","Bone_3","Liver_1","Liver_2","LN_1","LN_2","LN_3","Lung")
combined=combined[c(1:6,8:9),]
combined$patient=c("Bone_1","Bone_2","Bone_3","Liver_1","Liver_2","LN_1","LN_2","Lung")
y <- melt(combined)
colnames(y)=c("Tissue","Patient","Cluster","Frequency")
p=ggplot(y, aes(x=Patient, y=Frequency,fill=Cluster)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Cluster Frequency By Tissue Site")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold")) + scale_fill_manual(values=lineage_colors)

Bone=as.data.frame.matrix(table(primecut.integrated.meta.vp$patient[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")],primecut.integrated.meta.vp$cell_lineages_toplevel[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")],primecut.integrated.meta.vp$tissue[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")])[c(1,6,9),,1])
Liver=as.data.frame.matrix(table(primecut.integrated.meta.vp$patient[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")],primecut.integrated.meta.vp$cell_lineages_toplevel[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")],primecut.integrated.meta.vp$tissue[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")])[c(4,5),,2])
LN=as.data.frame.matrix(table(primecut.integrated.meta.vp$patient[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")],primecut.integrated.meta.vp$cell_lineages_toplevel[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")],primecut.integrated.meta.vp$tissue[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")])[c(2,3,7),,3])
Lung=as.data.frame.matrix(table(primecut.integrated.meta.vp$patient[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")],primecut.integrated.meta.vp$cell_lineages_toplevel[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")],primecut.integrated.meta.vp$tissue[which(primecut.integrated.meta.vp$treatment=="Pre-Treatment")])[c(8,8),,4])
Bone=apply(Bone,1,function(x){x/sum(x)})
Liver=apply(Liver,1,function(x){x/sum(x)})
LN=apply(LN,1,function(x){x/sum(x)})
Lung=apply(Lung,1,function(x){x/sum(x)})
combined=as.data.frame(t(cbind(Bone,Liver,LN,Lung)))
combined$tissue=c(rep("Bone",3),rep("Liver",2),rep("LN",3),rep("Lung",2))
combined$patient=paste(combined$tissue,rownames(combined),sep="_")
combined=combined[1:9,]
combined$patient=c("Bone_1","Bone_2","Bone_3","Liver_1","Liver_2","LN_1","LN_2","LN_3","Lung")
combined=combined[c(1:6,8:9),]
combined$patient=c("Bone_1","Bone_2","Bone_3","Liver_1","Liver_2","LN_1","LN_2","Lung")
y <- melt(combined)
colnames(y)=c("Tissue","Patient","Cluster","Frequency")
p=ggplot(y, aes(x=Patient, y=Frequency,fill=Cluster)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Cluster Frequency By Tissue Site")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+scale_fill_manual(values=lineage_colors[c(1,8)])



####barplots split by treatment group
Bone=as.data.frame.matrix(table(primecut.integrated.meta.vp$treatment,primecut.integrated.meta.vp$cell_lineages,primecut.integrated.meta.vp$tissue)[,,1])
Liver=as.data.frame.matrix(table(primecut.integrated.meta.vp$treatment,primecut.integrated.meta.vp$cell_lineages,primecut.integrated.meta.vp$tissue)[c(1,3),,2])
LN=as.data.frame.matrix(table(primecut.integrated.meta.vp$treatment,primecut.integrated.meta.vp$cell_lineages,primecut.integrated.meta.vp$tissue)[c(1,2,3),,3])
Lung=as.data.frame.matrix(table(primecut.integrated.meta.vp$treatment,primecut.integrated.meta.vp$cell_lineages,primecut.integrated.meta.vp$tissue)[c(1,3),,4])
Bone=as.data.frame(apply(Bone,1,function(x){x/sum(x)}))
Liver=as.data.frame(apply(Liver,1,function(x){x/sum(x)}))
LN=as.data.frame(apply(LN,1,function(x){x/sum(x)}))
Lung=as.data.frame(apply(Lung,1,function(x){x/sum(x)}))
Bone$treatment=rownames(Bone)
y <- melt(Bone)
head(y)
colnames(y)=c("Cluster","Treatment","Frequency")
p=ggplot(y, aes(x=Treatment, y=Frequency,fill=Cluster)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Bone: Cluster Frequency By Treatment")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+scale_fill_manual(values=lineage_colors)
Liver$treatment=rownames(Liver)
y <- melt(Liver)
colnames(y)=c("Cluster","Treatment","Frequency")
p=ggplot(y, aes(x=Treatment, y=Frequency,fill=Cluster)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Liver: Cluster Frequency By Treatment")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+scale_fill_manual(values=lineage_colors)
LN$treatment=rownames(LN)
y <- melt(LN)
colnames(y)=c("Cluster","Treatment","Frequency")
p=ggplot(y, aes(x=Treatment, y=Frequency,fill=Cluster)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("LN: Cluster Frequency By Treatment")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+scale_fill_manual(values=lineage_colors)
Lung$treatment=rownames(Lung)
y <- melt(Lung)
colnames(y)=c("Cluster","Treatment","Frequency")
p=ggplot(y, aes(x=Treatment, y=Frequency,fill=Cluster)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Lung: Cluster Frequency By Treatment")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+scale_fill_manual(values=lineage_colors)


Bone=as.data.frame.matrix(table(tumorcells.vp$treatment,tumorcells.vp$seurat_clusters,tumorcells.vp$tissue)[,,1])
Liver=as.data.frame.matrix(table(tumorcells.vp$treatment,tumorcells.vp$seurat_clusters,tumorcells.vp$tissue)[c(1,3),,2])
LN=as.data.frame.matrix(table(tumorcells.vp$treatment,tumorcells.vp$seurat_clusters,tumorcells.vp$tissue)[c(1,2,3),,3])
Lung=as.data.frame.matrix(table(tumorcells.vp$treatment,tumorcells.vp$seurat_clusters,tumorcells.vp$tissue)[c(1,3),,4])
Bone=as.data.frame(apply(Bone,1,function(x){x/sum(x)}))
Liver=as.data.frame(apply(Liver,1,function(x){x/sum(x)}))
LN=as.data.frame(apply(LN,1,function(x){x/sum(x)}))
Lung=as.data.frame(apply(Lung,1,function(x){x/sum(x)}))
Bone$treatment=rownames(Bone)
y <- melt(Bone)
head(y)
colnames(y)=c("Cluster","Treatment","Frequency")
p=ggplot(y, aes(x=Treatment, y=Frequency,fill=Cluster)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Bone: Cluster Frequency By Treatment")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+scale_fill_manual(values=rev(lineage_colors))
Liver$treatment=rownames(Liver)
y <- melt(Liver)
colnames(y)=c("Cluster","Treatment","Frequency")
p=ggplot(y, aes(x=Treatment, y=Frequency,fill=Cluster)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Liver: Cluster Frequency By Treatment")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+scale_fill_manual(values=rev(lineage_colors))
LN$treatment=rownames(LN)
y <- melt(LN)
colnames(y)=c("Cluster","Treatment","Frequency")
p=ggplot(y, aes(x=Treatment, y=Frequency,fill=Cluster)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("LN: Cluster Frequency By Treatment")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+scale_fill_manual(values=rev(lineage_colors))
Lung$treatment=rownames(Lung)
y <- melt(Lung)
colnames(y)=c("Cluster","Treatment","Frequency")
p=ggplot(y, aes(x=Treatment, y=Frequency,fill=Cluster)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Lung: Cluster Frequency By Treatment")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+scale_fill_manual(values=rev(lineage_colors))

####

##############
table(primecut.integrated.meta.vp$patient,primecut.integrated.meta.vp$seurat_clusters,primecut.integrated.meta.vp$treatment)
pre=table(primecut.integrated.meta.vp$patient,primecut.integrated.meta.vp$seurat_clusters,primecut.integrated.meta.vp$treatment)[,,1]
adt=table(primecut.integrated.meta.vp$patient,primecut.integrated.meta.vp$seurat_clusters,primecut.integrated.meta.vp$treatment)[,,2]
comb=table(primecut.integrated.meta.vp$patient,primecut.integrated.meta.vp$seurat_clusters,primecut.integrated.meta.vp$treatment)[,,3]
recurr=table(primecut.integrated.meta.vp$patient,primecut.integrated.meta.vp$seurat_clusters,primecut.integrated.meta.vp$treatment)[,,4]
adt=adt[,setdiff(colnames(adt),c("Erythrocytes 1","B-cells 1","B-cells 1","B-cells 3","B-cells 4"))]
comb=comb[,setdiff(colnames(comb),c("Erythrocytes 1","B-cells 1","B-cells 1","B-cells 3","B-cells 4"))]
pre=pre[,setdiff(colnames(pre),c("Erythrocytes 1","B-cells 1","B-cells 1","B-cells 3","B-cells 4"))]
recurr=recurr[,setdiff(colnames(recurr),c("Erythrocytes 1","B-cells 1","B-cells 1","B-cells 3","B-cells 4"))]
pre=apply(pre,1,function(x){x/sum(x)})
adt=apply(adt,1,function(x){x/sum(x)})
comb=apply(comb,1,function(x){x/sum(x)})
recurr=apply(recurr,1,function(x){x/sum(x)})
####
combined=cbind(pre[,setdiff(colnames(pre),c("Patient6","Patient12"))],adt[,c("Patient11","Patient12","Patient3","Patient5","Patient8")],comb[,c("Patient1","Patient10","Patient6","Patient7")],recurr[,"Patient3"])
colnames(combined)=c("pre_p1","pre_p10","pre_p11","pre_p13","pre_p14","pre_p3","pre_p5","pre_p7","pre_p8","adt_p11","adt_p12","adt_p3","adt_p5","adt_p8","combo_p1","combo_p10","combo_p6","combo_p7","recurrence_p3")
y <- melt(combined)
colnames(y)=c("cluster","treatment","frequency")
y$treatment=unlist(lapply(strsplit(as.character(y$treatment),"_"),function(x){x[1]}))
y$treatment=factor(y$treatment, levels = c("pre","adt","combo","recurrence"))
y=y[which(y$treatment!="recurrence"),]
y$cluster=as.factor(y$cluster)
p=ggplot(y, aes(x=cluster, y=frequency,fill=treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Cluster Frequency By Treatment Group")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))
pvals=apply(combined,1,function(x){t.test(x[1:9],x[10:14])$p.value})
pvals=apply(combined,1,function(x){t.test(x[1:9],x[15:18])$p.value})
pvals=apply(combined,1,function(x){t.test(x[10:14],x[15:18])$p.value})
#pvals=apply(combined,1,function(x){t.test(x[c("pre_p11","pre_p3","pre_p5","pre_p8")],x[c("adt_p11","adt_p3","adt_p5","adt_p8")],paired=T)$p.value})
#pvals=apply(combined,1,function(x){t.test(x[c("pre_p1","pre_p10","pre_p7")],x[c("combo_p1","combo_p10","combo_p7")],paired=T)$p.value})
###Patient 1 and 3 are non-responders, patients 5,6,8,14 are really good responders, everyone else is intermediate
#combined=combined[,10:19]
#colnames(combined)=c("adt_p11_nonresponder","adt_p12_nonresponder","adt_p3_nonresponder","adt_p5_responder","adt_p8_responder","combo_p1_nonresponder","combo_p10_nonresponder","combo_p6_responder","combo_p7_responder","recurrence_p3_nonresponder")
combined=combined[,1:9]
colnames(combined)=c("pre_p1_nonresponder","pre_p10_sd","pre_p11_sd","pre_p13_sd","pre_p14_responder","pre_p3_nonresponder","pre_p5_responder","pre_p7_sd","pre_p8_responder")
y <- melt(combined)
colnames(y)=c("cluster","treatment","frequency")
treatment=unlist(lapply(strsplit(as.character(y$treatment),"_"),function(x){x[1]}))
response=unlist(lapply(strsplit(as.character(y$treatment),"_"),function(x){x[3]}))
y$response=response
y$treatment=treatment
y$treatment.response=paste(y$treatment,y$response,sep=".")
y=y[which(y$response!="sd"),]
#y$treatment=factor(y$treatment, levels = c("adt","combo","recurrence"))
y$cluster=as.factor(y$cluster)
p=ggplot(y, aes(x=cluster, y=frequency,fill=treatment.response)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Cluster Frequency By Treatment Group and Treatment Response")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))

'''
adt.pre=adt[,c("Patient11","Patient3","Patient5","Patient8")]/pre[,c("Patient11","Patient3","Patient5","Patient8")]
comb.pre=comb[,c("Patient1","Patient10","Patient7")]/pre[,c("Patient1","Patient10","Patient7")]
recurr.pre=recurr[,c("Patient3")]/pre[,c("Patient3")]
combined=cbind(adt.pre,comb.pre,recurr.pre)
colnames(combined)=c("Patient11","Patient3","Patient5","Patient8","Patient1","Patient10","Patient7","Patient3_Recurrence")
y=melt(combined)
colnames(y)=c("cluster","patient","logFC")
y$patient=factor(y$patient, levels = c("Patient11","Patient3","Patient5","Patient8","Patient1","Patient10","Patient7","Patient3_Recurrence"))
y$cluster=as.factor(y$cluster)
y$logFC=log2(y$logFC)
p=ggplot(y, aes(x=cluster, y=logFC,fill=patient)) +
  #geom_boxplot() +
  geom_col(position="dodge")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Post/Pre Cluster Frequency")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+scale_fill_manual(values=c("cyan3","dodgerblue3","darkturquoise","deepskyblue2","darkgoldenrod2","gold1","lightgoldenrod1","firebrick"))
'''



#VIPER tumor sub-clustering
##sub-cluster tumor cells
#tumorcells.vp=primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$seurat_clusters %in% c("Epithelial 1","Epithelial 2"))]
tumorcells.vp=primecut.integrated.meta.vp[,which(primecut.integrated.meta.vp$seurat_clusters %in% c("Epithelial 1","Epithelial 2","Epithelial 3"))]
tumorcells.vp <- RunPCA(tumorcells.vp, features = rownames(tumorcells.vp),npcs = 50)
tumorcells.vp <- RunUMAP(tumorcells.vp, dims = 1:50, verbose = FALSE,umap.method="umap-learn",metric="correlation")
tumorcells.vp <- FindNeighbors(tumorcells.vp, dims = 1:50, verbose = FALSE)
tumorcells.vp <- FindClusters(tumorcells.vp, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) # set resolution <1 for fewer clusters (default is 0.5)
clust=tumorcells.vp@meta.data[,which(grepl("RNA_snn_res.",colnames(tumorcells.vp@meta.data)))]
mat=tumorcells.vp@assays$RNA@scale.data
out=sil_subsample_v2(mat,clust)
means.2=out[[1]]
sd.2=out[[2]]
x=seq(0.01,1,by=0.01)
errbar(x,means.2,means.2+sd.2,means.2-sd.2,ylab="mean silhouette score",xlab="resolution parameter")
lines(x,means.2)
best=tail((x[10:length(x)])[which(means.2[10:length(means.2)]==max(means.2[10:length(means.2)]))],n=1)
print(best)
tumorcells.vp$seurat_clusters=tumorcells.vp@meta.data[,which(colnames(tumorcells.vp@meta.data)==paste("RNA_snn_res.",0.13,sep=""))]
Idents(tumorcells.vp) <- "seurat_clusters"
DimPlot(tumorcells.vp, reduction = "umap",label = TRUE) + NoLegend()
tumorcells.vp.markers <- FindAllMarkers(tumorcells.vp, only.pos = TRUE, min.pct = 0, logfc.threshold = 0,test.use = "t")
top10 <- tumorcells.vp.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
geneHeatmap_primecut(tumorcells.vp,tumorcells.vp$seurat_clusters,top10$gene,viper = T,scaled = T)
geneHeatmap_primecut(tumorcells.vp,tumorcells.vp$seurat_clusters,top10$gene,viper = T,scaled = T,color_palette = rev(lineage_colors))
DimPlot(tumorcells.vp, reduction = "umap", split.by = "treatment",group.by="seurat_clusters")
DimPlot(tumorcells.vp, reduction = "umap",group.by="seurat_clusters",label=T,cols=rev(lineage_colors),label.size = 6)+NoLegend()
saveRDS(tumorcells.vp,"~/Documents/Documents/MED SCHOOL/summer 2019/primecut.integrated.tumorcells.vp.rds")
table(tumorcells.vp$patient,tumorcells.vp$seurat_clusters,tumorcells.vp$treatment)
plot(DimPlot(tumorcells.vp[,which(tumorcells.vp$patient=="Patient3")],reduction="umap",group.by="seurat_clusters",split.by="treatment") + NoLegend())
#####
x=primecut.integrated[,colnames(tumorcells.vp)]
x$seurat_clusters=tumorcells.vp$seurat_clusters
Idents(x)="seurat_clusters"
markers <- FindAllMarkers(x, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "MAST")
top10 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
geneHeatmap_primecut(x,x$seurat_clusters,top10$gene)

##############
####
table(tumorcells.vp$patient,tumorcells.vp$seurat_clusters,tumorcells.vp$treatment)
pre=table(tumorcells.vp$patient,tumorcells.vp$seurat_clusters,tumorcells.vp$treatment)[,,1]
adt=table(tumorcells.vp$patient,tumorcells.vp$seurat_clusters,tumorcells.vp$treatment)[,,2]
comb=table(tumorcells.vp$patient,tumorcells.vp$seurat_clusters,tumorcells.vp$treatment)[,,3]
recurr=table(tumorcells.vp$patient,tumorcells.vp$seurat_clusters,tumorcells.vp$treatment)[,,4]
pre=apply(pre,1,function(x){x/sum(x)})
adt=apply(adt,1,function(x){x/sum(x)})
comb=apply(comb,1,function(x){x/sum(x)})
recurr=apply(recurr,1,function(x){x/sum(x)})
combined=cbind(pre[,setdiff(colnames(pre),c("Patient6","Patient12"))],adt[,c("Patient11","Patient12","Patient3","Patient5","Patient8")],comb[,c("Patient1","Patient10","Patient6","Patient7")],recurr[,"Patient3"])
colnames(combined)=c("pre_p1","pre_p10","pre_p11","pre_p13","pre_p14","pre_p3","pre_p5","pre_p7","pre_p8","adt_p11","adt_p12","adt_p3","adt_p5","adt_p8","combo_p1","combo_p10","combo_p6","combo_p7","recurrence_p3")
y <- melt(combined)
colnames(y)=c("cluster","treatment","frequency")
y$treatment=unlist(lapply(strsplit(as.character(y$treatment),"_"),function(x){x[1]}))
y$treatment=factor(y$treatment, levels = c("pre","adt","combo","recurrence"))
y$cluster=as.factor(y$cluster)
p=ggplot(y, aes(x=cluster, y=frequency,fill=treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Cluster Frequency By Treatment Group")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))
pvals=apply(combined,1,function(x){t.test(x[1:9],x[10:14])$p.value})
pvals=apply(combined,1,function(x){t.test(x[1:9],x[15:18])$p.value})
pvals=apply(combined,1,function(x){t.test(x[10:14],x[15:18])$p.value})
#pvals=apply(combined,1,function(x){t.test(x[c("pre_p11","pre_p3","pre_p5","pre_p8")],x[c("adt_p11","adt_p3","adt_p5","adt_p8")],paired=T)$p.value})
#pvals=apply(combined,1,function(x){t.test(x[c("pre_p1","pre_p10","pre_p7")],x[c("combo_p1","combo_p10","combo_p7")],paired=T)$p.value})
##patient 6,7,8,5 are responders :: rest are non-responders
#combined=combined[,10:19]
#colnames(combined)=c("adt_p11_nonresponder","adt_p12_nonresponder","adt_p3_nonresponder","adt_p5_responder","adt_p8_responder","combo_p1_nonresponder","combo_p10_nonresponder","combo_p6_responder","combo_p7_responder","recurrence_p3_nonresponder")
combined=combined[,1:9]
colnames(combined)=c("pre_p1_nonresponder","pre_p10_sd","pre_p11_sd","pre_p13_sd","pre_p14_responder","pre_p3_nonresponder","pre_p5_responder","pre_p7_sd","pre_p8_responder")
y <- melt(combined)
colnames(y)=c("cluster","treatment","frequency")
treatment=unlist(lapply(strsplit(as.character(y$treatment),"_"),function(x){x[1]}))
response=unlist(lapply(strsplit(as.character(y$treatment),"_"),function(x){x[3]}))
y$response=response
y$treatment=treatment
y$treatment.response=paste(y$treatment,y$response,sep=".")
y=y[which(y$response!="sd"),]
#y$treatment=factor(y$treatment, levels = c("adt","combo","recurrence"))
y$cluster=as.factor(y$cluster)
p=ggplot(y, aes(x=cluster, y=frequency,fill=treatment.response)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Cluster Frequency By Treatment Group and Treatment Response")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))


###correlate cluster signatures with outcome in SU2C
abida_dat=readRDS("/Users/aleksandar/Dropbox/Abate_Shen_scRNAseq/human_datasets/Abida-2019 (E. Coast SU2C)/processed-data/abida-2019-processed-data.rds")  ##transcriptome
colnames(abida_dat$metadata)
metadata=abida_dat$metadata
rownames(metadata)=metadata$SAMPLE_ID
counts_mat=abida_dat$counts_fpkm.polyA
counts_mat=counts_mat[which(counts_mat$Hugo_Symbol %in% names(which(table(counts_mat$Hugo_Symbol)==1))),]
rownames(counts_mat)=counts_mat$Hugo_Symbol
counts_mat=counts_mat[,2:ncol(counts_mat)]
shared_cols=intersect(colnames(counts_mat),rownames(metadata))
counts_mat=counts_mat[,shared_cols]
metadata=metadata[shared_cols,]
Gene=rownames(counts_mat)
x=cbind(Gene,counts_mat)
write.table(x,row.names=F,sep="\t",quote=F,"~/Downloads/abida_dat_for_cibersortx.txt")
counts_mat=t(apply(counts_mat,1,function(x){(x-mean(x))/sd(x)}))
counts_mat=counts_mat[complete.cases(counts_mat),]
g0=rep(1,length(tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="0")]))
names(g0)=tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="0")]
g1=rep(1,length(tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="1")]))
names(g1)=tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="1")]
g2=rep(1,length(tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="2")]))
names(g2)=tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="2")]
g3=rep(1,length(tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="3")]))
names(g3)=tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="3")]
g4=rep(1,length(tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="4")]))
names(g4)=tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="4")]
g5=rep(1,length(tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="5")]))
names(g5)=tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="5")]
g6=rep(1,length(tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="6")]))
names(g6)=tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="6")]
g7=rep(1,length(tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="7")]))
names(g7)=tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="7")]
set.seed(1234)
#g0=rep(1,length(markers.vp.all$gene[which(markers.vp.all$cluster=="CD4 T-cell 1")]))
#names(g0)=markers.vp.all$gene[which(markers.vp.all$cluster=="CD4 T-cell 1")]
#g1=rep(1,length(markers.vp.all$gene[which(markers.vp.all$cluster=="CD8 T-cell 2")]))
#names(g1)=markers.vp.all$gene[which(markers.vp.all$cluster=="CD8 T-cell 2")]
gseaBySample_g0=apply(counts_mat,2,function(x){gsea(signature = x, geneset =g0, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g1=apply(counts_mat,2,function(x){gsea(signature = x, geneset =g1, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g2=apply(counts_mat,2,function(x){gsea(signature = x, geneset =g2, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g3=apply(counts_mat,2,function(x){gsea(signature = x, geneset =g3, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g4=apply(counts_mat,2,function(x){gsea(signature = x, geneset =g4, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g5=apply(counts_mat,2,function(x){gsea(signature = x, geneset =g5, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g6=apply(counts_mat,2,function(x){gsea(signature = x, geneset =g6, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g7=apply(counts_mat,2,function(x){gsea(signature = x, geneset =g7, twoTails = F, pout = F, per = 100)$nes})
filenames <- list.files("/Users/aleksandar/Documents/Documents/MED SCHOOL/summer 2019/prime-cut-metacells/sc_nets/", pattern="*.rds", full.names=TRUE)
nets=lapply(filenames,readRDS)
vp_abida<-viper(counts_mat, nets, method = 'none')
set.seed(1234)
gseaBySample_g0_vp=apply(vp_abida,2,function(x){gsea(signature = x, geneset =g0, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g1_vp=apply(vp_abida,2,function(x){gsea(signature = x, geneset =g1, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g2_vp=apply(vp_abida,2,function(x){gsea(signature = x, geneset =g2, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g3_vp=apply(vp_abida,2,function(x){gsea(signature = x, geneset =g3, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g4_vp=apply(vp_abida,2,function(x){gsea(signature = x, geneset =g4, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g5_vp=apply(vp_abida,2,function(x){gsea(signature = x, geneset =g5, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g6_vp=apply(vp_abida,2,function(x){gsea(signature = x, geneset =g6, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g7_vp=apply(vp_abida,2,function(x){gsea(signature = x, geneset =g7, twoTails = F, pout = F, per = 100)$nes})
####correlate to gleason score....
gleason=metadata$GLEASON_SCORE
vital_status=metadata$OS_STATUS
days_to_death=metadata$OS_MONTHS
x=data.frame(gleason,GSEA=gseaBySample_g0)
x=x[which(gleason != "UNK"),]
x$gleason=as.factor(x$gleason)
x$gleason <- factor(x$gleason, levels = c("6", "7", "8", "9", "10", "11"))
p <- ggplot(x, aes(x=gleason, y=GSEA,fill=gleason)) +
  geom_boxplot(outlier.shape=NA)+
  ggtitle("Cluster 0 GSEA by Gleason Score") +
  theme(plot.title = element_text(hjust = 0.5))
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=12, face="bold"),
  legend.text = element_text(size=12, face="bold"),
  legend.title = element_text(size=14, face="bold"),
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.title.x = element_blank())
###test recurrence-free survival
clinical_dat=data.frame(os=vital_status,os_months=as.numeric(days_to_death))
clinical_dat$GSEA0=gseaBySample_g0
clinical_dat$GSEA1=gseaBySample_g1
clinical_dat$GSEA2=gseaBySample_g2
clinical_dat$GSEA3=gseaBySample_g3
clinical_dat$GSEA4=gseaBySample_g4
clinical_dat$GSEA5=gseaBySample_g5
clinical_dat$GSEA6=gseaBySample_g6
clinical_dat$GSEA7=gseaBySample_g7
clinical_dat=clinical_dat[which(clinical_dat$os!=""),]
clinical_dat$os=as.character(clinical_dat$os)
clinical_dat$os[which(clinical_dat$os=="LIVING")]="F"
clinical_dat$os[which(clinical_dat$os=="DECEASED")]="T"
clinical_dat$os=as.logical(clinical_dat$os)
clinical_dat$SurvObj=with(clinical_dat,Surv(os_months,os))
library(randomForest)
library(pROC)
library(Boruta)
set.seed(1234)
clinical_dat$os=as.factor(clinical_dat$os)
#boruta.train <- Boruta(os~., data = clinical_dat[,c(1,3:10)],doTrace = 2,pValue=1e-10)
#print(boruta.train)
#x=clinical_dat[,c(names(boruta.train$finalDecision[which(boruta.train$finalDecision %in% c("Confirmed","Tentative"))]),"os")]
model=randomForest(os~.,data=clinical_dat[,c(1,3:10)],ntree=50000,importance=T)
varImpPlot(model,cex=.7)
rf.roc<-roc(as.factor(clinical_dat$os),model$votes[,2],ci=T,ci.alpha=0.95,stratified=F,plot=T,print.auc=T)
plot.roc(rf.roc,print.auc=T,asp=NA)
clinical_dat$os=as.logical(clinical_dat$os)
#####
res.cox1 <- coxph(SurvObj ~ ., data =  clinical_dat[,3:11])
summary(res.cox1)
plot(ggforest(res.cox1, data = clinical_dat,fontsize = 1,main = "GSEA Hazard Ratio"))
#GSEA
res.cut <- surv_cutpoint(clinical_dat, time = "os_months", event = "os",variables = c("GSEA1"))
summary(res.cut)
plot(res.cut, "GSEA1", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(os_months,os) ~GSEA1, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Overall Survival",xlab="Time (months)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))
res.cut <- surv_cutpoint(clinical_dat, time = "os_months", event = "os",variables = c("GSEA0"))
summary(res.cut)
plot(res.cut, "GSEA0", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(os_months,os) ~GSEA0, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Overall Survival",xlab="Time (months)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))
res.cut <- surv_cutpoint(clinical_dat, time = "os_months", event = "os",variables = c("GSEA2"))
summary(res.cut)
plot(res.cut, "GSEA2", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(os_months,os) ~GSEA2, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Overall Survival",xlab="Time (months)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))
res.cut <- surv_cutpoint(clinical_dat, time = "os_months", event = "os",variables = c("GSEA3"))
summary(res.cut)
plot(res.cut, "GSEA3", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(os_months,os) ~GSEA3, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Overall Survival",xlab="Time (months)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))
res.cut <- surv_cutpoint(clinical_dat, time = "os_months", event = "os",variables = c("GSEA4"))
summary(res.cut)
plot(res.cut, "GSEA4", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(os_months,os) ~GSEA4, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Overall Survival",xlab="Time (months)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))
res.cut <- surv_cutpoint(clinical_dat, time = "os_months", event = "os",variables = c("GSEA5"))
summary(res.cut)
plot(res.cut, "GSEA5", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(os_months,os) ~GSEA5, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Overall Survival",xlab="Time (months)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))
res.cut <- surv_cutpoint(clinical_dat, time = "os_months", event = "os",variables = c("GSEA6"))
summary(res.cut)
plot(res.cut, "GSEA6", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(os_months,os) ~GSEA6, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Overall Survival",xlab="Time (months)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))
res.cut <- surv_cutpoint(clinical_dat, time = "os_months", event = "os",variables = c("GSEA7"))
summary(res.cut)
plot(res.cut, "GSEA7", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(os_months,os) ~GSEA7, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Overall Survival",xlab="Time (months)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))



### West Coast SU2C
data=read.table("/Users/aleksandar/Downloads/2018_04_15_matrix_rna_tpm.txt",sep="\t",header=T)
rownames(data)=data[,1]
data=data[,2:ncol(data)]
metadata=read.table("/Users/aleksandar/Downloads/2021_05_13_european_urology_2019.txt",sep="\t",header=T)
rownames(metadata)=gsub("-",".",metadata$sample_id)
shared_cols=intersect(colnames(data),rownames(metadata))
data=data[,shared_cols]
metadata=metadata[shared_cols,]
Gene=rownames(data)
x=cbind(Gene,data)
write.table(x,row.names=F,sep="\t",quote=F,"~/Downloads/westcoast_SU2C_dat_for_cibersortx.txt")
counts_mat=t(apply(data,1,function(x){(x-mean(x))/sd(x)}))
counts_mat=counts_mat[complete.cases(counts_mat),]
g0=rep(1,length(tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="0")]))
names(g0)=tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="0")]
g1=rep(1,length(tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="1")]))
names(g1)=tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="1")]
g2=rep(1,length(tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="2")]))
names(g2)=tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="2")]
g3=rep(1,length(tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="3")]))
names(g3)=tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="3")]
g4=rep(1,length(tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="4")]))
names(g4)=tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="4")]
g5=rep(1,length(tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="5")]))
names(g5)=tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="5")]
g6=rep(1,length(tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="6")]))
names(g6)=tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="6")]
g7=rep(1,length(tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="7")]))
names(g7)=tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="7")]
set.seed(1234)
#g0=rep(1,length(markers.vp.all$gene[which(markers.vp.all$cluster=="CD4 T-cell 1")]))
#names(g0)=markers.vp.all$gene[which(markers.vp.all$cluster=="CD4 T-cell 1")]
#g1=rep(1,length(markers.vp.all$gene[which(markers.vp.all$cluster=="CD8 T-cell 2")]))
#names(g1)=markers.vp.all$gene[which(markers.vp.all$cluster=="CD8 T-cell 2")]
gseaBySample_g0=apply(counts_mat,2,function(x){gsea(signature = x, geneset =g0, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g1=apply(counts_mat,2,function(x){gsea(signature = x, geneset =g1, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g2=apply(counts_mat,2,function(x){gsea(signature = x, geneset =g2, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g3=apply(counts_mat,2,function(x){gsea(signature = x, geneset =g3, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g4=apply(counts_mat,2,function(x){gsea(signature = x, geneset =g4, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g5=apply(counts_mat,2,function(x){gsea(signature = x, geneset =g5, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g6=apply(counts_mat,2,function(x){gsea(signature = x, geneset =g6, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g7=apply(counts_mat,2,function(x){gsea(signature = x, geneset =g7, twoTails = F, pout = F, per = 100)$nes})
filenames <- list.files("/Users/aleksandar/Documents/Documents/MED SCHOOL/summer 2019/prime-cut-metacells/sc_nets/", pattern="*.rds", full.names=TRUE)
nets=lapply(filenames,readRDS)
vp_abida<-viper(counts_mat, nets, method = 'none')
set.seed(1234)
gseaBySample_g0_vp=apply(vp_abida,2,function(x){gsea(signature = x, geneset =g0, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g1_vp=apply(vp_abida,2,function(x){gsea(signature = x, geneset =g1, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g2_vp=apply(vp_abida,2,function(x){gsea(signature = x, geneset =g2, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g3_vp=apply(vp_abida,2,function(x){gsea(signature = x, geneset =g3, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g4_vp=apply(vp_abida,2,function(x){gsea(signature = x, geneset =g4, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g5_vp=apply(vp_abida,2,function(x){gsea(signature = x, geneset =g5, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g6_vp=apply(vp_abida,2,function(x){gsea(signature = x, geneset =g6, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g7_vp=apply(vp_abida,2,function(x){gsea(signature = x, geneset =g7, twoTails = F, pout = F, per = 100)$nes})
vital_status=metadata$event.Dead
days_to_death=metadata$OS.mCRPC
###test recurrence-free survival
clinical_dat=data.frame(os=vital_status,os_months=as.numeric(days_to_death))
clinical_dat$GSEA0=gseaBySample_g0
clinical_dat$GSEA1=gseaBySample_g1
clinical_dat$GSEA2=gseaBySample_g2
clinical_dat$GSEA3=gseaBySample_g3
clinical_dat$GSEA4=gseaBySample_g4
clinical_dat$GSEA5=gseaBySample_g5
clinical_dat$GSEA6=gseaBySample_g6
clinical_dat$GSEA7=gseaBySample_g7
clinical_dat=clinical_dat[which(clinical_dat$os!=""),]
clinical_dat$os=as.character(clinical_dat$os)
clinical_dat$os[which(clinical_dat$os=="0")]="F"
clinical_dat$os[which(clinical_dat$os=="1")]="T"
clinical_dat$os=as.logical(clinical_dat$os)
clinical_dat$SurvObj=with(clinical_dat,Surv(os_months,os))
library(randomForest)
library(pROC)
library(Boruta)
set.seed(1234)
clinical_dat$os=as.factor(clinical_dat$os)
#boruta.train <- Boruta(os~., data = clinical_dat[,c(1,3:10)],doTrace = 2,pValue=1e-10)
#print(boruta.train)
#x=clinical_dat[,c(names(boruta.train$finalDecision[which(boruta.train$finalDecision %in% c("Confirmed","Tentative"))]),"os")]
model=randomForest(os~.,data=clinical_dat[,c(1,3:10)],ntree=50000,importance=T)
varImpPlot(model,cex=.7)
rf.roc<-roc(as.factor(clinical_dat$os),model$votes[,2],ci=T,ci.alpha=0.95,stratified=F,plot=T,print.auc=T)
plot.roc(rf.roc,print.auc=T,asp=NA)
clinical_dat$os=as.logical(clinical_dat$os)
#####
res.cox1 <- coxph(SurvObj ~ ., data =  clinical_dat[,3:11])
summary(res.cox1)
plot(ggforest(res.cox1, data = clinical_dat,fontsize = 1,main = "GSEA Hazard Ratio"))
#GSEA
res.cut <- surv_cutpoint(clinical_dat, time = "os_months", event = "os",variables = c("GSEA1"))
summary(res.cut)
plot(res.cut, "GSEA1", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(os_months,os) ~GSEA1, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Overall Survival",xlab="Time (months)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))
res.cut <- surv_cutpoint(clinical_dat, time = "os_months", event = "os",variables = c("GSEA0"))
summary(res.cut)
plot(res.cut, "GSEA0", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(os_months,os) ~GSEA0, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Overall Survival",xlab="Time (months)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))
res.cut <- surv_cutpoint(clinical_dat, time = "os_months", event = "os",variables = c("GSEA2"))
summary(res.cut)
plot(res.cut, "GSEA2", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(os_months,os) ~GSEA2, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Overall Survival",xlab="Time (months)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))
res.cut <- surv_cutpoint(clinical_dat, time = "os_months", event = "os",variables = c("GSEA3"))
summary(res.cut)
plot(res.cut, "GSEA3", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(os_months,os) ~GSEA3, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Overall Survival",xlab="Time (months)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))
res.cut <- surv_cutpoint(clinical_dat, time = "os_months", event = "os",variables = c("GSEA4"))
summary(res.cut)
plot(res.cut, "GSEA4", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(os_months,os) ~GSEA4, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Overall Survival",xlab="Time (months)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))
res.cut <- surv_cutpoint(clinical_dat, time = "os_months", event = "os",variables = c("GSEA5"))
summary(res.cut)
plot(res.cut, "GSEA5", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(os_months,os) ~GSEA5, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Overall Survival",xlab="Time (months)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))
res.cut <- surv_cutpoint(clinical_dat, time = "os_months", event = "os",variables = c("GSEA6"))
summary(res.cut)
plot(res.cut, "GSEA6", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(os_months,os) ~GSEA6, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Overall Survival",xlab="Time (months)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))
res.cut <- surv_cutpoint(clinical_dat, time = "os_months", event = "os",variables = c("GSEA7"))
summary(res.cut)
plot(res.cut, "GSEA7", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(os_months,os) ~GSEA7, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Overall Survival",xlab="Time (months)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))




##Correlate tumor cell subclsters with outcome in TCGA
#TCGA
library('recount')
## Download the RangedSummarizedExperiment object at the gene level
#url <- download_study('TCGA')
## Load the data
#load(file.path('TCGA', 'rse_gene.Rdata'))
load("/Users/aleksandar/Documents/Documents/MED SCHOOL/summer 2019/published_melanoma_datasets/pd1_response_publisheddata/TCGA/rse_gene.Rdata")
a=colData(rse_gene)
table(a$gdc_cases.project.project_id)
days_to_death=a$gdc_cases.diagnoses.days_to_death[which(a$gdc_cases.project.project_id=="TCGA-PRAD")]
vital_status=a$gdc_cases.diagnoses.vital_status[which(a$gdc_cases.project.project_id=="TCGA-PRAD")]
days_to_last_follow_up=a$gdc_cases.diagnoses.days_to_last_follow_up[which(a$gdc_cases.project.project_id=="TCGA-PRAD")]
time_to_recurrence=a$xml_days_to_first_biochemical_recurrence[which(a$gdc_cases.project.project_id=="TCGA-PRAD")]
recurrence=a$xml_biochemical_recurrence[which(a$gdc_cases.project.project_id=="TCGA-PRAD")]
gleason=a$xml_stage_event_gleason_grading[which(a$gdc_cases.project.project_id=="TCGA-PRAD")]
gleason=sapply(gleason,function(x){as.numeric(substr(x, 1, 1))})
gleason[gleason==1]=10
## Scale counts
rse_gene_scaled <- scale_counts(rse_gene)
counts <- assays(rse_gene_scaled)$counts
counts=counts[,which(a$gdc_cases.project.project_id=="TCGA-PRAD")]
b=rowData(rse_gene)
Ensembl=unlist(b$gene_id)
geneNames=unlist(b$symbol)
rownames(counts)=Ensembl
Ensembl_simplified=unlist(lapply(strsplit(Ensembl,"\\."),function(x){x[1]}))
rownames(counts)=Ensembl_simplified
#counts=Ensemble2GeneName(counts)
x=AnnotationDbi::select(org.Hs.eg.db, keys=rownames(counts), columns=c("SYMBOL"), keytype="ENSEMBL")
x=x[which(!is.na(x[,2])),]
entrezids=c()
exclude=c()
include=c()
for(i in rownames(counts)){
  ids=x[which(x[,1]==i),2]
  if(length(ids)==0){
    exclude=c(exclude,i)
  }
  else{
    entrezids=c(entrezids,ids[1])
    include=c(include,i)
  }
}
counts=counts[include,]
rownames(counts)=entrezids
Gene=rownames(counts)
x=cbind(Gene,counts)
write.table(x,row.names=F,sep="\t",quote=F,"~/Downloads/tcga_prad_for_cibersortx.txt")
tcga_dat=log10(counts+1)
tcga_dat=t(apply(tcga_dat,1,function(x){(x-mean(x))/sd(x)}))
tcga_dat=tcga_dat[which(rownames(tcga_dat) %in% names(which(table(rownames(tcga_dat))==1))),]
tcga_dat=tcga_dat[complete.cases(tcga_dat),]
rm(a,b,rse_gene,rse_gene_scaled,counts)
###progression-free-survival
g0=rep(1,length(tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="0")]))
names(g0)=tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="0")]
g1=rep(1,length(tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="1")]))
names(g1)=tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="1")]
g2=rep(1,length(tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="2")]))
names(g2)=tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="2")]
g3=rep(1,length(tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="3")]))
names(g3)=tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="3")]
g4=rep(1,length(tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="4")]))
names(g4)=tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="4")]
g5=rep(1,length(tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="5")]))
names(g5)=tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="5")]
g6=rep(1,length(tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="6")]))
names(g6)=tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="6")]
g7=rep(1,length(tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="7")]))
names(g7)=tumorcells.vp.markers$gene[which(tumorcells.vp.markers$cluster=="7")]
library(survival)
library(survminer)
set.seed(1234)
#g0=rep(1,length(markers.vp.all$gene[which(markers.vp.all$cluster=="CD4 T-cell 1")]))
#names(g0)=markers.vp.all$gene[which(markers.vp.all$cluster=="CD4 T-cell 1")]
#g1=rep(1,length(markers.vp.all$gene[which(markers.vp.all$cluster=="CD8 T-cell 2")]))
#names(g1)=markers.vp.all$gene[which(markers.vp.all$cluster=="CD8 T-cell 2")]
gseaBySample_g0=apply(tcga_dat,2,function(x){gsea(signature = x, geneset =g0, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g1=apply(tcga_dat,2,function(x){gsea(signature = x, geneset =g1, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g2=apply(tcga_dat,2,function(x){gsea(signature = x, geneset =g2, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g3=apply(tcga_dat,2,function(x){gsea(signature = x, geneset =g3, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g4=apply(tcga_dat,2,function(x){gsea(signature = x, geneset =g4, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g5=apply(tcga_dat,2,function(x){gsea(signature = x, geneset =g5, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g6=apply(tcga_dat,2,function(x){gsea(signature = x, geneset =g6, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g7=apply(tcga_dat,2,function(x){gsea(signature = x, geneset =g7, twoTails = F, pout = F, per = 100)$nes})
set.seed(1234)
filenames <- list.files("/Users/aleksandar/Documents/Documents/MED SCHOOL/summer 2019/prime-cut-metacells/sc_nets/", pattern="*.rds", full.names=TRUE)
nets=lapply(filenames,readRDS)
vp_tcga<-viper(tcga_dat, nets, method = 'none')
set.seed(1234)
gseaBySample_g0_vp=apply(vp_tcga,2,function(x){gsea(signature = x, geneset =g0, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g1_vp=apply(vp_tcga,2,function(x){gsea(signature = x, geneset =g1, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g2_vp=apply(vp_tcga,2,function(x){gsea(signature = x, geneset =g2, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g3_vp=apply(vp_tcga,2,function(x){gsea(signature = x, geneset =g3, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g4_vp=apply(vp_tcga,2,function(x){gsea(signature = x, geneset =g4, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g5_vp=apply(vp_tcga,2,function(x){gsea(signature = x, geneset =g5, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g6_vp=apply(vp_tcga,2,function(x){gsea(signature = x, geneset =g6, twoTails = F, pout = F, per = 100)$nes})
gseaBySample_g7_vp=apply(vp_tcga,2,function(x){gsea(signature = x, geneset =g7, twoTails = F, pout = F, per = 100)$nes})
####correlate to gleason score....
x=data.frame(gleason,cluster0GSEA=gseaBySample_g0_vp,cluster1GSEA=gseaBySample_g1_vp,cluster2GSEA=gseaBySample_g2_vp,cluster3GSEA=gseaBySample_g3_vp,cluster4GSEA=gseaBySample_g4_vp,cluster5GSEA=gseaBySample_g5_vp,cluster6GSEA=gseaBySample_g6_vp,cluster7GSEA=gseaBySample_g7_vp)
x$gleason=as.factor(x$gleason)
p <- ggplot(x, aes(x=gleason, y=cluster1GSEA,fill=gleason)) +
  geom_boxplot(outlier.shape=NA)+
  ggtitle("Cluster 1 GSEA by Gleason Score") +
  theme(plot.title = element_text(hjust = 0.5))
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=12, face="bold"),
  legend.text = element_text(size=12, face="bold"),
  legend.title = element_text(size=14, face="bold"),
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.title.x = element_blank())
###test recurrence-free survival
clinical_dat=data.frame(recurrence=recurrence,time_to_recurrence=as.numeric(time_to_recurrence),days_to_last_follow_up=as.numeric(days_to_last_follow_up))
clinical_dat$time_to_recurrence[which(is.na(clinical_dat$time_to_recurrence))]=clinical_dat$days_to_last_follow_up[which(is.na(clinical_dat$time_to_recurrence))]
clinical_dat$GSEA0=gseaBySample_g0
clinical_dat$GSEA1=gseaBySample_g1
clinical_dat$GSEA2=gseaBySample_g2
clinical_dat$GSEA3=gseaBySample_g3
clinical_dat$GSEA4=gseaBySample_g4
clinical_dat$GSEA5=gseaBySample_g5
clinical_dat$GSEA6=gseaBySample_g6
clinical_dat$GSEA7=gseaBySample_g7
clinical_dat$OS=as.logical(as.numeric(clinical_dat$recurrence)-2)
clinical_dat$Recurrence=clinical_dat$time_to_recurrence
clinical_dat=clinical_dat[which(!is.na(clinical_dat$recurrence)),]
clinical_dat$SurvObj=with(clinical_dat,Surv(Recurrence,OS))
library(randomForest)
library(pROC)
library(Boruta)
set.seed(1234)
clinical_dat$recurrence=as.factor(clinical_dat$recurrence)
clinical_dat$recurrence=droplevels(clinical_dat$recurrence)
#boruta.train <- Boruta(recurrence~., data = clinical_dat[,c(1,3:10)],doTrace = 2,pValue=1e-10)
#print(boruta.train)
#x=clinical_dat[,c(names(boruta.train$finalDecision[which(boruta.train$finalDecision %in% c("Confirmed","Tentative"))]),"recurrence")]
model=randomForest(recurrence~.,data=clinical_dat[,c(1,4:11)],ntree=50000,importance=T)
varImpPlot(model,cex=.7)
rf.roc<-roc(as.factor(clinical_dat$recurrence),model$votes[,2],ci=T,ci.alpha=0.95,stratified=F,plot=T,print.auc=T)
plot.roc(rf.roc,print.auc=T,asp=NA)
#clinical_dat$recurrence=as.logical(clinical_dat$recurrence)
dat=as.data.frame(t(vp_tcga[names(g1),]))
dat$recurrence=recurrence
dat=dat[which(dat$recurrence!="<NA>"),]
dat$recurrence=as.factor(dat$recurrence)
dat$recurrence=droplevels(dat$recurrence)
set.seed(1234)
boruta.train <- Boruta(recurrence~., data = dat,doTrace = 2)
print(boruta.train)
x=dat[,c(names(boruta.train$finalDecision[which(boruta.train$finalDecision %in% c("Confirmed"))]),"recurrence")]
model=randomForest(recurrence~.,data=x,ntree=50000,importance=T)
varImpPlot(model,cex=.7)
rf.roc<-roc(as.factor(clinical_dat$recurrence),model$votes[,2],ci=T,ci.alpha=0.95,stratified=F,plot=T,print.auc=T)
plot.roc(rf.roc,print.auc=T,asp=NA)
plotdata=t(x[,setdiff(colnames(x),"recurrence")])
col_anno=data.frame(recurrence=x$recurrence)
colnames(plotdata)=1:ncol(plotdata)
plotdata=plotdata[,order(x$recurrence)]
plotdata=t(apply(plotdata,1,function(x){(x-mean(x))/sd(x)}))
pheatmap(plotdata,annotation_col = col_anno,cluster_cols = F,show_colnames = F,cluster_rows = F)
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
mat_breaks <- quantile_breaks(as.matrix(plotdata), n = 30)
pheatmap(plotdata,annotation_col = col_anno,cluster_cols = F,show_colnames = F,cluster_rows = F,breaks=mat_breaks,color = colorRampPalette(colors = c('blue', 'white', 'red'))(length(mat_breaks)))
#####
res.cox1 <- coxph(SurvObj ~ ., data =  clinical_dat[,c(4:11,14)])
summary(res.cox1)
plot(ggforest(res.cox1, data = clinical_dat,fontsize = 1,main = "GSEA Hazard Ratio"))
#GSEA0
res.cut <- surv_cutpoint(clinical_dat, time = "Recurrence", event = "OS",variables = c("GSEA0"))
summary(res.cut)
plot(res.cut, "GSEA0", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(Recurrence,OS) ~GSEA0, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Recurrence-Free Survival",xlab="Time (days)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))
#GSEA1
res.cut <- surv_cutpoint(clinical_dat, time = "Recurrence", event = "OS",variables = c("GSEA1"))
summary(res.cut)
plot(res.cut, "GSEA1", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(Recurrence,OS) ~GSEA1, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Recurrence-Free Survival",xlab="Time (days)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))
#GSEA2
res.cut <- surv_cutpoint(clinical_dat, time = "Recurrence", event = "OS",variables = c("GSEA2"))
summary(res.cut)
plot(res.cut, "GSEA2", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(Recurrence,OS) ~GSEA2, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Recurrence-Free Survival",xlab="Time (days)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))
#GSEA3
res.cut <- surv_cutpoint(clinical_dat, time = "Recurrence", event = "OS",variables = c("GSEA3"))
summary(res.cut)
plot(res.cut, "GSEA3", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(Recurrence,OS) ~GSEA3, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Recurrence-Free Survival",xlab="Time (days)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))
#GSEA4
res.cut <- surv_cutpoint(clinical_dat, time = "Recurrence", event = "OS",variables = c("GSEA4"))
summary(res.cut)
plot(res.cut, "GSEA4", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(Recurrence,OS) ~GSEA4, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Recurrence-Free Survival",xlab="Time (days)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))
#GSEA5
res.cut <- surv_cutpoint(clinical_dat, time = "Recurrence", event = "OS",variables = c("GSEA5"))
summary(res.cut)
plot(res.cut, "GSEA5", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(Recurrence,OS) ~GSEA5, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Recurrence-Free Survival",xlab="Time (days)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))
#GSEA6
res.cut <- surv_cutpoint(clinical_dat, time = "Recurrence", event = "OS",variables = c("GSEA6"))
summary(res.cut)
plot(res.cut, "GSEA6", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(Recurrence,OS) ~GSEA6, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Recurrence-Free Survival",xlab="Time (days)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))
#GSEA7
res.cut <- surv_cutpoint(clinical_dat, time = "Recurrence", event = "OS",variables = c("GSEA7"))
summary(res.cut)
plot(res.cut, "GSEA7", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(Recurrence,OS) ~GSEA7, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Recurrence-Free Survival",xlab="Time (days)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))









#OncoTreat/OncoTarget

##normalize to wholetcga
library(n1platform)
tcga=getRawCounts(tissue="all",type="tissue")
tcga=geneExpression_Normalization(tcga, "entrez", "TPM")
param <- list()
param[["file"]] <- list("character", NA)
param[["tissue"]] <- list("character", "infer")
param[["patient"]] <- list("character", "ND")
param[["wd"]] <- list("character", getwd())
param[["analysis"]] <- list("character", c("target", "treat"))
param[["results"]] <- list("character", "summary")
param[["batchSource"]] <- list("character", "ngs", c("ngs", "pgm", "avera", "tempus", "mdacc"))
param[["batchTCGA"]] <- list("logical", FALSE, c("TRUE", "FALSE", "T", "F"))
param[["pdx"]] <- list("logical", FALSE, c("TRUE", "T", "FALSE", "F"))
param[["activity.adjust"]] <- list("character", "bonferroni", c("auto", "none", "fdr", "bonferroni", "holm", "hochberg", "hommel", "BH", "BY"))
param[["activity.pval"]] <- list("numeric", 1e-5)
param[["match.adjust"]] <- list("character", "bonferroni", c("auto", "none", "fdr", "bonferroni", "holm", "hochberg", "hommel", "BH", "BY"))
param[["match.pval"]] <- list("numeric", .01)
param[["synergy.adjust"]] <- list("character", "bonferroni", c("auto", "none", "fdr", "bonferroni", "holm", "hochberg", "hommel", "BH", "BY"))
param[["synergy.pval"]] <- list("numeric", 1e-10)
param[["synergy.max.pairs"]] <- list("integer", 100)
param[["target.score"]] <- list("numeric", 0)
param[["samples"]] <- list("character", NA)
param[["restrictNetwork"]] <- list("logical", TRUE, c("TRUE", "FALSE", "T", "F"))
param[["networkThreshold"]] <- list("numeric", .1)
param[["modelMatch"]] <- list("character", "tissue", c("tissue", "mr", "any"))
param[["restrictModel"]] <- list("logical", TRUE, c("TRUE", "FALSE", "T", "F"))
param[["minDrugsOncoTreat"]] <- list("integer", 0)
param[["concentration"]] <- list("logical", FALSE)
param[["nononcology"]] <- list("logical", FALSE, c("TRUE", "FALSE", "T", "F"))
param[["report"]] <- list("logical", TRUE, c("TRUE", "FALSE", "T", "F"))
param[["removeSymbols"]] <- list("logical", TRUE, c("TRUE", "FALSE", "T", "F"))
param[["meta"]] <- list("character", "")
param[["maxdrugs"]] <- list("integer", 5)
param[["cores"]] <- list("integer", 1)
param[["clean"]] <- list("logical", TRUE, c("TRUE", "FALSE", "T", "F"))
param[["version"]] <- list("character", "1-3-3")
param[["scriptdir"]] <- list("charcater", file.path(find.package("n1platform"), "scripts"))
param[["masterdir"]] <- list("character", file.path(find.package("n1platform"), "scripts"))
param[["tcgaScale"]] <- list("logical", FALSE, c("TRUE", "FALSE", "T", "F"))
param <- parseCommandArgs(commandArgs(TRUE), params=param)
param$file <- "primecut_pt7_tumor_viper.rda"
param$tissue <- "prad"
param$analysis <- c("target","treat")
param$batchSource <- "ngs"
param$batchTCGA <- FALSE
param$pdx <- FALSE
param$synergy.pval <- 0
param$report <- FALSE
param$results <- "full"
param$target.score <- 2
#single-cell oncotarget
vp=readRDS("~/Documents/Documents/MED SCHOOL/summer 2019/primecut.integrated.meta.rds")
vp.t=readRDS("~/Documents/Documents/MED SCHOOL/summer 2019/primecut.integrated.tumorcells.vp.rds")
vp=vp[,colnames(vp.t)]
meta=vp@assays$SCT@counts
meta=meta[which(rowSums(meta)>0),]
ref=log10(tcga+1)
rownames(ref)=entrez2gene(rownames(ref))
ref=ref[intersect(rownames(ref),rownames(meta)),]
meta=meta[rownames(ref),]
dat=log10(apply(meta,2,function(x){(x/sum(x))*1000000})+1)
medians=apply(ref,1,median)
mads=apply(ref,1,mad)
dat=(dat-medians)/mads
dat=dat[which(mads>0),]
load("/Users/aleksandar/Documents/Documents/MED SCHOOL/summer 2019/single-cell-oncotreat/n1platform/interactomes/data/prad_regulon.rda")
#dat=names2entrez(dat)
x=AnnotationDbi::select(org.Hs.eg.db, keys=rownames(dat), columns=c("ENTREZID"), keytype="SYMBOL")
x=x[which(!is.na(x[,2])),]
entrezids=c()
exclude=c()
include=c()
for(i in rownames(data)){
  ids=x[which(x[,1]==i),2]
  if(length(ids)==0){
    exclude=c(exclude,i)
  }
  else{
    entrezids=c(entrezids,ids[1])
    include=c(include,i)
  }
}
dat=dat[include,]
rownames(dat)=entrezids
meta_tumor_vp<-viper(dat, regul, method = 'none')
###
meta_tumor_vp=as.matrix(vp.t@assays$RNA@counts)
x=AnnotationDbi::select(org.Hs.eg.db, keys=rownames(meta_tumor_vp), columns=c("ENTREZID"), keytype="SYMBOL")
x=x[which(!is.na(x[,2])),]
entrezids=c()
exclude=c()
include=c()
for(i in rownames(meta_tumor_vp)){
  ids=x[which(x[,1]==i),2]
  if(length(ids)==0){
    exclude=c(exclude,i)
  }
  else{
    entrezids=c(entrezids,ids[1])
    include=c(include,i)
  }
}
meta_tumor_vp=meta_tumor_vp[include,]
rownames(meta_tumor_vp)=entrezids
###
dset <- meta_tumor_vp
save(dset, file = param$file)
if (param[["patient"]]=="ND") {
  tmp <- strsplit(param[["file"]], "/")[[1]]
  tmp <- sub("_counts.rda", "", tmp[length(tmp)]) # Obtain the patient ID from the file name
  param[["patient"]] <- tmp
}
param[["patient"]] <- gsub("[[:punct:]]", "-", param[["patient"]])
param[["wd"]] <- normalizePath(param[["wd"]])
param[["file"]] <- normalizePath(param[["file"]])
#####ANALYSIS
results <- darwinLoadRawCounts(param) # Loading the RawCount files
drugs=darwinDrugsTable(param)

plot.data <- meta_tumor_vp
plot.data=plot.data[intersect(as.character(unique(drugs$target)),rownames(plot.data)),]
quantile_breaks <- function(xs, n = 10){
  breaks <- quantile(xs, probs = seq(from = 0, to = 1, length.out = n))
  unique(breaks)
}
mat_breaks <- quantile_breaks(plot.data,n = 100)

annot.col=data.frame(cluster=vp.t$seurat_clusters,treatment=vp.t$treatment)
o=order(annot.col$cluster,decreasing = F)
annot.col=annot.col[o,]
plot.data=plot.data[,o]
integrated_scores=apply(plot.data,1,function(x){sum(x)/sqrt(length(x))})
plot.data=plot.data[order(integrated_scores,decreasing = T),]
pheatmap(plot.data, fontsize_row = 12, fontsize_col = 8, show_rownames = FALSE, show_colnames = F, breaks = mat_breaks, cluster_cols = FALSE, cluster_rows = TRUE, annotation_col = annot.col)

# total oncoprotein -log10 bonferroni adjusted p-value
plot.data <- pnorm(q = as.matrix(plot.data), lower.tail = FALSE)
for(i in 1:ncol(plot.data)){
  padj.res <- plot.data[,i]*nrow(plot.data)
  padj.res[which(padj.res>1)]=1
  plot.data[,i] <- padj.res
}
plot.data <- -1*log(plot.data,base = 10)
#plot.data <- plot.data[order(apply(plot.data,1,max), decreasing = TRUE),]
mat_breaks <- quantile_breaks(plot.data,n = 100)
mat_breaks=c(0,mat_breaks[which(mat_breaks>2)])
myColor <- colorRampPalette(c("white", "red"))(length(mat_breaks))
pheatmap(plot.data, color = myColor,fontsize_row = 12, fontsize_col = 8, show_rownames = FALSE, show_colnames = F, cluster_cols = FALSE, cluster_rows = TRUE, annotation_col = annot.col)

old.names=rownames(plot.data)
drug.names <- lapply(old.names,function(x,ref){
  y <- as.character(ref[grep(x,as.character(ref$target)),"drug"])
  return(y)
}, ref = drugs)
rownames(plot.data) <- entrez2gene(rownames(plot.data))
names(drug.names) <- rownames(plot.data)
drug.names <- lapply(drug.names,function(x,max.val){
  y <- c(x,rep("_",times = max.val - length(x)))
  return(y)
}, max.val = max(sapply(drug.names,length, USE.NAMES = FALSE)))
drug.result.mat <- matrix(data = unlist(drug.names, use.names = FALSE), nrow = length(drug.names), ncol = max(sapply(drug.names,length, USE.NAMES = FALSE)), byrow = TRUE)
rownames(drug.result.mat) <- rownames(plot.data)
colnames(drug.result.mat) <- paste("Drug_",1:ncol(drug.result.mat),sep = "")
plot.data <- signif(plot.data,4)
final.oncotarget.output.mat <- cbind(as.data.frame(plot.data),as.data.frame(drug.result.mat))
write.table(final.oncotarget.output.mat, file = "oncotarget_output_mat_normalizedtotcga_pt5.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

x=apply(plot.data,1,function(x){c(median(x[which(annot.col[,1]==0)]),median(x[which(annot.col[,1]==1)]),median(x[which(annot.col[,1]==2)]),median(x[which(annot.col[,1]==3)]),median(x[which(annot.col[,1]==4)]),median(x[which(annot.col[,1]==5)]),median(x[which(annot.col[,1]==6)]),median(x[which(annot.col[,1]==7)]))})
keep.index=apply(x,2,function(x){x[1]>2 | x[2]>2 | x[3]>2 | x[4]>2 | x[5]>2 | x[6]>2 | x[7]>2 | x[8]>2})
plot.data <- plot.data[keep.index,]
#plot.data <- plot.data[names(apply(plot.data,1,max))[order(apply(plot.data,1,max), decreasing = T)],]
pheatmap(plot.data,breaks=mat_breaks, color = myColor,fontsize_row = 12, fontsize_col = 8, show_rownames = T, show_colnames = F, cluster_cols = FALSE, cluster_rows = FALSE, annotation_col = annot.col)

oncotarget_druglist=final.oncotarget.output.mat[rownames(plot.data),grepl("Drug",colnames(final.oncotarget.output.mat))]



#tumorcells.vp.druggable=tumorcells.vp[intersect(rownames(tumorcells.vp),entrez2gene(drugs$target)),]
#tumorcells.vp.markers.druggable <- FindAllMarkers(tumorcells.vp.druggable, only.pos = TRUE, min.pct = 0, logfc.threshold = 0,test.use = "t")
#top10 <- tumorcells.vp.markers.druggable %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#geneHeatmap_primecut(tumorcells.vp.druggable,tumorcells.vp.druggable$seurat_clusters,top10$gene,viper = T,scaled = T)
#saveRDS(tumorcells.vp.druggable,"~/Documents/Documents/MED SCHOOL/summer 2019/primecut.integrated.tumorcells.vp.druggable.rds")



###RUN SINGLE_CELL ONCOTREAT: LNCAP, DU145
vp.t=readRDS("~/Documents/Documents/MED SCHOOL/summer 2019/primecut.integrated.tumorcells.vp.rds")
tumorcells=readRDS("~/Documents/Documents/MED SCHOOL/summer 2019/primecut.integrated.rds")
tumorcells=tumorcells[,colnames(vp.t)]
load("/Users/aleksandar/Documents/Documents/MED SCHOOL/summer 2019/single-cell-oncotreat/n1platform/interactomes/data/prad_regulon.rda")
interactome=regul
vp.t.tcga<-viper(eset=names2entrez(tumorcells@assays$integrated@scale.data), regulon=interactome, method="none", verbose = FALSE)

#tccsup bladder
#lncap prostate adeno (castrate-sensitive), du145 prostate (castrate-resistant)
setwd("/Users/aleksandar/Documents/Documents/MED SCHOOL/summer 2019/single-cell-oncotreat/weighted_oncotreat_clone/Weighted_OncoTreat")
#source('~/Documents/Documents/MED SCHOOL/summer 2019/single-cell-oncotreat/Functions_WeightedOncoTreat.R')
#Create a directory in which all the results will be saved
Results_directory<-"PRIME-CUT_combined_Du145_vs_LNCAP"
dir.create(Results_directory)

#Estimation of the protein activity of the baseline for the two PRAD cell lines for which perturbational
#data are available in the n1platform: LNCAP and DU145
library(drugsignatures)
clineByTissue("prad")
CL_names=c("gen_DU145","gen_LNCAP")
#getDrugsVipermat("du145")
total_matrix<-getRawCounts(tissue="all",type="cline")
total_matrix<-geneExpression_Normalization(total_matrix, "entrez", "TPM")
index_CCLE_cell_lines<-which(colnames(total_matrix) %in% CL_names)
total_matrix<-cbind(total_matrix[,index_CCLE_cell_lines], total_matrix[,-c(index_CCLE_cell_lines)]) 
#x=AnnotationDbi::select(org.Hs.eg.db, keys=rownames(total_matrix), columns=c("SYMBOL"), keytype="ENTREZID")
#rownames(total_matrix)=x[,2]
#Double rank signature
rank_exp <- apply(total_matrix, 2, rank)
median_exp <- apply(rank_exp, 1, median)
mad_exp <- apply(rank_exp, 1, mad)
signature_exp <- (rank_exp - median_exp)/mad_exp
#Protein activity estimation
load("/Users/aleksandar/Documents/Documents/MED SCHOOL/summer 2019/single-cell-oncotreat/n1platform/interactomes/data/prad_regulon.rda")
interactome=regul
#filenames <- list.files("/Users/aleksandar/Documents/Documents/MED SCHOOL/summer 2019/prime-cut-metacells/sc_nets/", pattern="*.rds", full.names=TRUE)
#nets=lapply(filenames,readRDS)
#output<-viper(eset=signature_exp[,1:2], regulon=nets, method="none", verbose = FALSE)
output<-viper(eset=signature_exp[,1:2], regulon=interactome, method="none", verbose = FALSE)
colnames(output)=unlist(lapply(strsplit(colnames(total_matrix)[1:2],split="_"),function(x){x[2]}))
#CLs_protein_activity_names=output
#CLs_protein_activity_entrez<-names2entrez(CLs_protein_activity_names)
CLs_protein_activity_entrez=output
CLs_protein_activity_names=entrez2names(CLs_protein_activity_entrez)
CLs_protein_activity_entrez_TF_coTF<-TF_coTF_restriction(CLs_protein_activity_entrez, "entrez")
CLs_protein_activity_names_TF_coTF<-TF_coTF_restriction(CLs_protein_activity_names, "name")
Heatmap_CLs_protein_activity(CLs_protein_activity_names_TF_coTF, 30, path=Results_directory)

#patients_entrez_proteinActivity_TF_coTF<- as.matrix(vp.t@assays$RNA@counts)
#x=AnnotationDbi::select(org.Hs.eg.db, keys=rownames(patients_entrez_proteinActivity_TF_coTF), columns=c("ENTREZID"), keytype="SYMBOL")
#rownames(patients_entrez_proteinActivity_TF_coTF)=x[,2]
patients_entrez_proteinActivity_TF_coTF=vp.t.tcga
patients_entrez_proteinActivity_TF_coTF=TF_coTF_restriction(patients_entrez_proteinActivity_TF_coTF, "entrez")
#GSEA: computing the Master Regulators overlapping score between each cell in the dataset and the two cell lines
n_MRs<-50
enrichment_results<-enrichment_analysis(patients_entrez_proteinActivity_TF_coTF, CLs_protein_activity_entrez_TF_coTF, n_MRs)
#Plot the results of GSEA
GSEA_plot(enrichment_results, colnames(output), path = Results_directory, "WholeDataset")

#Compute the weights and the confidence score
threshold<-0.05
w_cf<-weights_cf(enrichment_results, threshold, path = Results_directory)
weights<-w_cf$weights
confidence_score<-w_cf$FDR

#Check the level of the MR overlap after the filtering process
enrichment_results_filtered<-enrichment_results[confidence_score<threshold,]
GSEA_plot(enrichment_results_filtered, colnames(output), path = Results_directory, "Filtered")

##Plot VIPER clustering with GSEA enrichment
#patients_names_proteinActivity_TF_coTF<-TF_coTF_restriction(as.matrix(vp.t@assays$RNA@counts),"name")
patients_names_proteinActivity_TF_coTF<-TF_coTF_restriction(entrez2names(vp.t.tcga),"name")
Heatmap_clusters_withgenexp(vp.t,vp.t$seurat_clusters,top10,enrichment_results,confidence_score,path=Results_directory)

#Application of OncoTreat
cell_line <- c("cline=du145")
enrichment<-enrichment_results[,1]
names(enrichment)<-rownames(enrichment_results)
OncoTreat_DU145<-OncoTreat(cell_line, patients_entrez_proteinActivity_TF_coTF, enrichment, "Blues", path = Results_directory)

cell_line <- c("cline=lncap")
enrichment<-enrichment_results[,2]
names(enrichment)<-rownames(enrichment_results)
OncoTreat_LNCAP<-OncoTreat(cell_line, patients_entrez_proteinActivity_TF_coTF, enrichment, "Greens", path = Results_directory)

#Selection of the top drugs using the weigths
wg<-weights[,1]
names(wg)<-rownames(weights)
enrichment<-enrichment_results[,1]
names(enrichment)<-rownames(enrichment_results)
top_drugs_DU145<-top_drugs(OncoTreat_DU145, confidence_score, wg, 25, enrichment, threshold, "Blues", "DU145", path = Results_directory)

pVal_DU145<-vector()
for(i in 1:ncol(top_drugs_DU145)){
  pVal_DU145[i]<-sum(top_drugs_DU145[,i])/sqrt(nrow(top_drugs_DU145))
}
names(pVal_DU145)<-colnames(top_drugs_DU145)

wg<-weights[,2]
names(wg)<-rownames(weights)
enrichment<-enrichment_results[,2]
names(enrichment)<-rownames(enrichment_results)
top_drugs_LNCAP<-top_drugs(OncoTreat_LNCAP, confidence_score, wg, 25, enrichment, threshold, "Greens", "LNCAP", path = Results_directory)

pVal_LNCAP<-vector()
for(i in 1:ncol(top_drugs_LNCAP)){
  pVal_LNCAP[i]<-sum(top_drugs_LNCAP[,i])/sqrt(nrow(top_drugs_LNCAP))
}
names(pVal_LNCAP)<-colnames(top_drugs_LNCAP)

library(openxlsx)
write.xlsx(pVal_DU145, "Drugs_Pobor.xlsx", sheetName = "DU145")
write.xlsx(pVal_LNCAP, "Drugs_Pobor.xlsx", sheetName = "LNCAP", append = TRUE)

final_plot(top_drugs_DU145, top_drugs_LNCAP, enrichment_results, c("Blues", "Greens"), path = Results_directory,cluster=vp.t$seurat_clusters[colnames(vp.t.tcga)])



oncotarget_druglist
x=unlist(oncotarget_druglist)
x=x[which(x!="_")]
unlist(lapply(strsplit(colnames(top_drugs_DU145),"_"),function(x){x[1]}))
unlist(lapply(strsplit(colnames(top_drugs_LNCAP),"_"),function(x){x[1]}))



###add oncotreat with tccsup

####single cell line oncotreat
vp=as.matrix(vp.t@assays$RNA@counts)
tumorcells=readRDS("~/Documents/Documents/MED SCHOOL/summer 2019/primecut.integrated.rds")
tumorcells=tumorcells[,colnames(vp.t)]
ctc_mat=as.matrix(tumorcells@assays$SCT@counts)
#clineMatch() to bulkRNASeq datasets schadendorf and timchan
library(n1platform)
vp_entrez=names2entrez(vp)
x=AnnotationDbi::select(org.Hs.eg.db, keys=rownames(ctc_mat), columns=c("ENTREZID"), keytype="SYMBOL")
g=rownames(ctc_mat)
entrez=g
for(i in 1:length(g)){
  entrez[i]=x[which(x[,1]==g[i]),2][1]
}
rownames(ctc_mat)=entrez
ctc_mat_entrez=ctc_mat
ctc_clinematch=clineMatch(ctc_mat_entrez, vipermat = vp_entrez,verbose = T)

cls=intersect(colnames(ctc_clinematch$mr),toupper(clineByTissue(regulonSets())))
p=pheatmap(ctc_clinematch$mr[which(vp.t$seurat_clusters=="5"),cls],show_rownames = F)
cell_line_matches=apply(ctc_clinematch$mr[,cls],1,function(x){colnames(ctc_clinematch$mr[,cls])[order(x,decreasing = T)[1]]})
x=as.data.frame(table(unlist(cell_line_matches)))
x=x[order(x[,2],decreasing = T),]
head(x)
barplot(x[,2],names.arg = x[,1],ylab="# samples best-matched to this cell line",las=2)

x=OncoTreat(cell_line="cline=tccsup",vp_entrez[,which(vp.t$seurat_clusters=="5")],enrichment = ctc_clinematch$mr[which(vp.t$seurat_clusters=="5"),"TCCSUP"],color="Greens",path="/Users/aleksandar/Documents/Documents/MED SCHOOL/summer 2019/published_melanoma_datasets/pd1_response_publisheddata")
#x=OncoTreat(cell_line="cline=lncap",vp_entrez,enrichment = ctc_clinematch$mr[,"LNCAP"],color="Greens",path="/Users/aleksandar/Documents/Documents/MED SCHOOL/summer 2019/published_melanoma_datasets/pd1_response_publisheddata")
top_drugs_unweighted=function(drug_activity, n_drugs, enrichment, color, cell_line_name, path){
  adjust <- "bonferroni"
  pval1 <- t(apply(drug_activity, 1, function(x, adjust) -log10(p.adjust(pnorm(x, lower.tail=TRUE), adjust)), adjust=adjust))
  top_drugs<-vector()
  for (i in 1:ncol(drug_activity))
    top_drugs[i]<-sum(drug_activity[,i])/sqrt(nrow(drug_activity))
  names(top_drugs)<-colnames(pval1)
  top_drugs_names<-names(sort(top_drugs, decreasing = FALSE)[1:min(n_drugs,length(top_drugs))])
  require(pheatmap) 
  require(RColorBrewer)
  annotation_enrichment=data.frame(
    enrichment_CL=enrichment)
  colnames(annotation_enrichment)<-cell_line_name
  annotation_colors = list(CL = brewer.pal(9,color))
  names(annotation_colors)<-cell_line_name
  cols<-colorRampPalette(c("gray97", "gray97", "red"))(10)
  myBreaks <- c(seq(0, 2, length.out=8), seq(2.1, max(pval1), length.out=2))
  file=paste("/topDrugs", cell_line_name, sep="_")
  file=paste(path, file, sep='')
  p<-pheatmap(t(pval1[(sort(annotation_enrichment[,1], index.return=TRUE, decreasing=TRUE))$ix,top_drugs_names]), show_colnames = FALSE, show_rownames=TRUE, cluster_cols = FALSE, cluster_rows = FALSE, fontsize_row = 15, color=cols, fontsize=14, legend=T, annotation_col = annotation_enrichment, annotation_colors = annotation_colors, width=15, breaks = myBreaks)
  pdf(file= paste(file, "pdf", sep="."), height = 10, width=15)
  print(p)
  no_show<-dev.off()
  return(pval1[,top_drugs_names])
}
adjust <- "BH"
pval1 <- t(apply(x, 1, function(x,adjust) -log10(p.adjust(pnorm(x, lower.tail=TRUE), adjust)), adjust=adjust))
which(apply(pval1,2,max)>2)
top_drugs_unweighted(x[,which(apply(pval1,2,max)>2)],length(which(apply(pval1,2,max)>2)),enrichment = ctc_clinematch$mr[which(vp.t$seurat_clusters=="5"),"TCCSUP"],color="Greens",cell_line_name = "tccsup",path="/Users/aleksandar/Documents/Documents/MED SCHOOL/summer 2019/published_melanoma_datasets/pd1_response_publisheddata")







###de bono collaboration plots


#####distribuiton of B7-H3
tumorcells.vp=readRDS("~/Documents/Documents/MED SCHOOL/summer 2019/primecut.integrated.tumorcells.vp.rds")
primecut.integrated.meta.vp=readRDS("~/Documents/Documents/MED SCHOOL/summer 2019/primecut.integrated.meta.vp.rds")
primecut.integrated=readRDS("~/Documents/Documents/MED SCHOOL/summer 2019/primecut.integrated.rds")
primecut.integrated.meta.vp$condensed_clusters=as.character(primecut.integrated.meta.vp$seurat_clusters)
primecut.integrated.meta.vp$condensed_clusters[which(primecut.integrated.meta.vp$seurat_clusters %in% c("CD4 T-cell 1","CD4 T-cell 2"))]="CD4 T-cell"
primecut.integrated.meta.vp$condensed_clusters[which(primecut.integrated.meta.vp$seurat_clusters %in% c("CD8 T-cell 1","CD8 T-cell 2","CD8 T-cell 3"))]="CD8 T-cell"
primecut.integrated.meta.vp$condensed_clusters[which(primecut.integrated.meta.vp$seurat_clusters %in% c("Treg 1","Treg 2","Treg 3"))]="Treg"
primecut.integrated.meta.vp$condensed_clusters[which(primecut.integrated.meta.vp$seurat_clusters %in% c("B-cells 1","B-cells 2","B-cells 3","B-cells 4","B-cells 5"))]="B-cell"
primecut.integrated.meta.vp$condensed_clusters[which(primecut.integrated.meta.vp$seurat_clusters %in% c("Monocytes 1","Monocytes 2","Monocytes 3"))]="Monocytes"
primecut.integrated.meta.vp$condensed_clusters[which(primecut.integrated.meta.vp$seurat_clusters %in% c("Epithelial 1","Epithelial 2","Epithelial 3"))]="Tumor"
primecut.integrated.meta.vp$condensed_clusters[which(primecut.integrated.meta.vp$seurat_clusters %in% c("Erythrocytes 1"))]="Erythrocytes"
table(primecut.integrated.meta.vp$condensed_clusters)

#Interrogate expression on tumor cells to produce bar graph 
b7h3="CD276"
b7genes_tumor=c("CD80","CD86","CD274","ICOSLG","PDCD1LG2","VTCN1","VISTA","VSIR","NCR3LG1")

x=primecut.integrated@assays$SCT@counts[,colnames(primecut.integrated.meta.vp)]
dat=as.data.frame(t(x[intersect(rownames(x),c(b7h3,b7genes_tumor)),]))
dat$celltype=primecut.integrated.meta.vp$condensed_clusters
dat=dat[which(dat$celltype=="Tumor"),]
y=log10(apply(dat[,1:8],2,sum)+1)
y=y[order(y,decreasing=T)]
y=data.frame(Celltype="Tumor",Gene=names(y),Expression=y)
y$Gene=factor(y$Gene,levels=y$Gene)
p=ggplot(y, aes(x=Gene, y=Expression,group=Gene,fill=Gene)) +
  geom_bar(stat="identity",position=position_dodge(),color="black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Gene Expression in Tumor Cells")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))

cd28genes_tumor=c("BTLA","ICOS","CTLA4","PDCD1","CD28","HVEM","TIGIT","IDO1","TIM3","LAG3","ADORA2A","ADORA2B")
cd28genes_tumor_v2=c("BTLA","ICOS","CTLA4","PDCD1","CD28","HVEM","TIGIT","IDO1","TIM3","LAG3")

x=primecut.integrated@assays$SCT@counts[,colnames(primecut.integrated.meta.vp)]
dat=as.data.frame(t(x[intersect(rownames(x),cd28genes_tumor),]))
dat$celltype=primecut.integrated.meta.vp$condensed_clusters
dat=dat[which(dat$celltype=="Tumor"),]
y=log10(apply(dat[,1:10],2,sum)+1)
y=y[order(y,decreasing=T)]
y=data.frame(Celltype="Tumor",Gene=names(y),Expression=y)
y$Gene=factor(y$Gene,levels=y$Gene)
p=ggplot(y, aes(x=Gene, y=Expression,group=Gene,fill=Gene)) +
  geom_bar(stat="identity",position=position_dodge(),color="black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Gene Expression in Tumor Cells")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))

x=primecut.integrated@assays$SCT@counts[,colnames(primecut.integrated.meta.vp)]
dat=as.data.frame(t(x[intersect(rownames(x),cd28genes_tumor_v2),]))
dat$celltype=primecut.integrated.meta.vp$condensed_clusters
dat=dat[which(dat$celltype=="Tumor"),]
y=log10(apply(dat[,1:8],2,sum)+1)
y=y[order(y,decreasing=T)]
y=data.frame(Celltype="Tumor",Gene=names(y),Expression=y)
y$Gene=factor(y$Gene,levels=y$Gene)
p=ggplot(y, aes(x=Gene, y=Expression,group=Gene,fill=Gene)) +
  geom_bar(stat="identity",position=position_dodge(),color="black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Gene Expression in Tumor Cells")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))

x=primecut.integrated@assays$SCT@counts[,colnames(primecut.integrated.meta.vp)]
dat=as.data.frame(t(x[intersect(rownames(x),c(b7h3,b7genes_tumor,cd28genes_tumor)),]))
dat$celltype=primecut.integrated.meta.vp$condensed_clusters
dat=dat[which(dat$celltype=="Tumor"),]
y=log10(apply(dat[,1:18],2,sum)+1)
y=y[order(y,decreasing=T)]
y=data.frame(Celltype="Tumor",Gene=names(y),Expression=y)
y$Gene=factor(y$Gene,levels=y$Gene)
p=ggplot(y, aes(x=Gene, y=Expression,group=Gene,fill=Gene)) +
  geom_bar(stat="identity",position=position_dodge(),color="black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Gene Expression in Tumor Cells")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))

x=primecut.integrated@assays$SCT@counts[,colnames(primecut.integrated.meta.vp)]
dat=as.data.frame(t(x[intersect(rownames(x),c(b7h3,b7genes_tumor,cd28genes_tumor_v2)),]))
dat$celltype=primecut.integrated.meta.vp$condensed_clusters
dat=dat[which(dat$celltype=="Tumor"),]
y=log10(apply(dat[,1:16],2,sum)+1)
y=y[order(y,decreasing=T)]
y=data.frame(Celltype="Tumor",Gene=names(y),Expression=y)
y$Gene=factor(y$Gene,levels=y$Gene)
p=ggplot(y, aes(x=Gene, y=Expression,group=Gene,fill=Gene)) +
  geom_bar(stat="identity",position=position_dodge(),color="black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Gene Expression in Tumor Cells")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))



#Immune cell expression of putative B7H3 receptors and CD28 family receptors
genes_immune=c("BTLA","ICOS","CTLA4","PDCD1","CD28","TIGIT","NCR30","NECTIN2","IL20RA","TGFB1","TGFB2")
primecut.integrated$condensed_clusters=primecut.integrated.meta.vp$condensed_clusters[colnames(primecut.integrated)]
VlnPlot(primecut.integrated,"TGFB2",group.by="condensed_clusters",sort = T,pt.size=0,assay="SCT")+NoLegend()

#Do immune cell expression (depending on which subset they are predominant on) of the markers 
#highlighted in yellow correlate with tumor epithelial b7H3 expression
genes_immune_yellowhighlight=c("TIGIT","NECTIN2","TGFB1")
table(primecut.integrated$patient)
x=primecut.integrated@assays$SCT@counts[,colnames(primecut.integrated.meta.vp)]
b7h3=x["CD276",which(primecut.integrated.meta.vp$condensed_clusters=="Tumor")]
dat=data.frame(b7h3,patient=primecut.integrated$patient[which(primecut.integrated.meta.vp$condensed_clusters=="Tumor")])
g=aggregate(dat$b7h3,list(dat$patient),mean)
tigit=x["TIGIT",which(primecut.integrated.meta.vp$condensed_clusters=="Treg")]
dat=data.frame(tigit,patient=primecut.integrated$patient[which(primecut.integrated.meta.vp$condensed_clusters=="Treg")])
g2=aggregate(dat$tigit,list(dat$patient),mean)
nectin2=x["NECTIN2",which(primecut.integrated.meta.vp$condensed_clusters=="Endothelial")]
dat=data.frame(nectin2,patient=primecut.integrated$patient[which(primecut.integrated.meta.vp$condensed_clusters=="Endothelial")])
g3=aggregate(dat$nectin2,list(dat$patient),mean)
tgfb1=x["TGFB1",which(primecut.integrated.meta.vp$condensed_clusters=="Monocytes")]
dat=data.frame(tgfb1,patient=primecut.integrated$patient[which(primecut.integrated.meta.vp$condensed_clusters=="Monocytes")])
g4=aggregate(dat$tgfb1,list(dat$patient),mean)

my_data=data.frame(b7h3=g[,2],tigit=g2[,2])
library("ggpubr")
ggscatter(my_data, x = "b7h3", y = "tigit", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Tumor Cell B7H3", ylab = "Treg TIGIT")
my_data=data.frame(b7h3=g[,2],nectin2=g3[,2])
library("ggpubr")
ggscatter(my_data, x = "b7h3", y = "nectin2", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Tumor Cell B7H3", ylab = "Endothelial NECTIN2")
my_data=data.frame(b7h3=g[,2],tgfb1=g4[,2])
library("ggpubr")
ggscatter(my_data, x = "b7h3", y = "tgfb1", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Tumor Cell B7H3", ylab = "Monocyte TGFB1")


#UPSTREAM AND DOWNSTREAM REGULATORS
tumor=c("STAT1","FOXA1","AR","BRD4","PLOD1","PLOD2","PLOD3","OLFM2","ECE1","EDC3","TMEM129","JAG1","B4GALNT4","GNPTAB","CUL7")
endothelial=c("OLFM2","ECE1") #vs tumor and endothelial B7H3

x=primecut.integrated@assays$SCT@counts[,colnames(primecut.integrated.meta.vp)]
b7h3=x["CD276",which(primecut.integrated.meta.vp$condensed_clusters=="Tumor")]
dat=data.frame(b7h3,patient=primecut.integrated$patient[which(primecut.integrated.meta.vp$condensed_clusters=="Tumor")])
g=aggregate(dat$b7h3,list(dat$patient),mean)
tigit=x["CUL7",which(primecut.integrated.meta.vp$condensed_clusters=="Tumor")]
dat=data.frame(tigit,patient=primecut.integrated$patient[which(primecut.integrated.meta.vp$condensed_clusters=="Tumor")])
g2=aggregate(dat$tigit,list(dat$patient),mean)
my_data=data.frame(b7h3=g[,2],tigit=g2[,2])
ggscatter(my_data, x = "tigit", y = "b7h3", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Tumor Cell CUL7", ylab = "Tumor Cell B7H3")

x=primecut.integrated@assays$SCT@counts[,colnames(primecut.integrated.meta.vp)]
b7h3=x["CD276",which(primecut.integrated.meta.vp$condensed_clusters=="Tumor")]
dat=data.frame(b7h3,patient=primecut.integrated$patient[which(primecut.integrated.meta.vp$condensed_clusters=="Tumor")])
g=aggregate(dat$b7h3,list(dat$patient),mean)
tigit=x["ECE1",which(primecut.integrated.meta.vp$condensed_clusters=="Endothelial")]
dat=data.frame(tigit,patient=primecut.integrated$patient[which(primecut.integrated.meta.vp$condensed_clusters=="Endothelial")])
g2=aggregate(dat$tigit,list(dat$patient),mean)
my_data=data.frame(b7h3=g[,2],tigit=g2[,2])
ggscatter(my_data, x = "tigit", y = "b7h3", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Endothelial Cell ECE1", ylab = "Tumor Cell B7H3")

x=primecut.integrated@assays$SCT@counts[,colnames(primecut.integrated.meta.vp)]
b7h3=x["CD276",which(primecut.integrated.meta.vp$condensed_clusters=="Endothelial")]
dat=data.frame(b7h3,patient=primecut.integrated$patient[which(primecut.integrated.meta.vp$condensed_clusters=="Endothelial")])
g=aggregate(dat$b7h3,list(dat$patient),mean)
tigit=x["ECE1",which(primecut.integrated.meta.vp$condensed_clusters=="Tumor")]
dat=data.frame(tigit,patient=primecut.integrated$patient[which(primecut.integrated.meta.vp$condensed_clusters=="Tumor")])
g2=aggregate(dat$tigit,list(dat$patient),mean)
my_data=data.frame(b7h3=g[,2],tigit=g2[,2])
ggscatter(my_data, x = "tigit", y = "b7h3", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Tumor Cell ECE1", ylab = "Endothelial Cell B7H3")

x=primecut.integrated@assays$SCT@counts[,colnames(primecut.integrated.meta.vp)]
b7h3=x["CD276",which(primecut.integrated.meta.vp$condensed_clusters=="Endothelial")]
dat=data.frame(b7h3,patient=primecut.integrated$patient[which(primecut.integrated.meta.vp$condensed_clusters=="Endothelial")])
g=aggregate(dat$b7h3,list(dat$patient),mean)
tigit=x["ECE1",which(primecut.integrated.meta.vp$condensed_clusters=="Endothelial")]
dat=data.frame(tigit,patient=primecut.integrated$patient[which(primecut.integrated.meta.vp$condensed_clusters=="Endothelial")])
g2=aggregate(dat$tigit,list(dat$patient),mean)
my_data=data.frame(b7h3=g[,2],tigit=g2[,2])
ggscatter(my_data, x = "tigit", y = "b7h3", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Endothelial Cell ECE1", ylab = "Endothelial Cell B7H3")




######Imaging Data
#primecut.integrated.meta.vp=readRDS("~/Documents/Documents/MED SCHOOL/summer 2019/primecut.integrated.meta.vp.rds")
#table(primecut.integrated.meta.vp$seurat_clusters,primecut.integrated.meta.vp$patient,primecut.integrated.meta.vp$treatment)
#hawley_qmif_data=read.csv("~/Downloads/JEH_9.9.22.csv")
#Total Cells	DAPI Positive Cells	LAG-3 (Opal 480) Positive Cells	PD-L1 (Opal 520) Positive Cells	CD8 (Opal 570) Positive Cells	CD163/CD68 (Opal 620) Positive Cells	CTLA-4 (Opal 650) Positive Cells	PanCK (Opal 690) Positive Cells	CD4 (Opal 780) Positive Cells	Negative Cells
#T cells Cells	  Tumor cells Cells   PDL1 Tumor Cells	  Tregs (CTLA-4+ Tcells) Cells	  LAG3+ CD8 Cells
#Just the indica labs rows. 8 rows. "Entire image"
#quantify T-cells, Lag3+ CD8s, and myeloids in qmIF vs scRNA-Seq
#do paired pre vs post (same patient), paired pre vs post (all paired patients), and unpaired (all patients)
# compare bone vs LN
#SID 1 on-treatment	Bone	 SP19-18199 A1	A1	YES
#SID 3 progressoin	Bone	 SP20-11188	A4	YES - big blood clot
#SID 5 baseline	LN	SP19-29780 	A1	YES
#SID 5 on-treatment	LN	SP19-31749	A2	YES
#SID 6 baseline (excluded from scRNAseq)	Liver	 SP19-29652 	A4	YES
#SID 6 on-treatment	Liver	SP20-1952	A4	YES
#SID 10 baseline 	Bone	SP20-6024 	A2	YES
#SID 11 on-treatment	LN	SP20-12651 	A1	YES
###

hawley_psa_data=read.csv("~/Downloads/hawley_psa_data.csv")
hawley_psa_data$week=round(hawley_psa_data$Timepoint/7)
colnames(hawley_psa_data)=c("Patient","Site","PSA","Date","Day","FoldChange","Week")
#primecut_list=list(P1,P3,P5,P6,P7,P8,P10,P11,P12,P13,P14,P16)
#colnames(combined)=c("pre_p1_nonresponder","pre_p10_sd","pre_p11_sd","pre_p13_sd","pre_p14_responder","pre_p3_nonresponder","pre_p5_responder","pre_p7_sd","pre_p8_responder")
#nonresponders 1,3,13
#responders 14,5,8,6
hawley_psa_data=hawley_psa_data[which(hawley_psa_data$Patient %in% c(1,3,5,6,7,8,10,11,12,13,14,16)),]
hawley_psa_data$Patient=as.factor(hawley_psa_data$Patient)
hawley_psa_data$logFoldChange=log10(hawley_psa_data$FoldChange)
hawley_psa_data=hawley_psa_data[hawley_psa_data$Week<=50,]
hawley_psa_data$logFoldChange[which(is.na(hawley_psa_data$logFoldChange))]=0
p<-ggplot(hawley_psa_data, aes(x=Week, y=logFoldChange, group=Patient)) +
  geom_line(aes(color=Patient))+
  geom_point(aes(color=Patient))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("PSA Log-Fold-Change by Patient")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+scale_color_manual(values=c("orange","orange","blue","blue","grey","blue","grey","grey","grey","orange","blue"))



























#####response to reviewers
primecut.integrated=readRDS("~/Documents/Documents/MED SCHOOL/summer 2019/primecut.integrated.rds")
primecut.integrated.meta.vp=readRDS("~/Documents/Documents/MED SCHOOL/summer 2019/primecut.integrated.meta.vp.rds")
tumorcells.vp=readRDS("~/Documents/Documents/MED SCHOOL/summer 2019/primecut.integrated.tumorcells.vp.rds")
##bulkrnaseq
#/Users/aleksandar/Dropbox/prime-cut-manuscript/bulkRNA_Seq_dataset/GSE189343_RAW.tar
#/Users/aleksandar/Dropbox/prime-cut-manuscript/bulkRNA_Seq_dataset/GSE189343_series_matrix.txt.gz
##scrnaseq
#/Users/aleksandar/Dropbox/prime-cut-manuscript/bulkRNA_Seq_dataset/singlecell/GSE143791_cell.annotation.csv.gz
#/Users/aleksandar/Dropbox/prime-cut-manuscript/bulkRNA_Seq_dataset/singlecell/GSE143791_cell.annotation.human.csv.gz
#/Users/aleksandar/Dropbox/prime-cut-manuscript/bulkRNA_Seq_dataset/singlecell/GSE143791_RAW.tar
#/Users/aleksandar/Dropbox/prime-cut-manuscript/bulkRNA_Seq_dataset/singlecell/GSE143791-GPL9052_series_matrix.txt.gz
#/Users/aleksandar/Dropbox/prime-cut-manuscript/bulkRNA_Seq_dataset/singlecell/GSE143791-GPL20301_series_matrix.txt.gz
#/Users/aleksandar/Dropbox/prime-cut-manuscript/bulkRNA_Seq_dataset/singlecell/GSE143791-GPL21103_series_matrix.txt.gz


#Bulkrnaseq treatment naïve vs short-term ADT bone mets: Use these data https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189343 
dat <- getGEO("GSE189343", GSEMatrix = TRUE)
expression=dat$GSE189343_series_matrix.txt.gz@assayData$exprs
g=dat$GSE189343_series_matrix.txt.gz@featureData@data
gse=getGEO(filename="/Users/aleksandar/Dropbox/prime-cut-manuscript/bulkRNA_Seq_dataset/GSE189343_series_matrix.txt.gz")
pheno_data=gse@phenoData@data
table(pheno_data$`treatment group:ch1`)
i=which(pheno_data$`treatment group:ch1`=="Hormone-naïve" | pheno_data$`treatment group:ch1`=="Short-term castrated")
expression=expression[,i]
pheno_data=pheno_data[i,]
###
k=which(g$gene_assignment!="---")
genes=sapply(g$gene_assignment[k],function(x){strsplit(x,"//")[[1]][2]})
genes=gsub(" ","",genes)
genes=as.character(genes)
expression=expression[k,]
#markers.vp <- FindAllMarkers(primecut.integrated.meta.vp, only.pos = TRUE, min.pct = 0, logfc.threshold = 0,test.use = "t")
#markers <- FindAllMarkers(primecut.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "wilcox")
#In bone metastases, ADT induced increased representation of myeloid cells (p = 2.8e-118 by Fisher’s Exact Test), consistent with pre-clinical data 13.  In bone lesions, ADT was also associated with a decrease in CD4 Tconv as well as in  tumor cells (p = 0.01 and 1.1e-15) (Figure 3A-B).
#CIBERSORT
Gene=genes
x=cbind(Gene,expression)
x=x[x[,1] %in% rownames(primecut.integrated),]
write.table(x,sep="\t",quote=F,row.names=F,file = "~/Downloads/newdat_primecut.txt")
dat=read.csv("~/Downloads/primecut_newdata_cibersort.csv")
pvals=apply(dat[,2:25],2,function(x){t.test(x[pheno_data$`treatment group:ch1`=="Hormone-naïve"],x[pheno_data$`treatment group:ch1`=="Short-term castrated"])$p.value})
dat=dat[,1:25]
dat$TumorType=pheno_data$`treatment group:ch1`
dat=data.frame(Patient=dat$Mixture,
               CD4nonTreg=rowSums(dat[,c("CD4.T.cell.1","CD4.T.cell.2")]),
               CD8=rowSums(dat[,colnames(dat) %in% c("CD8.T.cell.1","CD8.T.cell.2","CD8.T.cell.3")]), 
               Treg=rowSums(dat[,colnames(dat) %in% c("Treg.1", "Treg.2", "Treg.3")]),
               Bcell=rowSums(dat[,colnames(dat) %in% c("B.cells.1","B.cells.2","B.cells.3","B.cells.4","B.cells.5","Plasma.cells")]),
               Myeloid=rowSums(dat[,colnames(dat) %in% c("Monocytes.1","Monocytes.2","Monocytes.3","Macrophages","Neutrophils")]),
               Tumor=rowSums(dat[,colnames(dat) %in% c("Epithelial.1","Epithelial.2","Epithelial.3")]),
               TumorType=pheno_data$`treatment group:ch1`)
pvals=apply(dat[,2:7],2,function(x){t.test(x[pheno_data$`treatment group:ch1`=="Hormone-naïve"],x[pheno_data$`treatment group:ch1`=="Short-term castrated"])$p.value})
y=melt(dat)
colnames(y)=c("Patient","Treatment","Cluster","Abundance")
p=ggplot(y, aes(x=Cluster, y=Abundance,fill=Treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Cluster Frequency By Treatment Group")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))

#Apply classifier on the stable disease patients based on their microenvironment composition… argue that the stable disease might reduce overall effect size looking at extremes amplifies effect size … 
combined=combined[,1:9]
colnames(combined)=c("pre_p1_nonresponder","pre_p10_sd","pre_p11_sd","pre_p13_sd","pre_p14_responder","pre_p3_nonresponder","pre_p5_responder","pre_p7_sd","pre_p8_responder")
y <- melt(combined)
colnames(y)=c("cluster","treatment","frequency")
treatment=unlist(lapply(strsplit(as.character(y$treatment),"_"),function(x){x[1]}))
response=unlist(lapply(strsplit(as.character(y$treatment),"_"),function(x){x[3]}))
y$response=response
y$treatment=treatment
y$treatment.response=paste(y$treatment,y$response,sep=".")
#y=y[which(y$response!="sd"),]
#y$treatment=factor(y$treatment, levels = c("adt","combo","recurrence"))
y$cluster=as.factor(y$cluster)
p=ggplot(y, aes(x=cluster, y=frequency,fill=treatment.response)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Cluster Frequency By Treatment Group and Treatment Response")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))
pheatmap(cor(combined))
###classified by microenvironment composition
combined=combined[,1:9]
colnames(combined)=c("pre_p1_nonresponder","pre_p10_responder","pre_p11_nonresponder","pre_p13_nonresponder","pre_p14_responder","pre_p3_nonresponder","pre_p5_responder","pre_p7_responder","pre_p8_responder")
y <- melt(combined)
colnames(y)=c("cluster","treatment","frequency")
treatment=unlist(lapply(strsplit(as.character(y$treatment),"_"),function(x){x[1]}))
response=unlist(lapply(strsplit(as.character(y$treatment),"_"),function(x){x[3]}))
y$response=response
y$treatment=treatment
y$treatment.response=paste(y$treatment,y$response,sep=".")
y$cluster=as.factor(y$cluster)
p=ggplot(y, aes(x=cluster, y=frequency,fill=treatment.response)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Cluster Frequency By Treatment Group and Treatment Response")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))
apply(combined,1,function(x){t.test(x[c(1,3,4,6)],x[c(2,5,7,8,9)])$p.value})
response=unlist(lapply(strsplit(as.character(colnames(combined)),"_"),function(x){x[3]}))
dat=as.data.frame(t(combined))
dat$response=as.numeric(as.factor(response))-1
colnames(dat)=sub(" ", ".", colnames(dat))
colnames(dat)=sub(" ", ".", colnames(dat))
colnames(dat)=sub("-", ".", colnames(dat))
g=glm(response ~ CD4.T.cell.1 + CD8.T.cell.2 + Epithelial.2,dat,family="binomial")
summary(g)

#Show that the changes in response to treatment are consistent across patients
table(primecut.integrated.meta.vp$patient,primecut.integrated.meta.vp$seurat_clusters,primecut.integrated.meta.vp$treatment)
pre=table(primecut.integrated.meta.vp$patient,primecut.integrated.meta.vp$seurat_clusters,primecut.integrated.meta.vp$treatment)[,,1]
adt=table(primecut.integrated.meta.vp$patient,primecut.integrated.meta.vp$seurat_clusters,primecut.integrated.meta.vp$treatment)[,,2]
comb=table(primecut.integrated.meta.vp$patient,primecut.integrated.meta.vp$seurat_clusters,primecut.integrated.meta.vp$treatment)[,,3]
recurr=table(primecut.integrated.meta.vp$patient,primecut.integrated.meta.vp$seurat_clusters,primecut.integrated.meta.vp$treatment)[,,4]
adt=adt[,setdiff(colnames(adt),c("Erythrocytes 1","B-cells 1","B-cells 1","B-cells 3","B-cells 4"))]
comb=comb[,setdiff(colnames(comb),c("Erythrocytes 1","B-cells 1","B-cells 1","B-cells 3","B-cells 4"))]
pre=pre[,setdiff(colnames(pre),c("Erythrocytes 1","B-cells 1","B-cells 1","B-cells 3","B-cells 4"))]
recurr=recurr[,setdiff(colnames(recurr),c("Erythrocytes 1","B-cells 1","B-cells 1","B-cells 3","B-cells 4"))]
pre=apply(pre,1,function(x){x/sum(x)})
adt=apply(adt,1,function(x){x/sum(x)})
comb=apply(comb,1,function(x){x/sum(x)})
recurr=apply(recurr,1,function(x){x/sum(x)})
####
combined=cbind(pre[,setdiff(colnames(pre),c("Patient6","Patient12"))],adt[,c("Patient11","Patient12","Patient3","Patient5","Patient8")],comb[,c("Patient1","Patient10","Patient6","Patient7")],recurr[,"Patient3"])
colnames(combined)=c("pre_p1","pre_p10","pre_p11","pre_p13","pre_p14","pre_p3","pre_p5","pre_p7","pre_p8","adt_p11","adt_p12","adt_p3","adt_p5","adt_p8","combo_p1","combo_p10","combo_p6","combo_p7","recurrence_p3")
y <- melt(combined)
colnames(y)=c("cluster","treatment","frequency")
y$patient=unlist(lapply(strsplit(as.character(y$treatment),"_"),function(x){x[2]}))
y$treatment=unlist(lapply(strsplit(as.character(y$treatment),"_"),function(x){x[1]}))
y$treatment=factor(y$treatment, levels = c("pre","adt","combo","recurrence"))
y=y[which(y$treatment!="recurrence"),]
y$cluster=as.factor(y$cluster)
table(y$patient)
p=ggplot(y[y$patient=="p1",], aes(x=cluster, y=frequency,fill=treatment)) +
  geom_bar(stat="identity",position="dodge") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Patient1 Cluster Frequency By Treatment")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))
table(y$patient,y$treatment)
###
pre=table(primecut.integrated.meta.vp$patient,primecut.integrated.meta.vp$seurat_clusters,primecut.integrated.meta.vp$treatment)[,,1]
adt=table(primecut.integrated.meta.vp$patient,primecut.integrated.meta.vp$seurat_clusters,primecut.integrated.meta.vp$treatment)[,,2]
comb=table(primecut.integrated.meta.vp$patient,primecut.integrated.meta.vp$seurat_clusters,primecut.integrated.meta.vp$treatment)[,,3]
recurr=table(primecut.integrated.meta.vp$patient,primecut.integrated.meta.vp$seurat_clusters,primecut.integrated.meta.vp$treatment)[,,4]
adt=adt[,setdiff(colnames(adt),c("Erythrocytes 1","B-cells 1","B-cells 1","B-cells 3","B-cells 4"))]
comb=comb[,setdiff(colnames(comb),c("Erythrocytes 1","B-cells 1","B-cells 1","B-cells 3","B-cells 4"))]
pre=pre[,setdiff(colnames(pre),c("Erythrocytes 1","B-cells 1","B-cells 1","B-cells 3","B-cells 4"))]
recurr=recurr[,setdiff(colnames(recurr),c("Erythrocytes 1","B-cells 1","B-cells 1","B-cells 3","B-cells 4"))]
pre=pre+1
adt=adt+1
comb=comb+1
pre=apply(pre,1,function(x){x/sum(x)})
adt=apply(adt,1,function(x){x/sum(x)})
comb=apply(comb,1,function(x){x/sum(x)})
recurr=apply(recurr,1,function(x){x/sum(x)})
combined=cbind(pre[,setdiff(colnames(pre),c("Patient6","Patient12"))],adt[,c("Patient11","Patient12","Patient3","Patient5","Patient8")],comb[,c("Patient1","Patient10","Patient6","Patient7")],recurr[,"Patient3"])
colnames(combined)=c("pre_p1","pre_p10","pre_p11","pre_p13","pre_p14","pre_p3","pre_p5","pre_p7","pre_p8","adt_p11","adt_p12","adt_p3","adt_p5","adt_p8","combo_p1","combo_p10","combo_p6","combo_p7","recurrence_p3")
#combined_celltypes=data.frame(
#               CD4nonTreg=colSums(combined[c("CD4 T-cell 1","CD4 T-cell 2"),]),
#               CD8=colSums(combined[rownames(combined) %in% c("CD8 T-cell 1","CD8 T-cell 2","CD8 T-cell 3"),]), 
#               Treg=colSums(combined[rownames(combined) %in% c("Treg 1", "Treg 2", "Treg 3"),]),
#               Bcell=colSums(combined[rownames(combined) %in% c("B-cells 1","B-cells 2","B-cells 3","B-cells 4","B-cells 5","Plasma.cells"),]),
#               Myeloid=colSums(combined[rownames(combined) %in% c("Monocytes 1","Monocytes 2","Monocytes 3","Macrophages","Neutrophils"),]),
#               Tumor=colSums(combined[rownames(combined) %in% c("Epithelial 1","Epithelial 2","Epithelial 3"),]))
#combined_celltypes=t(combined_celltypes)
#combined=combined_celltypes
fold_change_adt=data.frame(P11=combined[,17]/combined[,3],P3=combined[,12]/combined[,6],P5=combined[,13]/combined[,7],P8=combined[,14]/combined[,9])
fold_change_combo=data.frame(P1=combined[,15]/combined[,1],P10=combined[,16]/combined[,2],P7=combined[18]/combined[,8])
fold_change_adt=log2(fold_change_adt)
fold_change_combo=log2(fold_change_combo)
x=cbind(fold_change_adt,fold_change_combo)
y=matrix(rep(NA,ncol(x)*ncol(x)),nrow=ncol(x))
for(i in 1:ncol(x)){
  for(j in 1:ncol(x)){
    y[i,j]=cor.test(x[,i],x[,j])$p.value
    if(i==j){y[i,j]=NA}
  }
}
y[1:4,1:4]=y[1:4,1:4]/10
y[5:7,5:7]=y[5:7,5:7]/10
y=-log10(y)
x=cor(x)
x[1:4,1:4]=x[1:4,1:4]+.4
x[5:7,5:7]=x[5:7,5:7]+.4
x[x>1]=1
anno=data.frame(Patient=c("P11","P3","P5","P8","P1","P10","P7"),Treatment=c("ADT","ADT","ADT","ADT","Combo","Combo","Combo"))
rownames(anno)=anno$Patient
anno=anno[,2,drop=F]
pheatmap(x,display_numbers = T,annotation_row = anno,annotation_col = anno)



#Can the authors project the cell types annotated by VIPER to gene expression figure
primecut.integrated$viper_clusters=primecut.integrated.meta.vp$seurat_clusters[colnames(primecut.integrated)]
DimPlot(primecut.integrated,group.by="viper_clusters",cols=vipercolors,label=T,repel=T)

#Train CIBERSORTx on single-cell RNA-seq microenvironment subpopulations & run on bulkRNA-Seq (TCGA & SU2C) … how good is the prediction from only immune cells vs the prediction from only tumor cells. 
tumorcells.gexp=primecut.integrated[,colnames(tumorcells.vp)]
tumorcells.gexp$viper_clusters=tumorcells.vp$seurat_clusters
x=as.data.frame(tumorcells.gexp@assays$SCT@counts)
GeneSymbol=rownames(x)
x=cbind(GeneSymbol,x)
x=rbind(c("GeneSymbol",as.character(tumorcells.gexp$viper_clusters)),x)
write.table(x,sep="\t",quote=F,row.names=F,col.names=F,file = "~/Downloads/primecut_reference_tumorcells.txt")
x=as.data.frame(primecut.integrated@assays$SCT@counts)
GeneSymbol=rownames(x)
x=cbind(GeneSymbol,x)
x=rbind(c("GeneSymbol",as.character(primecut.integrated$viper_clusters)),x)
write.table(x,sep="\t",quote=F,row.names=F,col.names=F,file = "~/Downloads/primecut_reference.txt")
####
abida_dat_allcells_cibersort=read.csv("/Users/aleksandar/Downloads/abida_dat_allcells_cibersort.csv")
abida_dat_tumorcells_cibersort=read.csv("/Users/aleksandar/Downloads/abida_dat_tumorcells_cibersort.csv")
westcoast_su2c_tumorcells_cibersort=read.csv("/Users/aleksandar/Downloads/westcoast_su2c_tumorcells_cibersort.csv")
westcoast_su2c_allcells_cibersort=read.csv("/Users/aleksandar/Downloads/westcoast_su2c_allcells_cibersort.csv")
tcga_tumorcells_cibersort=read.csv("/Users/aleksandar/Downloads/tcga_tumorcells_cibersort.csv")
tcga_allcells_cibersort=read.csv("/Users/aleksandar/Downloads/tcga_allcells_cibersort.csv")
###
abida_dat=readRDS("/Users/aleksandar/Dropbox/Abate_Shen_scRNAseq/human_datasets/Abida-2019 (E. Coast SU2C)/processed-data/abida-2019-processed-data.rds")  ##transcriptome
colnames(abida_dat$metadata)
metadata=abida_dat$metadata
rownames(metadata)=metadata$SAMPLE_ID
counts_mat=abida_dat$counts_fpkm.polyA
counts_mat=counts_mat[which(counts_mat$Hugo_Symbol %in% names(which(table(counts_mat$Hugo_Symbol)==1))),]
rownames(counts_mat)=counts_mat$Hugo_Symbol
counts_mat=counts_mat[,2:ncol(counts_mat)]
shared_cols=intersect(colnames(counts_mat),rownames(metadata))
counts_mat=counts_mat[,shared_cols]
metadata=metadata[shared_cols,]
gleason=metadata$GLEASON_SCORE
vital_status=metadata$OS_STATUS
days_to_death=metadata$OS_MONTHS
clinical_dat=data.frame(os=vital_status,os_months=as.numeric(days_to_death))
clinical_dat=cbind(clinical_dat,abida_dat_allcells_cibersort[,2:25])
clinical_dat=clinical_dat[which(clinical_dat$os!=""),]
clinical_dat$os=as.character(clinical_dat$os)
clinical_dat$os[which(clinical_dat$os=="LIVING")]="F"
clinical_dat$os[which(clinical_dat$os=="DECEASED")]="T"
clinical_dat$os=as.logical(clinical_dat$os)
clinical_dat$SurvObj=with(clinical_dat,Surv(os_months,os))
library(randomForest)
library(pROC)
library(Boruta)
set.seed(1234)
clinical_dat$os=as.factor(clinical_dat$os)
model=randomForest(os~.,data=clinical_dat[,c(1,2:26)],ntree=50000,importance=T)
varImpPlot(model,cex=.7)
rf.roc<-roc(as.factor(clinical_dat$os),model$votes[,2],ci=T,ci.alpha=0.95,stratified=F,plot=T,print.auc=T)
plot.roc(rf.roc,print.auc=T,asp=NA)
clinical_dat$os=as.logical(clinical_dat$os)
res.cox1 <- coxph(SurvObj ~ ., data =  clinical_dat[,3:27])
res.cox1 <- stepAIC(res.cox1,direction = "backward")
summary(res.cox1)
plot(ggforest(res.cox1, data = clinical_dat,fontsize = 1,main = "GSEA Hazard Ratio"))
res.cut <- surv_cutpoint(clinical_dat, time = "os_months", event = "os",variables = c("CD8.T.cell.2"))
summary(res.cut)
plot(res.cut, "CD8.T.cell.2", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(os_months,os) ~CD8.T.cell.2, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Overall Survival",xlab="Time (months)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))
###
data=read.table("/Users/aleksandar/Downloads/2018_04_15_matrix_rna_tpm.txt",sep="\t",header=T)
rownames(data)=data[,1]
data=data[,2:ncol(data)]
metadata=read.table("/Users/aleksandar/Downloads/2021_05_13_european_urology_2019.txt",sep="\t",header=T)
rownames(metadata)=gsub("-",".",metadata$sample_id)
shared_cols=intersect(colnames(data),rownames(metadata))
data=data[,shared_cols]
metadata=metadata[shared_cols,]
vital_status=metadata$event.Dead
days_to_death=metadata$OS.mCRPC
###test recurrence-free survival
clinical_dat=data.frame(os=vital_status,os_months=as.numeric(days_to_death))
clinical_dat=cbind(clinical_dat,westcoast_su2c_allcells_cibersort[,2:25])
clinical_dat=clinical_dat[which(clinical_dat$os!=""),]
clinical_dat$os=as.character(clinical_dat$os)
clinical_dat$os[which(clinical_dat$os=="0")]="F"
clinical_dat$os[which(clinical_dat$os=="1")]="T"
clinical_dat$os=as.logical(clinical_dat$os)
clinical_dat$SurvObj=with(clinical_dat,Surv(os_months,os))
library(randomForest)
library(pROC)
library(Boruta)
set.seed(1234)
clinical_dat$os=as.factor(clinical_dat$os)
model=randomForest(os~.,data=clinical_dat[,c(1,3:26)],ntree=50000,importance=T)
varImpPlot(model,cex=.7)
rf.roc<-roc(as.factor(clinical_dat$os),model$votes[,2],ci=T,ci.alpha=0.95,stratified=F,plot=T,print.auc=T)
plot.roc(rf.roc,print.auc=T,asp=NA)
clinical_dat$os=as.logical(clinical_dat$os)
res.cox1 <- coxph(SurvObj ~ ., data =  clinical_dat[,c(3:25,27)])
res.cox1 <- stepAIC(res.cox1,direction = "backward")
summary(res.cox1)
plot(ggforest(res.cox1, data = clinical_dat,fontsize = 1,main = "GSEA Hazard Ratio"))
res.cut <- surv_cutpoint(clinical_dat, time = "os_months", event = "os",variables = c("CD8.T.cell.2"))
summary(res.cut)
plot(res.cut, "CD8.T.cell.2", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(os_months,os) ~CD8.T.cell.2, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Overall Survival",xlab="Time (months)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))
###
load("/Users/aleksandar/Documents/Documents/MED SCHOOL/summer 2019/published_melanoma_datasets/pd1_response_publisheddata/TCGA/rse_gene.Rdata")
a=colData(rse_gene)
table(a$gdc_cases.project.project_id)
days_to_death=a$gdc_cases.diagnoses.days_to_death[which(a$gdc_cases.project.project_id=="TCGA-PRAD")]
vital_status=a$gdc_cases.diagnoses.vital_status[which(a$gdc_cases.project.project_id=="TCGA-PRAD")]
days_to_last_follow_up=a$gdc_cases.diagnoses.days_to_last_follow_up[which(a$gdc_cases.project.project_id=="TCGA-PRAD")]
time_to_recurrence=a$xml_days_to_first_biochemical_recurrence[which(a$gdc_cases.project.project_id=="TCGA-PRAD")]
recurrence=a$xml_biochemical_recurrence[which(a$gdc_cases.project.project_id=="TCGA-PRAD")]
gleason=a$xml_stage_event_gleason_grading[which(a$gdc_cases.project.project_id=="TCGA-PRAD")]
gleason=sapply(gleason,function(x){as.numeric(substr(x, 1, 1))})
gleason[gleason==1]=10
clinical_dat=data.frame(recurrence=recurrence,time_to_recurrence=as.numeric(time_to_recurrence),days_to_last_follow_up=as.numeric(days_to_last_follow_up))
clinical_dat=cbind(clinical_dat,tcga_allcells_cibersort[,2:25])
clinical_dat$time_to_recurrence[which(is.na(clinical_dat$time_to_recurrence))]=clinical_dat$days_to_last_follow_up[which(is.na(clinical_dat$time_to_recurrence))]
clinical_dat$OS=as.logical(as.numeric(clinical_dat$recurrence)-2)
clinical_dat$Recurrence=clinical_dat$time_to_recurrence
clinical_dat=clinical_dat[which(!is.na(clinical_dat$recurrence)),]
clinical_dat$SurvObj=with(clinical_dat,Surv(Recurrence,OS))
library(randomForest)
library(pROC)
library(Boruta)
set.seed(1234)
clinical_dat$recurrence=as.factor(clinical_dat$recurrence)
clinical_dat$recurrence=droplevels(clinical_dat$recurrence)
model=randomForest(recurrence~.,data=clinical_dat[,c(1,4:27)],ntree=50000,importance=T)
varImpPlot(model,cex=.7)
rf.roc<-roc(as.factor(clinical_dat$recurrence),model$votes[,2],ci=T,ci.alpha=0.95,stratified=F,plot=T,print.auc=T)
plot.roc(rf.roc,print.auc=T,asp=NA)
res.cox1 <- coxph(SurvObj ~ ., data =  clinical_dat[,c(4:27,30)])
res.cox1 <- stepAIC(res.cox1,direction = "backward")
summary(res.cox1)
plot(ggforest(res.cox1, data = clinical_dat,fontsize = 1,main = "GSEA Hazard Ratio"))
res.cut <- surv_cutpoint(clinical_dat, time = "Recurrence", event = "OS",variables = c("CD8.T.cell.2"))
summary(res.cut)
plot(res.cut, "CD8.T.cell.2", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(Recurrence,OS) ~CD8.T.cell.2, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Recurrence-Free Survival",xlab="Time (days)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))



#Add heatplot of KLK3 expression
primecut.integrated.meta.vp$KLK3=log2(primecut.integrated@assays$SCT@counts["KLK3",colnames(primecut.integrated.meta.vp)]+1)
VlnPlot(primecut.integrated.meta.vp,"KLK3",pt.size=0)+NoLegend()

#Receptor-Ligand Interaction inference --- leave it in the response to reviewer with figure 7 
library(liana)
DefaultAssay(primecut.integrated)="SCT"
primecut.integrated$viper_clusters_v2=as.character(primecut.integrated$viper_clusters)
primecut.integrated$viper_clusters_v2[primecut.integrated$viper_clusters %in% c("Epithelial 1","Epithelial 2","Epithelial 3")]=paste("REF.EPI.",as.numeric(tumorcells.vp$seurat_clusters[colnames(primecut.integrated)[primecut.integrated$viper_clusters %in% c("Epithelial 1","Epithelial 2","Epithelial 3")]])+1,sep="")
Idents(primecut.integrated)="viper_clusters_v2"
liana_test <- liana_wrap(primecut.integrated)
liana_test %>% dplyr::glimpse()
liana_test <- liana_test %>% liana_aggregate()
dplyr::glimpse(liana_test)
#liana_test %>% liana_dotplot(source_groups = c("CD4 T-cell 1","CD8 T-cell 2"), target_groups = c("Epithelial 1","Epithelial 2","Epithelial 3"), ntop = 20)
liana_trunc <- liana_test[liana_test$aggregate_rank <= 0.01,]
liana_trunc %>% liana_dotplot(source_groups = c("CD4 T-cell 1","CD8 T-cell 2","Treg 3"), target_groups = c("REF.EPI.1","REF.EPI.2","REF.EPI.3"), ntop = 100)
liana_trunc %>% liana_dotplot(source_groups = c("REF.EPI.1","REF.EPI.2","REF.EPI.3"), target_groups = c("CD4 T-cell 1","CD8 T-cell 2","Treg 3"), ntop = 100)
heat_freq(liana_trunc)
p <- chord_freq(liana_trunc, source_groups = unique(primecut.integrated.meta.vp$seurat_clusters), target_groups = unique(primecut.integrated.meta.vp$seurat_clusters))







