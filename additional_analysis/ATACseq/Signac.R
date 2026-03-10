library(Seurat)
library(scran)
library(Signac)
library(Matrix)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)


srt_RNA <- readRDS("/trigos_team/CASCADE/Analysis/230113_seurat_intron/CA0027/dura_base_skull_13/CA0027_dura_base_skull_13_srt_new_filtering_doublets_non_tumour_removed.Rds")
counts <- readMM("/trigos_team/CASCADE/SingleCell/CA0027/dura_base_skull_13/multiome/CA0027_dura_base_skull_13_multi/outs/raw_feature_bc_matrix/matrix.mtx.gz")
features <- read_tsv("/trigos_team/CASCADE/SingleCell/CA0027/dura_base_skull_13/multiome/CA0027_dura_base_skull_13_multi/outs/raw_feature_bc_matrix/features.tsv.gz", col_names = F)
cell_ids <- read_tsv("/trigos_team/CASCADE/SingleCell/CA0027/dura_base_skull_13/multiome/CA0027_dura_base_skull_13_multi/outs/raw_feature_bc_matrix/barcodes.tsv.gz", col_names = FALSE)$X1
View(counts)
rownames(counts) <- features
colnames(counts) <- cell_ids

counts <- counts[, colnames(srt)]
annotation.file = "/trigos_team/CASCADE/Analysis/reference_annotation/multiome/genes.gtf"

test_mat <- Read10X_h5("/trigos_team/CASCADE/SingleCell/CA0027/dura_base_skull_13/multiome/CA0027_dura_base_skull_13_multi/outs/filtered_feature_bc_matrix.h5")
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
fragpath <- read_tsv("/trigos_team/CASCADE/SingleCell/CA0027/dura_base_skull_13/multiome/CA0027_dura_base_skull_13_multi/outs/atac_fragments.tsv.gz")

srt <- CreateSeuratObject(
  counts = test_mat$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
srt[["ATAC"]] <- CreateChromatinAssay(
  counts = test_mat$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

DefaultAssay(srt) <- "ATAC"

srt <- NucleosomeSignal(srt)
srt <- TSSEnrichment(srt)

VlnPlot(
  object = srt,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

srt <- subset(
  x = srt,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
srt <- subset(srt, cells = colnames(srt_RNA))
srt@meta.data <- srt@meta.data[colnames(srt), ]
CA27_lineage <- absolute_cell_annotation$CA0027_dura_base_skull_13
srt$lineage <- NA
for (f in names(CA27_lineage)){
  srt@meta.data[CA27_lineage[[f]], "lineage"] <- f
}

peaks <- CallPeaks(srt)

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(srt),
  features = peaks,
  cells = colnames(srt)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
srt[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)


DefaultAssay(srt) <- "RNA"
lib.median <- median(srt$nCount_RNA)
target_pseudocount <- 1
srt <- srt %>% 
  NormalizeData(scale.factor = lib.median / target_pseudocount) %>% 
  ScaleData(do.scale = F)
set.seed(100)
poisson_fit3 <- modelGeneVarByPoisson(srt@assays$RNA@counts, dispersion = 0.2, size.factors = target_pseudocount * srt$nCount_RNA / lib.median)
residuals3 <- poisson_fit3$total - poisson_fit3@metadata$trend(poisson_fit3$mean)
names(residuals3) <- rownames(srt)
top_genes3 <- rownames(srt)[order(residuals3, decreasing = TRUE)]

set.seed(100)
srt <- srt %>% 
  RunPCA(features = top_genes3[1:2000], weight.by.var = F) %>% 
  FindNeighbors(dims = 1:20) %>% 
  FindClusters() %>% 
  RunUMAP(dims = 1:20)
FeaturePlot(srt, c("PCNA","EXO1", "CDK1", "TOP2A"))
srt <- AddModuleScore(srt, features = list(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes))
colnames(srt@meta.data)[grep("Cluster", colnames(srt@meta.data))] <- c("S.Score", "G2M.Score")
FeaturePlot(srt, c("S.Score", "G2M.Score"))


#processing peaks, dim reduction
DefaultAssay(srt) <- "peaks"
srt <- FindTopFeatures(srt, min.cutoff = "q0")
srt <- RunTFIDF(srt)
srt <- RunSVD(srt)

DepthCor(srt, assay = "peaks", reduction = "lsi")
srt <- RunUMAP(srt, reduction = "lsi", dims = 2:30) %>%
  FindNeighbors(reduction = "lsi", dims = 2:30) %>% 
  FindClusters()
DimPlot(srt, label = T) + DimPlot(srt, group.by = "lineage")
FeaturePlot(srt, c("S.Score", "G2M.Score"))

#gene activity: estimate RNA expression from peaks information 
gene_act <- GeneActivity(srt, assay = "peaks", extend.upstream = 5000)
srt[["Activity"]] <- CreateAssayObject(counts = gene_act)
DefaultAssay(srt) <- "Activity" 

lib.median <- median(srt$nCount_Activity)
target_pseudocount <- 1
srt <- srt %>% 
  NormalizeData(scale.factor = lib.median / target_pseudocount) %>% 
  ScaleData(do.scale = F)
set.seed(100)
poisson_fit3 <- modelGeneVarByPoisson(srt@assays$Activity@counts, dispersion = 0.2, size.factors = target_pseudocount * srt$nCount_RNA / lib.median)
residuals3 <- poisson_fit3$total - poisson_fit3@metadata$trend(poisson_fit3$mean)
names(residuals3) <- rownames(srt)
top_genes3 <- rownames(srt)[order(residuals3, decreasing = TRUE)]

set.seed(100)
srt <- srt %>% 
  RunPCA(features = top_genes3[1:2000], weight.by.var = F) %>% 
  FindNeighbors(dims = 1:20) %>% 
  FindClusters() %>% 
  RunUMAP(dims = 1:20)
DimPlot(srt, label = T) + DimPlot(srt, group.by = "lineage")
FeaturePlot(srt, c("S.Score", "G2M.Score", "AR", "ASCL1", "BACE2", "TGM2"))




# build a joint neighbor graph using both assays
srt <- FindMultiModalNeighbors(
  object = srt,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:20, 2:30),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
srt <- RunUMAP(
  object = srt,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

srt <- FindClusters(srt, graph.name = "wsnn", resolution = 1)
DimPlot(srt, label = T, reduction = "umap") + DimPlot(srt, reduction = "umap", group.by = "lineage")
DefaultAssay(srt) <- "RNA"
FeaturePlot(srt, reduction = "umap",c("S.Score", "G2M.Score", "AR", "ASCL1", "WSB1", "TGM2", "HSPH1", "NNAT"))



DefaultAssay(srt) <- "peaks"

# first compute the GC content for each peak
srt <- RegionStats(srt, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
srt <- LinkPeaks(
  object = srt,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = c("AR", "ASCL1", "WSB1", "TGM2", "HSPH1", "NNAT","LMX1A", "FOLH1")
)

CoveragePlot(
  object = srt,
  region = "OSMR",
  features = "OSMR",
  expression.assay = "RNA",
  group.by = "lineage",
  extend.upstream = 10000,
  extend.downstream = 500
)

Idents(srt) <- srt$lineage
DE_peak <- FindAllMarkers(srt, test.use = "LR", latent.vars = "nCount_peaks")
DE_peak$closest_gene <- ClosestFeature(srt, rownames(DE_peak))$gene_name
DE_peak$distance <- ClosestFeature(srt, rownames(DE_peak))$distance
saveRDS(DE_peak, "DEpeaks.Rds")

branch1_depeak <- DE_peak %>% dplyr::filter(cluster == "branch1") %>% dplyr::filter(distance == 0) %>% dplyr::filter(p_val_adj<0.05) %>% dplyr::arrange(avg_log2FC)

head(branch1_depeak$closest_gene[order(branch1_depeak$avg_log2FC, decreasing = T)])
DefaultAssay(srt) <- "RNA"
FeaturePlot(srt, reduction = "umap", features = (head(branch1_depeak$closest_gene[order(branch1_depeak$avg_log2FC, decreasing = T)])))
FeaturePlot(srt_RNA, reduction = "umap", features = c("TEAD3", "NFE2L2"))


## motif analysis 
library(JASPAR2020)
library(TFBSTools)

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE, species ='Homo sapiens')
)

DefaultAssay(srt) <- "peaks"
srt <- AddMotifs(srt, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
top.da.peak <- DE_peak$gene[DE_peak$p_val_adj < 0.05]
saveRDS(srt, "CA27_dura_multiome_srt.Rds")

# find peaks open in Pvalb or Sst cells
open.peaks <- AccessiblePeaks(srt)

branch1_depeak
top_branch1.peak <- branch1_depeak$gene[branch1_depeak$p_val_adj < 0.05&branch1_depeak$avg_log2FC>0]

# match the overall GC content in the peak set
meta.feature <- GetAssayData(srt, assay = "peaks", slot = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top_branch1.peak, ],
  features.match = "GC.percent",
  n = 50000, 
  verbose = T
)

enriched.motifs  <- FindMotifs(
  object = srt,
  features = top_branch1.peak, 
  background = peaks.matched
)

MotifPlot(
  object = srt,
  motifs = head(rownames(enriched.motifs))
)

srt <- RunChromVAR(
  object = srt,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

DefaultAssay(srt) <- "chromvar"
FeaturePlot(srt, "MA0808.1") # TEAD3

de.motif.activity <- FindMarkers(
  object = srt,
  ident.1 = 'branch1',
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

motifs_names <- sapply(row.names(de.motif.activity), function(x) {name(getMatrixByID(JASPAR2020, ID = x))})
de.motif.activity$motif.name <- motifs_names

top_motif <- de.motif.activity %>% dplyr::filter(p_val_adj < 0.05 & avg_diff > 1.25) %>% dplyr::filter(pct.1-pct.2 > 0.5)

DefaultAssay(srt) <- "chromvar"
FeaturePlot(srt, c("MA0090.3", "MA0841.1", "MA0089.2"))
DefaultAssay(srt) <- 'RNA'
FeaturePlot(srt, c("BCL6", "TGM2", "OSMR", "RBMS3"))
rna_mat <- as.matrix(srt@assays$RNA@data[c("BCL6", "TGM2", "OSMR", "RBMS3"), ])
motif_mat <- as.matrix(srt@assays$chromvar@data)

cormat <- cor(t(rbind(rna_mat, motif_mat)))
interested <- cormat[c("BCL6", "TGM2", "OSMR", "RBMS3"), rownames(motif_mat)]
heatmap(t(interested))
cm <- colMeans(interested)
summary(cm)
names(which(cm > 0.2))
DefaultAssay(srt) <- "chromvar"
FeaturePlot(srt,names(which(cm > 0.2)))
sapply(names(which(cm > 0.2)), function(x) {name(getMatrixByID(JASPAR2020, ID = x))})

saveRDS(srt, "CA27_dura_multiome_srt.Rds")


###test regulon peaks and motifs

regulonlist <- readRDS("~/regulon/PSMA/regulonlist.Rds")

srt <- LinkPeaks(
  object = srt,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = regulonlist$TEAD3.pos
)

tead3_link <- Links(srt)
tead3_link <- tead3_link[tead3_link$pvalue<0.05&tead3_link$score>0, ]

peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[tead3_link$peak, ],
  features.match = "GC.percent",
  n = 50000, 
  verbose = T
)

enriched.motifs.tead3  <- FindMotifs(
  object = srt,
  features = tead3_link$peak, 
  background = peaks.matched
)
enriched.motifs.tead3$motif.name

test_regulonlist <- function(srt, regulon){
  srt2 <- LinkPeaks(
    object = srt,
    peak.assay = "peaks",
    expression.assay = "RNA",
    genes.use = regulon
  )
  links <- Links(srt2)
  links <- links[links$pvalue<0.05&links$score>0, ]
  peaks.matched <- MatchRegionStats(
    meta.feature = meta.feature[open.peaks, ],
    query.feature = meta.feature[links$peak, ],
    features.match = "GC.percent",
    n = 50000, 
    verbose = T
  )
  enriched.motifs  <- FindMotifs(
    object = srt2,
    features = links$peak, 
    background = peaks.matched
  )
  return(list(link = links, motif = enriched.motifs))
}

regulon_link_motif <- test_regulonlist(srt, regulonlist$NFE2L2.pos)

absolute_module <- list(
  # metabolism1 = c("FOXN3.pos","RORA.pos", "GRHL2.pos","HNF4G.pos","PKNOX2.pos","SOX6.pos", "RFX3.pos"), 
  # metabolism2 = c( "MYC.pos", "TBX3.pos", "FOXA1.pos", "NKX3-1.pos"),#"CREB3.pos", "ATF4.pos", "JUN.pos", "JUND.pos",  too weak
  
  AR_absolute = c("AR.pos", "GATA2.pos"), #exclusive
  NE_absolute = c("ASCL1.pos", "NFATC2.pos"),
  # NE_transition = c("HOXB3.pos", "ATOH1.pos", "LMX1A.pos", "MEIS1.pos"),
  
  branch2.1 = c( "ARID3A.pos","RORC.pos","MEF2D.pos","MAFF.pos"), 
  branch2.2 = c("GRHL1.pos", "STAT2.pos"), 
  branch1 = c( "TEAD3.pos","NFE2L2.pos")
)
absolute_TF <- unlist(absolute_module, use.names = F)

regulon_link_motif <- lapply(regulonlist[absolute_TF], test_regulonlist, srt = srt)
