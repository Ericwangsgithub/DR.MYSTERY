# script to identify cluster identity -----------------
# Finding markers in every cluster
# Finding conserved markers 
# Finding markers DE between conditions
getwd()
setwd("/home/mingsong-wang/workspaces/gxr1-0/File Download")

set.seed(1234)

library(Seurat)
library(tidyverse)

# load data
ifnb_harmony <- readRDS('/home/mingsong-wang/workspaces/gxr1-0/File Download/ifnb_harmony.rds')
str(ifnb_harmony)
View(ifnb_harmony@meta.data)

# visualize data
clusters <- DimPlot(ifnb_harmony, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
condition <- DimPlot(ifnb_harmony, reduction = 'umap', group.by = 'stim')

condition|clusters

# Visualization of clusters with annotation name
Clusterplot<- DimPlot(ifnb_harmony, reduction = "umap", group.by = "seurat_annotations", label = TRUE)

clusters|Clusterplot
# findAll markers -----------------
FindAllMarkers(ifnb_harmony,
               logfc.threshold = 0.25,
               min.pct = 0.1,
               only.pos = TRUE,
               test.use = 'DESeq2',
               slot = 'counts')


# findConserved markers by FindConservedMarkers Function -------------
install.packages('metap')
library(metap)
install.packages('devtools')
devtools::install_github('immunogenomics/presto')
Idents(ifnb_harmony) <- "seurat_annotations"
nk.markers <- FindConservedMarkers(ifnb_harmony, 
                                   ident.1 = "NK", 
                                   grouping.var = "stim", 
                                   verbose = FALSE)
head(nk.markers)

# View the conserved cell type markers across conditions, showing both the expression level and the percentage of cells in a cluster expressing any given gene. 
Idents(ifnb_harmony) <- factor(Idents(ifnb_harmony), levels = c("pDC", "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono",
                                                "B Activated", "B", "CD8 T", "NK", "T activated", "CD4 Naive T", "CD4 Memory T"))

markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5",
                     "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1",
                     "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ", "PRSS57")
DotPlot(ifnb_harmony, features = markers.to.plot, cols = c("purple", "orange"), dot.scale = 8, split.by = "stim") +
  RotatedAxis()

# Notes:
# slot depends on the type of the test used, 
# default is data slot that stores normalized data
DefaultAssay(ifnb_harmony) <- 'RNA'
DefaultAssay(ifnb_harmony)
markers_cluster3 <- FindConservedMarkers(ifnb_harmony,
                                         ident.1 = 3,
                                         grouping.var = 'stim')

head(markers_cluster3)

# let's visualize top features
FeaturePlot(ifnb_harmony, features = c('FCGR3A'), min.cutoff = 'q10')


# min-cut off explanation:
seq(1,5)
SetQuantile('q50', seq(1,5))
SetQuantile('q10', seq(1,5))



# rename cluster 3 ident
Idents(ifnb_harmony)
ifnb_harmony <- RenameIdents(ifnb_harmony, `3` = 'CD16 Mono')

DimPlot(ifnb_harmony, reduction = 'umap', label = T)

# cells already have annotations provided in the metadata
View(ifnb_harmony@meta.data)

# Settings cluster identities is an iterative step
# multiple approaches could be taken - automatic/manual anotations (sometimes both)
# need to make sure each cell type forms a separate cluster

# setting Idents as Seurat annotations provided (also a sanity check!)
Idents(ifnb_harmony) <- ifnb_harmony@meta.data$seurat_annotations
Idents(ifnb_harmony)

DimPlot(ifnb_harmony, reduction = 'umap', label = TRUE)


# findMarkers between conditions ---------------------
ifnb_harmony$celltype.cnd <- paste0(ifnb_harmony$seurat_annotations,'_', ifnb_harmony$stim)
View(ifnb_harmony@meta.data)
Idents(ifnb_harmony) <- ifnb_harmony$celltype.cnd

DimPlot(ifnb_harmony, reduction = 'umap', label = TRUE)

# find markers
b.interferon.response <- FindMarkers(ifnb_harmony, ident.1 = 'CD16 Mono_STIM', ident.2 = 'CD16 Mono_CTRL')

head(b.interferon.response)

# plotting conserved features vs DE features between conditions
head(markers_cluster3)


FeaturePlot(ifnb_harmony, features = c('FCGR3A', 'AIF1', 'IFIT1'), split.by = 'stim', min.cutoff = 'q10')