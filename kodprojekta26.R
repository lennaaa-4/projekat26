library(Seurat)
library(tidyverse)
library(BiocManager)
View(seurat_obj@meta.data)
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
str(seurat_obj)

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 5000)
top10 <- head(VariableFeatures(seurat_obj), 10)
top10
plot1 <- VariableFeaturePlot(seurat_obj)
plot1
LabelPoints(plot = plot1, points = c("LINC00958", top10), repel = TRUE)
plot1

all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
str(seurat_obj)

seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(seurat_obj, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(seurat_obj)
Loadings(seurat_obj[["pca"]])["LINC00958", ]

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))
View(seurat_obj@meta.data)
DimPlot(seurat_obj, group.by = "RNA_snn_res.0.5", label = TRUE)
DimPlot(seurat_obj, group.by = "RNA_snn_res.0.1", label = TRUE)
DimPlot(seurat_obj, group.by = "RNA_snn_res.0.3", label = TRUE)
DimPlot(seurat_obj, group.by = "RNA_snn_res.0.7", label = TRUE)
DimPlot(seurat_obj, group.by = "RNA_snn_res.1", label = TRUE)
