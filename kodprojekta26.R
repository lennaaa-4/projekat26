library(Seurat)
library(tidyverse)
library(BiocManager) # Učitavanje paketa potrebnih za rad
View(seurat_obj@meta.data) # Učitavanje meta podataka
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm') # Izrada Violin Plota 

seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
str(seurat_obj) # Logaritamska normalizacija podataka

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 5000) 
top10 <- head(VariableFeatures(seurat_obj), 10)
top10
plot1 <- VariableFeaturePlot(seurat_obj)
plot1
LabelPoints(plot = plot1, points = c("LINC00958", top10), repel = TRUE)
plot1 # Identifikacija najvarijabilnijih gena, određivanje top10 najvarijabilnijih i izrada VariableFeaturesPlota
### izvrsena promena broja gena sa 2000 na 5000!

all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
str(seurat_obj) # skaliranje podataka
### izvrsena promena: skalirani samo prethodno definisani VariableFeatures, ne svi geni!

seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)
library(ggplot2)
DimHeatmap(seurat_obj, dims = 1, cells = 500, balanced = TRUE, fast = FALSE) + 
  labs(title = "DimHeatMap za PC_1", x="celije",y="geni") +
    theme(
    legend.position = "right",
  )
ElbowPlot(seurat_obj)
Loadings(seurat_obj[["pca"]])["LINC00958", ] #izvrsena PCA analiza, urađen DimHeatMap (za PC_1), napravljen ElbowPlot, proverena LINC00958 u PC
### izvrsene promene: DimHeatMapu dodani modifikovan naslov, nazivi x-ose i y-ose, provera LINC00958 u PC

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))
View(seurat_obj@meta.data) # klasterovanje
library(plotly)
pca <- Embeddings(seurat_obj, reduction = "pca")[, 1:3]
plot_ly(
  x = pca[,1],
  y = pca[,2],
  z = pca[,3],
  color = as.factor(Idents(seurat_obj)),
  type = "scatter3d",
  mode = "markers"
)  ### 3d prikaz

pca_df <- data.frame(
  PC1 = pca[, 1],
  PC2 = pca[, 2],
  PC3 = pca[, 3],
  cluster = Idents(seurat_obj)
)
ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(size = abs(PC3), color = cluster), alpha = 0.6) +
  theme_classic() +
  labs(
    title = "PCA: PC1 vs PC2 with PC3",
    size = "|PC3|",
    color = "Cluster"
  )  ### PCA bubble plot -> vizualizacija PC1 PC2 i PC3 u kojoj su PC1 i PC2 prikazane na osama a PC3 je prikazana kroz velicinu tacaka

plot_ly(
  x = pca[, 1],
  y = pca[, 2],
  z = pca[, 3],
  color = as.factor(Idents(seurat_obj)),
  type = "scatter3d",
  mode = "markers"
) %>%
  layout(
    scene = list(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2"),
      zaxis = list(title = "PC3")
    )
  )
 ### 3d prikaz 2. deo (isto kao i prethodno)

library(ggplot2)
library(grid)

# Klasteri na rezoluciji 0.1
Idents(seurat_obj) <- "RNA_snn_res.0.1"

# Prve 3 PCA dimenzije
pca <- Embeddings(seurat_obj, reduction = "pca")[, 1:3]

# Data frame
pca_df <- data.frame(
  PC1 = pca[, 1],
  PC2 = pca[, 2],
  PC3 = pca[, 3],
  Cluster = Idents(seurat_obj)
)

# -------------------------------------------------
# PROJEKCIJA PC3 POD 45°
# -------------------------------------------------

# Faktor koji određuje koliko je PC3 osa dugačka
z_scale <- 0.5

# PC3 se projektuje pod 45°
pca_df$plot_x <- pca_df$PC1 + pca_df$PC3 * z_scale
pca_df$plot_y <- pca_df$PC2 + pca_df$PC3 * z_scale

# Granice
x_min <- min(pca_df$plot_x)
x_max <- max(pca_df$plot_x)
y_min <- min(pca_df$plot_y)
y_max <- max(pca_df$plot_y)

# Početak koordinatnog sistema
origin_x <- min(pca_df$PC1)
origin_y <- min(pca_df$PC2)

# Dužina PC3 ose
z_length <- max(abs(pca_df$PC3)) * z_scale

# -------------------------------------------------
# PLOT
# -------------------------------------------------

ggplot(pca_df, aes(x = plot_x, y = plot_y)) +
  
  # Tačke
  geom_point(
    aes(color = Cluster),
    size = 1.5,
    alpha = 0.7
  ) +
  
  # X osa = PC1
  geom_segment(
    aes(
      x = origin_x,
      y = origin_y,
      xend = x_max,
      yend = origin_y
    ),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.2, "cm"))
  ) +
  
  # Y osa = PC2
  geom_segment(
    aes(
      x = origin_x,
      y = origin_y,
      xend = origin_x,
      yend = y_max
    ),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.2, "cm"))
  ) +
  
  # Z osa = PC3, pod 45°
  geom_segment(
    aes(
      x = origin_x,
      y = origin_y,
      xend = origin_x + z_length,
      yend = origin_y + z_length
    ),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.2, "cm"))
  ) +
  
  # Oznake osa
  annotate(
    "text",
    x = x_max,
    y = origin_y,
    label = "PC1",
    hjust = -0.2,
    vjust = 1
  ) +
  
  annotate(
    "text",
    x = origin_x,
    y = y_max,
    label = "PC2",
    hjust = 1.2,
    vjust = -0.2
  ) +
  
  annotate(
    "text",
    x = origin_x + z_length,
    y = origin_y + z_length,
    label = "PC3",
    hjust = -0.1,
    vjust = -0.5
  ) +
  
  theme_classic() +
  labs(
    title = "PCA plot: PC1, PC2 and PC3",
    color = "Cluster"
  ) +
  
  coord_fixed() 
### vizualizacija 3 PC ali u 2D 

Idents(seurat_obj) <- "RNA_snn_res.0.1"
table(Idents(seurat_obj))

unique(seurat_obj$RNA_snn_res.0.1)


# Sve clustering kolone u objektu
grep("snn_res", colnames(seurat_obj@meta.data), value = TRUE)
# Broj klastera za svaku rezoluciju
lapply(
  seurat_obj@meta.data[
    , grep("snn_res", colnames(seurat_obj@meta.data))
  ],
  function(x) length(unique(x))
)
# Trenutni aktivni identitet
Idents(seurat_obj)

# Koje vrednosti trenutno postoje u 0.1
unique(seurat_obj$RNA_snn_res.0.1)
# provera broja klastera u 0.1 rezoluciji


DimPlot(seurat_obj, group.by = "RNA_snn_res.0.5", label = TRUE)
DimPlot(seurat_obj, group.by = "RNA_snn_res.0.1", label = TRUE)
DimPlot(seurat_obj, group.by = "RNA_snn_res.0.3", label = TRUE)
DimPlot(seurat_obj, group.by = "RNA_snn_res.0.7", label = TRUE)
DimPlot(seurat_obj, group.by = "RNA_snn_res.1", label = TRUE)
# DimPlotovi za sve rezolucije

Idents(seurat_obj)
Idents(seurat_obj) <- "RNA_snn_res.0.1"
Idents(seurat_obj)
# postavka identiteta na 0.1 rezoluciju

reticulate::py_install(packages ='umap-learn')
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
DimPlot(seurat_obj, reduction = "umap", label = TRUE)
# UMAP analiza i vizualizacija UMAP DimPlota
