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


### II PSEUDOBULK


library(Seurat)
library(DESeq2)
library(tidyverse) # ucitavanje paketa

cell_plot <- DimPlot(seurat_obj, reduction = 'umap', group.by = 'cell_type', label = TRUE)
cond_plot <- DimPlot(seurat_obj, reduction = 'umap', group.by = 'pathology')
cell_plot|cond_plot # vizualizacija 

seurat_obj$samples <- paste0(seurat_obj$pathology, "_", seurat_obj$NBB_case)
DefaultAssay(seurat_obj)
cts <- AggregateExpression(seurat_obj, 
                           group.by = c("cell_type", "samples"),
                           assays = 'RNA',
                           slot = "counts",
                           return.seurat = FALSE) # pravljenje counts matrice

cts <- cts$RNA # definisanje counts matrice kao RNA eseja 

cts.t <- t(cts) # transpovanje
cts.t <- as.data.frame(cts.t) # konvertovanje u data frame
splitRows <- gsub('_.*', '', rownames(cts.t)) # grupisanje po prefiksu
cts.split <- split.data.frame(cts.t,
                              f = factor(splitRows)) # deljenje po grupama (sample-ovima)
cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
  
}) # sredjivanje napravljenih grupa

counts_astrocytes <- cts.split.modified$astrocytes 
dim(counts_astrocytes) 
View(counts_astrocytes) # definisanje matrice s astrocitima

colData <- data.frame(samples = colnames(counts_astrocytes), stringsAsFactors = FALSE) # definisanje metapodatka
colData <- colData %>%
  mutate(
    condition = case_when(
      grepl("chronic-active-MS-lesion-edge", samples) ~ "chronic_active",
      grepl("MS-periplaque-white-matter", samples) ~ "periplaque",
      TRUE ~ NA_character_
    )
  ) %>%
  column_to_rownames(var = "samples") # uredjivanje metapodataka tako da ostanu samo hronicno aktivne lezije i periplakna bela masa
keep <- !is.na(colData$condition) 
colData <- colData[keep, ] #odbacivanje podataka za koje nema nikakvih rezultata upisanih i ocitanih (NA)
counts_astrocytes <- counts_astrocytes[, rownames(colData)]
dim(counts_astrocytes)
colData <- colData[!is.na(colData$condition), , drop = FALSE]
nrow(colData)
counts_astrocytes <- counts_astrocytes[, rownames(colData), drop = FALSE]
dim(counts_astrocytes)
View(colData)
View(counts_astrocytes) # dodatna formatiranja i provere matrice astrocita i metapodataka

# III DESeq2
table(colData$condition) # prebrojavanje vrednosti u koloni condition
library(DESeq2)
dds <- DESeqDataSetFromMatrix(
  countData = counts_astrocytes,
  colData = colData,
  design = ~ condition
)                        # pravljenje DESeq2 objekta
keep <- rowSums(counts(dds) > 0) >= ceiling(0.20 * ncol(dds)) # definisanje prvog filtera gena (minimum 20% uzoraka)
dds <- dds[keep, ] # primena filtera
nrow(dds) # provera broja gena nakon filtracije
dds <- DESeq(dds) # pokretanje analize
res <- results(
  dds,
  name = "condition_periplaque_vs_chronic_active"
) # izvlacenje rezultata analize
res_sig <- res[
  !is.na(res$padj) &
    res$padj < 0.05 &
    abs(res$log2FoldChange) > 1,
] # filtriranje znacajnih gena 
head(res_sig) # prikaz prvih 10 gena koji su prosli analizu
res_sig["LINC00958",]
res_mn["LINC00958",] # provera postojanja LINC00958 u rezultatima (filtriranim znacajnim i nefiltriranim)
res_sig_df <- as.data.frame(res_sig)
res_sig_df <- as.data.frame(res_sig) # definisanje dataframe-a za znacajne rezultate
View(res_sig_df)
summary(res_sig) 
sum(!is.na(res$padj) & res$padj < 0.05)
table(colData$condition)
dim(counts_astrocytes)
head(colnames(counts_astrocytes))
sum(rowSums(counts_astrocytes) > 0) ## ponavljanje koda i provere zbog nelogicnog rezultata

rm(dds, keep, res, res_sig, res_sig_df) 
exists("dds") 


dds <- DESeqDataSetFromMatrix(
  countData = counts_astrocytes,
  colData = colData,
  design = ~ condition
)
dim(dds)

keep <- rowSums(counts(dds) > 0) >= ceiling(0.20 * ncol(dds))
sum(keep)
dds <- dds[keep, ]
nrow(dds)

dds <- DESeqDataSetFromMatrix(
  countData = counts_astrocytes,
  colData = colData,
  design = ~ condition
)
dim(dds)
keep <- rowSums(counts(dds) > 0) >= ceiling(0.20 * ncol(dds))
sum(keep)
dds <- dds[keep, ]
nrow(dds)
dds <- DESeq(dds)
res <- results(dds, name = "condition_periplaque_vs_chronic_active")
summary(res)

table(colData$NBB_case, colData$condition)
rownames(colData)
colData$NBB_case <- sub('.*-(\\d{2}-\\d{3})$', '\\1', rownames(colData))
colData$NBB_case
length(colData$NBB_case)
nrow(colData)
str(colData)

table(seurat_obj$pathology)
colData <- colData %>%
  mutate(
    condition = case_when(
      grepl("chronic_active_MS_lesion_edge", samples) ~ "chronic_active",
      grepl("MS_periplaque_white_matter", samples) ~ "periplaque",
      grepl("control_white_matter", samples) ~ "control",
      TRUE ~ NAcharacter
    )
  ) %>%
  column_to_rownames(var = "samples")
colnames(counts_astrocytes)

colData <- data.frame(samples = colnames(counts_astrocytes), stringsAsFactors = FALSE)

colData <- colData %>%
  mutate(
    condition = case_when(
      grepl("chronic-active-MS-lesion-edge", samples) ~ "chronic_active",
      grepl("MS-periplaque-white-matter", samples) ~ "periplaque",
      grepl("control-white-matter", samples) ~ "control",
      TRUE ~ NA_character_
    )
  ) %>%
  column_to_rownames(var = "samples")

table(colData$condition, useNA = "ifany")

counts_astrocytes <- cts.split.modified$astrocytes
dim(counts_astrocytes)
colnames(counts_astrocytes) 
colData <- colData %>%
  mutate(
    condition = case_when(
      grepl("chronic-active-MS-lesion-edge", samples) ~ "chronic_active",
      grepl("MS-periplaque-white-matter", samples) ~ "periplaque",
      grepl("control-white-matter", samples) ~ "control",
      TRUE ~ NAcharacter
    )
  ) %>%
  column_to_rownames(var = "samples")

table(colData$condition, useNA = "ifany")

colData <- colData %>%
  mutate(
    condition = case_when(
      grepl("chronic-active-MS-lesion-edge", samples) ~ "chronic_active",
      grepl("MS-periplaque-white-matter", samples) ~ "periplaque",
      grepl("control-white-matter", samples) ~ "control",
      TRUE ~ NA_character_
    )
  ) %>%
  column_to_rownames(var = "samples")

colData <- data.frame(samples = colnames(counts_astrocytes), stringsAsFactors = FALSE)

colData <- colData %>%
  mutate(
    condition = case_when(
      grepl("chronic-active-MS-lesion-edge", samples) ~ "chronic_active",
      grepl("MS-periplaque-white-matter", samples) ~ "periplaque",
      grepl("control-white-matter", samples) ~ "control",
      TRUE ~ NAcharacter
    )
  ) %>%
  column_to_rownames(var = "samples")

table(colData$condition, useNA = "ifany")

dim(counts_astrocytes)
colnames(counts_astrocytes)
colData <- data.frame(samples = colnames(counts_astrocytes), stringsAsFactors = FALSE)
colnames(colData)
str(colData)

colData <- data.frame(samples = colnames(counts_astrocytes), stringsAsFactors = FALSE)

colData <- colData %>%
  mutate(
    condition = case_when(
      grepl("chronic-active-MS-lesion-edge", samples) ~ "chronic_active",
      grepl("MS-periplaque-white-matter", samples) ~ "periplaque",
      grepl("control-white-matter", samples) ~ "control",
      TRUE ~ NAcharacter
    )
  ) %>%
  column_to_rownames(var = "samples")

table(colData$condition, useNA = "ifany")


colData <- data.frame(samples = colnames(counts_astrocytes), stringsAsFactors = FALSE)

condition <- rep(NA, nrow(colData))
condition[grepl("chronic-active-MS-lesion-edge", colData$samples)] <- "chronic_active"
condition[grepl("MS-periplaque-white-matter", colData$samples)] <- "periplaque"
condition[grepl("control-white-matter", colData$samples)] <- "control"

colData$condition <- condition
rownames(colData) <- colData$samples
colData$samples <- NULL

table(colData$condition, useNA = "ifany")
keep <- !is.na(colData$condition)
colData <- colData[keep, ]
counts_astrocytes_filtered <- counts_astrocytes[, rownames(colData)]

dim(counts_astrocytes_filtered)   
table(colData) 

colData <- data.frame(samples = colnames(counts_astrocytes), stringsAsFactors = FALSE)

condition <- rep(NA, nrow(colData))
condition[grepl("chronic-active-MS-lesion-edge", colData$samples)] <- "chronic_active"
condition[grepl("MS-periplaque-white-matter", colData$samples)] <- "periplaque"
condition[grepl("control-white-matter", colData$samples)] <- "control"

colData$condition <- condition
rownames(colData) <- colData$samples
colData$samples <- NULL

class(colData)   # provera — mora da piše "data.frame"

keep <- !is.na(colData$condition)
colData <- colData[keep, , drop = FALSE]      # <- dodato drop = FALSE
counts_astrocytes_filtered <- counts_astrocytes[, rownames(colData)]

class(colData)   # opet provera
dim(colData)

colData$condition <- factor(colData$condition)
colData$condition <- relevel(colData$condition, ref = "control")

dds <- DESeqDataSetFromMatrix(
  countData = counts_astrocytes_filtered,
  colData = colData,
  design = ~ condition
)

keep <- rowSums(counts(dds) > 0) >= ceiling(0.20 * ncol(dds))
dds <- dds[keep, ]
nrow(dds)

dds <- DESeq(dds)
resultsNames(dds)
res_chronic_vs_control <- results(dds, name = "condition_chronic_active_vs_control")
summary(res_chronic_vs_control)

res_periplaque_vs_control <- results(dds, name = "condition_periplaque_vs_control")
summary(res_periplaque_vs_control)

res_chronic_sig <- as.data.frame(res_chronic_vs_control) %>%
  filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
nrow(res_chronic_sig)

res_periplaque_sig <- as.data.frame(res_periplaque_vs_control) %>%
  filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
nrow(res_periplaque_sig)
head(res_chronic_sig[order(res_chronic_sig$padj), ], 10)

### ponavljanje cele analize ali uz kasnije dodavanje kontrolne grupe (zdrava bela masa)

library(ggplot2)
volcano_chronic <- as.data.frame(res_chronic_vs_control)
volcano_chronic$gene <- rownames(volcano_chronic)
volcano_chronic$sig <- case_when(
  !is.na(volcano_chronic$padj) & volcano_chronic$padj < 0.05 & volcano_chronic$log2FoldChange > 1 ~ "Up",
  !is.na(volcano_chronic$padj) & volcano_chronic$padj < 0.05 & volcano_chronic$log2FoldChange < -1 ~ "Down",
  TRUE ~ "NS"
)

ggplot(volcano_chronic, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_classic() +
  labs(
    title = "Volcano plot: Chronic active vs Control",
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value",
    color = "Regulation"
  )

volcano_chronic <- as.data.frame(res_chronic_vs_control)
volcano_chronic$gene <- rownames(volcano_chronic)
volcano_chronic$sig <- case_when(
  !is.na(volcano_chronic$padj) & volcano_chronic$padj < 0.05 & volcano_chronic$log2FoldChange > 1 ~ "Up",
  !is.na(volcano_chronic$padj) & volcano_chronic$padj < 0.05 & volcano_chronic$log2FoldChange < -1 ~ "Down",
  TRUE ~ "NS"
)

volcano_periplaque <- as.data.frame(res_periplaque_vs_control)
volcano_periplaque$gene <- rownames(volcano_periplaque)
volcano_periplaque$sig <- case_when(
  !is.na(volcano_periplaque$padj) & volcano_periplaque$padj < 0.05 & volcano_periplaque$log2FoldChange > 1 ~ "Up",
  !is.na(volcano_periplaque$padj) & volcano_periplaque$padj < 0.05 & volcano_periplaque$log2FoldChange < -1 ~ "Down",
  TRUE ~ "NS"
)

ggplot(volcano_periplaque, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_classic() +
  labs(
    title = "Volcano plot: Periplaque vs Control",
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value",
    color = "Regulation"
  )
# pravljenje VolcanoPlotova za vizuelizaciju rezultata


"LINC00958" %in% rownames(res_chronic_sig)
"LINC00958" %in% rownames(res_periplaque_sig)
"LINC00486" %in% rownames(res_chronic_sig)
"LINC00486" %in% rownames(res_periplaque_sig)






coldata_copy <- as.data.frame(colData(dds))
izbaceni <- coldata_copy[rownames(coldata_copy) %in% c("chronic-active-MS-lesion-edge-13-047", "MS-periplaque-white-matter-13-015", "control-white-matter-11-69", 
                                                      "control-white-matter-12-002","control-white-matter-14-043"), ]
coldata_copy <- coldata_copy[
  !rownames(coldata_copy) %in% rownames(izbaceni),
]
coldata_copy_df <- as.data.frame(coldata_copy)

ddsproba <- DESeqDataSetFromMatrix(
  countData = counts_astrocytes_copy,
  colData = coldata_copy,
  design = ~ condition
)
keep <- rowSums(counts(ddsproba) > 0) >= ceiling(0.20 * ncol(ddsproba))
ddsproba <- ddsproba[keep, ]
nrow(ddsproba)
ddsproba <- DESeq(ddsproba)
resultsNames(ddsproba)
res_chronic_vs_periplaque <- results(ddsproba, name = "condition_periplaque_vs_chronic_active")
summary(res_chronic_vs_periplaque)

ncol(counts_astrocytes)   
nrow(coldata_copy)        
counts_astrocytes_copy <- coldata_copy[colnames(coldata_copy), ]
counts_astrocytes_copy <- counts_astrocytes[, colnames(counts_astrocytes) %in% rownames(coldata_copy)]

res_periplchron_sig <- as.data.frame(res_chronic_vs_periplaque) %>%
  filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
nrow(res_periplchron_sig)
head(res_periplchron_sig[order(res_periplchron_sig$padj), ], 10)

library(ggplot2)
volcano_periplchron <- as.data.frame(res_chronic_vs_periplaque)
volcano_periplchron$gene <- rownames(volcano_periplchron)
volcano_periplchron$sig <- case_when(
  !is.na(volcano_periplchron$padj) & volcano_periplchron$padj < 0.05 & volcano_periplchron$log2FoldChange > 1 ~ "Up",
  !is.na(volcano_periplchron$padj) & volcano_periplchron$padj < 0.05 & volcano_periplchron$log2FoldChange < -1 ~ "Down",
  TRUE ~ "NS"
)

ggplot(volcano_periplchron, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_classic() +
  labs(
    title = "Volcano plot: Chronic active vs Periplaque",
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value",
    color = "Regulation"
  )
# ponovno isprobavanje  chronic active vs periplaque