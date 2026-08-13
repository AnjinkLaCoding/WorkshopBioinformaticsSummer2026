# ==============================================================================
# TBSB Week 2 Workshop: Single-Cell RNA-seq Analysis with Seurat
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Global Setup & Directory Configuration
# ------------------------------------------------------------------------------
# Set the project directory used by all code chunks
project_dir <- path.expand("~/TBSB_W2")

# Create output directory
dir.create(
  file.path(project_dir, "output"),
  showWarnings = FALSE,
  recursive = TRUE
)

# Set working directory
setwd(project_dir)


# ------------------------------------------------------------------------------
# 1. Set directory and load required packages
# ------------------------------------------------------------------------------
# Load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Set random seed for reproducibility
set.seed(1234)


# ------------------------------------------------------------------------------
# 2. Data loading
# ------------------------------------------------------------------------------
## 2.1 Load the PBMC 3K count matrix
pbmc.data <- Read10X(
  data.dir = "sc/filtered_gene_bc_matrices/hg19"
)

## 2.2 Inspect the raw matrix
class(pbmc.data)
dim(pbmc.data)
pbmc.data[1:10, 1:3]


# ------------------------------------------------------------------------------
# 3. Create a Seurat object
# ------------------------------------------------------------------------------
pbmc <- CreateSeuratObject(
  counts = pbmc.data,
  project = "PBMC3K",
  min.cells = 3,
  min.features = 200
)

print(pbmc)


# ------------------------------------------------------------------------------
# 4. Quality control (QC)
# ------------------------------------------------------------------------------
## 4.1 Inspect cell metadata
head(pbmc[[]])

## 4.2 Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(
  object = pbmc,
  pattern = "^MT-"
)

head(pbmc[[]])

## 4.3 Summarize QC metrics
summary(pbmc$nCount_RNA)
summary(pbmc$nFeature_RNA)
summary(pbmc$percent.mt)

## 4.4 Visualize QC metrics before filtering
p.qc.before <- VlnPlot(
  object = pbmc,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3,
  pt.size = 0.01
)

p.qc.before

## Optional: save the QC figure
ggsave(
  filename = "output/01_violin_pre_QC.png",
  plot = p.qc.before,
  width = 10,
  height = 5,
  dpi = 300
)

## 4.5 Filter low-quality cells
pbmc <- subset(
  x = pbmc,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 2500 &
    percent.mt < 5
)

pbmc

## 4.6 Visualize QC metrics after filtering
VlnPlot(
  object = pbmc,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3,
  pt.size = 0.1
)


# ------------------------------------------------------------------------------
# 5. Normalization
# ------------------------------------------------------------------------------
pbmc <- NormalizeData(
  object = pbmc,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

## 5.1 Inspect Seurat v5 layers
Layers(pbmc[["RNA"]])

## 5.2 Inspect normalized expression values
pbmc[["RNA"]]$data[
  c("CD3D", "MS4A1", "LYZ"),
  1:5
]


# ------------------------------------------------------------------------------
# 6. Identify highly variable genes
# ------------------------------------------------------------------------------
pbmc <- FindVariableFeatures(
  object = pbmc,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

length(VariableFeatures(pbmc))

## 6.1 Visualize highly variable genes
top10.variable <- head(
  VariableFeatures(pbmc),
  10
)

LabelPoints(
  plot = VariableFeaturePlot(pbmc),
  points = top10.variable,
  repel = TRUE
)


# ------------------------------------------------------------------------------
# 7. Scale the data
# ------------------------------------------------------------------------------
pbmc <- ScaleData(
  object = pbmc,
  features = rownames(pbmc),
  verbose = FALSE
)


# ------------------------------------------------------------------------------
# 8. Principal component analysis (PCA)
# ------------------------------------------------------------------------------
## 8.1 Run PCA
pbmc <- RunPCA(
  object = pbmc,
  features = VariableFeatures(pbmc),
  npcs = 30,
  verbose = FALSE
)

## 8.2 Inspect genes contributing to major PCs
print(
  pbmc[["pca"]],
  dims = 1:5,
  nfeatures = 5
)

## 8.3 Visualize cells in PCA space
DimPlot(
  object = pbmc,
  reduction = "pca"
) +
  NoLegend()

## 8.4 Determine how many PCs to retain (Elbow Plot)
ElbowPlot(
  object = pbmc,
  ndims = 30
)

## 8.5 PCA heatmap
DimHeatmap(
  object = pbmc,
  dims = 1:10,
  cells = 500,
  balanced = TRUE
)


# ------------------------------------------------------------------------------
# 9. Graph construction and clustering
# ------------------------------------------------------------------------------
## 9.1 Construct the neighbor graph
pbmc <- FindNeighbors(
  object = pbmc,
  dims = 1:10,
  verbose = FALSE
)

## 9.2 Identify clusters
pbmc <- FindClusters(
  object = pbmc,
  resolution = 0.5,
  random.seed = 1234,
  verbose = FALSE
)

## 9.3 Store and inspect cluster IDs
pbmc$cluster_id <- as.character(Idents(pbmc))

table(pbmc$cluster_id)


# ------------------------------------------------------------------------------
# 10. UMAP visualization
# ------------------------------------------------------------------------------
## 10.1 Run UMAP
pbmc <- RunUMAP(
  object = pbmc,
  dims = 1:10,
  seed.use = 1234,
  verbose = FALSE
)

## 10.2 Visualize clusters in UMAP space
DimPlot(
  object = pbmc,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.7
) +
  NoLegend()


# ------------------------------------------------------------------------------
# 11. Identify cluster marker genes
# ------------------------------------------------------------------------------
pbmc.markers <- FindAllMarkers(
  object = pbmc,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

## 11.1 Select the top 10 markers per cluster
top10.markers <- pbmc.markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10,
    with_ties = FALSE
  ) %>%
  ungroup()

top10.markers

## Optional: export marker results
write.csv(
  pbmc.markers,
  file = "output/PBMC3K_all_cluster_markers.csv",
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 12. Cell-type annotation
# ------------------------------------------------------------------------------
## 12.1 Define canonical marker genes
pbmc.marker.panel <- list(
  "Naive CD4+ T"     = c("IL7R", "CCR7"),
  "CD14+ Monocyte"   = c("CD14", "LYZ"),
  "Memory CD4+ T"    = c("IL7R", "S100A4"),
  "B cell"           = c("MS4A1"),
  "CD8+ T"           = c("CD8A"),
  "FCGR3A+ Monocyte" = c("FCGR3A", "MS4A7"),
  "NK cell"          = c("GNLY", "NKG7"),
  "Dendritic cell"   = c("FCER1A", "CST3"),
  "Platelet"         = c("PPBP")
)

## 12.2 Prepare marker genes for visualization
pbmc.marker.plot <- unique(
  unlist(pbmc.marker.panel)
)

pbmc.marker.plot

## 12.3 Visualize canonical marker expression (DotPlot)
DotPlot(
  object = pbmc,
  features = pbmc.marker.plot
) +
  RotatedAxis() +
  labs(
    x = "Marker genes",
    y = "Cluster",
    title = "PBMC3K canonical marker expression"
  )

## 12.4 Inspect selected markers with violin plots
VlnPlot(
  object = pbmc,
  features = c(
    "CD14",
    "CD8A",
    "MS4A7"
  )
)

## 12.5 Assign biological cell-type labels
cluster.annotation <- c(
  "0" = "Naive CD4+ T",
  "1" = "CD14+ Monocyte",
  "2" = "Memory CD4+ T",
  "3" = "B cell",
  "4" = "CD8+ T",
  "5" = "FCGR3A+ Monocyte",
  "6" = "NK cell",
  "7" = "Dendritic cell",
  "8" = "Platelet"
)

pbmc <- RenameIdents(
  object = pbmc,
  cluster.annotation
)

pbmc$celltype <- as.character(Idents(pbmc))

table(pbmc$celltype)


# ------------------------------------------------------------------------------
# 13. Visualize annotated cell types & Save output
# ------------------------------------------------------------------------------
## 13.1 Visualize annotated cell types in UMAP
DimPlot(
  object = pbmc,
  reduction = "umap",
  group.by = "celltype",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.7
) +
  NoLegend()

## 13.2 Save the final Seurat object
saveRDS(
  pbmc,
  file = "output/PBMC3K_annotated_seurat.rds"
)