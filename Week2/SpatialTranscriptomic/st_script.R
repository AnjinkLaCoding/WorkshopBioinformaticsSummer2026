# ==============================================================================
# TBSB Week 2 Workshop: Spatial Transcriptomics Analysis with 10x Visium & Seurat
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

# Load required packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Set random seed for reproducibility
set.seed(1234)


# ------------------------------------------------------------------------------
# 1. Load the spatial transcriptomics dataset
# ------------------------------------------------------------------------------
## Read the 10x Visium data
breast <- Load10X_Spatial(
  data.dir = "st",
  filename = "Visium_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "breast_ffpe",
  filter.matrix = TRUE
)

breast


# ------------------------------------------------------------------------------
# 2. Quality control
# ------------------------------------------------------------------------------
## Visualize spot-level QC metrics
VlnPlot(
  breast,
  features = c(
    "nCount_Spatial",
    "nFeature_Spatial"
  ),
  pt.size = 0.1,
  ncol = 2
)


# ------------------------------------------------------------------------------
# 3. Normalization with SCTransform
# ------------------------------------------------------------------------------
breast <- SCTransform(
  breast,
  assay = "Spatial",
  new.assay.name = "SCT",
  vst.flavor = "v2",
  verbose = FALSE
)

## Inspect available assays
Assays(breast)


# ------------------------------------------------------------------------------
# 4. Principal component analysis (PCA)
# ------------------------------------------------------------------------------
## Run PCA
breast <- RunPCA(
  breast,
  assay = "SCT",
  npcs = 50,
  verbose = FALSE
)

## Select informative principal components (Elbow Plot)
ElbowPlot(
  breast,
  ndims = 50
)


# ------------------------------------------------------------------------------
# 5. Graph construction and clustering
# ------------------------------------------------------------------------------
## Construct the neighbor graph
breast <- FindNeighbors(
  breast,
  reduction = "pca",
  dims = 1:20,
  verbose = FALSE
)

## Identify clusters
breast <- FindClusters(
  breast,
  resolution = 0.3,
  verbose = FALSE
)

## Inspect cluster sizes
table(Idents(breast))


# ------------------------------------------------------------------------------
# 6. UMAP visualization
# ------------------------------------------------------------------------------
## Run UMAP
breast <- RunUMAP(
  breast,
  reduction = "pca",
  dims = 1:20,
  seed.use = 1234,
  verbose = FALSE
)

## Visualize clusters in UMAP space
DimPlot(
  breast,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  label.size = 4
) +
  ggtitle("UMAP clusters")


# ------------------------------------------------------------------------------
# 7. Visualize clusters in tissue space
# ------------------------------------------------------------------------------
SpatialDimPlot(
  breast,
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  pt.size.factor = 1.6,
  image.alpha = 1
) +
  ggtitle("Spatial distribution of clusters")


# ------------------------------------------------------------------------------
# 8. Compare cluster sizes
# ------------------------------------------------------------------------------
cluster_size <- as.data.frame(
  table(
    cluster = Idents(breast)
  )
)

cluster_size

## Plot the number of spots per cluster
ggplot(
  cluster_size,
  aes(
    x = cluster,
    y = Freq
  )
) +
  geom_col() +
  xlab("Cluster") +
  ylab("Number of spots") +
  ggtitle("Number of spots per cluster") +
  theme_classic()


# ------------------------------------------------------------------------------
# 9. Marker genes handling
# ------------------------------------------------------------------------------
## NOTE: FindAllMarkers is skipped here as it takes significant time.
## If you need to recompute, uncomment the following block:
# breast_markers <- FindAllMarkers(
#     breast,
#     assay = "SCT",
#     only.pos = TRUE,
#     min.pct = 0.25,
#     logfc.threshold = 0.25,
#     test.use = "wilcox",
#     verbose = TRUE
# )
# write.csv(breast_markers, "./output/breast_markers.csv")

## Load precomputed marker genes
breast_markers <- read.csv(
  "./output/breast_markers.csv",
  row.names = 1
)

head(breast_markers)
table(breast_markers$cluster)

## Filter significant marker genes
breast_markers_filtered <- breast_markers %>%
  filter(
    p_val_adj < 0.05,
    avg_log2FC > 0.25,
    pct.1 > 0.25
  )

head(breast_markers_filtered)

## Select the top 10 marker genes per cluster
breast_topN <- breast_markers_filtered %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10,
    with_ties = FALSE
  ) %>%
  ungroup()

breast_topN

## Collect unique marker genes
topN_genes <- unique(
  breast_topN$gene
)

topN_genes

## Summarize markers by cluster
marker_summary <- breast_topN %>%
  group_by(cluster) %>%
  summarise(
    top_markers = paste(
      gene,
      collapse = ", "
    ),
    .groups = "drop"
  )

marker_summary


# ------------------------------------------------------------------------------
# 10. Visualize marker expression
# ------------------------------------------------------------------------------
## DotPlot of top cluster markers
DotPlot(
  breast,
  features = topN_genes,
  assay = "SCT"
) +
  RotatedAxis() +
  ggtitle(
    "Top marker genes of spatial clusters"
  )

## Visualize gene expression in tissue space (top 4 genes)
SpatialFeaturePlot(
  breast,
  features = topN_genes[1:4],
  slot = "data",
  ncol = 2,
  pt.size.factor = 1.6,
  min.cutoff = "q05",
  max.cutoff = "q95"
)


# ------------------------------------------------------------------------------
# 11. Annotate spatial domains
# ------------------------------------------------------------------------------
## Assign broad biological identities
breast <- RenameIdents(
  breast,
  `0` = "Fibroblast / CAF",
  `1` = "Epithelial / Tumor",
  `2` = "Macrophage / Myeloid",
  `3` = "Macrophage / Myeloid",
  `4` = "Macrophage / Myeloid",
  `5` = "Lymphoid",
  `6` = "Fibroblast / CAF"
)

breast$spatial_domain <- Idents(breast)

## Inspect annotated domain sizes
table(breast$spatial_domain)

## Visualize annotated spatial domains on tissue map
SpatialDimPlot(
  breast,
  group.by = "spatial_domain",
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  pt.size.factor = 1.6
)

## Visualize annotated domains in UMAP space
DimPlot(
  breast,
  reduction = "umap",
  group.by = "spatial_domain",
  label = TRUE,
  repel = TRUE,
  label.size = 4
) +
  ggtitle("Annotated spatial domains")


# ------------------------------------------------------------------------------
# 12. Save the analyzed Seurat object
# ------------------------------------------------------------------------------
saveRDS(
  breast,
  file = "./output/breast_visium_annotated.rds"
)
