library(data.table)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(GEOquery)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)


#Fetch metadata for GSE89408
gse <- getGEO("GSE89408", GSEMatrix = TRUE)
pheno <- pData(gse[[1]])
colnames(pheno)
head(pheno)

#Open the raw count dataset
count_file <- "C:/Users/Matthew/Downloads/GSE89408_raw_counts_GRCh38.p13_NCBI.tsv.gz"
counts <- fread(count_file)

#process count matrix
gene_id <- counts[[1]]
count_matrix <- as.data.frame(counts[, -1])
rownames(count_matrix) <- gene_id

#confirming data
head(count_matrix[, 1:3])
dim(count_matrix)

#Check that GSM IDs match your count matrix columns
colnames(count_matrix)[1:6]
rownames(pheno)[1:6]   

#Reorder metadata to match your count matrix column order
pheno_matched <- pheno[colnames(count_matrix), ]
# sanity check - should return TRUE
all(rownames(pheno_matched) == colnames(count_matrix))

#separate samples
sample_info <- data.frame(
  row.names = colnames(count_matrix),
  condition = pheno_matched[["disease:ch1"]]
)

#Check how many conditions is inside the dataset
table(sample_info$condition)

#Keep only RA early, RA established, Normal
keep_groups <- c("Normal", "Rheumatoid arthritis (early)", "Rheumatoid arthritis (established)")
sample_info <- sample_info[sample_info$condition %in% keep_groups, , drop = FALSE]
table(sample_info$condition)

#Subset count_matrix to match the filtered samples
count_matrix_filtered <- count_matrix[, rownames(sample_info)]
all(colnames(count_matrix_filtered) == rownames(sample_info))
dim(count_matrix_filtered)

#set explicit order
sample_info$condition <- factor(
  sample_info$condition,
  levels = c("Normal", "Rheumatoid arthritis (early)", "Rheumatoid arthritis (established)")
)

levels(sample_info$condition)   # confirm the 3 groups
table(sample_info$condition)    # confirm counts match filtering

#Make DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_filtered,
  colData = sample_info,
  design = ~ condition
)

# Filter low-count genes
#how many of the 173 samples have at least 1 read, and keeps only the genes where that count is 3 or more
keep <- rowSums(counts(dds) >= 1) >= 3
dds <- dds[keep, ]

#check row number
nrow(count_matrix_filtered)
nrow(dds)

#run DESeq2
dds <- DESeq(dds)

#-------------------------------------------------------------------------------
#get DESeq2 result
#As we comparing 3 group so we need to make 3 result getting
#Normal vs early, established vs normal, early vs established
#-------------------------------------------------------------------------------
#Early vs normal
res_early_vs_normal <- results(
  dds,
  contrast = c("condition", "Rheumatoid arthritis (early)", "Normal")
)
#make into data.frame
res_df_EarlyNormal <- as.data.frame(res_early_vs_normal)
res_df_EarlyNormal$gene_id <- rownames(res_df_EarlyNormal)
#in order of padj
res_df_EarlyNormal <- res_df_EarlyNormal[
  order(res_df_EarlyNormal$padj),
]
head(res_df_EarlyNormal, 20)
#-------------------------------------------------------------------------------
#Established vs normal
res_estab_vs_normal <- results(
  dds,
  contrast = c("condition", "Rheumatoid arthritis (established)", "Normal")
)
#make into data.frame
res_df_EstabNormal <- as.data.frame(res_estab_vs_normal)
res_df_EstabNormal$gene_id <- rownames(res_df_EstabNormal)
#in order of padj
res_df_EstabNormal <- res_df_EstabNormal[
  order(res_df_EstabNormal$padj),
]
head(res_df_EstabNormal, 20)
#-------------------------------------------------------------------------------
#Established vs early
res_estab_vs_early <- results(
  dds,
  contrast = c("condition", "Rheumatoid arthritis (established)", "Rheumatoid arthritis (early)")
)
#make into data.frame
res_df_EstabEarly <- as.data.frame(res_estab_vs_early)
res_df_EstabEarly$gene_id <- rownames(res_df_EstabEarly)
#in order of padj
res_df_EstabEarly <- res_df_EstabEarly[
  order(res_df_EstabEarly$padj),
]
head(res_df_EstabEarly, 20)
#-------------------------------------------------------------------------------
#get the DEG (different expressed genes)
#We will separate each group DEG, so later we can try to analyse each group
#res_early_vs_normal       (early RA-specific changes)
#res_established_vs_normal (established RA-specific changes)
#res_established_vs_early  (disease progression changes)
#-------------------------------------------------------------------------------
#Early vs normal
DEG_Early_Normal <- subset(
  res_df_EarlyNormal,
  !is.na(padj) &
    padj < 0.05 &
    abs(log2FoldChange) >= 1
)
up_DEG_Early_Normal <- subset(
  DEG_Early_Normal,
  log2FoldChange >= 1
)

down_DEG_Early_Normal <- subset(
  DEG_Early_Normal,
  log2FoldChange <= -1
)
cat("Up-regulated:", nrow(up_DEG_Early_Normal), "\n")
cat("Down-regulated:", nrow(down_DEG_Early_Normal), "\n")
cat("Total DEG:", nrow(DEG_Early_Normal), "\n")
#-------------------------------------------------------------------------------
#Normal vs Established
DEG_Estab_Normal <- subset(
  res_df_EstabNormal,
  !is.na(padj) &
    padj < 0.05 &
    abs(log2FoldChange) >= 1
)
up_DEG_Estab_Normal <- subset(
  DEG_Estab_Normal,
  log2FoldChange >= 1
)

down_DEG_Estab_Normal <- subset(
  DEG_Estab_Normal,
  log2FoldChange <= -1
)
cat("Up-regulated:", nrow(up_DEG_Estab_Normal), "\n")
cat("Down-regulated:", nrow(down_DEG_Estab_Normal), "\n")
cat("Total DEG:", nrow(DEG_Estab_Normal), "\n")
#-------------------------------------------------------------------------------
#Established vs Early
DEG_Estab_Early <- subset(
  res_df_EstabEarly,
  !is.na(padj) &
    padj < 0.05 &
    abs(log2FoldChange) >= 1
)
up_DEG_Estab_Early <- subset(
  DEG_Estab_Early,
  log2FoldChange >= 1
)

down_DEG_Estab_Early <- subset(
  DEG_Estab_Early,
  log2FoldChange <= -1
)
cat("Up-regulated:", nrow(up_DEG_Estab_Early), "\n")
cat("Down-regulated:", nrow(down_DEG_Estab_Early), "\n")
cat("Total DEG:", nrow(DEG_Estab_Early), "\n")
#-------------------------------------------------------------------------------
#We will use Volcano plot to visualize our DEG results per group
#-------------------------------------------------------------------------------
#Early vs Normal
res_df_EarlyNormal$significance <- "Not significant"

res_df_EarlyNormal$significance[
  res_df_EarlyNormal$padj < 0.05 &
    res_df_EarlyNormal$log2FoldChange >= 1
] <- "Up"

res_df_EarlyNormal$significance[
  res_df_EarlyNormal$padj < 0.05 &
    res_df_EarlyNormal$log2FoldChange <= -1
] <- "Down"

#Volcano plot - 畫圖
ggplot(
  res_df_EarlyNormal,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = significance
  )
) +
  scale_color_manual(
    values = c(
      "Not significant" = "grey",
      "Up" = "red",
      "Down" = "blue"
    )
  ) +
  geom_point(alpha = 0.6) +
  theme_classic() +
  labs(
    title = "Early vs Normal",
    x = "log2 Fold Change",
    y = "-log10(adjusted P-value)"
  )
#-------------------------------------------------------------------------------
#Established vs Normal
res_df_EstabNormal$significance <- "Not significant"

res_df_EstabNormal$significance[
  res_df_EstabNormal$padj < 0.05 &
    res_df_EstabNormal$log2FoldChange >= 1
] <- "Up"

res_df_EstabNormal$significance[
  res_df_EstabNormal$padj < 0.05 &
    res_df_EstabNormal$log2FoldChange <= -1
] <- "Down"

#Volcano plot - 畫圖
ggplot(
  res_df_EstabNormal,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = significance
  )
) +
  scale_color_manual(
    values = c(
      "Not significant" = "grey",
      "Up" = "red",
      "Down" = "blue"
    )
  ) +
  geom_point(alpha = 0.6) +
  theme_classic() +
  labs(
    title = "Established vs Normal",
    x = "log2 Fold Change",
    y = "-log10(adjusted P-value)"
  )
#-------------------------------------------------------------------------------
#Early vs Established
res_df_EstabEarly$significance <- "Not significant"

res_df_EstabEarly$significance[
  res_df_EstabEarly$padj < 0.05 &
    res_df_EstabEarly$log2FoldChange >= 1
] <- "Up"

res_df_EstabEarly$significance[
  res_df_EstabEarly$padj < 0.05 &
    res_df_EstabEarly$log2FoldChange <= -1
] <- "Down"

#Volcano plot - 畫圖
ggplot(
  res_df_EstabEarly,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = significance
  )
) +
  scale_color_manual(
    values = c(
      "Not significant" = "grey",
      "Up" = "red",
      "Down" = "blue"
    )
  ) +
  geom_point(alpha = 0.6) +
  theme_classic() +
  labs(
    title = "Established vs Early",
    x = "log2 Fold Change",
    y = "-log10(adjusted P-value)"
  )
#-------------------------------------------------------------------------------
#PCA Analysis
vsd <- vst(dds, blind = FALSE)
p <- plotPCA(
  vsd,
  intgroup = "condition"
)
p + 
  ylim(-20, 20)
#-------------------------------------------------------------------------------
#Write the DEG into a csv file
#-------------------------------------------------------------------------------
#Early vs Normal DEG
write.csv(
  DEG_Early_Normal,
  "GSE177038_DEG_Early_Normal.csv",
  row.names = FALSE
)

#Established vs Normal DEG
write.csv(
  DEG_Estab_Normal,
  "GSE177038_DEG_Established_Normal.csv",
  row.names = FALSE
)

#Early vs Established
write.csv(
  DEG_Estab_Early,
  "GSE177038_DEG_Established_Early.csv",
  row.names = FALSE
)
#-------------------------------------------------------------------------------
#GO enrichment
#-------------------------------------------------------------------------------
#confirm DEG ID is Entrez ID
#Entrez Gene IDs are plain numeric integers (e.g., 672, 7157, 1017) — no letters
#no prefixes like ENSG (Ensembl) or NM_/NP_ (RefSeq).
#We will check it just to be sure

ego <- enrichGO(
  gene = rownames(DEG),
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

head(ego)
#畫圖
dotplot(
  ego,
  showCategory = 20
) +
  ggtitle("GO Biological Process")

#####KEGG enrichment
ekegg <- enrichKEGG(
  gene = rownames(DEG),
  organism = "hsa",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)
#畫圖
dotplot(
  ekegg,
  showCategory = 20
) +
  ggtitle("KEGG Pathway")


#####GSEA

##建立 GSEA ranking
#使用stat
#因為 stat 是 DESeq2 的 Wald statistic，
#同時考慮了 effect size 和 variability，
#比單純按照 log2FoldChange 排序更適合 GSEA。
gene_rank <- res_df$stat
names(gene_rank) <- rownames(res_df)

gene_rank <- gene_rank[!is.na(gene_rank)]
gene_rank <- sort(gene_rank, decreasing = TRUE)

##GO Biological Process
gsea_go <- gseGO(
  geneList = gene_rank,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

#畫圖
dotplot(
  gsea_go,
  showCategory = 20,
  title = "GSEA - GO Biological Process"
)

#選特定GO
head(gsea_go@result[, c(
  "ID",
  "Description",
  "NES",
  "p.adjust"
)])

#畫單一 pathway 的 GSEA curve
gseaplot2(
  gsea_go,
  geneSetID = "GO:0033489",
  title = gsea_go$Description[1]
)

gseaplot2(
  gsea_go,
  geneSetID = "GO:0002483",
  title = gsea_go$Description[1]
)


###KEGG GSEA
gsea_kegg <- gseKEGG(
  geneList = gene_rank,
  organism = "hsa",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)
#畫圖
dotplot(
  gsea_kegg,
  showCategory = 20,
  title = "GSEA - KEGG"
)

#畫單一 pathway 的 GSEA curve
gseaplot2(
  gsea_kegg,
  geneSetID = "hsa04382",
  title = "Cornified envelope formation"
)

