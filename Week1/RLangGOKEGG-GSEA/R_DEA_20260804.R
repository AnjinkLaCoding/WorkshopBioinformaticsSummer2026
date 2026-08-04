#載入已安裝的套件
library(data.table)
library(DESeq2)
library(ggplot2)
library(pheatmap)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)


#讀取檔案
count_file <- "C:/Users/Matthew/Downloads/GSE177038_raw_counts_GRCh38.p13_NCBI.tsv.gz"
counts <- fread(count_file)
#確認資料
dim(counts) #查列數與欄數
head(counts) #列前10行
colnames(counts) #欄位名稱

#建立count matrix
gene_id <- counts[[1]]
count_matrix <- as.data.frame(counts[, -1])
rownames(count_matrix) <- gene_id

#確認資料
head(count_matrix[, 1:3])
dim(count_matrix)

#分samples
sample_info <- data.frame(
  row.names = colnames(count_matrix),
  condition = c(
    "A2780",
    "A2780",
    "A2780",
    "A2780_ADR",
    "A2780_ADR",
    "A2780_ADR"
  )
)

#A2780_ADR vs A2780
sample_info$condition <- factor(
  sample_info$condition,
  levels = c("A2780", "A2780_ADR")
)
#log2FoldChange > 0 == "A2780-ADR higher than A2780"

#建立DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ condition
)

# Filter low-count genes
#3 個 samples 有 count ≥ 1
keep <- rowSums(counts(dds) >= 1) >= 3
dds <- dds[keep, ]

#確認篩選資料
nrow(count_matrix)
nrow(dds)

#執行DESeq2
dds <- DESeq(dds)

#取得DESeq2結果
res <- results(
  dds,
  contrast = c("condition", "A2780_ADR", "A2780")
)

#整理成 data.frame
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)


#按照padj排序
res_df <- res_df[
  order(res_df$padj),
]
#顯示最顯著top20
head(res_df, 20)


#定義 DEG
DEG <- subset(
  res_df,
  !is.na(padj) &
    padj < 0.05 &
    abs(log2FoldChange) >= 1
)
up_DEG <- subset(
  DEG,
  log2FoldChange >= 1
)

down_DEG <- subset(
  DEG,
  log2FoldChange <= -1
)
#輸出數量
cat("Up-regulated:", nrow(up_DEG), "\n")
cat("Down-regulated:", nrow(down_DEG), "\n")
cat("Total DEG:", nrow(DEG), "\n")

#MA plot
plotMA(
  res,
  ylim = c(-5, 5)
)

#Volcano plot - 建立分類
res_df$significance <- "Not significant"

res_df$significance[
  res_df$padj < 0.05 &
    res_df$log2FoldChange >= 1
] <- "Up"

res_df$significance[
  res_df$padj < 0.05 &
    res_df$log2FoldChange <= -1
] <- "Down"

#Volcano plot - 畫圖
ggplot(
  res_df,
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
    title = "A2780-ADR vs A2780",
    x = "log2 Fold Change",
    y = "-log10(adjusted P-value)"
  )



#PCA分析
vsd <- vst(dds, blind = FALSE)
p <- plotPCA(
  vsd,
  intgroup = "condition"
)
p + 
  ylim(-20, 20)

#輸出DEG
write.csv(
  res_df,
  "GSE177038_DEG_all.csv",
  row.names = FALSE
)

#####GO enrichment

#確認 DEG ID 是Entrez ID

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

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

