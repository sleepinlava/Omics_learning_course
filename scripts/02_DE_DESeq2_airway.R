#!/usr/bin/env Rscript
###############################################################################
# 02_DE_DESeq2_airway.R
#
# 完整走一遍 bulk RNA-seq 差异表达分析（DESeq2）
# 数据集: airway (GSE52778) —— 人支气管平滑肌细胞, 地塞米松(dex)处理 vs 未处理
#   4 个细胞系, 每个细胞系都有 untreated / dexamethasone 两种处理 => 配对设计
#
# 学到的技能:
#   1. 认识 SummarizedExperiment / count matrix + metadata
#   2. 构建 DESeqDataSet, 设定 design 公式
#   3. 过滤低表达基因, 归一化(size factor), VST 变换
#   4. QC: PCA、样本距离热图
#   5. 运行差异分析, 提取结果, log2FC 收缩(shrinkage)
#   6. 可视化: MA plot、火山图、top gene 热图
#   7. 导出结果表
#
# 用法:  Rscript scripts/02_DE_DESeq2_airway.R
###############################################################################

suppressPackageStartupMessages({
  library(airway)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})
set.seed(1234)

# 自动定位项目根目录 (scripts/ 的上一级)
.args <- commandArgs(trailingOnly = FALSE)
.f <- sub("^--file=", "", .args[grep("^--file=", .args)])
ROOT <- if (length(.f)) dirname(dirname(normalizePath(.f))) else getwd()
OUT  <- file.path(ROOT, "results"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

###############################################################################
# 1. 读入数据 ---------------------------------------------------------------
###############################################################################
data("airway")                       # SummarizedExperiment
cat("== airway ==\n"); print(airway)

# 表达矩阵 (rows = 基因, cols = 样本)
counts <- assay(airway)
# 样本注释 (colData)
coldata <- as.data.frame(colData(airway))
cat("\n样本信息:\n"); print(coldata[, c("cell", "dex", "Run")])
cat("\n计数矩阵维度:", dim(counts), "\n")

###############################################################################
# 2. 构建 DESeqDataSet 并设定 design ------------------------------------------------
###############################################################################
# design 里放什么很重要: ~ cell + dex 表示"扣除细胞系差异后"看 dex 的效应
dds <- DESeqDataSet(airway, design = ~ cell + dex)

# 预过滤: 至少 4 个样本(共8个)中 count >= 10 的基因才保留 (降维、提速、增稳)
keep <- rowSums(counts(dds) >= 10) >= 4
cat("\n过滤前基因数:", nrow(dds), "-> 过滤后:", sum(keep), "\n")
dds <- dds[keep, ]

# 指定比较的参考水平 (baseline = untreated)
dds$dex <- relevel(dds$dex, ref = "untrt")

###############################################################################
# 3. 归一化 与 方差稳定变换 ----------------------------------------------------
###############################################################################
dds <- estimateSizeFactors(dds)           # 每个样本一个 size factor
cat("\nSize factors:\n"); print(sizeFactors(dds))

vsd <- vst(dds, blind = TRUE)             # 用于 QC / 可视化 (近似 log2)
# 也可以: rlog(dds, blind=TRUE)  —— 小样本时更稳但更慢

norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, file.path(OUT, "airway_normalized_counts.csv"))

###############################################################################
# 4. 质量评估 QC --------------------------------------------------------------
###############################################################################
# 4.1 PCA
pca <- prcomp(t(assay(vsd)))
pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                     cell = dds$cell, dex = dds$dex)
p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = dex, shape = cell)) +
  geom_point(size = 4) +
  labs(title = "PCA (VST-transformed counts)",
       x = paste0("PC1: ", pct[1], "% var"),
       y = paste0("PC2: ", pct[2], "% var")) +
  theme_bw(base_size = 13)
ggsave(file.path(OUT, "airway_PCA.png"), p_pca, width = 6, height = 4.5, dpi = 150)

# 4.2 样本间距离热图
sampleDists <- dist(t(assay(vsd)))
samDistMat  <- as.matrix(sampleDists)
rownames(samDistMat) <- paste(dds$cell, dds$dex, sep = "_")
png(file.path(OUT, "airway_sample_distance_heatmap.png"), width = 800, height = 700)
pheatmap(samDistMat, clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample-to-sample distances (VST)")
dev.off()

###############################################################################
# 5. 差异表达分析 -------------------------------------------------------------
###############################################################################
dds <- DESeq(dds)          # 一步到位: 估计 dispersion + 拟合 GLM + Wald test
res <- results(dds, contrast = c("dex", "trt", "untrt"), alpha = 0.05)
cat("\n=== 差异分析结果 summary ===\n"); print(summary(res))

# log2 fold change 收缩: 降
低低表达基因的假阳性(推荐用于排序/可视化)
resLFC <- lfcShrink(dds, coef = "dex_trt_vs_untrt", type = "apeglm", res = res)

# 按 padj 排序
resOrdered <- resLFC[order(resLFC$padj), ]

# 加上基因名(symbol)
if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  resOrdered$symbol <- mapIds(org.Hs.eg.db, keys = rownames(resOrdered),
                              keytype = "ENSEMBL", column = "SYMBOL")
}

# 导出全部结果 + 显著基因
write.csv(as.data.frame(resOrdered), file.path(OUT, "airway_DESeq2_all_genes.csv"))
sig <- subset(resOrdered, padj < 0.05 & abs(log2FoldChange) > 1)
cat("\n显著基因 (padj<0.05 & |log2FC|>1):", nrow(sig), "\n")
write.csv(as.data.frame(sig), file.path(OUT, "airway_DESeq2_significant.csv"))
cat("\nTop 5 上调/下调:\n"); print(head(resOrdered[order(resOrdered$log2FoldChange, decreasing = TRUE), ], 5))

###############################################################################
# 6. 可视化: MA plot / 火山图 / 热图 -------------------------------------------
###############################################################################
# 6.1 MA plot
png(file.path(OUT, "airway_MA_plot.png"), width = 800, height = 650)
plotMA(resLFC, ylim = c(-6, 6), main = "MA plot (apeglm-shrunk LFC)")
dev.off()

# 6.2 火山图 (EnhancedVolcano 若可用, 否则 ggplot 版本)
if (requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  library(EnhancedVolcano)
  png(file.path(OUT, "airway_volcano.png"), width = 1100, height = 900, res = 130)
  print(EnhancedVolcano(resOrdered, lab = as.character(resOrdered$symbol),
        x = "log2FoldChange", y = "padj",
        pCutoff = 0.05, FCcutoff = 1,
        title = "airway: dexamethasone vs untreated",
        subtitle = "DESeq2, apeglm shrinkage"))
  dev.off()
} else {
  df <- as.data.frame(resOrdered)
  df$sig <- ifelse(!is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) > 1, "sig", "ns")
  p_vol <- ggplot(df, aes(log2FoldChange, -log10(padj), colour = sig)) +
    geom_point(alpha = 0.6, size = 1.2) + theme_bw(base_size = 13) +
    labs(title = "Volcano: dexamethasone vs untreated")
  ggsave(file.path(OUT, "airway_volcano.png"), p_vol, width = 7, height = 5, dpi = 150)
}

# 6.3 显著基因表达热图
top_ids <- head(rownames(resOrdered)[!is.na(resOrdered$padj)], 40)
mat <- assay(vsd)[top_ids, ]
mat <- t(scale(t(mat)))                       # 按基因 z-score
ann <- data.frame(dex = dds$dex, cell = dds$cell, row.names = colnames(mat))
png(file.path(OUT, "airway_top40_heatmap.png"), width = 900, height = 900)
pheatmap(mat, annotation_col = ann, show_rownames = TRUE, fontsize_row = 7,
         main = "Top 40 DE genes (VST, z-scored)")
dev.off()

cat("\n✅ 完成。结果文件在:", OUT, "\n")
cat("   - airway_DESeq2_all_genes.csv / ..._significant.csv\n")
cat("   - airway_PCA.png / sample_distance_heatmap.png / MA_plot.png / volcano.png / top40_heatmap.png\n")
