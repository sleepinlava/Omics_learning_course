#!/usr/bin/env Rscript
###############################################################################
# 03_DE_edgeR_limma_pasilla.R
#
# 用同一个数据集对比三种主流差异表达方法: edgeR / limma-voom / DESeq2
# 数据集: pasilla —— 果蝇(Drosophila), 基因敲低(treated) vs 对照(untreated), 7 个样本
#   (经典 Pasilla 数据集, 注意样本类型 single/paired 不一致, 归一化需小心)
#
# 学到的技能:
#   - 从原始 count 表 + 注释表构建 DGEList
#   - filterByExpr 自动过滤低表达
#   - TMM 归一化, MDS 图
#   - edgeR: estimateDisp -> glmQLFit -> QLF test
#   - limma-voom: voom -> lmFit -> eBayes -> topTable
#   - 三种方法结果的一致性比较 (Venn / 相关性)
#
# 用法:  Rscript scripts/03_DE_edgeR_limma_pasilla.R
###############################################################################

suppressPackageStartupMessages({
  library(edgeR); library(limma); library(pasilla)
  library(ggplot2); library(pheatmap); library(RColorBrewer)
})

.args <- commandArgs(trailingOnly = FALSE)
.f <- sub("^--file=", "", .args[grep("^--file=", .args)])
ROOT <- if (length(.f)) dirname(dirname(normalizePath(.f))) else getwd()
OUT  <- file.path(ROOT, "results"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(1234)

###############################################################################
# 1. 读入数据 ---------------------------------------------------------------
###############################################################################
counts <- read.delim(system.file("extdata", "pasilla_gene_counts.tsv", package = "pasilla"),
                     row.names = 1, check.names = FALSE)
anno   <- read.csv(system.file("extdata", "pasilla_sample_annotation.csv", package = "pasilla"))
rownames(anno) <- sub("fb$", "", anno$file)
anno <- anno[colnames(counts), ]              # 与 count 列对齐
stopifnot(identical(rownames(anno), colnames(counts)))
anno$condition <- factor(anno$condition, levels = c("untreated", "treated"))
cat("counts:", dim(counts), " | 样本:\n"); print(anno[, c("condition", "type")])

###############################################################################
# 2. edgeR 分析 --------------------------------------------------------------
###############################################################################
dge <- DGEList(counts = counts, samples = anno, group = anno$condition)

# 过滤: 自动保留"能被有效检验"的基因
keep <- filterByExpr(dge, group = dge$samples$group)
cat("\nfilterByExpr:", nrow(dge), "->", sum(keep), "基因\n")
dge  <- dge[keep, , keep.lib.sizes = FALSE]

# TMM 归一化
dge <- calcNormFactors(dge, method = "TMM")
cat("归一化因子:\n"); print(dge$samples$norm.factors)

# MDS 图 (类似 PCA)
png(file.path(OUT, "pasilla_MDS.png"), width = 800, height = 650)
plotMDS(dge, col = as.numeric(dge$samples$group),
        labels = rownames(dge$samples), main = "pasilla MDS (edgeR, TMM)")
legend("topright", legend = levels(dge$samples$group),
       col = 1:2, pch = 16)
dev.off()

# 离散度估计 + 拟似然 F 检验 (QLF, 推荐用于小样本)
dge <- estimateDisp(dge, design = model.matrix(~ group, data = dge$samples))
fit <- glmQLFit(dge, design = model.matrix(~ group, data = dge$samples))
qlf <- glmQLFTest(fit, coef = "grouptreated")
et  <- topTags(qlf, n = Inf)$table
cat("\nedgeR 显著基因 (FDR<0.05):", sum(et$FDR < 0.05), "\n")
write.csv(et, file.path(OUT, "pasilla_edgeR_results.csv"))

###############################################################################
# 3. limma-voom 分析 (同一数据, 另一种方法) ------------------------------------
###############################################################################
png(file.path(OUT, "pasilla_voom_meanvar.png"), width = 800, height = 650)
v <- voom(dge, design = model.matrix(~ group, data = dge$samples), plot = TRUE)
dev.off()

fit2 <- lmFit(v)
fit2 <- eBayes(fit2)
lim <- topTable(fit2, coef = "grouptreated", n = Inf, sort.by = "P")
cat("limma-voom 显著基因 (adj.P.Val<0.05):", sum(lim$adj.P.Val < 0.05), "\n")
write.csv(lim, file.path(OUT, "pasilla_limma_voom_results.csv"))

###############################################################################
# 4. 与 DESeq2 比较 -----------------------------------------------------------
###############################################################################
if (requireNamespace("DESeq2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(DESeq2))
  dds <- DESeqDataSetFromMatrix(counts, colData = anno, design = ~ condition)
  dds <- dds[rowSums(counts(dds) >= 5) >= 3, ]
  dds <- DESeq(dds, quiet = TRUE)
  ds  <- as.data.frame(results(dds, contrast = c("condition", "treated", "untreated")))
  write.csv(ds, file.path(OUT, "pasilla_DESeq2_results.csv"))
  cat("DESeq2 显著基因 (padj<0.05):", sum(ds$padj < 0.05, na.rm = TRUE), "\n")

  # 方法一致性: 显著基因重叠 + log2FC 相关性
  sig_edge <- rownames(et)[et$FDR < 0.05]
  sig_lim  <- rownames(lim)[lim$adj.P.Val < 0.05]
  sig_ds   <- rownames(ds)[!is.na(ds$padj) & ds$padj < 0.05]
  common <- Reduce(intersect, list(sig_edge, sig_lim, sig_ds))
  cat("\n三者共同显著基因:", length(common), "\n")
  cat("edgeR:", length(sig_edge), "| limma:", length(sig_lim), "| DESeq2:", length(sig_ds), "\n")

  # log2FC 散点图 (显式对齐基因名, 避免行名丢失)
  cmp <- data.frame(
    edgeR  = et[, "logFC"],
    limma  = lim[rownames(et), "logFC"],
    DESeq2 = ds[rownames(et), "log2FoldChange"]
  )
  rownames(cmp) <- rownames(et)
  cmp <- cmp[complete.cases(cmp), , drop = FALSE]
  write.csv(round(cor(cmp), 3), file.path(OUT, "pasilla_method_logFC_correlation.csv"))
  p <- ggplot(cmp, aes(edgeR, DESeq2)) + geom_point(alpha = .3, size = .8) +
    geom_abline(slope = 1, colour = "red", linetype = 2) +
    labs(title = "log2FC: edgeR vs DESeq2", x = "edgeR logFC", y = "DESeq2 log2FC") +
    theme_bw(base_size = 13)
  ggsave(file.path(OUT, "pasilla_logFC_edgeR_vs_DESeq2.png"), p, width = 6, height = 5, dpi = 150)
}

###############################################################################
# 5. 结果热图 ---------------------------------------------------------------
###############################################################################
top <- head(rownames(lim)[lim$adj.P.Val < 0.05], 40)
if (length(top) > 5) {
  logcpm <- cpm(dge, log = TRUE)
  mat <- t(scale(t(logcpm[top, ])))
  ann <- data.frame(condition = anno$condition, row.names = colnames(mat))
  png(file.path(OUT, "pasilla_top40_heatmap.png"), width = 800, height = 900)
  pheatmap(mat, annotation_col = ann, fontsize_row = 7, main = "pasilla top 40 DE genes")
  dev.off()
}

cat("\n✅ 完成, 结果在:", OUT, "\n")
