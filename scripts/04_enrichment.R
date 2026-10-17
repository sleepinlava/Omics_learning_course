#!/usr/bin/env Rscript
###############################################################################
# 04_enrichment.R  --  差异基因的 GO / KEGG 功能富集
#
# 输入: results/airway_DESeq2_results.csv (由 02 脚本产生)
# 输出: results/enrichment_GO_*.png / enrichment_KEGG_*.png + csv 表
#
# 学到的技能: ID 转换(ENSEMBL->ENTREZ/SYMBOL)、enrichGO / enrichKEGG、富集气泡图
#
# 用法:  Rscript scripts/04_enrichment.R
###############################################################################

if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  stop("请先运行 scripts/01_setup_R.R 安装 clusterProfiler")
}
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db)
  library(ggplot2); library(enrichplot)
})

.args <- commandArgs(trailingOnly = FALSE)
.f <- sub("^--file=", "", .args[grep("^--file=", .args)])
ROOT <- if (length(.f)) dirname(dirname(normalizePath(.f))) else getwd()
OUT  <- file.path(ROOT, "results")
resfile <- file.path(OUT, "airway_DESeq2_all_genes.csv")
if (!file.exists(resfile)) stop("找不到 ", resfile, " —— 请先运行 02 脚本")

res <- read.csv(resfile, row.names = 1, check.names = FALSE)
res <- res[!is.na(res$padj), ]

# 1. 定义基因列表
sig_up   <- rownames(res)[res$padj < 0.05 & res$log2FoldChange >  1]
sig_down <- rownames(res)[res$padj < 0.05 & res$log2FoldChange < -1]
cat("上调基因:", length(sig_up), " 下调基因:", length(sig_down), "\n")

# 2. ENSEMBL -> ENTREZ
up_ent   <- bitr(sig_up,   fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
down_ent <- bitr(sig_down, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
# 背景基因 (所有被检验的基因)
bg_ent   <- bitr(rownames(res), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

# 3. GO 富集 (上调基因, Biological Process)
ego_up <- enrichGO(up_ent, universe = bg_ent, OrgDb = org.Hs.eg.db,
                   ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2, readable = TRUE)
if (!is.null(ego_up) && nrow(as.data.frame(ego_up)) > 0) {
  write.csv(as.data.frame(ego_up), file.path(OUT, "enrichment_GO_BP_up.csv"))
  p <- dotplot(ego_up, showCategory = 15) + ggtitle("GO BP enriched in UP genes (dex)")
  ggsave(file.path(OUT, "enrichment_GO_BP_up.png"), p, width = 8, height = 6, dpi = 150)
}

# 4. GO 富集 (下调基因)
ego_dn <- enrichGO(down_ent, universe = bg_ent, OrgDb = org.Hs.eg.db,
                   ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2, readable = TRUE)
if (!is.null(ego_dn) && nrow(as.data.frame(ego_dn)) > 0) {
  write.csv(as.data.frame(ego_dn), file.path(OUT, "enrichment_GO_BP_down.csv"))
  p <- dotplot(ego_dn, showCategory = 15) + ggtitle("GO BP enriched in DOWN genes (dex)")
  ggsave(file.path(OUT, "enrichment_GO_BP_down.png"), p, width = 8, height = 6, dpi = 150)
}

# 5. KEGG 通路富集 (需要联网, 若失败不影响其它结果)
kegg <- tryCatch(
  enrichKEGG(c(up_ent, down_ent), universe = bg_ent,
             organism = "hsa", pvalueCutoff = 0.05),
  error = function(e) { message("KEGG 富集失败(通常是网络问题): ", conditionMessage(e)); NULL })
if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
  kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  write.csv(as.data.frame(kegg), file.path(OUT, "enrichment_KEGG.csv"))
  ggsave(file.path(OUT, "enrichment_KEGG.png"),
         dotplot(kegg, showCategory = 15), width = 8, height = 6, dpi = 150)
}

cat("\n✅ 富集分析完成, 结果在:", OUT, "\n")
