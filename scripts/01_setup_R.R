#!/usr/bin/env Rscript
###############################################################################
# 01_setup_R.R  -- 检查 / 安装 RNA-seq 分析所需的 R / Bioconductor 包
#
# 用法:  Rscript scripts/01_setup_R.R
###############################################################################
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# 按学习阶段分组的包
groups <- list(
  "1_core-DE (差异表达)" = c(
    "DESeq2", "edgeR", "limma", "apeglm", "ashr", "IHW", "tximport", "vsn"
  ),
  "2_visualisation (可视化)" = c(
    "ggplot2", "pheatmap", "RColorBrewer", "ggrepel", "EnhancedVolcano",
    "patchwork", "PCAtools", "corrplot", "gplots", "reshape2"
  ),
  "3_data-IO (数据获取)" = c(
    "GEOquery", "biomaRt", "AnnotationDbi", "airway", "pasilla"
  ),
  "4_enrichment (富集分析)" = c(
    "clusterProfiler", "enrichplot", "DOSE", "org.Hs.eg.db", "org.Mm.eg.db"
  )
)

all_pkgs <- unique(unlist(groups, use.names = FALSE))
missing  <- all_pkgs[!vapply(all_pkgs, requireNamespace, logical(1), quietly = TRUE)]

cat("已安装:", length(all_pkgs) - length(missing), "/", length(all_pkgs), "\n")
if (length(missing)) {
  cat("需要安装:", paste(missing, collapse = ", "), "\n")
  BiocManager::install(missing, ask = FALSE, update = FALSE)
} else {
  cat("所有包均已就绪 ✅\n")
}

# 最终核对
ok <- vapply(all_pkgs, requireNamespace, logical(1), quietly = TRUE)
cat("\n=== 最终状态 ===\n")
for (g in names(groups)) {
  cat("\n[", g, "]\n", sep = "")
  for (p in groups[[g]]) cat(sprintf("  %-14s %s\n", p, if (ok[[p]]) "OK" else "缺失"))
}
if (all(ok)) cat("\n环境准备完成，可以运行 02_DE_DESeq2_airway.R\n")
