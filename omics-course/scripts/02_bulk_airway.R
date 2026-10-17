#!/usr/bin/env Rscript
# DESeq2 package workflow; see resources/SOURCES.md and the official DESeq2 vignette.
# Reads the existing local CSV export; does not require the airway package.
f <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])
source(file.path(dirname(normalizePath(f)), "_common.R"))
need(c("DESeq2", "SummarizedExperiment"))
suppressPackageStartupMessages(library(DESeq2))
set.seed(1234)
out <- output_dir("bulk")
paths <- file.path(REPO, "data", "counts", c("airway_counts.csv", "airway_coldata.csv"))
if (!all(file.exists(paths))) stop("Local airway CSV inputs are missing.")
k <- as.matrix(read.csv(paths[1], row.names = 1, check.names = FALSE))
meta <- read.csv(paths[2], row.names = 1, check.names = FALSE)
if (!is.numeric(k) || anyNA(k) || any(!is.finite(k)) || any(k < 0) || any(k != round(k))) {
  stop("Expected finite, nonnegative integer raw counts.")
}
if (anyDuplicated(colnames(k)) || anyDuplicated(rownames(k)) || anyDuplicated(rownames(meta))) {
  stop("Duplicated gene/sample identifiers.")
}
if (!setequal(colnames(k), rownames(meta))) {
  stop("Sample mapping mismatch. Missing metadata: ", paste(setdiff(colnames(k), rownames(meta)), collapse = ","),
       "; extra metadata: ", paste(setdiff(rownames(meta), colnames(k)), collapse = ","))
}
meta <- meta[colnames(k), , drop = FALSE]
if (!all(c("cell", "dex") %in% names(meta)) || anyNA(meta[, c("cell", "dex")])) stop("Missing design metadata.")
if (!setequal(unique(meta$dex), c("untrt", "trt"))) stop("Unexpected treatment labels.")
meta$dex <- factor(meta$dex, levels = c("untrt", "trt"))
meta$cell <- factor(meta$cell)
if (any(table(meta$cell, meta$dex) != 1)) stop("Expected one trt/untrt library per paired cell line.")
design <- model.matrix(~ cell + dex, meta)
if (qr(design)$rank != ncol(design)) stop("Design is not identifiable.")
keep <- rowSums(k >= 10) >= 4
write.csv(data.frame(gene = rownames(k), passes_filter = keep), file.path(out, "filtering.csv"), row.names = FALSE)
dds <- DESeqDataSetFromMatrix(k[keep, , drop = FALSE], meta, ~ cell + dex)
dds <- DESeq(dds, quiet = TRUE)
res <- as.data.frame(results(dds, contrast = c("dex", "trt", "untrt"), alpha = 0.05))
res$wald_ci_low <- res$log2FoldChange - qnorm(0.975) * res$lfcSE
res$wald_ci_high <- res$log2FoldChange + qnorm(0.975) * res$lfcSE
res$gene <- rownames(res)
res$status <- ifelse(is.na(res$pvalue), "no_valid_pvalue",
                     ifelse(is.na(res$padj), "no_adjusted_pvalue", "tested"))
write.csv(res, file.path(out, "tested_genes.csv"), row.names = FALSE)
full <- res[match(rownames(k), res$gene), , drop = FALSE]
full$gene <- rownames(k)
full$status[!keep] <- "prefiltered_not_tested"
write.csv(full, file.path(out, "all_genes.csv"), row.names = FALSE)
sig <- res[!is.na(res$padj) & res$padj < 0.05, , drop = FALSE]
sig <- sig[order(sig$padj), , drop = FALSE]
write.csv(sig, file.path(out, "significant_genes.csv"), row.names = FALSE)
vsd <- vst(dds, blind = TRUE)
z <- SummarizedExperiment::assay(vsd)
pc <- prcomp(t(z))
pct <- 100 * pc$sdev^2 / sum(pc$sdev^2)
png(file.path(out, "PCA.png"), width = 900, height = 650)
plot(pc$x[, 1:2], col = as.integer(meta$dex), pch = as.integer(meta$cell) + 14,
     xlab = sprintf("PC1 %.1f%%", pct[1]), ylab = sprintf("PC2 %.1f%%", pct[2]), main = "Airway VST QC")
text(pc$x[, 1:2], labels = meta$cell, pos = 3, cex = 0.7)
legend("topright", legend = levels(meta$dex), col = 1:2, pch = 16)
dev.off()
png(file.path(out, "MA.png"), width = 900, height = 650)
plotMA(results(dds, contrast = c("dex", "trt", "untrt"), alpha = 0.05), main = "trt / untrt, unshrunken LFC")
dev.off()
png(file.path(out, "sample_distances.png"), width = 900, height = 750)
heatmap(as.matrix(dist(t(z))), symm = TRUE, scale = "none", margins = c(10, 10))
dev.off()
qc <- data.frame(sample_id = colnames(k), cell = meta$cell, dex = meta$dex,
                 library_counts = colSums(k), detected_genes = colSums(k > 0),
                 size_factor = sizeFactors(dds))
write.csv(qc, file.path(out, "sample_qc.csv"), row.names = FALSE)
write.csv(counts(dds, normalized = TRUE), file.path(out, "normalized_counts.csv"))
write.csv(design, file.path(out, "design_matrix.csv"))
saveRDS(dds, file.path(out, "dds.rds"))
record_inputs(paths, out)
record_session(out, c("design=~cell+dex; contrast=trt/untrt; seed=1234; alpha=0.05",
                      "LFC and Wald intervals are unshrunken. VST is used only for QC.",
                      paste("input genes", nrow(k), "retained", sum(keep), "FDR<0.05", nrow(sig))))
cat("Bulk finished:", nrow(k), "genes ->", sum(keep), "tested input genes;", nrow(sig), "FDR<0.05\n")
