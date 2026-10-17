#!/usr/bin/env Rscript
# Feature-level continuous-intensity statistics using limma, not a count model.
# See 05_metabolomics/README.md for the deliberately limited QC assumptions.
f <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])
source(file.path(dirname(normalizePath(f)), "_common.R"))
need("limma")
suppressPackageStartupMessages(library(limma))
out <- output_dir("metabolomics")
paths <- file.path(out, c("feature_intensity.csv", "samples.csv", "feature_annotations.csv"))
if (!all(file.exists(paths))) stop("Run 05_metabolomics_xcms.R first. No real metabolomics result is available yet.")
x <- as.matrix(read.csv(paths[1], row.names = 1, check.names = FALSE))
meta <- read.csv(paths[2], stringsAsFactors = FALSE)
ann <- read.csv(paths[3], row.names = 1, check.names = FALSE)
if (!is.numeric(x) || any(is.infinite(x)) || any(x < 0, na.rm = TRUE)) stop("Expected nonnegative peak areas or NA.")
if (anyDuplicated(meta$sample_id) || anyDuplicated(colnames(x)) || anyDuplicated(rownames(x)) ||
    !setequal(colnames(x), meta$sample_id)) stop("Intensity/sample mapping is not one-to-one.")
meta <- meta[match(colnames(x), meta$sample_id), , drop = FALSE]
if (anyNA(meta$condition) || !setequal(meta$condition, c("KO", "WT"))) stop("Expected KO and WT labels.")
meta$condition <- factor(meta$condition, levels = c("WT", "KO"))
if (any(table(meta$condition) < 3)) stop("This exercise requires at least 3 independent samples per group.")
if (!setequal(rownames(x), rownames(ann))) stop("Feature annotation mapping mismatch.")
ann <- ann[rownames(x), , drop = FALSE]
missing <- is.na(x) | x <= 0
keep <- rowSums(missing) == 0
positive_sd <- rep(FALSE, nrow(x))
positive_sd[keep] <- apply(log2(x[keep, , drop = FALSE]), 1, sd) > 0
keep <- keep & positive_sd
filter <- data.frame(feature_id = rownames(x), missing_or_nonpositive = rowSums(missing),
                     included = keep, reason = ifelse(keep, "included",
                       ifelse(rowSums(missing) > 0, "incomplete_positive_values", "zero_variance")))
write.csv(filter, file.path(out, "feature_filtering.csv"), row.names = FALSE)
if (sum(keep) < 5) stop("Too few complete positive variable features. Review missingness; do not silently impute.")
logx <- log2(x[keep, , drop = FALSE])
medians <- apply(logx, 2, median)
normalized <- sweep(logx, 2, medians, "-") + median(medians)
write.csv(data.frame(sample_id = colnames(x), condition = meta$condition,
                     missing_fraction = colMeans(missing), log2_median = medians),
          file.path(out, "intensity_qc.csv"), row.names = FALSE)
design <- model.matrix(~ condition, meta)
stopifnot(qr(design)$rank == ncol(design))
fit_table <- function(z) {
  fit <- eBayes(lmFit(z, design), trend = TRUE)
  tab <- topTable(fit, coef = "conditionKO", number = Inf, adjust.method = "BH", confint = TRUE)
  tab$feature_id <- rownames(tab)
  cbind(tab, ann[tab$feature_id, , drop = FALSE])
}
main <- fit_table(normalized)
write.csv(main, file.path(out, "limma_KO_vs_WT.csv"), row.names = FALSE)
write.csv(fit_table(logx), file.path(out, "sensitivity_no_median_normalization.csv"), row.names = FALSE)
write.csv(normalized, file.path(out, "log2_median_normalized.csv"))
# Scale only features with nonzero post-normalization variance for PCA.
variable <- apply(normalized, 1, sd) > 0
if (sum(variable) < 2) stop("Insufficient post-normalization variation for PCA.")
pc <- prcomp(t(normalized[variable, , drop = FALSE]), center = TRUE, scale. = TRUE)
pct <- 100 * pc$sdev^2 / sum(pc$sdev^2)
png(file.path(out, "PCA.png"), width = 900, height = 650)
plot(pc$x[, 1:2], col = as.integer(meta$condition), pch = 16,
     xlab = sprintf("PC1 %.1f%%", pct[1]), ylab = sprintf("PC2 %.1f%%", pct[2]),
     main = "faahKO: log2, median normalization, feature scaling")
text(pc$x[, 1:2], labels = meta$sample_id, pos = 3, cex = 0.7)
legend("topright", levels(meta$condition), col = 1:2, pch = 16)
dev.off()
record_session(out, c("Feature-level limma: KO - WT; continuous log2 intensities; trend=TRUE.",
                      "Complete-positive-feature subset; median normalization assumes balanced overall changes.",
                      "No metabolite identities, pooled-QC drift correction, or causal effects claimed.",
                      paste("Total features", nrow(x), "analyzed", sum(keep))),
               filename = "statistics_sessionInfo.txt")
cat("Metabolomics statistics completed:", sum(keep), "features\n")
