#!/usr/bin/env Rscript
# For reviewed, normalized LOG2 expression (e.g. microarrays), not RNA-seq counts.
f <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])
source(file.path(dirname(normalizePath(f)), "_common.R"))
args <- commandArgs(TRUE)
if (length(args) != 6 || !args[5] %in% c("paired", "independent") || args[6] != "log2_verified") {
  stop("Usage: 11_geo_limma.R expression.csv[.gz] samples.csv reference treatment paired|independent log2_verified")
}
need("limma")
x <- as.matrix(read.csv(args[1], row.names = 1, check.names = FALSE))
m <- read.csv(args[2], stringsAsFactors = FALSE)
required <- c("sample_id", "condition", "donor_id")
if (!all(required %in% names(m)) || anyNA(m[, required]) ||
    anyDuplicated(m$sample_id) || anyDuplicated(colnames(x)) || anyDuplicated(rownames(x))) stop("Invalid/duplicate sample or feature IDs.")
if (!setequal(colnames(x), m$sample_id)) stop("Expression columns and samples.csv must map one-to-one, with no extras.")
if (!is.numeric(x) || anyNA(x) || any(!is.finite(x))) stop("Require finite reviewed log2 expression; handle missing values explicitly first.")
m <- m[match(colnames(x), m$sample_id), , drop = FALSE]
if (args[3] == args[4] || !setequal(m$condition, args[3:4])) stop("Exactly the specified two conditions are required.")
m$condition <- factor(m$condition, levels = args[3:4])
m$donor_id <- factor(m$donor_id)
if (args[5] == "paired") {
  pairs <- table(m$donor_id, m$condition)
  if (any(pairs != 1) || nrow(pairs) < 3) stop("This teaching wrapper requires at least 3 complete donor pairs.")
  formula <- ~ donor_id + condition
} else {
  if (anyDuplicated(m$donor_id) || any(table(m$condition) < 3)) stop("Independent design requires unique donors and at least 3 per group.")
  formula <- ~ condition
}
# Optional batch is modeled, never removed from the expression matrix for testing.
if ("batch" %in% names(m)) {
  if (anyNA(m$batch)) stop("Missing batch labels.")
  if (length(unique(m$batch)) > 1) {m$batch <- factor(m$batch); formula <- update(formula, ~ . + batch)}
}
design <- model.matrix(formula, m)
if (qr(design)$rank != ncol(design) || nrow(design) <= ncol(design)) stop("Nonidentifiable or saturated design.")
condition_col <- which(attr(design, "assign") == match("condition", attr(terms(formula), "term.labels")))
stopifnot(length(condition_col) == 1)
keep <- apply(x, 1, sd) > 0
if (sum(keep) < 2) stop("Need at least two variable features for this analysis and PCA.")
out <- output_dir("geo_limma")
write.csv(data.frame(feature_id = rownames(x), retained = keep), file.path(out, "filtering.csv"), row.names = FALSE)
fit <- limma::lmFit(x[keep, , drop = FALSE], design)
eb <- limma::eBayes(fit, trend = TRUE)
tab <- limma::topTable(eb, coef = condition_col, number = Inf, adjust.method = "BH", confint = TRUE)
idx <- match(rownames(tab), rownames(fit$coefficients))
tab$feature_id <- rownames(tab)
tab$SE_classical <- fit$stdev.unscaled[idx, condition_col] * fit$sigma[idx]
tab$SE_moderated <- fit$stdev.unscaled[idx, condition_col] * sqrt(eb$s2.post[idx])
write.csv(tab, file.path(out, "all_features.csv"), row.names = FALSE)
write.csv(design, file.path(out, "design_matrix.csv"), row.names = FALSE)
write.csv(m, file.path(out, "samples.csv"), row.names = FALSE)
png(file.path(out, "PCA.png"), width = 900, height = 650)
pca <- prcomp(t(x[keep, , drop = FALSE]))
plot(pca$x[, 1:2], col = as.integer(m$condition), pch = 16, main = "Reviewed log2 matrix: sample QC")
legend("topright", levels(m$condition), col = 1:2, pch = 16)
dev.off()
saveRDS(list(fit = fit, moderated = eb, design = design), file.path(out, "limma_fit.rds"))
record_inputs(args[1:2], out)
record_session(out, c(paste("Design:", paste(deparse(formula), collapse = " ")),
  paste("Direction:", args[4], "minus", args[3]), "Input explicitly declared log2-normalized; this cannot be proven from values alone.",
  "SE_classical is unmoderated sampling SE; p/CI from limma use its moderated inference.",
  "Feature IDs remain probes if the input uses probes. No automatic gene collapse."))
cat("GEO log2-matrix analysis completed:", out, "\n")
