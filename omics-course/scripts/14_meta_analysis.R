#!/usr/bin/env Rscript
# Mature REML/Hartung-Knapp meta-analysis via metafor, with feature-level BH.
f <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])
source(file.path(dirname(normalizePath(f)), "_common.R"))
args <- commandArgs(TRUE)
if (length(args) != 1) stop("Usage: 14_meta_analysis.R effects.csv | --demo")
demo <- identical(args, "--demo")
path <- if (demo) file.path(COURSE, "examples", "meta_demo_effects.csv") else args[1]
d <- read.csv(path, stringsAsFactors = FALSE)
required <- c("feature_id", "study_id", "cohort_id", "effect", "se", "effect_scale", "contrast", "organism", "independence_reviewed")
if (!all(required %in% names(d)) || anyNA(d[, required])) stop("Missing/NA meta-analysis schema fields.")
if (!nrow(d)) stop("No effect estimates supplied.")
for (field in c("feature_id", "study_id", "cohort_id", "effect_scale", "contrast", "organism")) {
  if (any(!nzchar(trimws(d[[field]])))) stop("Empty identifiers in ", field)
}
if (!is.numeric(d$effect) || !is.numeric(d$se) || any(!is.finite(d$effect)) ||
    any(!is.finite(d$se)) || any(d$se <= 0)) stop("Require finite effects and strictly positive sampling standard errors.")
if (!all(d$independence_reviewed %in% TRUE)) stop("Review participant overlap and mark independence_reviewed TRUE before combining studies.")
if (anyDuplicated(d[, c("feature_id", "study_id")]) || anyDuplicated(d[, c("feature_id", "cohort_id")])) {
  stop("Repeated study/cohort within feature: overlapping estimates require a covariance-aware model.")
}
if (anyDuplicated(unique(d[, c("study_id", "cohort_id")])$study_id)) stop("A study_id maps to multiple cohorts; clarify study units.")
if (length(unique(d$effect_scale)) != 1 || !unique(d$effect_scale) %in% c("log2FC", "SMD", "logRR", "logOR") ||
    length(unique(d$contrast)) != 1 || length(unique(d$organism)) != 1) stop("Require a common scale, contrast direction, and organism across this run.")
need("metafor")
out <- output_dir(if (demo) "meta_demo" else "meta_analysis")
if (demo) writeLines("SYNTHETIC teaching effects, not real GEO studies or biological evidence.", file.path(out, "SYNTHETIC_DATA.txt"))
summary <- list(); failures <- list(); loo <- list(); first <- TRUE
for (feature in sort(unique(d$feature_id))) {
  z <- d[d$feature_id == feature, , drop = FALSE]
  if (nrow(z) < 3) {
    failures[[feature]] <- data.frame(feature_id = feature, reason = "fewer_than_3_independent_studies")
    next
  }
  tryCatch({
    random <- metafor::rma.uni(yi = z$effect, sei = z$se, method = "REML", test = "knha", slab = z$study_id)
    fixed <- metafor::rma.uni(yi = z$effect, sei = z$se, method = "FE")
    pred <- predict(random)
    feature_summary <- data.frame(feature_id = feature, k = random$k,
      effect_scale = unique(z$effect_scale), contrast = unique(z$contrast),
      pooled_effect = as.numeric(coef(random)), se = random$se, ci_low = random$ci.lb,
      ci_high = random$ci.ub, pvalue = random$pval, tau2 = random$tau2, I2 = random$I2,
      Q = random$QE, Q_pvalue = random$QEp, prediction_low = pred$pi.lb, prediction_high = pred$pi.ub,
      fixed_effect = as.numeric(coef(fixed)), fixed_pvalue = fixed$pval)
    diagnostics <- as.data.frame(metafor::leave1out(random))
    diagnostics$omitted_study <- z$study_id
    diagnostics$feature_id <- feature
    loo[[feature]] <- diagnostics
    summary[[feature]] <- feature_summary
    if (first) {
      png(file.path(out, "forest_first_feature.png"), width = 1100, height = 800)
      tryCatch(metafor::forest(random, xlab = paste(unique(z$effect_scale), unique(z$contrast)),
                       main = paste(if (demo) "SIMULATION" else "Meta-analysis", feature)),
               finally = dev.off())
      writeLines(feature, file.path(out, "plotted_feature.txt"))
      first <- FALSE
    }
  }, error = function(e) {
    failures[[feature]] <<- data.frame(feature_id = feature, reason = conditionMessage(e))
  })
}
if (length(failures)) write.csv(do.call(rbind, failures), file.path(out, "skipped_or_failed.csv"), row.names = FALSE)
if (!length(summary)) stop("No feature successfully fitted; inspect failures. No meta result claimed.")
tab <- do.call(rbind, summary)
# Reserve the full submitted feature family, including features that could not be fitted.
tab$padj <- p.adjust(tab$pvalue, method = "BH", n = length(unique(d$feature_id)))
write.csv(tab, file.path(out, "all_meta_results.csv"), row.names = FALSE)
if (length(loo)) write.csv(do.call(rbind, loo), file.path(out, "leave_one_out.csv"), row.names = FALSE)
record_inputs(path, out)
record_session(out, c("Random effects: REML; Hartung-Knapp; fixed effect is a sensitivity comparison.",
  "BH family includes all submitted features. >=3 studies is a teaching floor, not adequate-power assurance.",
  "No automatic publication-bias test or funnel-based correction; overlap must be reviewed outside this script.",
  if (demo) "SYNTHETIC INPUT" else "User-reviewed independent study effect estimates."))
cat("Meta-analysis completed:", out, "\n")
