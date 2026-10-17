#!/usr/bin/env Rscript
# Calls xcms's established OnDisk workflow on the faahKO teaching data.
# Reference: https://github.com/sneumann/xcms/blob/master/vignettes/xcms.Rmd
# Not validated end-to-end here; see VALIDATION.md. Parameters are dataset-specific.
f <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])
source(file.path(dirname(normalizePath(f)), "_common.R"))
need(c("xcms", "faahKO", "MSnbase", "Biobase", "BiocParallel"))
suppressPackageStartupMessages({ library(xcms); library(MSnbase); library(Biobase) })
BiocParallel::register(BiocParallel::SerialParam())
set.seed(1234)
out <- output_dir("metabolomics")
files <- sort(list.files(system.file("cdf", package = "faahKO"), pattern = "\\.cdf$",
                         recursive = TRUE, full.names = TRUE, ignore.case = TRUE))
if (length(files) != 12) stop("Expected 12 faahKO CDF files; inspect the data package before proceeding.")
group <- toupper(basename(dirname(files)))
if (!setequal(group, c("KO", "WT")) || any(table(group) != 6)) stop("Unexpected faahKO KO/WT layout.")
pd <- data.frame(sample_name = tools::file_path_sans_ext(basename(files)),
                 sample_group = factor(group, levels = c("WT", "KO")))
if (anyDuplicated(pd$sample_name)) stop("Duplicate sample names.")
rownames(pd) <- pd$sample_name
write.csv(data.frame(sample_id = pd$sample_name, condition = pd$sample_group),
          file.path(out, "samples.csv"), row.names = FALSE)
record_inputs(files, out)
raw <- MSnbase::readMSData(files = files, pdata = new("NAnnotatedDataFrame", pd), mode = "onDisk")
png(file.path(out, "BPC.png"), width = 1100, height = 650)
plot(chromatogram(raw, aggregationFun = "max"), col = ifelse(group == "KO", "coral", "steelblue"))
dev.off()
# Parameters for the legacy faahKO exercise, not defaults for arbitrary instruments.
params <- list(peakwidth = c(20, 80), ppm = 25, noise = 5000,
               prefilter = c(6, 5000), grouping_bw_initial = 5,
               grouping_bw_final = 5, minFraction = 0.5)
dput(params, file.path(out, "xcms_parameters.R"))
xp <- findChromPeaks(raw, param = CentWaveParam(peakwidth = params$peakwidth,
                     ppm = params$ppm, noise = params$noise, prefilter = params$prefilter))
gp <- PeakDensityParam(sampleGroups = pd$sample_group, minFraction = params$minFraction, bw = 5)
xp <- groupChromPeaks(xp, param = gp)
xp <- adjustRtime(xp, param = PeakGroupsParam(minFraction = 0.5))
png(file.path(out, "retention_time_adjustment.png"), width = 1000, height = 600)
plotAdjustedRtime(xp, col = ifelse(group == "KO", "coral", "steelblue"))
dev.off()
xp <- groupChromPeaks(xp, param = gp)
before <- featureValues(xp, value = "into")
xp <- fillChromPeaks(xp)
after <- featureValues(xp, value = "into")
if (ncol(after) != nrow(pd)) stop("Unexpected intensity/sample dimensions.")
colnames(before) <- colnames(after) <- pd$sample_name
write.csv(before, file.path(out, "feature_intensity_before_fill.csv"))
write.csv(after, file.path(out, "feature_intensity.csv"))
defs <- as.data.frame(featureDefinitions(xp))
defs <- defs[, !vapply(defs, is.list, logical(1)), drop = FALSE]
write.csv(defs, file.path(out, "feature_annotations.csv"))
saveRDS(xp, file.path(out, "xcms_object.rds"))
record_session(out, c("faahKO: 6 KO / 6 WT, positive mode, truncated m/z and RT ranges.",
                      "OnDisk/MSnbase workflow; no chemical identities assigned.",
                      "Original CDF paths must remain accessible when reopening the saved object.",
                      "QC/blank/injection order are unavailable here; no drift correction claimed."),
               filename = "xcms_sessionInfo.txt")
cat("xcms features exported:", nrow(after), "\n")
