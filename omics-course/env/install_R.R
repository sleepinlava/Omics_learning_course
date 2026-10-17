#!/usr/bin/env Rscript
# Install only the selected module into a course-local library.
a <- commandArgs(trailingOnly = FALSE)
script <- normalizePath(sub("^--file=", "", a[grepl("^--file=", a)]))
course <- dirname(dirname(script))
lib <- file.path(course, "env", "R-library")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1 || !args %in% c("core", "metabolomics", "pseudobulk", "geo", "kegg", "meta")) {
  stop("Usage: Rscript omics-course/env/install_R.R core|metabolomics|pseudobulk|geo|kegg|meta")
}
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", lib = lib)
pkgs <- switch(args,
  core = c("DESeq2", "Seurat", "ggplot2", "limma"),
  metabolomics = c("xcms", "faahKO", "MSnbase", "Biobase", "BiocParallel", "limma"),
  pseudobulk = c("muscData", "scuttle", "SingleCellExperiment", "DESeq2"),
  geo = c("GEOquery", "limma", "Biobase", "SummarizedExperiment", "jsonlite"),
  kegg = c("KEGGREST", "clusterProfiler", "AnnotationDbi", "org.Hs.eg.db"),
  meta = c("metafor", "metadat"))
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  if (args == "meta") install.packages(missing, lib = lib)
  else BiocManager::install(missing, lib = lib, ask = FALSE, update = FALSE)
}
ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
print(ok)
writeLines(capture.output(sessionInfo()), file.path(course, "env", paste0(args, "-install-session.txt")))
if (!all(ok)) stop("Some dependencies are still unavailable; see installation messages.")
