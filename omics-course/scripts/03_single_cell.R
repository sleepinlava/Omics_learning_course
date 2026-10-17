#!/usr/bin/env Rscript
# Seurat PBMC teaching workflow; sources: resources/SOURCES.md.
# No donor-level testing, doublet removal, or automatic cell-type claims.
f <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])
source(file.path(dirname(normalizePath(f)), "_common.R"))
need(c("Seurat", "ggplot2"))
suppressPackageStartupMessages(library(Seurat))
set.seed(1234)
args <- commandArgs(TRUE)
allowed <- c("--demo", "--qc-only", "--min-features", "--max-features", "--max-mt")
if (any(grepl("^--", args) & !args %in% allowed)) stop("Unknown option.")
demo <- "--demo" %in% args
out <- output_dir(if (demo) "single_cell_demo" else "single_cell_pbmc3k")
if (demo) {
  data("pbmc_small", package = "SeuratObject")
  raw <- GetAssayData(pbmc_small, assay = "RNA", layer = "counts")
  obj <- CreateSeuratObject(raw, project = "PBMC_SMALL_DEMO_ONLY")
  d <- file.path(COURSE, "data", "pbmc_small")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  write.csv(as.matrix(raw), file.path(d, "counts.csv"))
  write.csv(obj[[]], file.path(d, "metadata.csv"))
  writeLines("SeuratObject::pbmc_small; reduced real teaching object, NOT full PBMC3k or synthetic data.",
             file.path(d, "SOURCE.txt"))
} else {
  path <- file.path(COURSE, "data", "pbmc3k", "filtered_gene_bc_matrices", "hg19")
  if (!dir.exists(path)) stop("Download PBMC3k first: python omics-course/scripts/04_fetch_data.py pbmc3k")
  raw <- Read10X(path)
  if (is.list(raw)) stop("Expected the PBMC3k gene-expression matrix, not a multimodal list.")
  obj <- CreateSeuratObject(raw, project = "PBMC3k", min.cells = 0, min.features = 0)
  if (!any(grepl("^MT-", rownames(obj)))) stop("Expected human gene symbols with MT- genes.")
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  record_inputs(list.files(path, full.names = TRUE), out)
}
qc <- obj[[]]
qc$barcode <- rownames(qc)
png(file.path(out, "qc_before.png"), width = 1200, height = 500)
par(mfrow = c(1, 3))
hist(qc$nFeature_RNA, breaks = 40, main = "Detected genes", xlab = "genes / cell")
plot(qc$nCount_RNA, qc$nFeature_RNA, pch = 16, cex = 0.3, xlab = "UMIs", ylab = "Detected genes")
if (!demo) hist(qc$percent.mt, breaks = 40, main = "Mitochondrial fraction", xlab = "percent") else {
  plot.new(); text(0.5, 0.5, "Reduced feature demo\nMT fraction not interpretable")
}
dev.off()
numeric_option <- function(name) {
  i <- match(name, args)
  if (is.na(i) || i == length(args)) stop("Provide ", name, " after reviewing QC, or use --qc-only.")
  val <- suppressWarnings(as.numeric(args[i + 1]))
  if (!is.finite(val)) stop("Invalid value for ", name)
  val
}
if ("--qc-only" %in% args || (!demo && !any(c("--min-features", "--max-features", "--max-mt") %in% args))) {
  qc$passes_qc <- NA
  write.csv(qc, file.path(out, "qc_cells.csv"), row.names = FALSE)
  record_session(out, "QC inspection only; no cells filtered and no clustering performed.")
  cat("QC written. Inspect plots before choosing thresholds.\n")
  quit(status = 0)
}
if (demo) {
  qc$passes_qc <- qc$nCount_RNA > 0
} else {
  lo <- numeric_option("--min-features")
  hi <- numeric_option("--max-features")
  mt <- numeric_option("--max-mt")
  if (lo < 0 || hi <= lo || mt <= 0 || mt > 100) stop("Invalid QC threshold range.")
  qc$passes_qc <- qc$nFeature_RNA > lo & qc$nFeature_RNA < hi & qc$percent.mt < mt
}
write.csv(qc, file.path(out, "qc_cells.csv"), row.names = FALSE)
if (sum(qc$passes_qc) < 20) stop("Too few cells after QC for this teaching workflow.")
obj <- subset(obj, cells = qc$barcode[qc$passes_qc])
# Gene presence filtering is performed after preserving all barcode QC records.
count_layer <- GetAssayData(obj, layer = "counts")
obj <- subset(obj, features = rownames(obj)[Matrix::rowSums(count_layer > 0) >= 3])
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = min(2000L, nrow(obj)), verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
npcs <- min(if (demo) 10L else 30L, ncol(obj) - 1L, length(VariableFeatures(obj)) - 1L)
if (npcs < 2) stop("Insufficient varying features.")
obj <- RunPCA(obj, npcs = npcs, seed.use = 1234, verbose = FALSE)
dims <- seq_len(min(10L, npcs))
obj <- FindNeighbors(obj, dims = dims, k.param = min(20L, ncol(obj) - 1L), verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, random.seed = 1234, verbose = FALSE)
ggplot2::ggsave(file.path(out, "PCA.png"), DimPlot(obj, reduction = "pca", label = TRUE),
                width = 7, height = 5, dpi = 150, bg = "white")
if (!demo) {
  obj <- RunUMAP(obj, dims = dims, seed.use = 1234, verbose = FALSE)
  ggplot2::ggsave(file.path(out, "UMAP.png"), DimPlot(obj, reduction = "umap", label = TRUE),
                  width = 7, height = 5, dpi = 150, bg = "white")
}
markers <- intersect(c("CD3D", "CD3E", "MS4A1", "CD79A", "LST1", "LYZ", "NKG7", "GNLY"), rownames(obj))
if (length(markers)) {
  ggplot2::ggsave(file.path(out, "marker_panel.png"), DotPlot(obj, features = markers) + RotatedAxis(),
                  width = 8, height = 5, dpi = 150, bg = "white")
}
write.csv(obj[[]], file.path(out, "cell_metadata.csv"))
saveRDS(obj, file.path(out, "seurat.rds"))
record_session(out, c(paste("Arguments:", paste(args, collapse = " ")),
                      paste("Cells retained", ncol(obj), "genes", nrow(obj), "PC dimensions", length(dims)),
                      "No cell types assigned automatically. No disease/group DE. No doublet/ambient correction.",
                      if (demo) "Reduced PBMC demo: QC metrics are not representative of full data." else
                        "Single-donor PBMC3k; cell clusters are not biological replicates."))
cat("Single-cell workflow completed:", out, "\n")
