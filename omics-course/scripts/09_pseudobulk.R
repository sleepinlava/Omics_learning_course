#!/usr/bin/env Rscript
# Thin scuttle + DESeq2 wrapper for a reviewed, complete paired SCE only.
# Teaching reference: hbctraining/Pseudobulk-for-scRNAseq (see resources/SOURCES.md).
f <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])
source(file.path(dirname(normalizePath(f)), "_common.R"))
args <- commandArgs(TRUE)
if (length(args) != 5) stop("Usage: 09_pseudobulk.R reviewed.rds cell_type reference treatment min_cells")
need(c("SingleCellExperiment", "scuttle", "DESeq2"))
suppressPackageStartupMessages({ library(SingleCellExperiment); library(DESeq2) })
sce <- readRDS(args[1])
if (!is(sce, "SingleCellExperiment") || !"counts" %in% assayNames(sce)) stop("Expected an SCE with a raw counts assay.")
required <- c("sample_id", "donor_id", "condition", "cell_type", "passes_qc")
if (!all(required %in% names(colData(sce)))) stop("Missing required cell metadata: ", paste(required, collapse = ", "))
m <- as.data.frame(colData(sce))
if (anyNA(m[, required]) || !is.logical(m$passes_qc)) stop("Complete metadata and an explicitly reviewed logical passes_qc column are required.")
for (field in required[1:4]) {
  m[[field]] <- as.character(m[[field]])
  if (any(!nzchar(m[[field]]))) stop("Empty metadata values in ", field)
}
if (anyDuplicated(colnames(sce)) || anyDuplicated(rownames(sce))) stop("Require unique cell and feature identifiers.")
min_cells <- suppressWarnings(as.integer(args[5]))
if (is.na(min_cells) || min_cells < 1) stop("min_cells must be a positive integer.")
if (args[3] == args[4] || !setequal(unique(m$condition), args[3:4])) stop("This wrapper requires exactly the two specified conditions.")
sample_meta <- unique(m[, c("sample_id", "donor_id", "condition")])
if (anyDuplicated(sample_meta$sample_id)) stop("A sample ID maps to more than one donor/condition.")
if (anyDuplicated(sample_meta[, c("donor_id", "condition")])) stop("Multiple samples per donor/condition need an explicit design beyond this wrapper.")
take <- m$passes_qc & m$cell_type == args[2]
if (!any(take)) stop("Requested cell type has no QC-passing cells; inspect labels.")
cells_per_sample <- table(factor(m$sample_id[take], levels = sample_meta$sample_id))
out <- output_dir("pseudobulk")
sample_meta$cells <- as.integer(cells_per_sample)
sample_meta$passes_cell_number <- sample_meta$cells >= min_cells
write.csv(sample_meta, file.path(out, "sample_cell_numbers.csv"), row.names = FALSE)
if (!all(sample_meta$passes_cell_number)) stop("Some samples have too few cells; inspect sample_cell_numbers.csv and review the paired design. No DE performed.")
pairing <- table(sample_meta$donor_id, sample_meta$condition)
if (nrow(pairing) < 3 || any(pairing != 1)) stop("Require at least three complete donor pairs; this is not a power guarantee.")
sce <- sce[, take]
count_values <- assay(sce, "counts")
if (anyNA(count_values) || any(!is.finite(count_values)) || any(count_values < 0) ||
    any(count_values != round(count_values))) stop("Expected nonnegative integer raw counts, not log/integrated values.")
agg <- scuttle::aggregateAcrossCells(sce, ids = S4Vectors::DataFrame(sample_id = m$sample_id[take]),
                                    use.assay.type = "counts")
k <- as.matrix(assay(agg, "counts"))
sample_ids <- as.character(colData(agg)$sample_id)
sample_meta <- sample_meta[match(sample_ids, sample_meta$sample_id), , drop = FALSE]
if (anyNA(sample_meta$sample_id)) stop("Aggregation metadata mismatch.")
colnames(k) <- rownames(sample_meta) <- sample_ids
sample_meta$condition <- factor(sample_meta$condition, levels = args[3:4])
sample_meta$donor_id <- factor(sample_meta$donor_id)
design <- model.matrix(~ donor_id + condition, sample_meta)
if (qr(design)$rank != ncol(design)) stop("Paired design is not full rank.")
keep <- rowSums(k >= 10) >= nrow(pairing)
if (sum(keep) < 2) stop("Too few genes after prefiltering.")
write.csv(k, file.path(out, "pseudobulk_raw_counts.csv"))
write.csv(data.frame(gene = rownames(k), retained = keep), file.path(out, "gene_filtering.csv"), row.names = FALSE)
write.csv(sample_meta, file.path(out, "samples.csv"), row.names = FALSE)
write.csv(design, file.path(out, "design_matrix.csv"))
dds <- DESeqDataSetFromMatrix(k[keep, , drop = FALSE], sample_meta, ~ donor_id + condition)
dds <- DESeq(dds, quiet = TRUE)
res <- as.data.frame(results(dds, contrast = c("condition", args[4], args[3]), alpha = 0.05))
write.csv(res, file.path(out, "all_tested_genes.csv"))
write.csv(counts(dds, normalized = TRUE), file.path(out, "normalized_counts.csv"))
saveRDS(dds, file.path(out, "dds.rds"))
record_inputs(args[1], out)
record_session(out, c(paste("Arguments:", paste(args, collapse = " ")),
                      "Complete paired design only; raw counts summed by sample within selected cell type.",
                      "BH correction within the selected cell type; cross-cell-type multiplicity not handled."))
cat("Paired pseudobulk completed:", out, "\n")
