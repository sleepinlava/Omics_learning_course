#!/usr/bin/env Rscript
# Human ENSEMBL DESeq2 results -> reviewed ID mapping -> cached KEGG ORA or GSEA.
f <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])
source(file.path(dirname(normalizePath(f)), "_common.R"))
args <- commandArgs(TRUE)
if (length(args) != 3 || !args[3] %in% c("ora", "gsea")) stop("Usage: 13_kegg_enrichment.R all_DESeq2_results.csv hsa_snapshot.rds ora|gsea")
need(c("clusterProfiler", "AnnotationDbi", "org.Hs.eg.db"))
de <- read.csv(args[1], stringsAsFactors = FALSE)
cache <- readRDS(args[2])
if (!identical(cache$organism, "hsa") || isTRUE(cache$synthetic)) stop("This wrapper requires a real human KEGG snapshot; synthetic data are not KEGG.")
if (!all(c("term", "gene") %in% names(cache$term2gene)) ||
    !all(c("term", "name") %in% names(cache$term2name)) ||
    anyNA(cache$term2gene) || !length(cache$retrieved_at)) stop("Incomplete KEGG snapshot schema.")
if (!all(c("gene", "log2FoldChange", "pvalue", "padj", "stat") %in% names(de))) stop("Expected complete DESeq2 results with gene, log2FoldChange, pvalue, padj, stat.")
if (anyNA(de$gene) || anyDuplicated(de$gene)) stop("Missing/duplicate input gene IDs; resolve before enrichment.")
de$ENSEMBL <- sub("\\.[0-9]+$", "", de$gene)
if (anyDuplicated(de$ENSEMBL)) stop("Stripping Ensembl versions created duplicates; review mapping.")
de <- de[is.finite(de$pvalue) & is.finite(de$stat) & is.finite(de$log2FoldChange), , drop = FALSE]
if (!nrow(de) || any(!grepl("^ENSG[0-9]+$", de$ENSEMBL))) stop("Require tested human Ensembl genes, not symbols or probes.")
map <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = de$ENSEMBL,
                              keytype = "ENSEMBL", columns = "ENTREZID")
map <- unique(map[!is.na(map$ENTREZID), , drop = FALSE])
# A conservative one-to-one subset avoids inflating significance by duplicated IDs.
ok <- !duplicated(map$ENSEMBL) & !duplicated(map$ENSEMBL, fromLast = TRUE) &
      !duplicated(map$ENTREZID) & !duplicated(map$ENTREZID, fromLast = TRUE)
map$one_to_one <- ok
out <- output_dir(paste0("kegg_", args[3]))
write.csv(map, file.path(out, "id_mapping.csv"), row.names = FALSE)
write.csv(data.frame(ENSEMBL = setdiff(de$ENSEMBL, map$ENSEMBL)), file.path(out, "unmapped.csv"), row.names = FALSE)
one <- map[ok, , drop = FALSE]
de <- de[match(one$ENSEMBL, de$ENSEMBL), , drop = FALSE]
de$ENTREZID <- one$ENTREZID
in_pathway <- de$ENTREZID %in% cache$term2gene$gene
write.csv(data.frame(gene = de$gene, ENTREZID = de$ENTREZID, in_snapshot = in_pathway),
          file.path(out, "universe_mapping.csv"), row.names = FALSE)
de <- de[in_pathway, , drop = FALSE]
if (nrow(de) < 10) stop("Too few uniquely mapped pathway-annotated genes.")
universe <- unique(de$ENTREZID)
sizes <- table(unique(cache$term2gene[cache$term2gene$gene %in% universe, c("term", "gene")])$term)
eligible <- names(sizes)[sizes >= 10 & sizes <= 500]
write.csv(data.frame(term = names(sizes), measured_genes = as.integer(sizes),
  eligible = names(sizes) %in% eligible), file.path(out, "pathway_test_family.csv"), row.names = FALSE)
set.seed(1234)
if (args[3] == "ora") {
  sig <- !is.na(de$padj) & de$padj < 0.05
  combined <- list()
  for (direction in c("up", "down")) {
    choose <- sig & if (direction == "up") de$log2FoldChange > 0 else de$log2FoldChange < 0
    ids <- unique(de$ENTREZID[choose])
    write.csv(data.frame(ENTREZID = ids), file.path(out, paste0(direction, "_input.csv")), row.names = FALSE)
    if (!length(ids)) next
    result <- clusterProfiler::enricher(ids, universe = universe, TERM2GENE = cache$term2gene,
      TERM2NAME = cache$term2name, minGSSize = 10, maxGSSize = 500,
      pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
    if (is.null(result)) next
    tab <- as.data.frame(result)
    if (nrow(tab)) {tab$direction <- direction; combined[[direction]] <- tab}
  }
  tab <- if (length(combined)) do.call(rbind, combined) else data.frame()
  # Reserve both directions for every eligible pathway, including zero-overlap p=1 tests.
  if (nrow(tab)) tab$BH_across_directions <- p.adjust(tab$pvalue, "BH", n = 2 * length(eligible))
} else {
  ranking <- sort(setNames(de$stat, de$ENTREZID), decreasing = TRUE)
  write.csv(data.frame(ENTREZID = names(ranking), Wald_statistic = unname(ranking)), file.path(out, "ranking.csv"), row.names = FALSE)
  result <- clusterProfiler::GSEA(ranking, TERM2GENE = cache$term2gene, TERM2NAME = cache$term2name,
    minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", verbose = FALSE, seed = TRUE)
  tab <- as.data.frame(result)
}
write.csv(tab, file.path(out, "all_pathways.csv"), row.names = FALSE)
if (!nrow(tab)) writeLines("No eligible pathways returned. This is not evidence that pathways are unchanged.", file.path(out, "EMPTY_RESULT.txt"))
record_inputs(args[1:2], out)
record_session(out, c(paste("Mode:", args[3], "snapshot:", cache$retrieved_at),
  paste("Effective universe:", length(universe), "uniquely mapped, tested, KEGG-annotated genes"),
  "ORA: BH_across_directions reserves both directions across ALL size-eligible pathways, including zero overlaps.",
  "GSEA: full available signed Wald ranking, not just significant genes; no pathway activation claims."))
cat("KEGG analysis completed:", out, "\n")
