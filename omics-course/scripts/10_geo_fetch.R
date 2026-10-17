#!/usr/bin/env Rscript
# GEOquery acquisition only; expression scale and sample groups require review.
f <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])
source(file.path(dirname(normalizePath(f)), "_common.R"))
args <- commandArgs(TRUE)
if (length(args) != 2 || !grepl("^GSE[0-9]+$", args[1]) ||
    !args[2] %in% c("matrix", "counts", "supp-list")) {
  stop("Usage: 10_geo_fetch.R GSE_ID matrix|counts|supp-list")
}
need(c("GEOquery", "Biobase", "SummarizedExperiment", "jsonlite"))
dest <- file.path(COURSE, "data", "geo", args[1], args[2])
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
options(timeout = 180)
stamp <- format(Sys.time(), tz = "UTC", usetz = TRUE)
tryCatch({
  if (args[2] == "supp-list") {
    files <- GEOquery::getGEOSuppFiles(args[1], fetch_files = FALSE)
    write.csv(files, file.path(dest, "supplementary_files.csv"), row.names = FALSE)
  } else if (args[2] == "matrix") {
    series <- GEOquery::getGEO(args[1], GSEMatrix = TRUE, getGPL = FALSE, destdir = dest)
    if (!length(series)) stop("No Series Matrix returned; inspect supplementary files or counts availability.")
    for (i in seq_along(series)) {
      es <- series[[i]]
      subdir <- file.path(dest, paste0("platform_", i, "_", Biobase::annotation(es)))
      dir.create(subdir, showWarnings = FALSE)
      x <- Biobase::exprs(es)
      if (!nrow(x) || !ncol(x)) stop("Empty expression matrix; no analyzable matrix claimed.")
      con <- gzfile(file.path(subdir, "expression.csv.gz"), "wt")
      tryCatch(write.csv(x, con), finally = close(con))
      write.csv(Biobase::pData(es), file.path(subdir, "metadata.csv"))
      write.csv(Biobase::fData(es), file.path(subdir, "feature_metadata.csv"))
      saveRDS(es, file.path(subdir, "ExpressionSet.rds"))
      writeLines(c("Scale: UNKNOWN; review processing protocol before analysis.",
                    "Do not infer control/case from column order or title heuristics.",
                    paste("Dimensions:", paste(dim(x), collapse = " x "))), file.path(subdir, "REVIEW_REQUIRED.txt"))
    }
  } else {
    api <- getNamespaceExports("GEOquery")
    if (!all(c("hasRNASeqQuantifications", "getRNASeqData") %in% api)) stop("This GEOquery version lacks RNA-seq quantification APIs.")
    if (!isTRUE(GEOquery::hasRNASeqQuantifications(args[1]))) stop("No GEO-calculated RNA-seq counts available for this accession.")
    se <- GEOquery::getRNASeqData(args[1])
    if (!is(se, "SummarizedExperiment")) stop("Unexpected response type; inspect before exporting counts.")
    saveRDS(se, file.path(dest, "RNASeqData.rds"))
    write.csv(as.data.frame(SummarizedExperiment::colData(se)), file.path(dest, "metadata.csv"))
    write.csv(as.data.frame(SummarizedExperiment::rowData(se)), file.path(dest, "feature_metadata.csv"))
    writeLines(c("Check assay names, annotation/genome build, counts representation and sample mapping before DE.",
                  paste("Assays:", paste(SummarizedExperiment::assayNames(se), collapse = ", "))),
               file.path(dest, "REVIEW_REQUIRED.txt"))
  }
  jsonlite::write_json(list(accession = args[1], mode = args[2], time_utc = stamp,
    source = paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", args[1]),
    status = "downloaded_not_scientifically_reviewed", GEOquery = as.character(packageVersion("GEOquery"))),
    file.path(dest, "receipt.json"), auto_unbox = TRUE, pretty = TRUE)
  record_session(dest, paste("GEO retrieval only:", paste(args, collapse = " ")))
  cat("Saved GEO materials:", dest, "\n")
}, error = function(e) {
  jsonlite::write_json(list(accession = args[1], mode = args[2], time_utc = stamp,
    status = "failed", error = conditionMessage(e)), file.path(dest, "last_failure.json"), auto_unbox = TRUE, pretty = TRUE)
  stop(conditionMessage(e), call. = FALSE)
})
