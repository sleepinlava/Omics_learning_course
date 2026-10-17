# Shared local paths and provenance; no package installation or downloads.
a <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("^--file=", "", a[grepl("^--file=", a)]))
COURSE <- dirname(dirname(script_path))
REPO <- dirname(COURSE)
local_lib <- file.path(COURSE, "env", "R-library")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))
need <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "),
                           ". See omics-course/00_setup/README.md")
}
output_dir <- function(name) {
  out <- file.path(COURSE, "results", name)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out
}
record_session <- function(out, details, filename = "sessionInfo.txt") {
  writeLines(c(paste("Run time:", Sys.time()), details, "", capture.output(sessionInfo())),
             file.path(out, filename))
}
record_inputs <- function(paths, out) {
  write.csv(data.frame(path = normalizePath(paths), bytes = file.info(paths)$size,
                       md5 = unname(tools::md5sum(paths))),
            file.path(out, "input_files.csv"), row.names = FALSE)
}
