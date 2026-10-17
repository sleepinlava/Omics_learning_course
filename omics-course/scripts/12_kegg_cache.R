#!/usr/bin/env Rscript
# Fetch a dated KEGG pathway snapshot via KEGGREST; keep locally, out of Git.
f <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])
source(file.path(dirname(normalizePath(f)), "_common.R"))
args <- commandArgs(TRUE)
if (length(args) != 1 || !grepl("^[a-z][a-z0-9]+$", args[1])) stop("Usage: 12_kegg_cache.R hsa|mmu|other_KEGG_organism")
need("KEGGREST")
organism <- args[1]
dest <- file.path(COURSE, "data", "kegg", paste0(organism, "_", format(Sys.Date(), "%Y%m%d"), ".rds"))
if (file.exists(dest)) stop("Today's snapshot already exists; retaining it: ", dest)
options(timeout = 120)
links <- KEGGREST::keggLink("pathway", organism)
Sys.sleep(0.4)
terms <- KEGGREST::keggList("pathway", organism)
if (!length(links) || !length(terms) || !all(startsWith(names(links), paste0(organism, ":")))) stop("Unexpected KEGG response.")
t2g <- unique(data.frame(term = sub("^path:", "", unname(links)),
                         gene = sub(paste0("^", organism, ":"), "", names(links))))
t2n <- data.frame(term = sub("^path:", "", names(terms)), name = unname(terms))
cache <- list(organism = organism, retrieved_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
              source = "https://rest.kegg.jp", term2gene = t2g, term2name = t2n,
              KEGGREST = as.character(packageVersion("KEGGREST")), synthetic = FALSE)
dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
saveRDS(cache, dest)
cat("KEGG snapshot saved:", dest, "\n")
