# Repository Guidelines

## Project Structure & Module Organization

This is a bulk RNA-seq learning workspace. Maintain executable examples in `scripts/`, numbered `00`–`07` in learning order. `data/` holds count matrices, sample annotations, FASTQ files, and references. `results/` contains generated tables, plots, and pipeline outputs. `tutorials/` contains downloaded documentation and upstream repositories; preserve their attribution and avoid unrelated edits. `README.md` describes the learning workflow; `MANIFEST.txt` inventories workspace files.

## Build, Test, and Development Commands

Run commands from the repository root. There is no top-level build step.

- `Rscript scripts/01_setup_R.R`: check and install R/Bioconductor dependencies.
- `Rscript scripts/02_DE_DESeq2_airway.R`: run airway differential expression and generate plots/tables.
- `Rscript scripts/03_DE_edgeR_limma_pasilla.R`: compare pasilla differential-expression methods.
- `Rscript scripts/04_enrichment.R`: run enrichment after script `02`; KEGG requires network access.
- `bash scripts/06_install_tools.sh`: install Miniforge and the `rnaseq` environment; then use `conda activate rnaseq`.
- `bash scripts/05_raw_pipeline.sh`: run FASTQ QC, quantification, alignment, and counting.
- `bash scripts/07_run_nfcore_rnaseq.sh`: run the pinned nf-core test workflow using Nextflow and a running Docker daemon.

Before using `bash scripts/00_download.sh`, adjust its hard-coded `ROOT="$HOME/rnaseq-learning"` to this checkout. It downloads resources and can replace tutorial directories lacking Git metadata.

## Coding Style & Naming Conventions

Use two-space indentation, descriptive names, and numbered filenames such as `08_new_analysis.R`. Follow existing R conventions: `<-` assignment, `file.path()` for paths, and explicit package loading. Quote Bash path variables and follow the pipeline scripts' `set -euo pipefail` convention. Derive paths from script locations. Name outputs by dataset and analysis, such as `pasilla_edgeR_results.csv`. No top-level formatter or linter is configured.

## Testing Guidelines

There is no workspace-wide test suite or coverage threshold; tests inside downloaded tutorials belong to those projects. Check changed shell scripts with `bash -n scripts/<file>.sh` and R syntax with `Rscript -e 'parse(file="scripts/<file>.R")'`. Run the affected example and inspect generated tables and plots. Preserve sample-order assertions, explicit contrasts, and random seeds. Use raw counts for count-based differential expression.

## Commit & Pull Request Guidelines

Git history is unavailable in this checkout, so no existing commit convention can be verified. Use concise imperative subjects, such as `Fix pasilla sample alignment`. Describe the change, affected scripts, validation commands, and output differences in pull requests. Link relevant issues and include plots when visualization changes. Avoid unrelated generated-data changes.
