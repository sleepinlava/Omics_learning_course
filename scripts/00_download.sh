#!/usr/bin/env bash
###############################################################################
# RNA-seq learning workspace bootstrap
#   - practice datasets (processed counts + raw FASTQ + small reference)
#   - tutorials (markdown repos, HTML/PDF vignettes)
# Idempotent: re-running only fetches what is missing.
###############################################################################
set -uo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
DATA="$ROOT/data"
TUT="$ROOT/tutorials"
RAW="$DATA/raw_fastq"
mkdir -p "$DATA"/{airway,pasilla,counts,reference} "$RAW" "$TUT" "$ROOT/results" "$ROOT/scripts"

CURL=(curl -L --fail --retry 8 --retry-all-errors --retry-delay 3 \
      --connect-timeout 20 --max-time 1800 --progress-bar -sS)
log() { printf '\n\033[1;34m==> %s\033[0m\n' "$*"; }

###############################################################################
log "1/5  GEO processed data: GSE52778 (human airway, dexamethasone vs untreated)"
###############################################################################
A="$DATA/airway"
get() { [ -s "$2" ] || "${CURL[@]}" -o "$2" "$1"; }   # get <url> <file>
get "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE52nnn/GSE52778/suppl/GSE52778_All_Sample_FPKM_Matrix.txt.gz" \
    "$A/GSE52778_All_Sample_FPKM_Matrix.txt.gz"
get "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE52nnn/GSE52778/suppl/GSE52778_Dex_vs_Untreated_gene_exp.diff.gz" \
    "$A/GSE52778_Dex_vs_Untreated_gene_exp.diff.gz"

###############################################################################
log "2/5  Raw FASTQ + matched small reference (nf-core/rnaseq test dataset)"
###############################################################################
# Self-contained kit: testdata/*.fastq.gz + reference/genome.fa+genes.gtf
# and prebuilt salmon / hisat2 / bowtie2 indexes -> full pipeline runnable.
NF="$RAW/nf-core-testdata"
if [ ! -d "$NF/reference" ]; then
  TARBALL=$(mktemp /tmp/rnaseq-testdata.XXXX.tar.gz)
  echo "   downloading branch tarball (~250 MB)..."
  "${CURL[@]}" -o "$TARBALL" \
    "https://codeload.github.com/nf-core/test-datasets/tar.gz/refs/heads/rnaseq" \
  && mkdir -p "$NF" \
  && tar -xzf "$TARBALL" -C "$NF" --strip-components=1 \
  && rm -f "$TARBALL"
  # prune heavy files we don't need (BAMs, STAR/RSEM indexes, uncompressed fastq)
  rm -rf "$NF"/testdata/uncompressed "$NF"/reference/star.tar.gz "$NF"/reference/rsem.tar.gz
  find "$NF" -name '*.bam' -delete 2>/dev/null || true
else
  echo "   already present, skipping"
fi

# Convenience symlinks so the flat data layout is easy to browse
ln -sfn "$NF/reference" "$DATA/reference/nf-core-small-reference" 2>/dev/null || true

###############################################################################
log "3/5  Count matrices for DE analysis (exported from Bioc packages)"
###############################################################################
if Rscript -e 'suppressMessages(library(airway)); data(airway)' 2>/dev/null; then
  Rscript -e 'suppressMessages(library(airway)); data(airway); d<-assay(airway)
    write.csv(as.data.frame(d), "'"$DATA"'/counts/airway_counts.csv")
    write.csv(as.data.frame(colData(airway)), "'"$DATA"'/counts/airway_coldata.csv")
    cat("   airway counts exported (", nrow(d), "genes x", ncol(d), "samples )\n")' 2>/dev/null || true
else echo "   airway not installed yet - run scripts/01_setup_R.R"; fi

if Rscript -e 'suppressMessages(library(pasilla))' 2>/dev/null; then
  Rscript -e 'suppressMessages(library(pasilla))
    file.copy(system.file("extdata","pasilla_gene_counts.tsv",package="pasilla"),
              "'"$DATA"'/counts/pasilla_gene_counts.tsv", overwrite=TRUE)
    file.copy(system.file("extdata","pasilla_sample_annotation.csv",package="pasilla"),
              "'"$DATA"'/counts/pasilla_sample_annotation.csv", overwrite=TRUE)
    cat("   pasilla counts exported\n")' 2>/dev/null || true
else echo "   pasilla not installed yet - run scripts/01_setup_R.R"; fi

###############################################################################
log "4/5  Tutorials: HTML / PDF vignettes and books"
###############################################################################
T="$TUT"
dl() { [ -s "$2" ] && { echo "   skip $(basename "$2")"; return; }
       echo "   get  $(basename "$2")"; "${CURL[@]}" -o "$2" "$1" || echo "   !! failed: $1"; }

dl "https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html" \
   "$T/01_rnaseqGene-workflow-Love2016.html"
dl "https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html" \
   "$T/02_limma-RNAseq123-workflow.html"
dl "https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html" \
   "$T/03_DESeq2-vignette.html"
dl "https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf" \
   "$T/04_edgeR-UsersGuide.pdf"
dl "https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf" \
   "$T/05_limma-UsersGuide.pdf"
dl "https://bioconductor.org/packages/release/bioc/vignettes/apeglm/inst/doc/apeglm.html" \
   "$T/06_apeglm-shrinkage.html"
dl "https://bioconductor.org/books/release/OSCA/" \
   "$T/07_OSCA-single-cell-book.html"
dl "https://bioconductor.org/books/release/csawBook/" \
   "$T/08_csawBook-differential-analysis.html"
dl "https://nf-co.re/rnaseq/latest/" \
   "$T/09_nf-core-rnaseq-docs.html"
dl "https://training.galaxyproject.org/training-material/topics/transcriptomics/" \
   "$T/10_galaxy-transcriptomics.html"
dl "https://docs.seqera.io/multiqc/" \
   "$T/11_multiqc-docs.html"
dl "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/" \
   "$T/12_fastqc-home.html"

###############################################################################
log "5/5  Tutorials: git repositories (markdown lessons + code)"
###############################################################################
clone() {
  local dst="$T/$2"
  [ -d "$dst/.git" ] && { echo "   skip $2"; return; }
  rm -rf "$dst"
  echo "   clone $2"
  # 优先直连 github (绕过 ghproxy 镜像); 失败再回退到默认配置
  GIT_CONFIG_GLOBAL=/dev/null git clone --depth 1 --quiet "$1" "$dst" 2>/dev/null \
    || git clone --depth 1 --quiet "$1" "$dst" 2>/dev/null \
    || { echo "   !! clone failed: $1"; rm -rf "$dst"; }
}
clone https://github.com/hbctraining/Intro-to-Rnaseq-HPC-O2.git    hbc-intro-to-rnaseq
clone https://github.com/hbctraining/Intro-to-DGE.git              hbc-intro-to-dge
clone https://github.com/COMBINE-Australia/RNAseq-R.git            combine-rnaseq-in-R
clone https://github.com/MonashBioinformaticsPlatform/RNAseq-DE-analysis-with-R.git monash-rnaseq-DE
clone https://github.com/nf-core/rnaseq.git                        nf-core-rnaseq-pipeline
clone https://github.com/csoneson/conquer_comparison.git           conquer-DE-benchmark

log "DONE.  Workspace: $ROOT"
du -sh "$DATA" "$TUT" 2>/dev/null
