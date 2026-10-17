#!/usr/bin/env bash
###############################################################################
# 07_run_nfcore_rnaseq.sh  --  用工业级流程 nf-core/rnaseq 跑一遍完整分析
#
# nf-core/rnaseq 封装了 FastQC/TrimGalore/STAR/Salmon/RSEQC/MultiQC/DESeq2 等，
# 是当前最主流的可复现 RNA-seq 流程。用官方 test profile 可在几分钟内跑通。
#
# 前置:  nextflow (scripts/06_install_tools.sh 已安装) + docker
#         docker 守护进程需运行:  sudo systemctl start docker
#
# 用法:  bash scripts/07_run_nfcore_rnaseq.sh
###############################################################################
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$ROOT/results/nf-core-rnaseq-test"
mkdir -p "$OUT"

command -v nextflow >/dev/null 2>&1 || {
  echo "❌ 未找到 nextflow。请先运行: bash scripts/06_install_tools.sh"; exit 1; }
docker info >/dev/null 2>&1 || {
  echo "❌ docker 守护进程未运行。请执行:  sudo systemctl start docker"; exit 1; }

# Docker Hub 上已有所有工具镜像, 无需本地安装任何生信软件
export NXF_VER="${NXF_VER:-24.10.0}"

echo "==> 运行 nf-core/rnaseq (test profile, docker)..."
nextflow run nf-core/rnaseq \
  -r 3.14.0 \
  -profile test,docker \
  --outdir "$OUT" \
  -resume

cat <<EOF

✅ 完成! 打开报告: $OUT/multiqc/multiqc_report.html
   星号(*)的样本关系图、比对率、基因计数、差异分析一应俱全。

学习建议: 用 --help 查看参数; 用你自己的 samplesheet.csv 替换 test:
   nextflow run nf-core/rnaseq -profile docker \\
     --input samplesheet.csv --fasta genome.fa --gtf genes.gtf --outdir results/my_run
   官方文档: tutorials/09_nf-core-rnaseq-docs.html
EOF
