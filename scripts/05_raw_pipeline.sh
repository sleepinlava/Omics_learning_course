#!/usr/bin/env bash
###############################################################################
# 05_raw_pipeline.sh
#
# 从原始 FASTQ 到 count 矩阵的完整命令行流程 (可练手)
#   输入: data/raw_fastq/nf-core-testdata/{testdata,reference}
#   步骤: FastQC -> fastp 修剪 -> Salmon 定量 -> (可选) tximport/DESeq2
#        以及 HISAT2 比对 -> samtools -> featureCounts 的 BAM 路线
#
# 依赖工具(见 scripts/06_install_tools.sh 或 nf-core 的 docker 方式):
#   fastqc fastp salmon hisat2 samtools subread(featureCounts) multiqc
#
# 用法:  bash scripts/05_raw_pipeline.sh
###############################################################################
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NF="$ROOT/data/raw_fastq/nf-core-testdata"
REF="$NF/reference"
FQ="$NF/testdata/GSE110004"
WORK="$ROOT/results/raw_pipeline"
mkdir -p "$WORK"/{qc,trimmed,salmon,hisat2}

need() { command -v "$1" >/dev/null 2>&1 || { echo "❌ 缺少工具: $1  ($2)"; MISSING=1; }; }
MISSING=0
need fastqc   "conda install -c bioconda fastqc"
need fastp    "conda install -c bioconda fastp"
need salmon   "conda install -c bioconda salmon"
need hisat2   "conda install -c bioconda hisat2"
need samtools "conda install -c bioconda samtools"
need featureCounts "conda install -c bioconda subread"
[ "$MISSING" = 1 ] && { echo; echo "请先安装依赖:  bash scripts/06_install_tools.sh"; exit 1; }

if [ ! -d "$NF/reference" ]; then
  echo "❌ 找不到测试数据, 请先运行: bash scripts/00_download.sh"; exit 1
fi

# 选取 4 个样本练手 (2 个对照 + 1 个处理, 及 1 个单端)
declare -A PAIRED=( [SRR6357070]=1 [SRR6357072]=1 [SRR6357076]=1 )
SINGLES=( SRR6357074 )

echo "==> 1/4  FastQC (原始 read 质量)"
fastqc -t 4 -o "$WORK/qc" "$FQ"/*.fastq.gz

echo "==> 2/4  fastp 修剪 (去接头 + 质量过滤)"
for s in "${!PAIRED[@]}"; do
  fastp -i "$FQ/${s}_1.fastq.gz" -I "$FQ/${s}_2.fastq.gz" \
        -o "$WORK/trimmed/${s}_1.fq.gz" -O "$WORK/trimmed/${s}_2.fq.gz" \
        --detect_adapter_for_pe -q 20 -l 30 \
        -h "$WORK/qc/${s}.fastp.html" -j "$WORK/qc/${s}.fastp.json" 2>/dev/null
done
for s in "${SINGLES[@]}"; do
  fastp -i "$FQ/${s}_1.fastq.gz" -o "$WORK/trimmed/${s}.fq.gz" \
        -q 20 -l 30 -h "$WORK/qc/${s}.fastp.html" -j "$WORK/qc/${s}.fastp.json" 2>/dev/null
done

echo "==> 3/4  Salmon 定量 (选择性比对, 推荐快速路线)"
# 优先使用预构建索引 (reference/salmon.tar.gz), 否则现场构建
IDX="$WORK/salmon_index"
if [ ! -d "$WORK/salmon_index/salmon" ]; then
  if [ -f "$REF/salmon.tar.gz" ]; then
    echo "    解压预构建 salmon 索引"
    mkdir -p "$IDX"
    tar -xzf "$REF/salmon.tar.gz" -C "$IDX"
  else
    echo "    现场构建索引 (transcriptome.fasta)"
    salmon index -t "$REF/transcriptome.fasta" -i "$IDX" -k 31
  fi
fi
# 找到含 versionInfo.json 的索引目录
IDX=$(dirname "$(find "$WORK/salmon_index" -name versionInfo.json | head -1)")
echo "    salmon index = $IDX"
for s in "${!PAIRED[@]}"; do
  salmon quant -i "$IDX" -l A \
    -1 "$WORK/trimmed/${s}_1.fq.gz" -2 "$WORK/trimmed/${s}_2.fq.gz" \
    --validateMappings -p 4 -o "$WORK/salmon/$s"
done
for s in "${SINGLES[@]}"; do
  salmon quant -i "$IDX" -l A -r "$WORK/trimmed/${s}.fq.gz" \
    --validateMappings -p 4 -o "$WORK/salmon/$s"
done

echo "==> 4/4  HISAT2 比对 + featureCounts (BAM 路线)"
if [ ! -d "$WORK/hisat2_index" ]; then
  mkdir -p "$WORK/hisat2_index"
  if [ -f "$REF/hisat2.tar.gz" ]; then
    tar -xzf "$REF/hisat2.tar.gz" -C "$WORK/hisat2_index"
  else
    hisat2-build "$REF/genome.fa" "$WORK/hisat2_index/genome"
  fi
fi
# hisat2 -x 需要的是索引前缀(不含 .1.ht2)
HIDX=$(find "$WORK/hisat2_index" -name '*.1.ht2' | head -1 | sed 's/\.1\.ht2$//')
echo "    hisat2 index = $HIDX"

for s in "${!PAIRED[@]}"; do
  hisat2 -p 4 -x "$HIDX" \
    -1 "$WORK/trimmed/${s}_1.fq.gz" -2 "$WORK/trimmed/${s}_2.fq.gz" 2>"$WORK/hisat2/${s}.log" \
  | samtools sort -@ 2 -o "$WORK/hisat2/${s}.sorted.bam" -
  samtools index "$WORK/hisat2/${s}.sorted.bam"
done

featureCounts -T 4 -p -a "$REF/genes.gtf" -o "$WORK/counts.txt" \
  "$WORK"/hisat2/*.sorted.bam

if command -v multiqc >/dev/null 2>&1; then
  multiqc -o "$WORK/multiqc" "$WORK" >/dev/null 2>&1 || true
fi

cat <<EOF

✅ 原始数据流程完成! 输出在: $WORK
   - QC 报告:      $WORK/qc/ , $WORK/multiqc/multiqc_report.html
   - Salmon 定量:  $WORK/salmon/<sample>/quant.sf
   - count 矩阵:   $WORK/counts.txt   (可直接喂给 DESeq2/edgeR)

下一步: 用 tximport 把 salmon 结果合并成基因水平矩阵, 或直接读 counts.txt,
        参考教程 hbc-intro-to-rnaseq / nf-core。
EOF
