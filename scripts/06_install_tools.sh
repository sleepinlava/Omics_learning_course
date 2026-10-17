#!/usr/bin/env bash
###############################################################################
# 06_install_tools.sh  --  安装 RNA-seq 命令行工具 (Miniforge/mamba + bioconda)
#
# 会安装: FastQC, fastp, MultiQC, Salmon, HISAT2, STAR, samtools,
#         subread(featureCounts), SeqKit, trim-galore, nextflow
# 安装位置: ~/miniforge3   (环境名: rnaseq)
#
# 用法:  bash scripts/06_install_tools.sh
# 注意:  需要联网, 首次安装约 1-3 GB, 视网速需几分钟到十几分钟。
###############################################################################
set -euo pipefail

PREFIX="$HOME/miniforge3"
INSTALLER="/tmp/Miniforge3-Linux-x86_64.sh"

if [ ! -x "$PREFIX/bin/mamba" ]; then
  echo "==> 下载并安装 Miniforge 到 $PREFIX"
  curl -L --retry 5 --retry-all-errors -o "$INSTALLER" \
    "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
  bash "$INSTALLER" -b -p "$PREFIX"
else
  echo "==> Miniforge 已安装, 跳过"
fi

export PATH="$PREFIX/bin:$PATH"
conda init bash >/dev/null 2>&1 || true

echo "==> 创建/更新 conda 环境: rnaseq"
if mamba env list | grep -qE '^rnaseq\s'; then
  mamba install -y -n rnaseq -c conda-forge -c bioconda \
    fastqc fastp multiqc salmon hisat2 star samtools subread seqkit trim-galore nextflow
else
  mamba create -y -n rnaseq -c conda-forge -c bioconda \
    fastqc fastp multiqc salmon hisat2 star samtools subread seqkit trim-galore nextflow
fi

cat <<'EOF'

✅ 安装完成。启用方式(任选):
   1) 新开终端会自动加载; 或执行:
        conda activate rnaseq
   2) 不激活, 直接用绝对路径:
        ~/miniforge3/envs/rnaseq/bin/salmon --version

验证:
   conda activate rnaseq && salmon --version && hisat2 --version && fastqc --version

随后即可运行原始数据流程:
   bash scripts/05_raw_pipeline.sh
EOF
