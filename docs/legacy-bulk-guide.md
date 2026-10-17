# RNA-seq 组学分析 · 练手工作区

一个自带**数据 + 教程 + 可运行脚本**的 bulk RNA-seq 学习环境。
从「已经数好的 count 矩阵」一路练到「原始 FASTQ → 比对 → 定量 → 差异表达 → 富集」。

> 位置: `~/rnaseq-learning`　|　系统已装: R 4.6.1 + 大量 Bioconductor 包、Python、Docker

---

## 1. 快速开始 (3 步)

```bash
cd ~/rnaseq-learning

# ① 数据/教程若还没下全，补下 (幂等，可重复运行)
bash scripts/00_download.sh

# ② 检查 R 依赖 (airway/pasilla/DESeq2/edgeR/limma/clusterProfiler...)
Rscript scripts/01_setup_R.R

# ③ 跑第一个完整差异表达分析，结果输出到 results/
Rscript scripts/02_DE_DESeq2_airway.R
```

跑完打开 `results/` 里的 PCA、火山图、MA 图、热图，你就完成了第一次 RNA-seq 分析 🎉

---

## 2. 目录结构

```
rnaseq-learning/
├── data/
│   ├── airway/         GSE52778 处理后的表达矩阵 (FPKM / 差异表)
│   ├── counts/         airway & pasilla 的 count 矩阵 (csv/tsv，可直接读)
│   ├── raw_fastq/
│   │   └── nf-core-testdata/    原始 FASTQ + 小型参考基因组 + 预建索引
│   │       ├── testdata/GSE110004/   GSE110004 的 fastq.gz (7 样本)
│   │       ├── reference/            genome.fa, genes.gtf, transcriptome.fasta
│   │       │                         salmon.tar.gz / hisat2.tar.gz (预建索引!)
│   │       └── samplesheet/          官方 samplesheet.csv 示例
│   ├── reference/      参考基因组 (软链到上面的小参考)
│   └── pasilla/        (占位)
├── tutorials/          离线教程：网页版 vignette + clone 下来的课程仓库
├── scripts/            所有可运行脚本 (见下表)
└── results/            所有输出 (图、表)
```

---

## 3. 数据清单

| 数据集 | 物种 / 设计 | 大小 | 适合练什么 | 如何加载 |
|---|---|---|---|---|
| **airway** (GSE52778) | 人气道平滑肌细胞；地塞米松处理 vs 未处理；4 细胞系×2 处理，**配对** | 小 | 从 count 矩阵做 DESeq2 全套（含配对 design） | `library(airway); data(airway)` 或 `data/counts/airway_counts.csv` |
| **pasilla** | 果蝇；基因敲低 vs 对照；7 样本（single/paired 混合） | 小 | edgeR / limma-voom，多方法对比 | `data/counts/pasilla_gene_counts.tsv` |
| **GSE110004** | 酿酒酵母 *S. cerevisiae*；WT vs RAP1 突变；101bp 双端链特异 RNA-seq (Wu et al. 2018, Mol Cell)，nf-core 官方测试集，已下采样到 5万 reads/样本 | ~50 MB | 原始 FASTQ → QC → 比对 → 定量 | `data/raw_fastq/nf-core-testdata/testdata/GSE110004/` |
| **小型参考基因组** | *S. cerevisiae* R64-1-1 (Ensembl) **下采样子集**：genome.fa + genes.gtf + transcriptome.fasta | < 1 MB | 建/用索引、比对、featureCounts | `data/raw_fastq/nf-core-testdata/reference/` |
| **原核测试集** | *S. enterica* Typhimurium SL1344 模拟 reads + 参考 | 小 | 原核 RNA-seq 练习 | `.../reference/prokaryotic/`, `.../testdata/prokaryotic/` |

**为什么要配这个小参考？** 它自带 `salmon.tar.gz`、`hisat2.tar.gz`、`bowtie2_rrna_index.tar.gz` **预建索引**，你不用花几十分钟建索引、也不用下载几十 GB 的人类基因组，就能把完整比对流程跑通。
（注意：这是**下采样**的参考，所以比对率不会很高，这是测试流程用的小数据，属正常现象。）

> 想练真实物种（人/小鼠/酵母）时，到 **Ensembl** 下 genome+GTF、**GENCODE** 下注释，
> 再用 STAR/HISAT2 建索引即可（README 末有链接）。

---

## 4. 脚本清单

| 脚本 | 作用 | 命令 |
|---|---|---|
| `00_download.sh` | 下载全部数据 + 教程（幂等） | `bash scripts/00_download.sh` |
| `01_setup_R.R` | 检查/安装 R/Bioconductor 依赖 | `Rscript scripts/01_setup_R.R` |
| `02_DE_DESeq2_airway.R` | **DESeq2 完整流程**：QC→归一化→DE→可视化+导出 | `Rscript scripts/02_DE_DESeq2_airway.R` |
| `03_DE_edgeR_limma_pasilla.R` | **edgeR + limma-voom + 方法对比** | `Rscript scripts/03_DE_edgeR_limma_pasilla.R` |
| `04_enrichment.R` | 差异基因的 GO/KEGG 富集 | `Rscript scripts/04_enrichment.R` |
| `05_raw_pipeline.sh` | 原始 FASTQ 全流程（FastQC→fastp→Salmon/HISAT2→featureCounts） | `bash scripts/05_raw_pipeline.sh` |
| `06_install_tools.sh` | 用 Miniforge/mamba 安装生信命令行工具 | `bash scripts/06_install_tools.sh` |
| `07_run_nfcore_rnaseq.sh` | 用 nf-core/rnaseq（Docker）跑工业级流程 | `bash scripts/07_run_nfcore_rnaseq.sh` |

建议顺序：**02 → 03 → 04 → (安装工具) 05 → 07**。

---

## 5. 推荐学习路线

### 阶段 0 · 打基础（可选，1–2 天）
- R 基础 + tidyverse；`tutorials/nf-core-rnaseq-pipeline` 里的文档不是入门，先看下面 HBC。
- 看 StatQuest 的 RNA-seq / DESeq2 / PCA 视频（链接见文末）。

### 阶段 1 · 差异表达分析（核心，2–4 天）
从 count 矩阵出发，这是 90% 课题的日常。
1. 读 **`tutorials/hbc-intro-to-dge/lessons/`**（Harvard Chan 的 DGE 课程，最通俗）
2. 跑 **`scripts/02_DE_DESeq2_airway.R`**，逐行对照理解
3. 读官方 vignette：`tutorials/03_DESeq2-vignette.html`、`02_limma-RNAseq123-workflow.html`
4. 跑 **`scripts/03_DE_edgeR_limma_pasilla.R`** 体会三种方法差异
5. 概念重点：**归一化 / size factor / 离散度(dispersion) / 设计矩阵 / 配对设计 / LFC 收缩 / 多重检验校正**

### 阶段 2 · 上游原始数据处理（3–5 天）
1. 安装工具：`bash scripts/06_install_tools.sh`（或直接用 Docker）
2. 读 **`tutorials/hbc-intro-to-rnaseq/`**（HPC 版，含 FASTQ/比对/计数全过程）
3. 跑 **`scripts/05_raw_pipeline.sh`**：FastQC → fastp → Salmon → featureCounts
4. 理解：**read QC / 接头修剪 / 比对 vs 选择性比对 / 计数 / 链特异性(strandedness) / rRNA 污染**

### 阶段 3 · 下游功能解释（2–3 天）
1. 跑 **`scripts/04_enrichment.R`**（GO / KEGG）
2. 读 `tutorials/nf-core-rnaseq-pipeline` 与 clusterProfiler 相关章节
3. 进阶：GSEA、WGCNA（共表达网络）、可变剪接（rMATS/DEXSeq）

### 阶段 4 · 可复现流程（进阶）
1. 跑 **`scripts/07_run_nfcore_rnaseq.sh`**，看 MultiQC 报告
2. 学 Nextflow / Snakemake 写自己的流程

---

## 6. 离线教程清单 (`tutorials/`)

**网页 / PDF vignette**
| 文件 | 内容 |
|---|---|
| `01_rnaseqGene-workflow-Love2016.html` | Bioconductor 官方 RNA-seq 工作流（DESeq2 作者亲写，**首推**） |
| `02_limma-RNAseq123-workflow.html` | limma/edgeR/voom 官方完整工作流 |
| `03_DESeq2-vignette.html` | DESeq2 包官方手册（最权威，当字典用） |
| `04_edgeR-UsersGuide.pdf` | edgeR 官方手册 |
| `05_limma-UsersGuide.pdf` | limma 官方手册 |
| `06_apeglm-shrinkage.html` | LFC 收缩原理 |
| `07_OSCA-single-cell-book.html` | 单细胞分析书（进阶，扩展视野） |
| `08_csawBook-differential-analysis.html` | csaw 差异分析书（计数数据统计建模，很推荐） |
| `09_nf-core-rnaseq-docs.html` | nf-core/rnaseq 流程文档 |
| `10_galaxy-transcriptomics.html` | Galaxy 转录组教程（点鼠标版） |
| `11_multiqc-docs.html` | MultiQC 汇总报告 |
| `12_fastqc-home.html` | FastQC 质量报告 |

**代码课程仓库（markdown 逐课）**
| 目录 | 内容 |
|---|---|
| `hbc-intro-to-rnaseq/` | Harvard Chan：RNA-seq 全流程（HPC 版，含课程数据） |
| `hbc-intro-to-dge/` | Harvard Chan：差异表达 + 可视化（Mov10 数据集） |
| `combine-rnaseq-in-R/` | COMBINE Australia：R 做 RNA-seq（bookdown 源） |
| `monash-rnaseq-DE/` | Monash：RNA-seq DE 分析 |
| `nf-core-rnaseq-pipeline/` | nf-core/rnaseq 流程源码与文档 |
| `conquer-DE-benchmark/` | 差异表达方法基准测评（进阶） |

> 用 VS Code 或 Typora 打开 `.md`/`.html` 阅读；`.html` 直接双击用浏览器打开。

---

## 7. 命令行工具安装

生信命令行工具用 **conda/bioconda** 管理最省心：

```bash
bash scripts/06_install_tools.sh      # 安装 Miniforge + 名为 rnaseq 的环境
conda activate rnaseq                 # 之后即可用 fastqc/salmon/hisat2/star/...
```

或直接用 **Docker**（无需本地装任何软件，见 `scripts/07_run_nfcore_rnaseq.sh`）：

```bash
sudo systemctl start docker          # 启动 docker 守护进程
bash scripts/07_run_nfcore_rnaseq.sh
```

---

## 8. 常见坑（血泪总结）

1. **别用 FPKM/TPM 做 DESeq2/edgeR 的输入** —— 它们要的是**原始 count**。FPKM 仅供你自己看表达量。
2. **count 必须是整数**，别四舍五入成小数。
3. **design 公式决定一切**：配对/批次要写进去，如 `~ cell + dex`、`~ batch + condition`。
4. **参考水平(reference level)** 决定了 log2FC 的方向，用 `relevel()` 明确指定。
5. **样本名要对齐** —— count 矩阵列顺序必须和 metadata 行顺序一致（`03` 脚本里有 `stopifnot` 检查）。
6. **低表达基因先过滤**（`filterByExpr` / `rowSums(counts>=10)>=n`），能显著提升检验功效。
7. **链特异性(strandedness)** 选错会导致计数几乎为 0 —— 不确定就先用 `reverse`/`auto` 试。
8. **p 值看 `padj`/`FDR`，不要看 `pvalue`**。
9. 生物学重复 **n≥3**；n=2 只能做探索性分析。
10. **KEGG 富集需要联网** 访问 `rest.kegg.jp`；国内网络可能超时（`04_enrichment.R` 已做容错，GO 可离线跑）。可改用 `clusterProfiler::enrichKEGG` 的本地注释或 `ReactomePA`。

---

## 9. 优质在线资源

- **StatQuest**（B站/YouTube，最好的入门视频）: DESeq2、PCA、RNA-seq 系列
- Bioconductor 官方 workflows: https://bioconductor.org/packages/release/workflows/
- Harvard Chan Bioinformatics Core: https://hbctraining.github.io/
- nf-core: https://nf-co.re/ 　|　Galaxy Training: https://training.galaxyproject.org/
- Biostars（问答）: https://www.biostars.org/ 　|　SEQanswers 论坛
- 参考基因组: Ensemble https://ftp.ensembl.org/ ｜ GENCODE https://www.gencodegenes.org/ ｜ UCSC https://genome.ucsc.edu/
- 公共数据: GEO https://www.ncbi.nlm.nih.gov/geo/ ｜ SRA https://www.ncbi.nlm.nih.gov/sra ｜ ArrayExpress/BioStudies ｜ recount3 https://rna.recount.bio/

---

## 10. 许可与出处

数据均来自公开资源：`airway`/`pasilla`（Bioconductor）、GSE52778（GEO）、
GSE110004 与配套参考（[nf-core/test-datasets](https://github.com/nf-core/test-datasets)，MIT/CC）。
教程版权归各自作者所有，仅作个人学习使用。
