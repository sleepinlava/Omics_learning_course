# 02 · bulk RNA-seq：从计数到可解释的处理效应

## 学习问题和上游路线

问题：在同一细胞系内部比较地塞米松处理与未处理，哪些基因呈现表达变化？独立配对单位是 4 个细胞系，8 个测序文库不是 8 个互不相关的受试者。

先读本地 `tutorials/hbc-intro-to-dge/lessons/` 中实验设计、count normalization、DESeq2 与结果解释相关课程，再查 [DESeq2 官方文档](https://bioconductor.org/packages/DESeq2/)。本课程调用同一个成熟包，额外加入输入校验与证据记录。

## 实验 A：已有 count 矩阵的完整分析

```bash
Rscript omics-course/scripts/02_bulk_airway.R
```

输入是仓库现成 `data/counts/airway_counts.csv` 和 `airway_coldata.csv`，不需要再次下载 FASTQ 或安装 airway。脚本检查非负整数、唯一 ID、样本一对一映射、配对结构和设计满秩；不允许默默丢失样本后继续。

| 步骤 | 做什么 | 先问自己的问题 |
|---|---|---|
| 读入 | 基因×样本整数表，匹配元数据 | 列顺序真的和样本行顺序一致吗？ |
| 过滤 | 至少 4 个样本 count≥10 | 是否在看组间 p 值之前规定了规则？ |
| 建模 | `~ cell + dex`，对比 trt/untrt | 配对与参考组是否正确？ |
| 归一化 | DESeq2 size factors | 是否有极端文库组成或异常样本？ |
| QC | VST、PCA、样本距离 | 处理、细胞系还是异常样本主导变化？ |
| 统计 | NB GLM、Wald、BH | 效应、标准误、padj 是否一起报告？ |
| 检查 | 完整表、MA、候选基因逐样本图 | 结果是否由一个样本驱动？ |

输出 `results/bulk/`：`all_genes.csv`、`tested_genes.csv`、`significant_genes.csv`、`filtering.csv`、`sample_qc.csv`、PCA/MA/样本距离图、`dds.rds` 和会话记录。完整基因表含预过滤状态；未检验不等于无效应。

`log2FoldChange > 0` 在本脚本中表示处理高于未处理。`significant_genes.csv` 使用 `padj < 0.05`，不额外用事后选定的 fold change 截断；研究中的实际意义门槛应预先讨论。

## 背后的模型与归一化

本章模型为：

$$\log\mu_{gi}=\log s_i+\beta_{g0}+\beta_{g,\mathrm{cell}}+\beta_{g,\mathrm{dex}}.$$

DESeq2 的经典 median-of-ratios 思路先对每个基因跨样本构造几何均值，再取样本相对该值的比值中位数。零值、参考基因选择和多数基因的整体变化都会影响适用性；标准处理交给包实现，不自己用总 reads 数强行替代。

VST 把均值依赖的方差变得更平稳，用于距离/PCA；DE 的拟合仍基于原始整数计数。p 值与 padj 可能为 NA，原因包括无法检验、异常值或独立过滤。脚本保留这些情况，不替换为零。

普通 Wald 置信区间是模型近似区间，不是收缩后的区间。低计数 fold change 的排名可进一步学习 `lfcShrink()`；课程基础脚本不把未收缩估计冒充稳定效应排名。原理参见 [DESeq2 方法论文](https://doi.org/10.1186/s13059-014-0550-8)。

## 实验 B：从原始 reads 回到矩阵

本仓库已有 `scripts/05_raw_pipeline.sh` 和 `07_run_nfcore_rnaseq.sh`，学习时沿用这些流程及 `tutorials/nf-core-rnaseq-pipeline/` 文档。先安装相应工具、检查输入与链特异性，再运行；它们不是本课程默认验证的一部分。

顺序：FASTQ/样本表 → FastQC/汇总 → 必要的接头处理 → Salmon 或 HISAT2/featureCounts → 样本 QC → 合格计数矩阵。不同路线解决问题不同，不要求同时跑所有比对工具。

- FASTA、GTF、转录本与索引必须是匹配版本；小参考仅用于执行练习。
- Salmon 的估计 `NumReads` 可能为小数；采用 tximport 与相应接口，不随意四舍五入后声称是直接观测 counts。
- 链特异性、paired-end 参数和计数软件版本需要按协议核对；完成比对不代表定量正确。
- 下采样小参考会影响比对率；不能套用真实全参考的经验阈值。

## 实验 C：方法比较与富集

已有 pasilla 脚本可用来学 edgeR 的 TMM、NB QL 与 limma-voom 的均值—方差权重。运行前查看 `condition × type` 列联表，评估是否需 `~ type + condition` 及是否可识别。现有教学脚本主要使用 condition，不能自动当作适合所有实验的最终模型。

比较方法时统一样本、对比方向、过滤集合和多重校正，再比较效应量与排序；不要只数“谁发现的基因多”。

富集是“候选基因是否相对背景过度集中于某类功能”。超几何检验中：$N$ 为可检验背景基因数，$K$ 为背景中通路基因数，$n$ 为候选数，$k$ 为候选落在通路中的数。背景应是测到、可检验且能映射的基因集合；用全基因组背景可能引入偏差。GSEA 使用排序而非先二分显著名单，但仍依赖有效的排序和基因集。

仓库 `scripts/04_enrichment.R` 是可读示例。它读取根目录 `results/` 的历史文件；本课程输出位于 `omics-course/results/`，不要不检查输入就混读旧结果。作为扩展练习，在副本中改为本课程路径并记录新的背景集合。

## 常见错误与作业

错误包括：拿 TPM/FPKM 输入 DESeq2；按列位置拼元数据；把配对当独立；拿 VST 当 raw counts；只展示显著基因、不交完整表；把 GO 富集当作机制已被证明。

作业：选择一个候选基因，画出每个细胞系处理前后连线；比较有/无配对模型的标准误；解释差异来自哪里。报告必须包含设计、对比、过滤数、QC、完整表、版本和局限。过关标准是“能解释设计”，不是复现一个固定显著基因数。
