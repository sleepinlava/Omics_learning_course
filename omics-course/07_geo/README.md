# 07 · GEO 公共数据获取与重分析

## 目标与上游材料

GEO 是公共数据仓库，不是某一种测序技术或统计模型。先识别平台、测量尺度和设计，再选择分析方法。主教程为 [GEOquery 官方 vignette](https://bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.html)，差异分析继续使用前面的 DESeq2/limma 课程。

- GSE：一个研究系列，可能含多个平台或子研究。
- GSM：提交的样本记录，不自动等于独立供者。
- GPL：平台/探针注释；芯片 probe ID 不自动等于基因 ID。
- SRA Run：测序运行，多个 run 可能对应同一个生物样本。

## 实验一：获取并保存原始元数据

```bash
Rscript omics-course/env/install_R.R geo
Rscript omics-course/scripts/10_geo_fetch.R GSE2553 matrix
Rscript omics-course/scripts/10_geo_fetch.R GSE52778 supp-list
Rscript omics-course/scripts/10_geo_fetch.R GSE164073 counts
```

GSE2553 和 GSE164073 是 GEOquery 文档中的获取示例，不构成可以直接合并的研究队列。`matrix` 按平台分别保存表达矩阵、元数据、特征表和 ExpressionSet；`supp-list` 仅列附加文件；`counts` 先查询 NCBI 是否提供预计算 RNA-seq 数据。所有结果在 `data/geo/<GSE>/<mode>/`。首次使用需联网；失败会留下 `last_failure.json`，不会把失败标成完成。

## 实验二：建立可审阅的样本表

阅读论文、GEO 的 characteristics 和处理流程，填写 `templates/geo_samples.csv`。必须明确：sample_id、donor_id、condition、是否配对、组织、平台、批次、排除理由和原文依据。逐一比对矩阵列与 GSM，不能按文件名正则随便猜疾病/对照。

若只选择部分样本，先显式导出对应子矩阵，写出排除清单，再交给分析脚本。脚本要求输入列与样本表严格一对一，不静默丢弃额外样本。发现供者复用时采用配对/重复测量模型；GSM 数不是自动的独立样本数。

## 实验三：依据输入尺度选择模型

| 输入 | 正确入口 | 需避免 |
|---|---|---|
| 整数 RNA-seq counts | DESeq2 / edgeR；参照课程 02 | 把 log、TPM、FPKM 四舍五入成 counts |
| Salmon 转录本估计 | tximport 与对应模型接口 | 把估计数当作直接观测整数 |
| 已归一化 log2 芯片表达 | `11_geo_limma.R` | 再次 log2；把 probe 直接当基因 |
| CEL 等原始芯片文件 | 先按平台进行背景校正、归一化与 QC | 与其他平台原始值直接拼接 |
| 未知尺度 Series Matrix | 回查处理协议，暂不拟合 | 根据最大值自动猜是否 log |

示例命令（路径与分组必须替换成已核对的实际内容）：

```bash
Rscript omics-course/scripts/11_geo_limma.R expression.csv.gz samples.csv control treated paired log2_verified
```

可选 `independent`；若样本表有多个 batch，模型会纳入 batch 并检查满秩。程序不能从数值本身证明数据已正确归一化，`log2_verified` 表示分析者已核对协议。默认输出到 `results/geo_limma/`，换研究前归档该目录。

## 输出、注释与 meta 的衔接

输出完整 feature 表、模型、设计矩阵、PCA、过滤记录和输入摘要。`SE_classical` 是普通线性模型的采样标准误；`SE_moderated`、p 值和区间来自 limma 的经验贝叶斯推断。跨研究合并前，应预定可比的效应和标准误定义，不拿调整后的 p 值反推 SE。

probe→gene 映射先保留一对多关系，再预定去歧义/汇总规则；不能每个研究都挑 p 最小的 probe。基因版本、物种、比较方向一致后，才能进入 meta 分析。研究内 batch 校正并不保证不同研究的效应天然可比。

## 作业

为一个 GSE 提交研究设计摘要、GSM 与供者对应表、平台和尺度证据、排除清单、QC 与完整结果。列出一个“这个数据下载得到，但无法回答当前问题”的例子。跨研究练习必须先检查样本是否来自同一原始队列。
