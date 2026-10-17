# 多供者 pseudobulk 扩展：从细胞图谱走到处理效应

## 上游课程与数据

按 [HBC pseudobulk 自学路线](https://github.com/hbctraining/Pseudobulk-for-scRNAseq/blob/main/schedule/self-learning.md) 学习“聚合 counts → DESeq2 → 差异可视化”，完整材料可用 `08_get_tutorials.py fetch pseudobulk` 获取。

本扩展建议使用 [muscData 的 Kang18_8vs8](https://bioconductor.org/packages/release/data/experiment/vignettes/muscData/inst/doc/muscData.html)，来源 GSE96583。它具有 8 位供者的对照/IFN-β 刺激配对信息，适合学习处理效应；不要写成健康人与狼疮患者的疾病对比。首次加载通过 ExperimentHub 下载，需要联网及额外磁盘空间。

```bash
Rscript omics-course/env/install_R.R pseudobulk
```

在仓库根目录启动 R，先加载而不直接开始检验：

```r
.libPaths(c("omics-course/env/R-library", .libPaths()))
library(SingleCellExperiment)
sce <- muscData::Kang18_8vs8()
assayNames(sce)
names(colData(sce))
table(sce$ind, sce$stim)
table(sce$cell, useNA = "ifany")
table(sce$multiplets, useNA = "ifany")
```

上游当前示例字段包括 `ind`、`stim`、`cell` 等；未来对象变化时依据文档核实。保留原字段，再创建标准字段：`donor_id <- ind`、`condition <- stim`、`sample_id <- donor_id × condition`、`cell_type <- cell`。拼接后的 ID 要保证唯一映射；保留原 barcode，多样本的 barcode 可能重复。

## 必须完成的前置实验

1. 查看 donor/condition 的 QC 分布、原作者 multiplet 信息及注释；沿用 HBC/OSCA 的 QC 思路，不机械复制 PBMC3k 阈值。
2. 为每个细胞记录逻辑型 `passes_qc`，不能默认所有细胞合格；保存排除原因和支持图。
3. 选择一种细胞类型，检查每个 donor×condition 的细胞数，预定最低细胞数。低细胞数导致缺少半对样本时应先处理设计，不悄悄补零。
4. 保存审核后的 SingleCellExperiment 为 `omics-course/data/kang_reviewed.rds`。它必须具有 counts assay 和 `sample_id/donor_id/condition/cell_type/passes_qc` 五列。

这是科学判断环节，不是程序报错后的随意“修参数”。本课程不替尚未下载的数据伪造 QC 通过标记。

## 汇总与配对检验

脚本 `09_pseudobulk.R` 调用 `scuttle::aggregateAcrossCells` 求和，再调用 DESeq2。它只支持这里的简单完整配对设计，避免把复杂纵向或批次设计默默当成通用情况。使用实际读取到的细胞类型、condition 拼写；下面标签是格式示例：

```bash
Rscript omics-course/scripts/09_pseudobulk.R \
  omics-course/data/kang_reviewed.rds 'B cells' ctrl stim 20
```

最后的 20 是示例最低细胞数，须在看结果之前根据设计确定；不是通用阈值。若实际标签不是 `B cells/ctrl/stim`，先用 `table()` 确认后替换。脚本遇到不完整配对、少于 3 对、混乱样本映射或错误数值尺度时停止，并输出细胞数供检查。这个教学下限不是充分功效保证。

模型为 `~ donor_id + condition`，输出方向是最后一个条件相对于参考条件。结果保存在 `results/pseudobulk/`，一次只分析一个类型；换类型前先归档输出。多类型联合声明应预先定义跨类型检验族，不能把每个类型内的 FDR 当作所有类型共同 FDR。

## 解释与练习

对候选基因同时画供者级配对表达与细胞级分布。若某个供者贡献一半细胞，pseudobulk 的列仍是该供者的一个样本，不会变成更多生物学重复。但细胞数差异、细胞状态组成变化和 QC 仍然需要检查。

本扩展本次未在真实 Kang 数据上执行。过关目标是能说明聚合单位、完整配对和测试范围，并解释为什么 cluster marker 与供者处理效应不是同一个问题。
