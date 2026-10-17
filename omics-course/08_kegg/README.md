# 08 · KEGG 通路分析：ORA 与 GSEA

## 来源与用途

本章复用 [KEGGREST](https://bioconductor.org/packages/release/bioc/vignettes/KEGGREST/inst/doc/KEGGREST-vignette.html)、[clusterProfiler](https://github.com/YuLab-SMU/clusterProfiler) 和 [KEGG API](https://www.kegg.jp/kegg/rest/keggapi.html)。通路注释会更新，分析时保存日期与版本；本仓库发布获取代码和引用，不打包上传 KEGG 数据库缓存。

## 获取日期快照

```bash
Rscript omics-course/env/install_R.R kegg
Rscript omics-course/scripts/12_kegg_cache.R hsa
```

结果是 `data/kegg/hsa_YYYYMMDD.rds`，包含 pathway→gene 和 pathway→name、来源与时间。`hsa` 为人，`mmu` 为小鼠。获取脚本可以缓存不同物种，但本章分析封装明确只接受**人类 ENSEMBL→ENTREZ**，不会把人类注释强行用于其他物种。

## 从已完成的 airway 分析开始

将下面 `YYYYMMDD` 替换为实际快照日期：

```bash
Rscript omics-course/scripts/13_kegg_enrichment.R omics-course/results/bulk/tested_genes.csv omics-course/data/kegg/hsa_YYYYMMDD.rds ora
Rscript omics-course/scripts/13_kegg_enrichment.R omics-course/results/bulk/tested_genes.csv omics-course/data/kegg/hsa_YYYYMMDD.rds gsea
```

脚本记录 unmapped、歧义映射和有效背景；保守地使用一对一 ENSEMBL/ENTREZ 子集。去除 ID 版本后出现重复会停止。该选择会漏掉部分基因，应报告映射率与潜在偏倚；不要把未映射当作“不参与通路”。

## ORA 的数学与生物学意义

在有效背景 $N$ 中，有 $K$ 个通路基因；从中选择 $n$ 个候选，观察到 $k$ 个属于该通路。超几何右尾概率为：

$$P(X\ge k)=\sum_{j=k}^{\min(K,n)}\frac{\binom Kj\binom{N-K}{n-j}}{\binom Nn}.$$

它问的是候选是否比背景更集中于某通路。背景应是确实可检验、能映射且有对应注释的基因；本脚本先取有有效 p/stat 的 tested genes，再与快照相交。不能拿全人类基因组当任何实验的默认背景。

上调和下调列表分开分析，候选规则为 padj<0.05，方向来自 log2FC；报告集合大小和比例。`BH_across_directions` 对所有满足背景集大小限制的通路×两个方向合并校正，包含未返回的零交集检验（p=1）；每个通路共享基因、通路之间相关，因此多个类似命中不能被当作多份独立机制证据。

## GSEA 的数学与输入

GSEA 沿完整排序表累积集合成员的加权贡献，与非成员的背景步进比较，得到富集分数；标准化后得到 NES。脚本使用完整可用的、有符号 Wald statistic 排序，而不是仅使用显著基因，也不把 padj 当作有方向的效应量。

排序定义、ID 去重、基因集大小、软件版本和随机种子都会影响结果。NES 正/负表示在当前对比的排序顶端/底端富集；并不自动意味着通路被激活/抑制。若用其他软件的排序分数，需要重新说明其方向和统计含义。

## 结果与练习

输出位于 `results/kegg_ora/` 或 `results/kegg_gsea/`：完整 pathway 表、ID 映射、有效 universe/ranking、输入摘要与会话。没有返回 eligible pathway 时留下说明，不用空表宣称“没有生物学变化”。

练习：选择同一通路比较 ORA 和 GSEA，列出实际驱动信号的基因，查看逐样本表达；改变预先定义的候选阈值后做敏感性分析。说明“富集”“活性”“因果机制”的区别，并写出所需验证实验。
