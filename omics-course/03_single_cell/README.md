# 03 · scRNA-seq：学会识别细胞，也学会承认不确定

## 学习目标与阅读顺序

目标是从过滤后的 UMI 矩阵建立一个可审阅的细胞图谱，理解 QC、归一化、降维、聚类和注释。随后才学习多供者差异分析。

先按 [Seurat PBMC3k](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) 理解对象和函数，再将步骤映射到 [Scanpy tutorials](https://scanpy.readthedocs.io/en/stable/tutorials/index.html)。需要判断统计方法时查 [OSCA](https://bioconductor.org/books/release/OSCA/) 与 [Single-cell Best Practices](https://www.sc-best-practices.org/)。R 与 Python 是两条实现路线，不必初学时同时安装两套。

## 实验 A：离线认识对象

```bash
Rscript omics-course/scripts/03_single_cell.R --demo
```

数据来自 SeuratObject 的 `pbmc_small`：80 个细胞，基因已被缩减。脚本重新从其 counts 构建教学对象，导出矩阵/元数据，并执行归一化、PCA 和邻接聚类。它不能代表完整测序数据的基因检测率或线粒体比例，也不用于生物学 QC 阈值训练。这个模式不做差异显著性推断。

打开 `results/single_cell_demo/`：核对细胞数、特征数、PCA 与聚类。`demo` 和完整 PBMC 输出分开，防止误把小对象结果当作正式分析。

## 实验 B：完整 PBMC3k

```bash
python omics-course/scripts/04_fetch_data.py pbmc3k
Rscript omics-course/scripts/03_single_cell.R --qc-only
```

先看 `results/single_cell_pbmc3k/qc_before.png` 与 `qc_cells.csv`。逐细胞计算检测基因数、UMI 总量和线粒体比例；看分布和它们之间的关系，再决定过滤规则。

为复现这个历史 PBMC 教学案例，可显式运行：

```bash
Rscript omics-course/scripts/03_single_cell.R --min-features 200 --max-features 2500 --max-mt 5
```

这些数字只用于该案例的比较练习，不移植到脑组织、肿瘤、单核或其他平台。脚本保留每个原始 barcode 的 `passes_qc`，并记录参数；不能仅保存最终“好看”的细胞。

## 每一步到底做了什么

| 阶段 | Seurat | Scanpy 对应概念 | 解释 |
|---|---|---|---|
| 计数存储 | RNA counts layer | `layers['counts']` | 不可丢失的原始 UMI；AnnData `.raw` 名字不保证内容是 raw |
| QC | nFeature/nCount/PercentageFeatureSet | calculate_qc_metrics | 对异常信号做观察；不同批次分别看 |
| 归一化 | NormalizeData | normalize_total + log1p | 调整总量并压缩动态范围 |
| 选特征 | FindVariableFeatures | highly_variable_genes | 找高于噪声背景的变异特征 |
| 表达空间 | ScaleData + RunPCA | scale + pca | 压缩冗余变化；避免全矩阵无必要稠密化 |
| 图和聚类 | FindNeighbors/FindClusters | neighbors/leiden | 在表达近邻图中找社区 |
| 展示 | RunUMAP | umap | 展示局部结构，不证明谱系 |
| 注释证据 | DotPlot/marker summaries | dotplot/rank_genes_groups | 多标记、多证据核查，不让 cluster 编号充当细胞类型 |

本课程用成熟 Seurat 函数完成分析。默认 Louvain 是本案例的实现选择，非声称其全面优于 Leiden。不同算法、resolution、PC 数都会改变社区，真实细胞类型也不总是离散簇。

## 注释：从“簇 3”到有证据的细胞类型

先看粗谱系，再细分。PBMC 可检查 T 细胞的 CD3D/CD3E，B 细胞的 MS4A1/CD79A，髓系的 LST1/LYZ，NK 相关 NKG7/GNLY；单个 marker 不是唯一判据，缺失一个 marker 也不立即否定身份。

脚本输出 marker 面板图；请自行填写：

| cluster | 候选类型 | 支持 marker | 冲突 marker/证据 | 参考是否匹配 | 最终标签/置信度 |
|---|---|---|---|---|---|
| 自填 | T / B / NK / 髓系 / unknown | 至少两个 | 必须检查 | 物种、组织、平台 | 高/中/低 |

混合 marker 可能是 doublet、ambient RNA、真实过渡状态或低分辨率聚类。完整分析要按真实 capture/channel 评估双细胞；可进一步按 scDblFinder 官方教程学习。过滤后的 PBMC3k 缺少空液滴背景，不能声称已全面纠正 ambient RNA。本课程脚本不自动完成这些扩展步骤。

细胞级 marker 检验描述本数据中的簇差异，不能取代疾病效应检验；在同一矩阵上先聚类再检验还涉及选择性推断。不要将其 p 值解释为跨人群的证据。

## 从单供者探索到多供者推断

PBMC3k 只有一个供者，不能比较健康/疾病。多供者数据应保留 `sample_id`、`donor_id`、`condition`、`batch` 和 `cell_type`。对某一细胞类型 $t$：

$$K_{gdt}=\sum_{c\in(d,t)}Y_{gc}.$$

这就是按样本/供者×细胞类型的 pseudobulk 求和，得到基因×生物学样本的整数矩阵。配对实验先按 `sample_id` 汇总，不把同一供者不同 condition 合并。之后使用 DESeq2/edgeR 等适当模型，并保留供者配对关系。

需要事先规定每个 sample/type 的最低细胞数，报告排除情况。某类型完全没采到，不能机械补零并当作“所有基因都不表达”。重复不足或完全混杂时停止组间检验。也可使用经过设计的细胞级混合模型，但不是把所有细胞直接 t 检验。

细胞比例变化与细胞内表达变化是两个问题。每个供者某细胞类型的比例，分母是纳入分析的哪些细胞？解离、捕获和 QC 都影响比例，不能直接等同于原组织绝对丰度。

## 练习与验收

准备进入有生物学重复的分析时，继续 [多供者 pseudobulk 实验](pseudobulk.md)，使用另一个明确具备配对设计的数据源；不要改造 PBMC3k 来假装重复。

1. 使用 10/20 PCs 或不同 resolution 重跑，记录不稳定的簇，并说明是否影响粗谱系标签。
2. 写出一张注释证据表，至少保留一个可接受的 `unknown` 情形。
3. 设计“6 名患者、6 名对照”的 pseudobulk 样本表；说明为什么不能把每人 500 个细胞写成 n=3000。
4. 列出 PBMC3k 无法支持的三个结论，并说明需要增加什么数据。

验收重点：counts 与变换矩阵分开、QC 可追溯、注释有证据、没有伪造生物学重复。
