# 上游教程、代码与引用

核对日期：2026-10-17。下面是阅读与复用映射，不是本地下载完成声明。网页检索可访问的资源与本机终端下载能力是两件事。以官方教程/原作者代码为主；不按热度排名、不镜像整本书。

| 模块 | 首选教程/代码 | 本课程如何使用 |
|---|---|---|
| bulk | [HBC Intro-to-DGE](https://github.com/hbctraining/Intro-to-DGE) | 从实验设计、计数到差异分析的阅读路线；本仓库已带课程副本 |
| bulk | [DESeq2 官方文档](https://bioconductor.org/packages/DESeq2/)；[方法论文](https://doi.org/10.1186/s13059-014-0550-8) | 调用 DESeq2；补充 NB、设计矩阵和效应解释 |
| 单细胞 R | [Seurat PBMC3k](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)；[代码仓库](https://github.com/satijalab/seurat) | 复用成熟对象与标准函数序列；新增本地路径、QC 导出与限制说明 |
| 单细胞 Python | [Scanpy 教程仓库](https://github.com/scverse/scanpy-tutorials)；[教程索引](https://scanpy.readthedocs.io/en/stable/tutorials/index.html) | R/Python 对照阅读；PBMC3k 页面已被官方标为 legacy，不将它说成最新最佳流程 |
| 单细胞统计 | [OSCA](https://bioconductor.org/books/release/OSCA/)；[Best Practices](https://www.sc-best-practices.org/)；[代码](https://github.com/theislab/single-cell-best-practices) | QC、注释、重复结构与下游推断的进阶路线 |
| 多样本单细胞 | [HBC Pseudobulk](https://github.com/hbctraining/Pseudobulk-for-scRNAseq)；[自学路线](https://github.com/hbctraining/Pseudobulk-for-scRNAseq/blob/main/schedule/self-learning.md)；[muscData 数据文档](https://bioconductor.org/packages/release/data/experiment/vignettes/muscData/inst/doc/muscData.html) | 供者级汇总、配对分析与 GSE96583 的进一步练习 |
| 单细胞推断 | [Squair 等，2021](https://www.nature.com/articles/s41467-021-25960-2)；[Zimmerman 等，2021](https://www.nature.com/articles/s41467-021-21038-1) | 区分细胞内相关、pseudobulk 与混合模型，解释伪重复 |
| 空间 | [Squidpy Visium H&E](https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_hne.html)；[教程代码](https://github.com/scverse/squidpy-tutorials) | 读取预处理练习对象，调用空间图与 Moran 分析 |
| 空间统计 | [空间邻接图](https://squidpy.readthedocs.io/en/stable/notebooks/examples/graph/compute_spatial_neighbors.html)；[Moran 示例](https://squidpy.readthedocs.io/en/stable/notebooks/examples/graph/compute_moran.html) | 邻接定义与检验族的解释，参数敏感性练习 |
| HD/成像 | [SpatialData tutorials](https://github.com/scverse/spatialdata-tutorials)；[10x 数据入口](https://www.10xgenomics.com/datasets) | 坐标系统、图像、分割与分辨率进阶；按具体数据页面保存技术版本 |
| 代谢组原始处理 | [xcms](https://bioconductor.org/packages/xcms/)；[vignette 源代码](https://github.com/sneumann/xcms/blob/master/vignettes/xcms.Rmd) | 通过包函数实现峰检测/对齐/分组/补峰；参数仅面向 faahKO 练习 |
| 现代代谢组 | [Metabonaut](https://rformassspectrometry.github.io/Metabonaut/)；[代码](https://github.com/rformassspectrometry/Metabonaut) | 进阶主线：原始信号到 QC、统计及 MS/MS 注释 |
| 代谢组补充 | [MetaboAnalystR](https://github.com/xia-lab/MetaboAnalystR) | 可视化/统计/通路分析补充；使用其自带 vignette，不混抄不同版本 API |
| 数据库 | [MetaboLights](https://www.ebi.ac.uk/metabolights/)；[EBI 数据访问课程](https://www.ebi.ac.uk/training/online/courses/chemical-biology-collection/explore-biological-interactions/metabolights/) | 找到原始数据、样本表与实验协议；不得只下载强度表丢弃实验设计 |

## 复用与许可

中文原理说明、练习、路径封装和验证记录是本课程新增内容。算法实现由所调用的软件包提供；本课程不是这些包的替代实现。上游教程的全文、图和代码保持原作者许可。脚本注释标明主要参考入口，下载上游仓库时保留 LICENSE/CITATION。不要将第三方数据或代码统一声明为本课程所有；发布或分发数据前查看相应数据页的条款。

正式报告引用实际使用的软件与数据论文：在 R 中运行 `citation("DESeq2")`、`citation("Seurat")`、`citation("xcms")`、`citation("faahKO")`、`citation("limma")`；Python 工具按其文档提供的论文引用。

## 已核对但需要注意的变化

- Squidpy 教程仓库旧名 `squidpy_notebooks` 已重定向到 `squidpy-tutorials`；SpatialData 的 `spatialdata-notebooks` 也重定向到 `spatialdata-tutorials`。
- Squidpy 的 `visium_hne_adata` 是预处理 AnnData，不能把 `.X` 自动当作原始 UMI。
- 当前 [grid 邻接 API](https://squidpy.readthedocs.io/en/stable/api/squidpy.gr.spatial_neighbors_grid.html) 为 `spatial_neighbors_grid`；课程采用它，避开计划移除的通用 `spatial_neighbors` 入口。
- xcms 新旧容器 API 并存。课程 `05` 使用官方示例中的 MSnbase/OnDisk 路线；现代大数据分析转向 Metabonaut 的当前工作流，不混用对象方法。
- 网页核对不等于每个 URL 的二进制下载已验证。细节与失败状态见 `../VALIDATION.md`。

## GEO、KEGG、meta 扩展来源

- [GEOquery 官方教程](https://bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.html)：复用 series matrix、多平台拆分、supplement 列表与 RNA-seq quantification 获取接口；不自动猜测分组或量纲。
- [KEGG REST API](https://www.kegg.jp/kegg/rest/keggapi.html)、[KEGGREST](https://bioconductor.org/packages/release/bioc/vignettes/KEGGREST/inst/doc/KEGGREST-vignette.html)、[clusterProfiler 源码](https://github.com/YuLab-SMU/clusterProfiler)：数据库快照和成熟富集函数；课程补充测试背景、ID 丢失与方向解释。
- [metafor rma.uni](https://wviechtb.github.io/metafor/reference/rma.uni.html)、[leave1out](https://wviechtb.github.io/metafor/reference/leave1out.html)、[predict](https://wviechtb.github.io/metafor/reference/predict.rma.html)、[作者代码](https://github.com/wviechtb/metafor)：采用 REML、Hartung–Knapp、预测区间与逐研究剔除；不重新实现估计器。

报告补充 `citation("GEOquery")`、`citation("clusterProfiler")`、`citation("metafor")`；记录 KEGG 快照日期。
