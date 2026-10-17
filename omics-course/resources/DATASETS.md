# 数据目录：来源、实验单位与实际用途

## 入门必做数据

| ID | 来源与获取 | 数据状态/实验单位 | 能回答什么 | 本课程位置 |
|---|---|---|---|---|
| airway / GSE52778 | [Bioconductor airway](https://bioconductor.org/packages/airway/)；仓库已有导出 CSV | 人气道平滑肌，4 个细胞系各处理/未处理；8 个文库不是 8 个独立供者 | 配对设计下的处理相关表达变化 | `../../data/counts/airway_counts.csv`、`airway_coldata.csv` |
| pasilla | [Bioconductor pasilla](https://bioconductor.org/packages/pasilla/)；已有 TSV/CSV | 果蝇，condition 和文库类型需检查 | 比较 edgeR、limma-voom、DESeq2 | `../../data/counts/pasilla_*` |
| PBMC small | SeuratObject 自带 `pbmc_small` | 80 细胞、特征已缩减的教学对象 | 对象结构、函数与图是否能运行；不能评估完整 QC | `03_single_cell.R --demo` 导出到课程 `data/pbmc_small/` |
| PBMC3k | [官方压缩包](https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz)；[原教程](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) | 约 2700 个 PBMC；单供者，hg19，过滤后的 UMI 矩阵 | 细胞类型探索；无健康/疾病独立组别 | `data/pbmc3k/filtered_gene_bc_matrices/hg19/` |
| Kang18_8vs8 / GSE96583（扩展） | [muscData](https://bioconductor.org/packages/release/data/experiment/vignettes/muscData/inst/doc/muscData.html)：`muscData::Kang18_8vs8()` | 8 位供者的对照/IFN-β 刺激配对；实际字段和 QC 需读取后核实 | 同一细胞类型内的配对处理效应，不是患者/健康比较 | 见 `03_single_cell/pseudobulk.md`，本次未下载 |
| Visium H&E | [官方 API](https://squidpy.readthedocs.io/en/stable/api/squidpy.datasets.visium_hne_adata.html)：`sq.datasets.visium_hne_adata(path=...)` | 小鼠脑冠状切片，2688 spot × 18078 genes 的预处理对象/坐标；高分辨率图像另由 `visium_hne_image()` 获取 | 单切片空间模式；不是供者层面的处理比较 | `data/spatial/visium_hne.h5ad` |
| faahKO | [Bioconductor 数据包](https://bioconductor.org/packages/release/data/experiment/html/faahKO.html)：`BiocManager::install("faahKO")` | 6 KO 与 6 WT 小鼠脊髓；正离子模式 NetCDF；m/z 200–600、RT 2500–4500 秒子集 | LC–MS 预处理与特征差异；不自动得到确定化合物身份 | 数据包 `cdf/`；`05` 导出位置和 MD5 清单 |

样本量来自对应教程/数据说明；下载后必须再核对维度、样本表和数据版本。不要用文件名自行推断其他数据集的组别。本课程 faahKO 文件夹 KO/WT 的映射仅适用于该官方包，并在脚本中检查。

## 进阶数据的选择规则

1. **空间 HD / Xenium**：从 [10x Datasets](https://www.10xgenomics.com/datasets) 或 [Xenium Explorer 演示入口](https://www.10xgenomics.com/support/software/xenium-explorer/latest/resources/xenium-explorer-demos) 选一个带完整输出和图像的数据集。记录平台版本、panel、组织处理、分割算法、bin 大小、坐标单位和样本数。先下载处理后矩阵与必要图像，再决定是否获取 FASTQ/分子坐标。
2. **供者层面单细胞比较**：从 OSCA/Best Practices 的多样本案例继续，要求显式 `donor_id`、`sample_id`、`condition`，并证明其关系。PBMC3k 不可通过随机分组“制造”健康与疾病供者。
3. **真实代谢组重分析**：在 MetaboLights 按物种、组织和 LC–MS 检索；保存 accession、ISA 样本表、原始文件清单、QC/blank、离子模式、色谱方法和 MS/MS。只有符合设计的项目才进入差异比较，不因为能下载就适合统计检验。
4. **现代完整 LC–MS 教程数据**：用 [Metabonaut](https://rformassspectrometry.github.io/Metabonaut/) 的数据入口和端到端 vignette，保持其数据版本与运行环境一致。

## 下载与复现记录

`04_fetch_data.py` 将实际下载字节数、来源和 SHA256 写入 `data/download_receipts.json`。文件尚未下载时没有摘要，不凭空填写大小或哈希。现有 bulk 输入的真实摘要可由 `00_check.py` 生成。原始文件保持不变；下游选择/排除写入新的表。

入门案例避免把大文件捆绑到课程中。联网终端可按 setup 中的逐模块命令获取；当前受限终端不能解析外网域名，因此本次不声称 PBMC3k、Visium 或 faahKO 原始文件已全部到位。

## 新增数据入口

| 数据 | 获取方式 | 用途和边界 |
|---|---|---|
| GSE2553 | `Rscript omics-course/scripts/10_geo_fetch.R GSE2553 matrix` | GEOquery 教程公开 series matrix 案例；平台/分组/尺度必须人工核对，不自动当作病例对照 |
| GSE52778 | `Rscript omics-course/scripts/10_geo_fetch.R GSE52778 supp-list` | 查看 airway 研究补充文件；FPKM 不是整数 counts |
| GSE164073 | `Rscript omics-course/scripts/10_geo_fetch.R GSE164073 counts` | GEOquery 官方 RNA-seq quantification 示例；先检查返回 assay 和注释 |
| KEGG hsa | `Rscript omics-course/scripts/12_kegg_cache.R hsa` | 在线获取物种通路映射并保存日期快照，数据库不是独立实验数据 |
| meta 模拟效应 | `examples/meta_demo_effects.csv` 已随源码提供 | 3 个虚构特征×6 项虚构研究，只练习模型与图 |

前三个 GEO accession 是不同获取模式的示例，**不能因为可以下载就合并为 meta 分析**。真实 meta 项目先完成 `templates/meta_studies.csv` 的系统筛选和队列重叠审阅，再填效应表。当前未声称下载了上述新增 GEO 数据或 KEGG 快照。

真实 meta 方法练习：[metadat::dat.bcg](https://wviechtb.github.io/metadat/reference/dat.bcg.html)，13 项历史疫苗研究的汇总数据，安装 metafor/metadat 后可直接载入；导出与审核步骤见第 09 章。此例不与基因 log2FC 合并。
