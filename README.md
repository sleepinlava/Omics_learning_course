# ximcs · 组学分析中文教学工作区

从可直接使用的练习数据出发，学习 bulk RNA-seq、单细胞、空间转录组、代谢组，以及 GEO 数据再分析、KEGG 富集和跨研究 meta 分析。每章结合现有开源教程、数学原理、执行步骤、生物学解释与练习。

## 学习入口

从 [完整课程与学习路线](omics-course/README.md) 开始；先阅读 [环境准备](omics-course/00_setup/README.md) 和 [数学与生物统计](omics-course/01_foundations/README.md)。

| 分析任务 | 教程入口 |
|---|---|
| 配对 bulk 差异表达 | [DESeq2 与 airway](omics-course/02_bulk_rnaseq/README.md) |
| 单细胞及多供者推断 | [scRNA-seq](omics-course/03_single_cell/README.md) · [pseudobulk](omics-course/03_single_cell/pseudobulk.md) |
| 空间坐标与组织模式 | [空间转录组](omics-course/04_spatial/README.md) |
| LC–MS 峰处理与差异特征 | [代谢组学](omics-course/05_metabolomics/README.md) |
| 公共研究下载与再分析 | [GEO](omics-course/07_geo/README.md) |
| 通路过度代表与排序富集 | [KEGG ORA/GSEA](omics-course/08_kegg/README.md) |
| 多研究效应合并与异质性 | [meta 分析](omics-course/09_meta_analysis/README.md) |

## 最小运行示例

在仓库根目录执行，需要 R；安装包和下载完整数据需要联网。

```bash
python omics-course/scripts/00_check.py
Rscript omics-course/scripts/01_statistics_lab.R
Rscript omics-course/env/install_R.R core
Rscript omics-course/scripts/02_bulk_airway.R
Rscript omics-course/scripts/03_single_cell.R --demo
```

GEO、KEGG 和 meta 模块分别安装 `geo`、`kegg`、`meta` 依赖，完整命令见对应章节。meta 演示输入明确为模拟效应，不能当作真实研究结论。下载失败会报错；查阅 [实际验证记录](omics-course/VALIDATION.md) 区分已运行案例与待联网步骤。

## 目录与复现

- `omics-course/`：主课程、脚本、模板、来源清单与模拟示例。
- `data/counts/`：已有 airway、pasilla 小型计数与样本表；[来源说明](data/counts/README.md)。
- `scripts/`：保留原 bulk 课程 `00`–`07` 脚本。
- `data/`、`results/`、`tutorials/`：本地下载、生成结果、上游教程；大文件不进入 Git。
- `docs/`：目录整理规则、历史 bulk 指南与原清单。

原始数据、历史输出和上游教程保留原位置。Git 发布包含课程源码与小型练习输入，不包含本地安装环境、FASTQ、完整第三方仓库或生成的 RDS。详见 [目录管理](docs/WORKSPACE.md)、[文件清单](MANIFEST.txt) 和 [第三方来源](THIRD_PARTY.md)。

## 贡献与检查

遵循 [AGENTS.md](AGENTS.md)。新增分析需说明数据来源、样本单位、设计矩阵、效应方向、验证方法与局限。可运行 `python -m unittest discover -s omics-course/tests -v` 检查 GEO 输入契约和 meta 输入保护；它不替代真实数据端到端验证。
