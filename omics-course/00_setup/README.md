# 00 · 环境、文件和数据准备

## 学习目标

能分清终端命令与 R/Python 代码，能定位一个输入文件，能说明矩阵每行每列代表什么，能在新环境复现同一个分析。

## 开始前的最小知识

终端的 `pwd` 显示当前位置，`ls` 查看文件，`cd` 切换目录。`Rscript file.R` 执行 R 脚本，`python file.py` 执行 Python 脚本。不要把 Markdown 的三个反引号复制到终端。路径含空格时加引号。

R 中 `x <- ...` 是赋值；`dim(x)` 看维度；`head(x)` 看前几行；`str(x)` 看类型。`x[行, 列]` 的顺序不能反过来。常见 bulk 表是基因×样本，AnnData 通常是细胞/spot×基因。CSV 第一列经常是行名，读入时务必检查。

先阅读现有 [仓库 README](../../README.md)；原 bulk 下载脚本 `scripts/00_download.sh` 已按脚本位置定位仓库；新课程请使用分模块下载入口。

## 环境分开建立

### R：bulk 与单细胞

```bash
Rscript omics-course/env/install_R.R core
```

该脚本使用课程内 `env/R-library/`，不会主动更新系统包。首次安装需要联网；Linux 编译依赖可能需要由系统管理员提供。R 与 Bioconductor 必须使用官方支持的对应版本，不混装 release/devel。安装后 `BiocManager::valid()` 检查兼容性。`sessionInfo.txt` 才是实际运行版本记录；安装脚本不是完整依赖锁。

### Python：空间转录组

```bash
python3.11 -m venv omics-course/env/spatial-venv
source omics-course/env/spatial-venv/bin/activate
python -m pip install -r omics-course/env/spatial-requirements.txt
python -m pip freeze > omics-course/env/spatial-resolved.txt
```

建议 Python 3.11 或 3.12；若安装的是 3.12，把命令中的 `python3.11` 换成 `python3.12`。当前终端的默认 Python 是 3.15，不应直接用它建立这组科学计算依赖环境；没有对应解释器时先按自己的系统安装受支持版本。版本范围是兼容目标，当前机器尚未实际安装验证。先做小数据，遇到依赖冲突保留报错，回到官方安装说明，不随机升级所有包。

### R：代谢组

```bash
Rscript omics-course/env/install_R.R metabolomics
```

xcms 的原始数据处理依赖较多。替代路线：按 [Metabonaut](https://rformassspectrometry.github.io/Metabonaut/) 的安装/Workshop 入口使用其预装环境，完成后下载自己的结果。上游环境版本与本课程脚本不同，仍需记录会话信息。

## 数据获取：从现成计数开始，逐步回溯原始数据

| 数据 | 准备方式 | 用途与边界 |
|---|---|---|
| airway | 仓库 `data/counts/airway_*.csv` 已存在 | 8 个文库、4 个细胞系的配对处理实验 |
| pasilla | 仓库 `data/counts/pasilla_*` 已存在 | 方法比较；建模前检查 single/paired 文库类型 |
| PBMC small | `03_single_cell.R --demo` 从 SeuratObject 内置数据导出 | 80 细胞的缩小对象，仅用于流程连通性 |
| PBMC3k | 下方下载命令 | 官方 10x 计数矩阵，单供者细胞类型练习 |
| Visium H&E | 下方 Squidpy 获取命令 | 已预处理对象；不冒充原始 counts |
| faahKO | 安装数据包，运行脚本 `05` | 12 个截取的 NetCDF 文件，KO/WT 练习 |

完整来源与分析单位见 [数据目录](../resources/DATASETS.md)。不要将教学小参考基因组用于真实全基因组实验。

```bash
python omics-course/scripts/04_fetch_data.py pbmc3k
Rscript omics-course/scripts/03_single_cell.R --qc-only
# 先审阅 QC，再复现该历史 PBMC 案例的阈值（不是通用阈值）
Rscript omics-course/scripts/03_single_cell.R --min-features 200 --max-features 2500 --max-mt 5

python omics-course/scripts/04_fetch_data.py spatial
python omics-course/scripts/06_spatial_visium.py

Rscript omics-course/scripts/05_metabolomics_xcms.R
Rscript omics-course/scripts/07_metabolomics_stats.R
```

下载工具失败时明确退出，不生成假数据。下载成功会记录来源、大小和 SHA256；这个摘要用于后续完整性比较，不是数据提供方签名。只有数据包的版本/文件摘要记录齐全，才有可追溯的输入。

## 上游教程如何保存

```bash
# 默认仅列出配置；显式指定 fetch 才下载
python omics-course/scripts/08_get_tutorials.py
python omics-course/scripts/08_get_tutorials.py fetch scanpy
python omics-course/scripts/08_get_tutorials.py fetch metabonaut
```

工具下载到 `resources/upstream/`，保留原仓库 LICENSE，并记录提交号。已有目录会跳过，不覆盖你的笔记。下载的是当时 HEAD，而不是伪称已固定的版本；得到 `upstream_receipts.json` 后应连同分析报告归档。

## 算力规划与问题排查

入门矩阵案例可先预留 8–16 GB RAM；空间带图像或原始 LC–MS 建议从 16 GB 起评估。HD 全切片可能远超这些配置，先选 ROI/较粗 bin。这里是规划建议，不是已测量的最低要求。下载前查官方文件大小；不要默认获取所有 FASTQ、图像和中间 BAM。

- `No module named ...`：确认激活了课程 venv，再安装对应模块依赖。
- R 提示缺包：使用课程安装脚本；脚本不会联网暗装依赖。
- 下载返回网页而非压缩包：检查代理、URL 和下载错误，不能改后缀后强行读取。
- 样本名不匹配：先修复元数据的一对一映射，不依靠位置拼接。
- 图存在但结果异常：先查输入尺度、对照水平、QC 和日志，不能只重画图。

作业：独立写出一个样本表、一个软件版本清单和一条完整分析命令。让另一人只根据它们找到对应输入。

## GEO、KEGG 与 meta 扩展依赖

```bash
Rscript omics-course/env/install_R.R geo
Rscript omics-course/env/install_R.R kegg
Rscript omics-course/env/install_R.R meta
Rscript omics-course/scripts/14_meta_analysis.R --demo
```

meta 使用 CRAN 的 metafor；其他模块使用对应 Bioconductor 依赖。只安装要学的模块；网络失败时不要把“无法访问包索引”误判为 R 版本一定不兼容。
