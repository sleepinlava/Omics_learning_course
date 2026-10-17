# 工作区整理与发布范围

## 本次整理

根 README 改为多组学总入口；原文逐字保存在 `docs/legacy-bulk-guide.md`，原清单保存在 `docs/legacy-manifest.txt`。它们是历史记录，内部 `~/rnaseq-learning` 路径和运行状态不代表当前 checkout。已有 `AGENTS.md` 不作修改。

原 bulk 脚本、数据、结果、上游教程保持路径，避免破坏已有分析。`scripts/00_download.sh` 的 ROOT 已改为从脚本位置定位仓库；AGENTS 中旧的硬编码提示保留为历史指导。下载脚本会处理第三方目录，运行前阅读其行为。

## 新内容放在哪里

中文课程放入 `omics-course/00_setup` 至 `09_meta_analysis`，可执行入口按学习顺序编号；课程数据与输出分别写到其 `data/` 和 `results/`。GEO 原下载矩阵与人工核对样本表分开存放，不覆盖下载元数据。KEGG 缓存按物种和日期保存。meta 输入每行是一项研究的一个特征效应，并注明队列独立性。

新研究建议在本地建立 `omics-course/data/<project>/`，按 `raw/`、`metadata/`、`derived/` 分类；报告列出输入 SHA256/版本和完整命令。参数变化前归档旧输出，不在原始矩阵内直接修改数值。

## Git 管理

`.gitignore` 排除下载文件、生成结果、本地依赖、嵌套上游仓库和凭据目录；忽略不会删除本地文件。Git 只保留课程、脚本、模板、模拟输入、来源说明及 `data/counts/` 的小型练习数据。根 `MANIFEST.txt` 与课程 `resources/FILE_MANIFEST.tsv` 是发布源码清单，不是声称所有数据已下载。

新克隆的 `results/` 为空是正常现象。先运行章节示例，再打开结果导读里的文件。完整上游教程通过 `08_get_tutorials.py fetch <name>` 获取，保留原 LICENSE/CITATION。
