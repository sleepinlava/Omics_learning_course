# 已生成结果：建议怎样阅读

> 这些文件由本地脚本生成，不随 Git 分发。新克隆仓库需先运行对应章节，下面路径才会出现；已运行状态以 `../VALIDATION.md` 为准。

这些是本机实际运行得到的文件。空间和代谢组的真实结果目录尚未生成，不能用下列结果代替它们。

## 1. 先看统计直觉

- [Poisson 与负二项分布图](statistics/count_distributions.png)：均值相近时，离散程度仍可不同。
- [伪重复模拟图](statistics/pseudoreplication.png)：比较忽略供者与供者级检验的零假设 p 值分布。
- [假阳性比例表](statistics/pseudoreplication_summary.csv)：这是特定模拟模型的表现，不是所有真实数据的固定错误率。
- [缺失机制敏感性](statistics/missingness_sensitivity.csv)：同一模拟真值，删失/填补策略可改变方向和 p 值。
- [模拟环境](statistics/sessionInfo.txt)。

## 2. 再看 airway 配对案例

- [PCA](bulk/PCA.png)：处理沿主要变化方向分布，细胞系也贡献变异；不能仅凭此图宣布机制。
- [样本距离](bulk/sample_distances.png)：检查异常样本与配对结构。
- [MA 图](bulk/MA.png)：低表达端的效应估计尤其需要谨慎；图中是未收缩 LFC。
- [完整基因表](bulk/all_genes.csv)：包括未通过预过滤的基因及状态。
- [显著基因表](bulk/significant_genes.csv)：本次 4081 行，条件为 padj<0.05；不是独立验证后的机制名单。
- [样本 QC](bulk/sample_qc.csv)、[设计矩阵](bulk/design_matrix.csv)、[输入摘要](bulk/input_files.csv)、[会话记录](bulk/sessionInfo.txt)。

先按 ID 选择一个基因，再查看 `normalized_counts.csv` 中的所有 8 个样本，按细胞系配对；不要只看组均值或 p 值最小的一行。

## 3. 小型单细胞演示

- [PCA](single_cell_demo/PCA.png)、[marker 面板](single_cell_demo/marker_panel.png)。
- [细胞元数据](single_cell_demo/cell_metadata.csv)、[QC 记录](single_cell_demo/qc_cells.csv)、[会话记录](single_cell_demo/sessionInfo.txt)。

这是 230 特征×80 细胞的缩小对象。它验证软件流程，不替代完整 PBMC3k，不用于制定线粒体过滤阈值或声明跨供者差异。marker 面板仅帮助练习读图。
