# 04 · 空间转录组：位置如何改变统计问题

## 为什么从传统 Visium 入门，再到 HD/Xenium

空间分析同时处理表达矩阵、空间位置、图像和实验设计。先在小型 Visium 对象学习这些对应关系，再进入高分辨率或分子级数据，可减少把软件复杂度误当作科学原理的风险。

主教程：[Squidpy Visium H&E](https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_hne.html) 与 [代码仓库](https://github.com/scverse/squidpy-tutorials)。进阶使用 [SpatialData tutorials](https://github.com/scverse/spatialdata-tutorials)，学习图像、点、形状、表格在坐标系中的关联。

## 数据与目标

调用 `sq.datasets.visium_hne_adata()` 获取 2688 spot×18078 genes 的**预处理**对象。目标是检查组织中的表达空间结构，不能将 `.X` 当原始 counts 再无条件 normalize/log；也不能据单切片推断处理效应。

```bash
python omics-course/scripts/04_fetch_data.py spatial
python omics-course/scripts/06_spatial_visium.py
```

脚本检查坐标维度、有限数值与 observation ID，记录 `.layers`、`.raw` 和处理尺度声明；对给定练习对象调用 Squidpy 的空间邻接和 Moran 分析，不用 UMAP 构建空间图。

输出：`results/spatial/coordinates.csv`、`spatial_expression.png`、`graph_qc.json`、`moran_panel.csv`、`analysis.h5ad` 和版本记录。脚本以按全局方差选出的有限基因面板进行计算并明确记录，默认不是全基因组空间筛选。

## 先理解数据结构

| 类型 | 观测单位 | 核心问题 |
|---|---|---|
| 传统 Visium spot | 捕获点，通常混合多个细胞 | 细胞组成与细胞内表达如何区分？ |
| Visium HD bin | 指定大小的空间 bin | bin 不必等于细胞；细分提高分辨率也增加稀疏性 |
| Xenium/MERFISH 等成像 | 分子位置和分割后的细胞 | panel 限制、背景、细胞边界和转录本归属是否可靠？ |

同一物理位置在原始图像像素、缩略图像素和微米坐标中数值不同。叠图前核对单位、变换、尺度、原点及轴方向；不要为了“看着对齐”随意平移旋转。脚本中的散点图仅使用保存的空间坐标，不冒充经过验证的 H&E 图像叠加。

## 邻接图：先定义什么叫邻居

空间权重矩阵 $W=(w_{ij})$ 指定 spot/细胞之间的邻接。

- Grid：适合规则 spot 排布，通常取最近一圈。
- Radius：固定物理尺度内相邻；半径必须注明单位。
- kNN：每个点固定邻居数量，局部稀疏处可能跨越不应连接的组织间隙。
- Delaunay：几何三角化；需要处理过长边和组织空洞。

表达相似性图与物理邻接图不能替换。多切片必须按 library/section 分开建图，不把重叠坐标的不同切片连在一起。对应 API 与例子见 [空间邻接教程](https://squidpy.readthedocs.io/en/stable/notebooks/examples/graph/compute_spatial_neighbors.html)。

检查孤立点数、节点度、最长边、边在组织图上的位置。图不合理时，后续统计再精确也只是在分析错误的邻接关系。

## Moran's I：空间自相关不等于空间机制

对一个基因在 $n$ 个 spot 上的数值 $x_i$：

$$I=\frac{n}{\sum_{ij}w_{ij}}\frac{\sum_{ij}w_{ij}(x_i-\bar{x})(x_j-\bar{x})}{\sum_i(x_i-\bar{x})^2}.$$

邻居的高低趋势一致通常产生正值；相邻值相反可能为负。随机置换零假设下期望通常为 $-1/(n-1)$，不是恒为 0；具体范围与权重结构有关，不能机械声称 I 永远在 [−1,1]。

置换检验打乱表达和位置的对应关系，与观测统计量比较。$B$ 次置换下使用加一规则的最小可报告 p 约为 $1/(B+1)$；99 次不足以支撑极小尾部概率。本课程默认 999 次用于教学，运行时间取决于机器。检验 100 个预选基因就只在这 100 个内校正；不能把这个 q 值叫全转录组 FDR。

脚本同时保留 Squidpy 返回的检验列。不同版本对置换 p 值的实现可能不同，解释时查 [Moran 示例](https://squidpy.readthedocs.io/en/stable/notebooks/examples/graph/compute_moran.html) 和当前 API。无约束置换检验的是随机空间标签假设；有组织区域、供者层级或其他非交换结构时，需要更合适的限制置换/模型。

组织细胞组成、测序深度和组织边缘都可能造成空间信号。显著 I 只是支持“该图和该尺度下有空间结构”，不证明该基因驱动组织形成。

## 去卷积、邻域富集和通讯

混合 spot 的表达可粗略表示为 $Y_{sg}\approx\sum_t a_{st}R_{tg}$。$R$ 是参考细胞类型表达，$a$ 是待估计贡献；参考不匹配或两类表达相近，会使解不稳定。cell2location/RCTD 等是进阶工具，必须按各自教程使用原始数据尺度、匹配参考和诊断，不把结果必然解释成精确细胞数。

邻域富集统计关注类型 A/B 的邻接是否超过某个随机背景；细胞总数、组织分区和图定义都影响结果。配体—受体共表达/共定位只是候选通讯证据，不能证明实际结合或因果传递。

## 进阶实验：HD 与成像

从 [10x 官方数据](https://www.10xgenomics.com/datasets) 选择一个公开案例，记录下载页面和输出版本。先完成：

1. 分辨 raw molecules、binned counts、segmented cells 和图像分别在哪里。
2. 同一 ROI 比较两种 bin 大小，记录每 bin 分子数、检测基因数和邻接图变化。
3. 成像数据查看背景/阴性对照、分割失败、边界细胞、面积异常与分子归属。
4. 明确 panel 限制：未测基因不能解释为不表达，富集背景不能用全转录组。
5. 多切片保留供者层级；切片是供者内重复时不能直接当独立人群样本。

## 作业与验收

对三个候选基因画空间表达，给出 I、检验族和校正值；换一种邻接或尺度重算，解释不稳定性。写出组织边缘效应、细胞组成、技术质量三种替代解释。需要进一步的病理标注、独立切片/供者或原位验证时明确提出。
