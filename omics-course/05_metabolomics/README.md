# 05 · 代谢组学：从离子峰到有证据的分子解释

## 学习范围与上游教程

本章聚焦非靶向 LC–MS；靶向定量、GC–MS、NMR、脂质组和空间代谢成像有各自的数据特性，不能完全套用。先用历史小数据理解信号处理，再进入现代 QC/MS/MS 工作流。

主线使用 [xcms 官方教程代码](https://github.com/sneumann/xcms/blob/master/vignettes/xcms.Rmd) 与 [Metabonaut 端到端课程](https://rformassspectrometry.github.io/Metabonaut/)。补充阅读 [MetaboAnalystR](https://github.com/xia-lab/MetaboAnalystR) 的统计/通路 vignettes。不要跨不同版本复制函数调用而忽略输入对象类型。

## 数据与问题

[faahKO](https://bioconductor.org/packages/release/data/experiment/html/faahKO.html) 提供 6 个 FAAH knockout 和 6 个 WT 小鼠脊髓的正离子 LC–MS 截取数据。文件只覆盖部分 m/z 和时间范围，不是完整实验，也没有足够的现代 pooled-QC/blank/进样顺序信息供本课程完成全面质控。

问题限定为：“这些数据中哪些离子特征的强度在 KO/WT 间不同？”它与“哪个确定代谢物的浓度变化”“哪个酶导致何种通量变化”不是同一个结论。

```bash
Rscript omics-course/env/install_R.R metabolomics
Rscript omics-course/scripts/05_metabolomics_xcms.R
Rscript omics-course/scripts/07_metabolomics_stats.R
```

本课程 `05` 是官方 OnDisk/MSnbase 路线的轻量运行封装。现代大规模项目优先在 Metabonaut 对应环境中跟随其当前流程。此次本机没有 xcms/faahKO，脚本只做语法检查，未声称得到真实 LC–MS 结果。

## 1. 先看信号，再做峰提取

仪器记录的是 $I(m/z,t)$：质荷比、保留时间和强度。一个 chromatographic peak 是某 m/z 窗口内强度随时间升高再下降的结构，面积近似 $A=\int I(t)\,dt$。

查看 TIC/BPC、单个特征的 EIC、内标和空白，问：有没有仪器漂移、强度突变、峰分裂、明显污染、峰宽变化？只看最后的 PCA 会漏掉信号层面的问题。

`ppm` 描述相对质量偏差：

$$\mathrm{ppm}=10^6\frac{|m_{\mathrm{obs}}-m_{\mathrm{ref}}|}{m_{\mathrm{ref}}}.$$

例如 m/z 500 相差 0.005，为 10 ppm。峰宽参数单位通常是秒；不能把分钟直接填进去。数据为 profile 还是 centroid、分辨率和离子模式都影响算法选择。

## 2. xcms 的处理链

| 步骤 | 解决的问题 | 要检查的输出 |
|---|---|---|
| 数据导入 | 文件与样本一对一对应 | 文件摘要、样本表、RT 范围 |
| 峰检测 | 找每个样本中的信号峰 | 峰数、峰宽、EIC、参数 |
| 初始分组 | 找跨样本相似峰 | m/z 与 RT 的对应关系 |
| RT 对齐 | 校正保留时间漂移 | 对齐前后漂移与标准物质 |
| 重新分组 | 在校正时间上对应特征 | 组内检测比例、错配风险 |
| 补峰积分 | 在预期峰区间重新积分 | 补峰前后缺失矩阵、低信号 |

补峰不是创造观测证据：没有信号的区域积分可能得到噪声。脚本同时输出 `feature_intensity_before_fill.csv` 和 `feature_intensity.csv`，保留“原来没检出、后来积分”的区分。

示例参数来自 faahKO 教学条件，仅用于复现练习。自己的仪器必须以 EIC、内标和 QC 调参，不能认为 `ppm=25` 或 `peakwidth=20–80` 适用于全部 LC–MS。对齐也不能通过扭曲时间把本来不同的化合物硬合并。

## 3. 先分清 feature 与 metabolite

一个 feature 常由 m/z、RT、离子模式定义。一个分子可产生不同加合物、同位素峰和碎片；多个分子也可能共享近似质量。因而“100 个显著 feature”不等于“100 个已鉴定代谢物”。

保存独立的 feature 表：ID、m/z、RT、离子模式、可能的 adduct、MS/MS 谱、库匹配、标准品信息、鉴定置信级别。只有 m/z 匹配时写“候选注释”，不直接替换成确定的化合物名称。

通常需要标准品在相同实验条件下的保留行为与谱图等正交证据，才能给出高置信身份。库谱支持的推测注释、化合物类别、未知 feature 应分开报告；不同规范的等级定义应引用实际采用标准，不混用编号。

## 4. 缺失、归一化与 QC

现代实验至少要了解 pooled QC、blank、内标和进样顺序的作用。QC 的相对标准偏差为 $100s/\bar{x}$，但适用阈值取决于方法与目标；先看空白比例、检测率、漂移和稳定性，而不是机械套一个百分数。

脚本 `07` 为这个小练习采用**完整正值特征**的保守子集：任一生物样本缺失/非正值的特征暂不进入统计，并记录排除原因。这会偏向稳定检出的特征，可能漏掉真实组特异信号，所以必须报告该选择偏倚。不要把该策略当作真实项目的最终缺失处理方案。

后续练习比较：预设检测率过滤后保留 NA 的适当模型、基于检出限的方案、以及其他与缺失机制匹配的方法。不要无条件补零或用半最小值包装成真实观测。填补若依赖组别，会进一步影响组间差异估计。

主分析对保留强度做 log2 转换，并进行中位数中心化作为教学归一化；另输出“不做中位数归一化”的敏感性分析。中位数策略假设大部分特征整体不发生同向偏移，没有内标/完整 QC 时这个假设无法得到充分验证。

## 5. 差异统计与 PCA

对 feature $f$ 的 log2 强度建模：

$$z_{fi}=\beta_{f0}+\beta_{f1}\mathrm{KO}_i+\epsilon_{fi}.$$

本课程调用 [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) 的线性模型与经验贝叶斯方差估计，输出 logFC、置信区间、原始 p 和 BH 调整值。这里是连续强度模型，不需要 voom，也不使用 DESeq2 的整数计数模型。

PCA 只用于探索，脚本使用按特征标准化的 log 强度并报告处理方式。PCA 分开不意味着因果差异，PCA 重叠也不表示所有特征都无差异。主表保留全部进入检验的 feature，并回连 m/z/RT。

PLS-DA/OPLS-DA 能在噪声里制造看似分开的投影。若做预测，需要供者分组的嵌套验证、训练折内预处理、标签置换和独立验证集。R² 高、VIP 大或训练集分类完美都不是足够证据。

## 6. 机制与通路解释

先解决化合物身份和映射歧义，再做通路富集；背景是本平台可测且可映射的对象。一个代谢物参与多个通路、多个 feature 属于同一分子，重复计数会夸大证据。

稳态丰度不等于代谢通量。简化动态关系 $dM/dt=v_{\mathrm{in}}-v_{\mathrm{out}}$ 表明：同样的池大小可以对应完全不同的合成/消耗速度。RNA 上调、蛋白增加、代谢物积累之间也没有必然的一一对应。需要机制结论时考虑稳定同位素示踪、酶活性、标准品及扰动实验。

## 作业与验收

选 3 个显著 feature：回看 EIC、缺失情况、组内分布、归一化敏感性和鉴定证据。分别写“统计上观察到什么”“身份有多确定”“尚不能说明什么”。不得只提交火山图或通路图。
