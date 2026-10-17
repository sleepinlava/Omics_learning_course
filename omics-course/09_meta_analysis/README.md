# 09 · Meta 分析：合并独立研究，而不是堆叠表达矩阵

## 范围与上游教程

本章采用 [metafor 官方文档](https://wviechtb.github.io/metafor/reference/rma.uni.html) 和 [原作者代码](https://github.com/wviechtb/metafor)，完成研究级效应量合并。它不是把不同组学、物种或不同处理随便混合，也不等同于单细胞批次整合。

先完成 GEO/研究内分析，保存每项研究的可比效应 $y_i$ 和采样标准误 $s_i$。相同数据被不同论文/GSE 重复提交时，仍是同一队列，不能算两个独立研究。现有 airway 与 pasilla 是不同物种和不同扰动，禁止作为同一生物学问题的 meta 输入。

## 运行入口

```bash
Rscript omics-course/env/install_R.R meta
# 明确标记为模拟的 3 个 feature × 6 个独立研究
Rscript omics-course/scripts/14_meta_analysis.R --demo
# 审核后替换为真实研究级输入
Rscript omics-course/scripts/14_meta_analysis.R effects.csv
```

模拟表在 `examples/meta_demo_effects.csv`，只用于教学和检查函数，不代表 GEO 真实结果。真实表的格式参照 `templates/meta_effects.csv`，研究筛选记录参照 `templates/meta_studies.csv`。

## 真实数据的入门复现：BCG 文献研究

[metadat 的 dat.bcg](https://wviechtb.github.io/metadat/reference/dat.bcg.html) 提供 13 项 BCG 疫苗研究的汇总事件数。它适合学习效应和抽样方差；这是历史文献方法案例，不是组学数据，也不是本课程的临床建议。包安装后，在仓库根目录的 R 会话中运行：

```r
.libPaths(c("omics-course/env/R-library", .libPaths()))
library(metafor)
data("dat.bcg", package = "metadat")
x <- escalc(measure = "RR", ai = tpos, bi = tneg,
            ci = cpos, di = cneg, data = dat.bcg)
effects <- data.frame(
  feature_id = "TB_incidence", study_id = paste0("BCG_", x$trial),
  cohort_id = paste0("BCG_", x$trial), effect = x$yi, se = sqrt(x$vi),
  effect_scale = "logRR", contrast = "vaccinated_vs_control",
  organism = "human", independence_reviewed = FALSE)
dir.create("omics-course/data/meta_bcg", recursive = TRUE, showWarnings = FALSE)
write.csv(effects, "omics-course/data/meta_bcg/effects.csv", row.names = FALSE)
```

先阅读数据页与原分析，记录采用独立研究假设的依据，再把确认过的 `independence_reviewed` 改为 TRUE，传给脚本 `14`。不能只为绕过检查修改字段。原教程默认推断与本课程 Hartung–Knapp 区间可能不同；区分方法差异和错误。本例在当前环境尚未运行，原因是缺少 metafor。

## 系统性检索与纳入规则

在看到效应方向之前写出：研究问题、组织/物种/人群、暴露/处理、比较组、终点、时间范围、技术平台、纳入排除条件、文献与 GEO 检索语句。保存每个排除理由，核对供者重叠，尽量获取完整效应表而不只收录“显著基因”。

同一个研究同一个基因只贡献一个适当效应；多个时间点、多个相关对照、重叠供者需要协方差或多层模型，本简化脚本会拒绝明显的重复研究/队列。`independence_reviewed=TRUE` 表示人为核对过独立性，不是算法推断的事实。

## 固定效应与随机效应

共同效应模型假设所有研究估计同一个真值，权重 $w_i=1/s_i^2$，合并值 $\hat\mu=\sum_iw_iy_i/\sum_iw_i$。

随机效应模型允许研究真效应不同：

$$y_i\sim N(\theta_i,s_i^2),\quad\theta_i\sim N(\mu,\tau^2),\quad
w_i=1/(s_i^2+\tau^2).$$

主分析用 REML 估计 $\tau^2$，采用 Hartung–Knapp 推断，另输出共同效应模型作为敏感性比较。随机效应模型不会自动消除偏倚、平台差异或不合理的纳入标准。

## 异质性、不确定性与逐研究剔除

- Q 衡量相对抽样误差的研究间不一致；研究少时检验能力有限。
- $\tau^2$ 是真效应间的方差，有实际效应尺度。
- $I^2$ 是相对异质性指标，受研究精度影响，不是“百分之多少研究是假的”。软件给出的模型估计与简单 Q 公式可能不同。
- 合并效应的置信区间描述平均效应的不确定性；预测区间描述模型下新的可比研究真效应范围，往往更宽。两者不能互换。
- [leave-one-out](https://wviechtb.github.io/metafor/reference/leave1out.html) 检查结论是否被单项研究主导；不能因为删掉某研究更显著就事后排除它。

脚本至少要求每 feature 有 3 项独立研究，这仅是教学下限，不是功效或异质性可靠性保证。少研究时应重点呈现宽区间与不确定性；不自动运行 Egger 检验或以漏斗图“证明无发表偏倚”。

## 如何得到可以合并的 effect 与 SE

| 研究内结果 | 候选效应 | 注意事项 |
|---|---|---|
| 同方向 RNA-seq 对比 | 未收缩 log2FC 与对应 lfcSE | 基因 ID/协变量/对比应可比；不要把收缩 LFC 配未收缩 SE |
| log2 芯片表达 | limma 的 logFC 与预定 SE 定义 | probe→gene 规则预先固定；本课程另导出 SE_classical |
| 连续指标且量纲不同 | 适当的标准化均差及方差 | 使用成熟 `escalc` 计算，不能直接与 log2FC 混合 |
| 二分类结局 | logRR/logOR 与 SE | log 是自然对数；RR 与 OR 不混用，不与组学 log2FC 混用 |

显著性 p 值不是效应大小，padj 不能直接转成 SE。Fisher/Stouffer 合并 p 值回答另一个问题，通常不能保留效应尺度与全部方向信息；本模块明确选择效应量合并。

## 输出和多重检验

`results/meta_analysis/`（模拟为 `meta_demo/`）包含完整合并表、k、Q、I²、tau²、置信/预测区间、逐研究剔除结果和第一项 feature 的森林图。跨基因 BH 家族包括所有提交 feature，拟合失败/研究不足另记，不通过丢掉失败基因缩小校正家族。基因间共享受试者是跨 feature 的相关性，不等同于同一 feature 内重复研究；仍需解释 FDR 条件与选择偏倚。

## 作业

为三个独立、同物种同对比的研究建立筛选表；展示一个方向一致、一个明显异质、一个证据不足的 feature。解释为什么森林图合并点不等于临床/生物学重要性。将可比分析后的完整排序送入后续通路分析时，重新声明排序规则和背景；不能拼接每项研究的显著列表后声称完成 meta 分析。
