# 明确标注的模拟输入

`meta_demo_effects.csv` 手工指定 3 个虚构特征、每个 6 个虚构独立研究。三组分别模拟方向一致、明显异质、接近零的效应。标准误依次设为 0.20、0.24、0.18、0.30、0.22、0.26；没有抽样数据、真实受试者或 GEO accession。

`independence_reviewed=TRUE` 只描述这个模拟设定，不表示审核过真实队列。运行 `Rscript omics-course/scripts/14_meta_analysis.R --demo` 需要 metafor。观察随机效应区间、预测区间和逐研究剔除；不能据此宣称任何基因或疾病结论。
