# 已有小型教学输入

- `airway_counts.csv`、`airway_coldata.csv`：从 Bioconductor airway 数据对象导出的 8 文库、4 细胞系配对实验整数计数与样本信息；[airway 数据包](https://bioconductor.org/packages/airway/)、[GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778)。这不是 GEO 的 FPKM 矩阵，不能混用。
- `pasilla_gene_counts.tsv`、`pasilla_sample_annotation.csv`：pasilla 练习计数与文库注释；[pasilla 数据包](https://bioconductor.org/packages/pasilla/)。分析需核对文库类型与分组。

这些文件在本次整理之前已经存在；没有把当前时间伪记为原始下载时间。使用时按数据包和原论文引用，保留样本说明；可用 `scripts/00_download.sh` 检查既有获取流程。文件摘要见根 MANIFEST.txt。
