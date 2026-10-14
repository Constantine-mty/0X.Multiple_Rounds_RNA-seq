# GO BP 分类与 GSEA Summary

新代码统一放在 `gobp_pipeline/`，按数据处理阶段命名：

| 模块 | 职责 |
| --- | --- |
| `go_ontology.py` | 读取项目根目录的 `go_bp_ancestors.pkl`，提供 GO `is_a` 祖先关系与类别锚点。 |
| `go_bp_classification.py` | 底层 GO BP 分类。输入通路名称与可选 GO ID；输出细分类、合并类、背景、信号和主题。 |
| `gsea_classification.py` | 调用底层分类器，为 GSEA 结果增加 `absNES`、方向及按比较/主题统计。 |
| `gsea_summary.py` | 调用 GSEA 层，生成按 `Category` 或 `Category_Merged` 汇总的 Excel Summary。 |

依赖方向固定为 `go_ontology -> go_bp_classification -> gsea_classification -> gsea_summary`。分类规则只在 `go_bp_classification.py` 和 GO 锚点模块中维护；Summary 不另有一套通路分类规则。

## Python 调用

```python
from gobp_pipeline.go_bp_classification import classify_go_bp_term, classify_go_bp_table
from gobp_pipeline.gsea_classification import classify_gsea_results
from gobp_pipeline.gsea_summary import generate_gsea_summary

one = classify_go_bp_term("response to virus", "GO:0009615")
classified = classify_go_bp_table(go_bp_df)  # 必须有 ID、Description 两列
gsea_classified = classify_gsea_results(gsea_df)
```

## 命令行

从项目根目录运行：

```bash
python -m gobp_pipeline.gsea_summary input.csv summary.xlsx
python -m gobp_pipeline.gsea_summary input.csv summary_fine.xlsx --category-col Category
python gobp_pipeline/gsea_summary.py input.csv summary.xlsx
```

默认按 `Category_Merged` 汇总。原来的 `go_hierarchy.py`、`go_bp_classifier.py`、`categorize_fine_en.py`、`gsea_summary_table_generator_en.py` 现在是兼容入口，原有导入与旧命令仍可使用；未修改更早的 `categorize_fine.py` 和 `gsea_summary_table_generator.py`。


实际项目运行：
在src路径下
```bash
# 分类
python 05_run_category_gsea.py

# summary汇总
python -m gobp_category.gsea_summary \
  ../table/enrichment_results/gsea_sig_group.merged.classified.csv \
  ../table/enrichment_results/gsea_sig_group.merged.summary.xlsx
```

就会在./table/enrichment_results下生成分类结果和summary结果

