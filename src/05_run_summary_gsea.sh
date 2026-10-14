cd /media/desk16/MaTianyu/Project/0X.Multiple_Rounds_RNA-seq/src

python -m gobp_category.gsea_comparison_summary \
  ../table/enrichment_results/gsea_sig_group.merged.classified.csv \
  ../table/enrichment_results/gsea_category_comparison_by_group \
  --category-col Category

python -m gobp_category.gsea_summary \
  ../table/enrichment_results/gsea_sig_group.merged.classified.csv \
  ../table/enrichment_results/gsea_sig_group.merged.summary.xlsx
