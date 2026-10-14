



setwd("/media/desk16/MaTianyu/Project/0X.Multiple_Rounds_RNA-seq")
source("src/F.gene_expression_trend.R")


# function(gene, project_dir = ".", groups = NULL,
#          group_colors = NULL, line_type = c("straight", "curve"),
#          annotation = c("tpm_fc", "deseq2_log2fc", "none"),
#          cycle_order = NULL, pseudocount = 0,
#          padj_cutoff = 0.05, log2fc_cutoff = 1,
#          significant_rule = c("padj", "padj_and_log2fc"),
#          digits = 2, point_size = 2.7, label_size = 2.8,
#          ggplot2_lib_loc = NULL) {
# }


p <- plot_gene_expression_trend(
  gene = "AKT3",
  project_dir = ".",
  groups = c("CQ", "N396"),
  group_colors = c(CQ = "#005596", N396 = "#CA6E28"),
  line_type = "straight",       # "curve" "straight" 即连接各轮次均值
  annotation = "tpm_fc"     # 也可用 "deseq2_log2fc" 或 "none"
)

print(p)



plot_gene_expression_trend(
  gene = "AKT3",
  project_dir = ".",
  groups = c("CQ", "N396"),
  group_colors = c(CQ = "#005596", N396 = "#CA6E28"),
  line_type = "curve",       # "curve" "straight" 即连接各轮次均值
  annotation = "deseq2_log2fc"     # 也可用 "deseq2_log2fc" 或 "none"
)



plot_gene_expression_trend(
  gene = "TP53",
  project_dir = ".",
  groups = c("CQ"),
  group_colors = c(CQ = "#005596"),
  line_type = "curve",       # "curve" "straight" 即连接各轮次均值
  annotation = "deseq2_log2fc"     # 也可用 "deseq2_log2fc" 或 "none"
)



plot_gene_expression_trend_grid(
  c("AKT3", "DNA2", "ENDOG"),
  ncol = 3,
  groups = c("CQ", "N396")
)







p <- plot_gene_expression_trend(
  gene = "AKT3", project_dir = ".",
  groups = c("CQ", "N396"),
  group_colors = c(CQ = "#005596", N396 = "#CA6E28"),
  group_shapes = c(CQ = 16, N396 = 17), show_legend = TRUE,
  line_type = "curve", annotation = "tpm_fc"
)

p_de <- plot_gene_expression_trend(
  "AKT3", project_dir = ".", line_type = "straight",
  annotation = "deseq2_log2fc", significant_rule = "padj"
)

p_grid <- plot_gene_expression_trend_grid(
  genes = c("AKT3", "DNA2", "ENDOG", "TP53", "SELL"), project_dir = ".",
  ncol = 5, groups = c("CQ", "N396"),
  group_colors = c(CQ = "#005596", N396 = "#CA6E28"),
  group_shapes = c(CQ = 16, N396 = 17), show_legend = FALSE,
  line_type = "straight", annotation = "tpm_fc",
  sort = "descending"
)
















