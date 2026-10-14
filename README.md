# 多轮次 RNA-seq 分析流程

这是一个以**单个组内、不同轮次之间比较**为单位的 RNA-seq 下游分析模板。它把 DESeq2 差异分析、相邻（或自定义）比较的趋势统计，以及 GO-BP ORA/GSEA 拆成可复用函数，并通过配置文件批量运行。

## 目录约定

```text
.
├── config.example.yml          # 复制为 config.yml 后按项目修改
├── data/
│   ├── counts.tsv              # 输入：原始整数 counts 矩阵（不提交到模板）
│   ├── tpm.tsv                 # 输入：GSVA 打分用 TPM 矩阵（不提交到模板）
│   ├── metadata.tsv            # 输入：样本信息（不提交到模板）
│   ├── contrasts.tsv           # 输入：组内 DE 比较清单
│   ├── trend_sets.tsv          # 输入：趋势统计清单
│   ├── GO_BP_TERM2GENE.tsv     # 输入：冻结的 GO-BP term-gene 对照表
│   ├── GO_BP_TERM2NAME.tsv     # 输入：冻结的 GO-BP term-name 对照表
│   ├── de_results/             # 自动生成：DESeq2 结果
│   ├── gsva/                   # 自动生成：全样本 GSVA 分数、RDS 与通路使用记录
│   ├── gsva_results/           # 自动生成：组内 GSVA-limma 结果
│   └── pipeline/               # 自动生成：运行清单和日志
├── figures/                    # 自动生成：PCA、MA 图等
├── notebooks/                  # 用于项目解释性分析
├── scripts/                    # 运行入口与旧脚本
└── src/                        # 可复用分析函数
```

`data/` 下的输入表均为 TSV。输出始终保留在本项目内，因而不同项目只需复制这套目录并替换数据和配置，不会相互覆盖。

如果服务器有多个 `ggplot2` 版本，可在配置顶层填写 `runtime.ggplot2_lib_loc`，例如 `"/opt/R/4.4.1/lib/R/library"`。所有入口会在加载 `patchwork` 或其他绘图依赖前，从该目录加载 `ggplot2`；若当前 R 会话已经加载了另一版本，会直接报错并要求使用新的会话，避免版本混用。

## 输入格式

`counts.tsv` 的第一列是基因 ID，其他列是样本名；必须为未标准化的、非负整数原始 count。样本名必须与 `metadata.tsv` 的 `Sample` 列一一对应。

`metadata.tsv` 至少包含：

```text
Sample  Group  Cycle  Batch
S01     A      R0     B1
S02     A      R0     B2
```

- `Group` 用于拆分多个独立实验组；模型只在同一 `Group` 内拟合。
- `Cycle` 是轮次/时间点，如 `R0`、`R1`、`R2`。
- `Batch` 是可选协变量。若没有批次列，配置中的 `design_covariates` 设为 `[]`。

GSVA 另需一个 `tpm.tsv`：第一列仍为 `Gene`，其余列与 counts/metadata 使用完全相同的样本名。流程以 raw counts 做表达过滤（默认 `> 10` 的样本数至少为 3），再以 `log2(TPM + 1)` 做 GSVA；TPM 不参与 DESeq2。

富集分析不依赖 `OrgDb`、`enrichGO()` 或 `gseGO()`。将预先保存的 GO-BP 对照表放入 `data/`：

```text
# GO_BP_TERM2GENE.tsv
term        gene
GO:0000002  AKT3

# GO_BP_TERM2NAME.tsv
term        name
GO:0000002  mitochondrial genome maintenance
```

流程用同一对文件调用 `clusterProfiler::enricher()` 和 `clusterProfiler::GSEA()`，因此后续任何使用这两份表的模块都会采用同一套通路边界。`Gene` 列必须与 `TERM2GENE$gene` 使用同一 ID 类型；流程不会自动转换 SYMBOL、ENSEMBL 或 ENTREZID。每次富集还会在 `data/pipeline/pathway_reference.tsv` 记录对照表的路径、行数和 MD5。

每个成功完成的比较除了阈值筛选后的 TSV，还会保存原始对象：`ora_up.rds`、`ora_down.rds`、`ora_up_down.rds`（三个 `enrichResult`）以及 `gsea_go_bp.rds`（一个 `gseaResult`）。RDS 在 TSV 导出前保存；ORA 的 `@result` 与 GSEA 的 `@result` 都保留全部通过基因集大小筛选的通路，当前阈值只影响 TSV。可用 `readRDS()` 读取后直接绘图，或调整对象阈值后重新导出：

```r
gsea <- readRDS("data/enrichment_results/CQ/R1_vs_R0/gsea_go_bp.rds")
gsea_table_relaxed <- subset(as.data.frame(gsea), pvalue <= 0.1 & p.adjust <= 0.1)

ora_up <- readRDS("data/enrichment_results/CQ/R1_vs_R0/ora_up.rds")
ora_up@pvalueCutoff <- 0.1
ora_up@qvalueCutoff <- 1
ora_up_table_relaxed <- as.data.frame(ora_up)
```

在 `contrasts.tsv` 中，每行定义一个组内比较：

```text
group_id  comparison_id  numerator  denominator  enabled
A         R1_vs_R0       R1         R0           TRUE
A         R2_vs_R1       R2         R1           TRUE
B         R1_vs_R0       R1         R0           TRUE
```

`numerator - denominator` 决定 `log2FoldChange` 的方向。流程会对每项比较单独检查每组重复数、设计矩阵是否满秩及低表达基因数；无法检验的比较会写入 manifest，而不是悄悄丢失。

`trend_sets.tsv` 让每个组选择用于趋势定义的比较，`comparison_ids` 以分号连接：

```text
group_id  trend_id       comparison_ids
A         adjacent       R1_vs_R0;R2_vs_R1;R3_vs_R2
B         adjacent       R1_vs_R0;R2_vs_R1
```

## 快速开始

```bash
cp config.example.yml config.yml
# 填入 data/counts.tsv、data/metadata.tsv、data/contrasts.tsv、data/trend_sets.tsv
Rscript scripts/check_environment.R --config config.yml
Rscript scripts/run_pipeline.R --config config.yml
```

也可以单独重跑某一步：

```bash
Rscript scripts/run_pipeline.R --config config.yml --stage de
Rscript scripts/run_pipeline.R --config config.yml --stage enrichment
Rscript scripts/run_pipeline.R --config config.yml --stage gsva
```

核心 pipeline 的顺序是 DESeq2、富集和 GSVA；趋势统计不在这里提前运行。完成富集、GSEA 分类后，再用 `run_post_category_trends.R` 同时计算基因和分类后通路趋势。富集步骤会读取 `data/pipeline/de_manifest.tsv`，所以单独运行前必须已经完成 DE 步骤。使用 `--dry-run` 可只检查核心步骤的配置和输入路径，不写结果。

源码编号现按实际依赖关系组织：`01_deseq2.R` → `02_enrichment.R`；`03_gsva.R` 为独立分支；手动分类工具位于 `src/05_*` 与 `src/gobp_category/`；分类完成后才运行 `06_post_category_trends.R`；随后为 `07_gobp_gsva_plots.R` 与 `08_result_plots.R`。因此不存在容易误解为“第二步”的 `02_trends.R`。

GSVA 与 DESeq2 使用同一份 `contrasts.tsv`：每个 `group_id` 内的 `numerator - denominator` 以 limma 检验 GSVA 分数。默认同时运行两条彼此独立的分支：全样本统一打分后按组内检验，以及每个 `Group` 仅用其自身样本重新打分后再做相同的组内检验。第二条分支由 `gsva.run_groupwise_scoring` 控制。两条分支使用相同的 counts 过滤和通路集，以保证通路定义可比；差异仅来自 GSVA 打分所用样本集合。它还会按组生成可选的 `specific_results.tsv`（某轮次相对于其余轮次均值）和严格峰值筛选结果。主要产物为：

- `data/gsva/gsva_scores.tsv` 与 `gsva_scores.rds`：全样本通路分数；
- `data/gsva_results/<group>/<comparison>/limma_results.tsv`：每个组内比较的完整检验表；
- `data/gsva/by_group/<group>/gsva_scores.tsv`：仅以该组样本计算的 GSVA 分数；
- `data/gsva_results/by_group/<group>/<comparison>/limma_results.tsv`：基于组别独立打分的检验表；
- `limma_significant.tsv`：按 GSVA 配置阈值筛选后的表；
- `data/pipeline/gsva_manifest.tsv`：每个比较的状态、显著通路数与结果路径；
- `data/pipeline/gsva_groupwise_score_manifest.tsv` 与 `gsva_groupwise_manifest.tsv`：分组独立打分与检验的状态和路径；
- `figures/gsva_results/<group>/pairwise_significant_pathway_counts.pdf`：上下调通路数概览。

## GO-BP 分类后趋势与 GSVA 绘图

`src/gobp_category/` 保持为独立的分类工具，不由主 pipeline 改写。完成 GSEA 后，在 [Mapping_GOBP_Category.Rmd](notebooks/Mapping_GOBP_Category.Rmd) 依次完成：导出 `gsea_sig_group.merged.csv`、运行既有的 `05_run_category_gsea.py` 与 `05_run_summary_gsea.sh`、读取分类结果并计算 Category 的时间趋势。

趋势计算不再从样本名解析轮次。Rmd 直接用 `metadata.tsv` 中配置的样本、组别和轮次列匹配每个 Group 的 `data/gsva/by_group/<group>/gsva_scores.rds`；`gobp_temporal_trend.cycle_order` 可固定任意项目的实际时间顺序。输出保存在：

```text
table/enrichment_results/gobp_temporal_trends/<group>/
├── temporal_pathway_analysis.rds
├── category_trend_summary.tsv
├── category_activity_score.tsv
├── direction_metrics.tsv
└── pathway_trajectory.tsv
```

分类及汇总表生成后，可执行独立绘图阶段：

```bash
Rscript scripts/run_gobp_gsva_plots.R --config config.yml
```

它为每个 `Group × Category` 输出该类别全部已分类 GSEA 通路的 GSVA 热图和时间轨迹；同一 Group 的所有类别轨迹共用纵轴范围。它还从 `gsea_sig_group.merged.summary.xlsx` 的每个 Group sheet 读取代表通路，并生成保留配置顺序的总热图和从上到下同序的纵向趋势图：

```text
figures/gobp_gsva/
├── plot_manifest.tsv
└── <group>/
    ├── categories/<category>/
    │   ├── gsva_heatmap.pdf
    │   ├── gsva_trajectory.pdf
    │   ├── gsva_trajectory.tsv
    │   └── pathway_mapping.tsv
    └── representative_pathways/
        ├── representative_pathways_heatmap.pdf
        ├── representative_pathways_trajectory.pdf
        ├── representative_pathways_trajectory.tsv
        └── representative_pathway_mapping.tsv
```

`gobp_gsva_plots.category_order` 和 `category_colors` 用于指定代表通路总热图/趋势图的类别顺序与颜色；`cycle_order`、`cycle_colors` 和每个 Group 的 `line_plot.y_limits` 也在同一配置节中控制。代表通路趋势图按类别纵向拼接、共用自动计算的对称 Y 轴，且不重复标注类别文字；类别颜色和上下顺序与总热图保持一致。`heatmap.representative_extra_width` 只会拓宽代表通路总热图的外围画布，不改变 heatmap cell 的长宽。未指定类别或颜色时会自动补全，不会被静默丢弃。

## DE、趋势与 GSEA 结果汇总绘图

DESeq2、趋势统计和富集完成后，可独立重跑汇总绘图，不会重新进行统计检验：

```bash
Rscript scripts/run_result_plots.R --config config.yml
```

该阶段生成四类图：每项 DE 比较的火山图及按 comparison 坐标拼接的总火山图；指定 `trend_id` 的趋势模式计数柱状图；将 DEG 上下调数和显著 GSEA 正/负富集通路数放入同一 comparison 布局的柱状图；以及逐 comparison 的富集图。富集图包括每个方向各取 Top 10/20/50 或全部显著结果的 GSEA 双向 NES 柱图、ORA 的 up/down/up_down Top N 或全部显著通路气泡图、分类后 GSEA 的类别通路数双向柱图和每个类别最显著通路气泡图。默认使用 `trend.fc_cutoff` / `trend.padj_cutoff`，因此与 ORA 的 DE 阈值一致；`result_plots.fc_cutoff` 和 `padj_cutoff` 可同时覆盖火山图与 DEG 计数阈值。

`result_plots.comparison_layout.positions` 用 `[行, 列]` 指定每个 comparison 的位置，空位保持为空，可表示 R0-R3 的三角布局，也适用于任意自定义时间点/比较。`volcano.axis_scope: within_group` 保证同组总图坐标一致，改为 `global` 则所有组完全一致；X 轴始终以 0 对称。火山图的手动数字位置、三类颜色、各图尺寸、PNG DPI、趋势柱的组别颜色和可选 Y 轴断轴均可在 `result_plots` 配置节调整。使用 `trend_bar.y_breaks` 时需额外安装 `ggbreak`；断轴趋势图会走专用单页 PDF 渲染，并由 `trend_bar.label_top_padding` 保留柱顶数字空间。

`result_plots.enrichment` 控制富集图。GSEA/ORA 优先读取原始 RDS、按绘图配置重新筛选；这不会修改原有 TSV。`gsea_bar.top_n` 和 `ora_dot.top_n` 均可写多个值；写入字符串 `"all"` 会另输出全部显著通路图。各图的 `save.width` / `save.height` 独立解析，因而 `{width: 4, height: null}` 会固定宽度而继续自动计算高度；需要按图中通路数分别设置宽度时，使用 `save.by_plot.top_10`、`save.by_plot.top_20`、`save.by_plot.top_50` 或 `save.by_plot.all`。`all` 图默认只保存 PDF：自动高度可能使 300 dpi 的 PNG 超过 ragg 的 50,000 px 单边限制。若确实需要 PNG，可在 `save.by_plot.all.formats` 显式覆盖，同时降低高度、宽度或 DPI。分类图读取 `classified_gsea.file`，通常就是分类脚本产生的 `gsea_sig_group.merged.classified.csv`；类别列可由 `category_column` 在 `Category`、`Category_CN`、`Category_Merged` 等列间切换。

`classified_gsea.category_sankey` 按 `trend_sets.tsv` 中每个 Group × trend_id（例如 `adjacent`、`ref0`）统计各比较的显著分类 GSEA 通路组成，并输出 `count` 与 `proportion` 两个直线过渡的 alluvial 图。它追踪的是类别组成变化，不是单通路显著性趋势。默认自动继承 `gobp_gsva_plots.category_order` 和 `category_colors`，保证类别顺序与颜色和 GSVA 热图/轨迹图一致；也可以在该配置节单独覆盖。

输出目录默认如下：

```text
figures/result_plots/
├── plot_manifest.tsv
├── volcano/<group>/<comparison>/volcano.{pdf,png}
├── volcano/<group>/volcano_combined.{pdf,png}
├── trend_bars/<group>/<trend_id>/trend_bar.{pdf,png}
├── deg_gsea_counts/<group>/deg_gsea_count_combined.{pdf,png}
└── enrichment/
    ├── gsea_bidirectional_bar/<group>/<comparison>/top_<N>.{pdf,png}；all.pdf
    ├── ora_dot/<group>/<comparison>/ora_<set>_top_<N>.{pdf,png}；ora_<set>_all.pdf
    ├── classified_gsea/<group>/<comparison>/
    │   ├── category_bidirectional_bar.{pdf,png}
    │   └── category_representative_bubble.{pdf,png}
    └── classified_gsea_category_sankey/<group>/<trend_id>/
        ├── category_composition_count.{pdf,png}
        └── category_composition_proportion.{pdf,png}
```

## 分类后基因与 GSEA 通路趋势追踪

在 [Mapping_GOBP_Category.Rmd](notebooks/Mapping_GOBP_Category.Rmd) 完成分类、生成 `gsea_sig_group.merged.classified.csv` 后，运行：

```bash
Rscript scripts/run_post_category_trends.R --config config.yml
```

这个阶段先按 `trend_sets.tsv` 重新计算基因 `trend_table.tsv`，再为分类后 GSEA 通路建立 `Positive / NS / Negative` 状态表；两条分支共享完全相同的 `Group × trend_id × comparison_ids` 定义。`adjacent` 和 `ref0` 只是 `trend_sets.tsv` 中不同的 `trend_id`，因此可在同一运行中各自处理任意数量的比较，例如：

```text
group_id  trend_id  comparison_ids
CQ        adjacent  R1_vs_R0;R2_vs_R1;R3_vs_R2
CQ        ref0      R1_vs_R0;R2_vs_R0;R3_vs_R0
```

每个 `Group × trend_id` 会输出基因和 GSEA pathway 的每一种非全稳定模式 TSV，并分别汇总为一个 Excel 工作簿（每种模式一个 sheet）；GSEA 的完整趋势表另存为 `data/gsea_trend_results/`。原有通路趋势及桑基图仍以“该组所有比较中曾显著的通路”为池子，其中 NS 表示该通路在当前比较未显著。桑基图同时输出带 strata 内标签和无标签版本，且对每个实际出现的非全稳定子趋势另绘制“仅高亮该簇 flow”的图。高亮版的带标签图只显示各状态及其总数，普通版保留相邻比较的转移计数。高亮 flow 带可配置深色描边（`highlight_outline_color`、`highlight_outline_width`），使数量较少的分支更容易看到。绘图时按完全相同的趋势路线汇总权重，重复使用汇总后的数据，保留柱高和流量统计，减少大量重复 flow 的渲染时间。`post_category_trends.sankey.stratum_width` 默认为 `0.48`，比参考函数中的窄柱更接近，可按审美继续调整。

另外，`post_category_trends.sankey.full_universe.enabled: true` 会额外绘制原始 GSEA 通路全集版：从每个比较的 `gsea_go_bp.rds` 的未过滤 `@result$ID` 取并集，再用分类后的显著结果赋予 Positive/Negative；其余通路（包括当前比较未受检者）统一显示为 NotSig，不单列 NotTested。`plot_state_counts.tsv` 的 `n_nonsig` 因此是“全集数量减当前比较显著数量”，不等同于仅在该比较受检的非显著数量。此分支只生成绘图及计数表，不改动前述趋势表、Excel 和显著通路池桑基图；默认只绘制带/不带标签的总图。完整池绘图失败也会单独记录在 `data/pipeline/post_category_gsea_full_universe_plot_manifest.tsv`，不影响原有趋势统计的状态。

原有基因与显著通路池的绘图数量仍随实际出现的趋势模式增加：每个 `Group × trend_id × 特征类型` 会生成 `2 × (1 + 模式数) × 格式数` 张图。完整通路池版默认只增加两种标签版本的总图，不生成逐模式高亮图。默认 PDF 与 PNG 都输出；只需要一种版本时可设置 `sankey.variants: ["with_labels"]`，只需要 PDF 时可设置 `sankey.save.formats: ["pdf"]`，进一步缩短导出时间。

```text
data/gsea_trend_results/<group>/<trend_id>/
├── trend_table.tsv
└── trend_summary.tsv
table/trend_feature_subsets/<group>/<trend_id>/<genes|gsea_pathways>/
├── <feature_type>_trend_patterns.xlsx
└── patterns/<Full_Trend>.tsv
figures/trend_sankey/<group>/<trend_id>/<genes|gsea_pathways>/
├── sankey_with_labels.{pdf,png}
├── sankey_without_labels.{pdf,png}
└── highlighted/<Full_Trend>_<variant>.{pdf,png}
figures/trend_sankey/<group>/<trend_id>/gsea_pathways_full_universe/
├── sankey_with_labels.{pdf,png}
├── sankey_without_labels.{pdf,png}
└── plot_state_counts.tsv
```

该步骤的 Excel 和桑基图分别需要 `openxlsx` 与 `ggalluvial`；`result_plots` 的分类通路气泡图需要 `ggrepel`。

## 单基因与多基因表达趋势图

[F.gene_expression_trend.R](src/F.gene_expression_trend.R) 直接读取当前项目的 `data/metadata.tsv` 和 `data/tpm.tsv`，以 metadata 的 `Sample`、`Group`、`Cycle` 匹配数据。`R0`、`R1`、`R2` 等轮次自动按数字排序和定位；其他轮次名沿 metadata 首次出现顺序，也可传入 `cycle_order`。某个 Group 缺少一个轮次时不补造观测点，TPM FC 按该组相邻的实际观测轮次计算（例如 `R1→R3`）。`line_type="curve"` 按 demo1 的做法对组内样本 TPM 使用 LOESS 局部平滑并绘制置信带；不足 3 个已观测轮次时用直线。图例只保留各 Group 的点形状与颜色，不显示线条。

```r
source("src/F.gene_expression_trend.R")

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
  genes = c("AKT3", "DNA2", "ENDOG"), project_dir = ".",
  ncol = 2, groups = c("CQ", "N396"),
  group_colors = c(CQ = "#005596", N396 = "#CA6E28"),
  group_shapes = c(CQ = 16, N396 = 17), show_legend = FALSE,
  line_type = "straight", annotation = "tpm_fc",
  sort = "descending"
)
```

多基因的 `sort` 可设为 `"none"`（默认，保持输入顺序）、`"descending"` 或 `"ascending"`；排序值是本次所选 Group 内全部样本的平均 TPM，按图中的从左到右、从上到下顺序排列。均值相同时保持输入顺序。
`group_shapes` 接受包含所有所选 Group 的命名整数向量（ggplot2 点形状编号 0–25）；不指定时沿用原来的自动形状。`show_legend` 对单图和多图均默认为 `TRUE`，设为 `FALSE` 可隐藏图例。

`annotation="deseq2_log2fc"` 依据 `data/contrasts.tsv` 的 `numerator - denominator` 查找相邻观测轮次的同方向比较，读取 `data/de_results/<group>/<comparison>/deseq2_results.tsv`；显著值后加 `*`。若存在 `data/pipeline/de_manifest.tsv`，仅使用其中 `status=success` 的比较。缺少比较或结果时给出警告并省略该区间标注，不会用 TPM FC 代替。两个函数返回绘图对象，由调用者用 `ggplot2::ggsave()` 保存；如项目 `config.yml` 设置了 `runtime.ggplot2_lib_loc`，函数会先从指定库加载 ggplot2。

## 统计约定

- DESeq2 仅接受原始 count；TPM、FPKM、CPM 或已做 log 转换的矩阵不能直接作为输入。
- 默认每个比较的两个轮次至少各有 2 个生物学重复；过滤要求至少 3 个样本的 count ≥ 10，可在配置中调整。
- 批次等协变量在组内模型中以 `~ Batch + Cycle` 形式建模；如果与轮次完全混杂，流程会跳过该比较并说明原因。
- ORA 的背景基因使用该比较中实际进入检验、且可映射到 `TERM2GENE` 的基因；GSEA 默认按 DESeq2 Wald statistic 排序。结果 manifest 会记录每项比较实际匹配到通路表的基因数。

旧的三个脚本保留在 `scripts/` 以便追溯，但新项目请从 `scripts/run_pipeline.R` 运行。
