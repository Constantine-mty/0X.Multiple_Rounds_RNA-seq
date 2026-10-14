# 多时间点 Pathway / GO-BP Category 动态趋势定量分析

## 1. 目的

该流程用于把“多个 pathway 在多个有序时间点的样本级评分”压缩为可解释的 **生物学类别（Category）时间动态**。输入分数可以是 GSVA、ssGSEA、module score、signature score 或其他连续型 pathway/activity score。

核心思想不是只看热图，而是同时描述四件事：

1. **中心轨迹**：一个 Category 在每个时间点总体处于什么水平（CAS）。
2. **变化幅度**：相邻时间点以及全程改变了多少（ΔCAS、Amplitude、Net Change）。
3. **内部一致性**：Category 中有多少 pathway 支持总体方向（NDI、DCS）。
4. **轨迹形态**：持续激活、持续抑制、瞬时变化、晚期变化、波动或内部异质等。

> 若 pathway 集合来自前置 GSEA 显著结果筛选，本流程的定位是 **对已筛选 pathway 的 temporal characterization（时间动态刻画）**，而不是一次独立的显著性发现。避免把后续 GSVA 趋势描述成对同一假设的再次独立验证。

---

## 2. 输入数据

### 2.1 score matrix

行为 pathway ID，列为样本：

```text
             Sample_R0_1 Sample_R0_2 Sample_R1_1 ...
GO:0000002      0.31         0.44        0.05
GO:0000012      0.21         0.41        0.13
```

### 2.2 annotation

至少包含：

- `Group`
- `ID`
- `Description`
- `Category`

可以额外包含 `Category_CN`、`Category_Merged`。

同一个 pathway ID 原则上应只属于一个用于本分析的 Category；函数会检查多重分类并给出 warning。

### 2.3 time_map

用正则表达式从样本名识别时间：

```r
time_map <- c(
  "_Non_" = "R0",
  "_R1_"  = "R1",
  "_R2_"  = "R2",
  "_R3_"  = "R3"
)
```

向量顺序就是时间顺序。

---

## 3. 最简使用方法

```r
source("temporal_pathway_trend.R")

time_map <- c(
  "_Non_" = "R0",
  "_R1_"  = "R1",
  "_R2_"  = "R2",
  "_R3_"  = "R3"
)

cq.res <- run_temporal_pathway_analysis(
  score_matrix = gsva.cq.res,
  annotation   = gsea.category.df,
  group_name   = "CQ",
  time_map     = time_map,
  zscore       = TRUE
)

cq.res
cq.res$summary
cq.res$diagnostics
```

N396 同理：

```r
n396.res <- run_temporal_pathway_analysis(
  score_matrix = gsva.n396.res,
  annotation   = gsea.category.df,
  group_name   = "N396",
  time_map     = time_map,
  zscore       = TRUE
)
```

---

# 4. 是否进行 Z-score

参数：

```r
zscore = TRUE
```

开启时，对**每一条 pathway 自己**跨全部样本标准化：

\[
Z_{ij}=\frac{S_{ij}-\bar S_i}{SD(S_i)}
\]

其中 i 是 pathway，j 是样本。

## 为什么 GSVA 仍然可以 Z-score？

GSVA score 已经是样本层面的 pathway enrichment/activity score，因此并不是“必须再次归一化”。Z-score 的目的不同：它把不同 pathway 转换到共同的相对变化尺度，使 Category 聚合时不会主要由数值范围更大的 pathway 决定。

### `zscore = TRUE` 更适合

- 主要问题是“每条 pathway 相对于自身基线如何变化”；
- 不同 pathway score 范围差异明显；
- 希望把 GSVA、module score 等统一放入同一套算法；
- Category 是多个异质 pathway 的聚合。

### `zscore = FALSE` 更适合

- 原始 score 的绝对尺度本身有明确意义；
- 所有 pathway 的 score 本来就严格可比；
- 希望保留原始 GSVA effect size。

建议在正式分析中固定一种方案，并在 Methods 中说明。可以用 TRUE/FALSE 做 sensitivity analysis，但不要根据哪个结果“更漂亮”来选择。

---

# 5. Pathway-level trajectory

每个时间点存在多个生物学重复。先计算每条 pathway 在每个时间点的平均分析分数：

\[
\bar S_{i,t}=\frac{1}{n_t}\sum_j S_{i,t,j}
\]

输出：

```r
res$pathway_trajectory
```

其中：

- `mean_raw_score`：原始 score 的时间点均值；
- `mean_score`：真正进入趋势算法的值（Z-score 开启时为 Z-score）；
- `sd_score`：该时间点重复间 SD；
- `n_sample`：重复数。

---

# 6. Category Activity Score（CAS）

对于 Category C，在时间 t：

\[
CAS_{C,t}=median_{i\in C}(\bar S_{i,t})
\]

默认使用 **median** 而非 mean，以降低少数极端 pathway 对类别中心趋势的影响。

同时输出 `CAS_mean`、Q25、Q75 和 pathway 数量。

```r
res$CAS
```

CAS 是后续 Category-level temporal trajectory 的核心变量。

### CAS 的实际意义与判读

CAS 回答的是：**“在某个时间点，这一生物学类别中典型 pathway 的相对活性处在什么位置？”**

- `CAS > 0`：若使用 pathway-wise Z-score，说明该类别中典型 pathway 在该时间点高于其自身跨时间总体平均水平；
- `CAS < 0`：说明典型 pathway 低于其自身跨时间总体平均水平；
- `CAS ≈ 0`：说明整体接近各 pathway 自身的平均状态。

需要特别注意：**CAS 的正负不是“激活/抑制”的充分条件。** 动态趋势应看 CAS 随时间的相对变化。例如 `[-0.8, -0.4, -0.1, 0.2]` 虽然前三个时间点 CAS 都为负，但这是清楚的持续增强轨迹。相反 `[0.8, 0.5, 0.2, -0.1]` 是持续减弱。

当 `zscore = FALSE` 时，CAS 保留原 score 的尺度，此时其绝对值是否能解释为“高/低活性”取决于原始评分方法本身，不能套用 Z-score 的零点解释。

---

# 7. Δ score 与 ΔCAS

Pathway 相邻时间变化：

\[
\Delta S_{i,t}=\bar S_{i,t+1}-\bar S_{i,t}
\]

Category 相邻时间变化：

\[
\Delta CAS_{C,t}=CAS_{C,t+1}-CAS_{C,t}
\]

输出：

```r
res$pathway_delta
res$category_delta
```

例如四个时间点会产生三个 interval：R0→R1、R1→R2、R2→R3。

### Δ score / ΔCAS 的实际意义

它们回答的是“**相邻两个时间点之间发生了什么**”，而不是整段实验的最终趋势。

- `ΔCAS > 0`：Category 中心活性向上移动；
- `ΔCAS < 0`：Category 中心活性向下移动；
- `|ΔCAS|` 越大：该 interval 的 Category-level shift 越明显；
- pathway-level `Δ score` 则用于判断 Category 内部各 pathway 是否共同支持这一移动。

例如 CAS 为 `[-0.4, 0.2, 0.5, 0.1]`，则三个 ΔCAS 为 `[+0.6, +0.3, -0.4]`：前两个阶段增强、最后阶段回落。单看终点与起点会漏掉这种中期峰值，因此 ΔCAS 是 Pattern classification 的主要形态信息。

---

# 8. epsilon：什么程度才算真正发生方向变化？

不能把任意微小的浮点变化都定义为 Up/Down，因此引入 epsilon：

```text
Δ > +epsilon  -> Up
Δ < -epsilon  -> Down
其余          -> Stable
```

本模块将其拆成两个阈值：

- `epsilon_pathway`：判断单 pathway 的 Δ score；
- `epsilon_category`：判断 Category 的 ΔCAS。

这是必要的，因为 pathway Δ 与聚合后的 ΔCAS 通常具有不同分布。

## 自动模式

默认：

```r
thresholds = list(
  mode = "auto",
  pathway_quantile = 0.25,
  category_quantile = 0.25,
  amplitude_quantile = 0.25,
  dcs_cutoff = 0.60
)
```

自动建议：

```text
epsilon_pathway  = Q25(|pathway Δ|)
epsilon_category = Q25(|category ΔCAS|)
```

这里的 Q25 **不是统计显著性 cutoff**，只是一个 scale-adaptive operational threshold：把数据集中最低幅度的一部分变化视为接近稳定，从而避免对不同尺度数据硬编码 0.05。

## 如何校准 epsilon

查看：

```r
res$diagnostics$abs_pathway_delta
res$diagnostics$abs_category_delta
```

会返回 Q10/Q25/Q50/Q75/Q90。

推荐把 Q25 作为探索性起点，再做 Q10、Q25、Q50 sensitivity analysis：如果 Category pattern 在合理阈值范围内稳定，结论更可靠；若分类高度依赖 epsilon，应报告其不稳定性，而不是寻找能产生预期结论的 cutoff。

---

# 9. NDI：Net Direction Index

对于一个 Category 的某个时间区间：

\[
NDI=\frac{N_{up}-N_{down}}{N_{total}}
\]

范围 [-1,1]。

- 接近 +1：绝大多数 pathway 上升；
- 接近 -1：绝大多数下降；
- 接近 0：方向互相抵消，或大量 pathway 稳定。

同时输出：

- `P_up`
- `P_down`
- `P_stable`

NDI 表达的是**净方向偏向**，不是显著性。

### NDI 的实际意义与判读

NDI 回答：**“这个 Category 内，上升 pathway 和下降 pathway 哪一方在数量上占优势？”**

例如 20 条 pathway 中 14 条 Up、4 条 Down、2 条 Stable，则 `NDI=(14-4)/20=0.50`，说明该 interval 明显偏向上升。若 9 条 Up、9 条 Down、2 条 Stable，则 NDI=0，表示内部方向相互抵消，即使个别 pathway 的变化幅度可能很大。

因此：

- 高正 NDI：广泛的正向重塑；
- 高负 NDI：广泛的负向重塑；
- NDI 接近 0：可能是真正稳定，也可能是“强烈但方向相反”的异质变化。

所以 **NDI≈0 不能单独解释为 Stable**，必须结合 Amplitude、P_stable 和 DCS。

---

# 10. DCS：Direction Consistency Score

模块输出两个版本。

### DCS_all

三分类（Up/Down/Stable）完全一致的 pathway 比例。

### DCS_directional

当 Category 为 Up 时：

\[
DCS=P(pathway=Up)
\]

Category 为 Down 时：

\[
DCS=P(pathway=Down)
\]

Category 为 Stable 时：

\[
DCS=P(pathway=Stable)
\]

趋势分类默认使用 `DCS_directional`，因为它更直接回答：**有多少 pathway 真正支持 Category 的总体方向？**

### DCS 的实际意义

DCS 描述的是 **coherence / coordination（协同性）**，而不是变化强度。

例如两个 Category 的 ΔCAS 都是 +0.30：

- Category A：80% pathway 同时 Up，`DCS_directional=0.80`，说明这是较协调的功能程序增强；
- Category B：52% Up、45% Down，虽然中位数仍可能上升，但 `DCS_directional≈0.52`，说明类别内部存在明显分化。

因此，同样的 CAS 轨迹可以有完全不同的生物学含义：高 DCS 更接近“该功能模块整体协同变化”，低 DCS 则提示该人工分类下包含多个不同响应亚程序。

`DCS_all` 更严格地要求 pathway 与 Category 在 Up/Down/Stable 三分类上完全一致，适合质量控制；`DCS_directional` 更适合回答生物学上的“有多少 pathway 支持总体方向”。

输出：

```r
res$direction_metrics
```

---

# 11. DCS cutoff 如何选择

默认：

```r
dcs_cutoff = 0.60
```

这是一个可解释的经验性一致性标准，而不是统计学天然阈值。0.60 表示至少约 60% pathway 支持总体方向，否则优先标记为 `Heterogeneous`。

查看实际分布：

```r
res$diagnostics$dcs
```

建议同时检查 0.55、0.60、0.67、0.70 等阈值下结果是否稳定。

Category pathway 数量很少时，DCS 是离散的。例如只有 3 条 pathway 时，DCS 只能以 1/3 为步长变化，因此不应过度解释 0.60 与 0.67 的细微区别。正式分析可额外设置最小 pathway 数，例如 n < 5 的 Category 单独标记为 low-support。

---

# 12. Global trend metrics

```r
res$global_metrics
```

包含：

### Trend_rho

Spearman correlation：

\[
\rho=cor_{Spearman}(time,CAS)
\]

反映单调上升/下降倾向。它是辅助描述指标，不负责识别瞬时峰值。

### Amplitude

\[
Amplitude=max(CAS)-min(CAS)
\]

表示全时间序列动态范围。

### Peak_Time / Trough_Time

分别是 CAS 最大/最小的时间点。

### Net_Change

\[
CAS_{last}-CAS_{first}
\]

表示终点相对于起点的净改变。

### 这些 global metrics 分别回答什么？

- **Trend_rho**：回答“整体是否近似单调”。接近 +1 表示时间越晚 CAS 越高，接近 -1 表示越晚越低；接近 0 不等于无变化，也可能存在明显峰谷。
- **Amplitude**：回答“整个过程中最多拉开了多大距离”，只看幅度、不看方向。
- **Peak_Time / Trough_Time**：回答“最高点/最低点发生在什么时候”，适合区分 early、transient、late pattern。
- **Net_Change**：回答“实验终点相对起点最后留下了多少改变”。它会忽略中间过程。

例如 CAS `[0, 1.0, 0.9, 0.1]` 的 Amplitude 很大，但 Net_Change 很小；这不是 Stable，而是典型的 transient activation。CAS `[0, 0.3, 0.6, 0.9]` 则同时具有正 Trend_rho、较大 Amplitude 和正 Net_Change，更符合 sustained activation。

---

# 13. amplitude cutoff

即使三个 ΔCAS 有微小正负变化，如果整个轨迹动态范围极小，也不应该称为 biological activation/suppression。

因此分类首先检查：

```text
Amplitude <= amplitude_cutoff -> Stable
```

自动模式默认：

```text
amplitude_cutoff = Q25(Category Amplitude)
```

同样，这是数据驱动的 operational threshold，而不是显著性阈值。

查看：

```r
res$diagnostics$amplitude
```

推荐进行 Q10/Q25 等 sensitivity analysis。如果有独立的技术噪声、阴性对照或重复实验误差估计，优先用这些信息定义 minimal meaningful change，而不是分位数。

---

# 14. Pattern classification（趋势模式分类）

Pattern classification 的目标不是再次检验显著性，而是把连续的 CAS 时间轨迹转换成**可重复、可解释的动态形态标签**。分类综合：

1. **Amplitude**：先判断整条轨迹是否具有足够动态范围；
2. **Category ΔCAS 的方向序列**：决定轨迹的主要形态；
3. **DCS_directional**：判断 Category 内 pathway 是否足够协同；
4. **NDI、Trend_rho、Peak/Trough、Net Change**：作为解释和质控信息。

## 14.1 Pattern 标签中英文定义

| Pattern | 中文名称 | 典型方向示意 | 实际含义 | 推荐描述 |
|---|---|---|---|---|
| `Stable` | 稳定型 / 低动态型 | → → → | 全程 CAS 动态范围低于 amplitude cutoff；没有足够证据把小幅波动解释为明确动态程序 | “该类别整体保持相对稳定/仅表现有限波动” |
| `Heterogeneous` | 异质型 / 内部不一致型 | pathway 方向分散 | Category 中心轨迹可能发生变化，但支持总体方向的 pathway 比例不足 | “该类别内部通路响应方向不一致，提示功能异质性” |
| `Sustained activation` | 持续增强型 / 持续激活型 | ↑ ↑ ↑ | 多个连续 interval 的 CAS 总体持续向上 | “该功能类别随时间持续增强” |
| `Sustained suppression` | 持续减弱型 / 持续抑制型 | ↓ ↓ ↓ | 多个连续 interval 的 CAS 总体持续向下 | “该功能类别随时间持续减弱” |
| `Early activation` | 早期激活后回落型 | ↑ ↓ ↓ 或 ↑ ↓ → | 早期迅速升高，随后持续回落或趋稳；Peak 通常靠前 | “该类别表现为早期激活，随后逐渐回落” |
| `Transient activation` | 瞬时/中期激活型 | ↑ ↑ ↓ 或 → ↑ ↓ | 中间阶段形成明显高峰，随后下降 | “该类别出现暂时性增强，在中期达到峰值后回落” |
| `Late activation` | 晚期激活型 | →/↓ →/↑ ↑ | 前期变化有限或下降，主要增强发生在后期 | “该类别的增强主要出现在晚期” |
| `Early suppression` | 早期抑制后恢复型 | ↓ ↑ ↑ 或 ↓ ↑ → | 早期快速下降，随后恢复；Trough 通常靠前 | “该类别早期受到抑制，随后逐渐恢复” |
| `Transient suppression` | 瞬时/中期抑制型 | ↓ ↓ ↑ 或 → ↓ ↑ | 中间阶段形成低谷，随后回升 | “该类别出现暂时性抑制，在中期达到低谷后恢复” |
| `Late suppression` | 晚期抑制型 | →/↑ →/↓ ↓ | 前期相对稳定或增强，主要下降发生在后期 | “该类别的抑制主要发生于晚期” |
| `Fluctuating / biphasic` | 波动型 / 双相型 | ↑ ↓ ↑ 或 ↓ ↑ ↓ | 至少出现多次方向反转，无法用单一激活/抑制阶段概括 | “该类别呈双相/波动性动态，提示阶段依赖性调控” |

这里的“activation / suppression”指的是**该 Category 的 pathway activity score 相对增强/减弱**。它不自动等价于细胞功能真正被激活或被抑制，更不能自动推出因果机制。尤其对于名称包含 negative regulation 的 GO term，score 增强表示该 GO biological process 对应的基因集活性增强，而不能简单翻译成其下游生物学对象一定下降。

## 14.2 Stable、Heterogeneous 与 Fluctuating 的区别

这三个标签最容易混淆：

- **Stable**：重点是“幅度小”。Category 从头到尾没有足够大的总体动态变化。
- **Heterogeneous**：重点是“内部不一致”。Category 可以变化很大，但不同 pathway 并不共同支持同一方向。
- **Fluctuating / biphasic**：重点是“时间方向反转”。Category 内部可以很一致，但整个类别在不同阶段先升后降、再升，形成真实的阶段性动态。

因此，`Heterogeneous ≠ Stable`，`Fluctuating ≠ noise`。前者是 pathway 组成层面的异质性，后者是时间维度上的非单调性。

## 14.3 四时间点的方向序列如何理解

四个时间点产生三个 ΔCAS。例如：

```text
+,+,+ -> 持续增强
-,-,- -> 持续减弱
+,-,- -> 早期增强后持续回落
+,+,- -> 中期/偏晚峰值后回落
-,+,+ -> 早期抑制后恢复
-,-,+ -> 中期/偏晚低谷后恢复
+,-,+ -> 双相/波动
-,+,- -> 双相/波动
```

实际分类不能只看符号，因为非常小的 ΔCAS 会先被 epsilon 归为 Stable，而且低 Amplitude 或低 DCS 可能覆盖单纯的符号形态。

## 14.4 Pattern 不是“好/坏”或“强/弱”排序

Pattern 描述的是**形状**。变化强度主要看 Amplitude / ΔCAS，协同性主要看 DCS/NDI。例如两个 Category 都可属于 `Sustained activation`，但一个 Amplitude=0.30、DCS=0.62，另一个 Amplitude=1.10、DCS=0.95；后者是幅度更大、内部更协调的持续增强程序。

---

# 15. 自动阈值与手动阈值

## 自动

```r
res <- run_temporal_pathway_analysis(
  ...,
  thresholds = list(
    mode = "auto",
    pathway_quantile = 0.25,
    category_quantile = 0.25,
    amplitude_quantile = 0.25,
    dcs_cutoff = 0.60
  )
)
```

适合新数据集探索和不同 score 类型迁移。

## 手动锁定

完成诊断和 sensitivity analysis 后，正式分析建议锁定参数：

```r
res <- run_temporal_pathway_analysis(
  ...,
  thresholds = list(
    mode = "manual",
    epsilon_pathway = 0.18,
    epsilon_category = 0.12,
    amplitude_cutoff = 0.35,
    dcs_cutoff = 0.60
  )
)
```

数值仅为代码示例，不是推荐固定值。

**原则：先依据整体数据分布/技术噪声确定规则，再统一应用于所有 Category；不要逐 Category 调阈值。**

如果 CQ 与 N396 使用相同评分方法、相同实验设计，并且目标是直接比较两组 pattern，最好考虑基于两组联合分布确定一套统一 cutoff；如果两组数据尺度本身不可比，则应分别校准，并避免把分类频率做简单直接比较。

---

# 16. 结果对象

`run_temporal_pathway_analysis()` 返回：

```text
settings              实际参数与阈值
annotation            匹配后的 pathway annotation
score_subset           GO ID 行名的 score matrix
score_readable         Description 行名的可读矩阵
long                   样本级长表
pathway_trajectory     pathway × time 轨迹
CAS                    category × time 中心轨迹
pathway_delta          pathway 相邻时间变化
category_delta         CAS 相邻时间变化
direction_metrics      NDI / DCS / P_up / P_down / P_stable
global_metrics         rho / amplitude / peak / trough / net change
summary                最终 Category trend summary
diagnostics            阈值校准所需分布
```

最常使用：

```r
res$summary
res$CAS
res$direction_metrics
res$diagnostics
```

---

# 17. 推荐的分析顺序

1. 检查 pathway ID 与 Category mapping。
2. 检查 score matrix 与 annotation 匹配率。
3. 决定是否 pathway-wise Z-score。
4. 计算 pathway trajectory。
5. 计算 CAS。
6. 检查 |Δ pathway|、|ΔCAS|、Amplitude、DCS 分布。
7. 设定/校准 epsilon、amplitude cutoff、DCS cutoff。
8. 做 sensitivity analysis。
9. 锁定参数。
10. 生成最终 Pattern。
11. 结合 CAS heatmap、trajectory plot 和 pathway-level heatmap 解释生物学结果。

---

# 18. 如何实际阅读一个 Category 的结果

推荐按固定顺序判读，而不是只看 `Pattern` 一列。

### 第一步：看 n_pathway

确认这个 Category 有多少 pathway 支撑。pathway 数过少时，median、NDI 和 DCS 都会更加离散，分类可信度应降低。

### 第二步：看各时间点 CAS

先直接阅读中心轨迹的形状。例如 `[-0.60, -0.10, 0.35, 0.70]` 显示持续向上；不要因为前两个值为负就称其为“抑制”。

### 第三步：看 ΔCAS 与 Amplitude

ΔCAS 告诉你变化发生在哪个 interval；Amplitude 告诉你整段实验的动态范围是否足够大。

### 第四步：看 DCS_directional 和 NDI

判断这个中心轨迹是否代表 Category 中多数 pathway 的协同行为。高 DCS + 同方向 NDI 比单独 CAS 更支持“coordinated biological program”的描述。

### 第五步：看 Peak/Trough、Trend_rho 和 Net_Change

用它们补充时间位置、单调性以及最终是否保留变化。

### 第六步：最后读取 Pattern

Pattern 是前面多个信息的压缩标签，不应替代原始指标。正式结果描述最好同时报告 Pattern + 关键数值。

## 18.1 完整数值示例

假设某 Category：

```text
CAS:             R0=-0.55, R1=-0.12, R2=0.31, R3=0.62
ΔCAS:            +0.43, +0.43, +0.31
Amplitude:       1.17
Trend_rho:       1.00
Mean DCS:        0.82
Mean |NDI|:      0.71
Pattern:         Sustained activation
```

合理解释为：**“该功能类别从 R0 至 R3 呈持续增强趋势（CAS -0.55→0.62；Amplitude=1.17），三个连续时间区间均表现为正向变化，并且多数成员 pathway 与类别总体方向一致（平均 DCS=0.82），提示较为协调的持续性功能增强。”**

这里不能仅因为 R0 CAS=-0.55 就写“R0 时该功能被显著抑制”；CAS 是相对评分，而且本描述框架本身不提供显著性检验。

另一个 Category：

```text
CAS:             R0=0.05, R1=0.72, R2=0.81, R3=0.08
Amplitude:       0.76
Peak_Time:       R2
Net_Change:      0.03
Pattern:         Transient activation
```

其终点几乎回到起点，因此 Net Change 很小，但不能称为 Stable；较大的 Amplitude 和 R2 峰值说明它经历了明显的暂时性增强。

---

# 19. 统计学边界

这套规则是 **descriptive temporal pattern framework**，不是 time-effect hypothesis test。

如果需要声明“该 Category 随时间显著变化”，应另外在样本层面建立 inferential model，例如：

- linear/mixed-effects model；
- spline model；
- GAM；
- repeated-measures framework（若样本确实为纵向重复测量）。

随后把两层结论分开：

```text
Inferential model -> 是否存在可靠的 time effect？
本框架          -> time effect 呈现什么动态形态？
```

这样最适合论文 Methods 与 Results 的严谨表述。

---

# 20. GSEA 预筛选数据的解释

若输入 pathway 是“所有 pairwise GSEA 显著 pathway 的并集”，建议 Methods 表述为：

> GO-BP terms showing significant enrichment in at least one pairwise comparison were first collected as a union set. Sample-level pathway activity scores were subsequently used to characterize the temporal trajectories of these preselected pathways and their predefined biological categories.

避免把后续趋势分类再次描述为独立 pathway discovery，因为 pathway 集合已经受到前置 GSEA selection 的条件化。

---

# 21. 推荐报告的核心字段

最终结果表至少保留：

- Category
- n_pathway
- 各时间点 CAS
- 各 interval ΔCAS
- 各 interval NDI
- 各 interval DCS_directional
- Trend_rho
- Amplitude
- Peak_Time
- Trough_Time
- Net_Change
- Pattern
- 实际使用的 epsilon / amplitude / DCS cutoff

这样任何 Category 的文字描述都可以追溯到明确的数值规则，而不是主观读图。


---

# 22. 指标速查表

| 指标 | 核心问题 | 高/正值通常意味着 | 低/负值通常意味着 | 不能单独说明什么 |
|---|---|---|---|---|
| CAS | 当前 Category 中心处于什么水平？ | Z-score 模式下高于 pathway 自身平均状态 | 低于 pathway 自身平均状态 | 不能单独定义动态方向或显著性 |
| ΔCAS | 这个 interval Category 往哪变？ | Category 中心增强 | Category 中心减弱 | 不能说明内部 pathway 是否一致 |
| Amplitude | 全程变化有多大？ | 动态范围大 | 动态范围小 | 不提供方向 |
| Net Change | 终点相对起点留下多少变化？ | 最终高于起点 | 最终低于起点 | 会遗漏中间峰谷 |
| Trend_rho | 是否近似单调？ | 随时间总体上升 | 随时间总体下降 | rho≈0 不等于无变化 |
| NDI | Up/Down pathway 哪边数量占优？ | 上升 pathway 占优 | 下降 pathway 占优 | 接近 0 不一定稳定 |
| DCS_directional | 有多少 pathway 支持 Category 总方向？ | 内部协同性高 | 内部异质性高 | 不表示变化幅度 |
| P_stable | 有多少 pathway 未越过 epsilon？ | 多数 pathway 局部变化有限 | 更多 pathway 有方向变化 | 依赖 epsilon 设置 |
| Peak_Time | 最高活性何时出现？ | — | — | 不表示峰值是否显著 |
| Trough_Time | 最低活性何时出现？ | — | — | 不表示低谷是否显著 |
| Pattern | 整体轨迹属于什么形态？ | — | — | 是描述性标签，不是统计检验 |

**一句话记忆：CAS 看水平，ΔCAS 看阶段方向，Amplitude 看幅度，rho 看单调性，NDI 看净方向偏向，DCS 看内部协同性，Peak/Trough 看时间位置，Pattern 负责把这些信息压缩成可读的动态形态。**
