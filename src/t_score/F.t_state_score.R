# Interactive T-cell state GSVA analysis; functions only, no automatic execution.

t_score_read_tpm <- function(path, gene_column = "Gene") {
  if (!file.exists(path)) stop("TPM 文件不存在：", path, call. = FALSE)
  tab <- readr::read_tsv(path, show_col_types = FALSE, name_repair = "minimal")
  if (sum(names(tab) == gene_column) != 1L || ncol(tab) < 2L) stop("TPM 需要 Gene 和样本列。", call. = FALSE)
  genes <- as.character(tab[[gene_column]])
  samples <- names(tab)[names(tab) != gene_column]
  if (anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes) ||
      anyNA(samples) || any(!nzchar(samples)) || anyDuplicated(samples)) {
    stop("TPM 基因/样本 ID 不能为空或重复。", call. = FALSE)
  }
  mat <- as.matrix(tab[, samples, drop = FALSE])
  suppressWarnings(storage.mode(mat) <- "numeric")
  if (any(!is.finite(mat)) || any(mat < 0)) stop("TPM 必须是有限非负数。", call. = FALSE)
  rownames(mat) <- genes
  mat
}

t_score_read_counts <- function(path, gene_column = "Gene") {
  mat <- t_score_read_tpm(path, gene_column)
  if (any(mat != round(mat))) stop("counts 必须是原始整数计数。", call. = FALSE)
  mat
}

t_score_prepare_input <- function(tpm, metadata, gene_sets, groups = NULL,
                                  cycle_order = NULL, counts = NULL,
                                  count_min = 10, count_min_samples = 3L,
                                  min_tpm = 0,
                                  min_samples = 1L, min_set_size = 2L,
                                  max_set_size = 300L,
                                  sample_column = "Sample", group_column = "Group",
                                  cycle_column = "Cycle") {
  required <- c(sample_column, group_column, cycle_column)
  if (!is.data.frame(metadata) || !all(required %in% names(metadata))) {
    stop("metadata 需要 Sample、Group、Cycle 列（或相应的列名参数）。", call. = FALSE)
  }
  meta <- as.data.frame(metadata, stringsAsFactors = FALSE)
  for (column in required) {
    meta[[column]] <- as.character(meta[[column]])
    if (anyNA(meta[[column]]) || any(!nzchar(meta[[column]]))) stop("metadata ", column, " 有空值。", call. = FALSE)
  }
  if (anyDuplicated(meta[[sample_column]])) stop("metadata 样本重复。", call. = FALSE)
  if (!is.null(groups)) meta <- meta[meta[[group_column]] %in% groups, , drop = FALSE]
  if (!nrow(meta)) stop("所选 Group 没有样本。", call. = FALSE)
  group_order <- if (is.null(groups)) unique(meta[[group_column]]) else groups[groups %in% meta[[group_column]]]
  if (is.null(cycle_order)) cycle_order <- unique(meta[[cycle_column]])
  if (anyDuplicated(cycle_order) || !setequal(unique(meta[[cycle_column]]), cycle_order)) {
    stop("cycle_order 必须恰好包含所选样本中的所有 Cycle。", call. = FALSE)
  }
  meta <- meta[order(match(meta[[group_column]], group_order),
                     match(meta[[cycle_column]], cycle_order), meta[[sample_column]]), , drop = FALSE]
  rownames(meta) <- meta[[sample_column]]
  tpm <- as.matrix(tpm)
  if (is.null(rownames(tpm)) || is.null(colnames(tpm)) || anyNA(rownames(tpm)) ||
      anyNA(colnames(tpm)) || any(!nzchar(rownames(tpm))) || any(!nzchar(colnames(tpm))) ||
      anyDuplicated(rownames(tpm)) || anyDuplicated(colnames(tpm))) {
    stop("TPM 必须有非空且唯一的基因和样本名。", call. = FALSE)
  }
  suppressWarnings(storage.mode(tpm) <- "numeric")
  if (any(!is.finite(tpm)) || any(tpm < 0)) stop("TPM 必须是有限非负数。", call. = FALSE)
  missing <- setdiff(meta[[sample_column]], colnames(tpm))
  if (length(missing)) stop("TPM 缺少样本：", paste(missing, collapse = ", "), call. = FALSE)
  tpm <- tpm[, meta[[sample_column]], drop = FALSE]
  if (!is.finite(min_tpm) || min_tpm < 0 || min_samples < 1L ||
      !is.finite(count_min) || count_min < 0 || count_min_samples < 1L ||
      min_set_size < 1L || max_set_size < min_set_size) stop("过滤参数无效。", call. = FALSE)
  keep <- rowSums(tpm > min_tpm) >= min_samples
  if (!is.null(counts)) {
    counts <- as.matrix(counts)
    if (is.null(rownames(counts)) || is.null(colnames(counts)) ||
        anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) {
      stop("counts 必须有唯一的基因和样本名。", call. = FALSE)
    }
    suppressWarnings(storage.mode(counts) <- "numeric")
    if (any(!is.finite(counts)) || any(counts < 0) || any(counts != round(counts))) {
      stop("counts 必须是有限的非负整数。", call. = FALSE)
    }
    missing <- setdiff(meta[[sample_column]], colnames(counts))
    if (length(missing)) stop("counts 缺少样本：", paste(missing, collapse = ", "), call. = FALSE)
    count_genes <- intersect(rownames(tpm), rownames(counts))
    keep_count <- rowSums(counts[count_genes, meta[[sample_column]], drop = FALSE] > count_min) >= count_min_samples
    keep <- keep & rownames(tpm) %in% count_genes[keep_count]
  }
  tpm <- tpm[keep, , drop = FALSE]
  if (nrow(tpm) < 2L) stop("表达过滤后少于 2 个基因。", call. = FALSE)
  if (!is.list(gene_sets) || is.null(names(gene_sets)) || anyNA(names(gene_sets)) ||
      any(!nzchar(names(gene_sets))) || anyDuplicated(names(gene_sets))) {
    stop("gene_sets 必须是具备唯一非空状态名的基因向量 list。", call. = FALSE)
  }
  sets <- lapply(gene_sets, function(x) intersect(unique(as.character(x[!is.na(x) & nzchar(as.character(x))])), rownames(tpm)))
  coverage <- data.frame(State = names(sets), n_input = lengths(gene_sets),
                         n_mapped = lengths(sets), stringsAsFactors = FALSE)
  eligible <- coverage$n_mapped >= min_set_size & coverage$n_mapped <= max_set_size
  if (any(!eligible)) warning("以下状态基因集未达到映射大小要求：",
                              paste(coverage$State[!eligible], collapse = ", "), call. = FALSE)
  sets <- sets[eligible]
  if (length(sets) < 3L) stop("雷达图至少需要 3 个有效状态基因集；请检查基因 ID 和 min_set_size。", call. = FALSE)
  list(expression_log = log2(tpm + 1), metadata = meta, gene_sets = sets,
       coverage = coverage, group_order = group_order, cycle_order = cycle_order,
       columns = list(sample = sample_column, group = group_column, cycle = cycle_column),
       min_set_size = min_set_size, max_set_size = max_set_size)
}

t_score_gsva_matrix <- function(expression_log, gene_sets, min_set_size = 2L,
                                max_set_size = 300L, kcdf = "Gaussian") {
  if (!requireNamespace("GSVA", quietly = TRUE)) stop("T 状态打分需要安装 GSVA。", call. = FALSE)
  if (ncol(expression_log) < 2L) stop("GSVA 至少需要 2 个样本。", call. = FALSE)
  if ("gsvaParam" %in% getNamespaceExports("GSVA")) {
    pars <- GSVA::gsvaParam(exprData = expression_log, geneSets = gene_sets,
                            minSize = min_set_size, maxSize = max_set_size, kcdf = kcdf)
    scores <- GSVA::gsva(pars)
  } else {
    scores <- GSVA::gsva(expr = expression_log, gset.idx.list = gene_sets,
                         method = "gsva", min.sz = min_set_size,
                         max.sz = max_set_size, kcdf = kcdf, verbose = FALSE)
  }
  scores <- as.matrix(scores)
  if (!all(names(gene_sets) %in% rownames(scores)) ||
      !identical(colnames(scores), colnames(expression_log)) || any(!is.finite(scores))) {
    stop("GSVA 输出的状态/样本或分数无效。", call. = FALSE)
  }
  scores[names(gene_sets), , drop = FALSE]
}

t_score_run <- function(prepared, kcdf = "Gaussian", score_function = t_score_gsva_matrix) {
  expr <- prepared$expression_log
  meta <- prepared$metadata
  col <- prepared$columns
  args <- list(gene_sets = prepared$gene_sets, min_set_size = prepared$min_set_size,
               max_set_size = prepared$max_set_size, kcdf = kcdf)
  run <- function(x) do.call(score_function, c(list(expression_log = x), args))
  all_samples <- run(expr)
  by_group <- setNames(lapply(prepared$group_order, function(group) {
    samples <- meta[[col$sample]][meta[[col$group]] == group]
    if (length(samples) < 2L) stop("Group ", group, " 少于 2 个样本，无法单独 GSVA。", call. = FALSE)
    run(expr[, samples, drop = FALSE])
  }), prepared$group_order)
  list(all_samples = all_samples, by_group = by_group, metadata = meta,
       gene_sets = prepared$gene_sets, coverage = prepared$coverage,
       group_order = prepared$group_order, cycle_order = prepared$cycle_order,
       columns = col)
}

t_score_long <- function(scores, result) {
  scores <- as.matrix(scores)
  sample_col <- result$columns$sample
  if (!all(colnames(scores) %in% result$metadata[[sample_col]])) stop("score matrix 中有 metadata 未知样本。", call. = FALSE)
  df <- data.frame(Sample = rep(colnames(scores), each = nrow(scores)),
                   State = rep(rownames(scores), times = ncol(scores)),
                   score = as.vector(scores), stringsAsFactors = FALSE)
  meta <- result$metadata
  meta$Sample <- meta[[sample_col]]
  meta <- meta[match(df$Sample, meta$Sample), , drop = FALSE]
  df$Group <- meta[[result$columns$group]]
  df$Cycle <- meta[[result$columns$cycle]]
  df
}

t_score_cycle_means <- function(scores, result) {
  df <- t_score_long(scores, result)
  aggregate(score ~ Group + Cycle + State, data = df, FUN = mean)
}

.t_score_colors <- function(values, levels, palette) {
  if (is.null(values)) values <- setNames(grDevices::hcl.colors(length(levels), palette), levels)
  if (is.null(names(values)) || !all(levels %in% names(values))) {
    stop("颜色向量必须按 Group/Cycle 名称命名并覆盖全部所需水平。", call. = FALSE)
  }
  values <- values[levels]
  tryCatch(grDevices::col2rgb(values), error = function(e) {
    stop("颜色向量包含无效颜色：", conditionMessage(e), call. = FALSE)
  })
  values
}

.t_score_radar <- function(means, states, line_var, facet_var = NULL,
                           line_order, colors, limits = NULL, title = NULL) {
  if (length(states) < 3L) stop("雷达图至少需要 3 个状态。", call. = FALSE)
  means <- means[means$State %in% states, , drop = FALSE]
  if (!nrow(means)) stop("雷达图无可绘制数据。", call. = FALSE)
  keys <- if (is.null(facet_var)) line_var else c(facet_var, line_var)
  splits <- split(means, interaction(means[keys], drop = TRUE))
  closed <- do.call(rbind, lapply(splits, function(df) {
    missing <- setdiff(states, df$State)
    if (length(missing)) stop("雷达图缺少状态：", paste(missing, collapse = ", "), call. = FALSE)
    df <- df[match(states, df$State), , drop = FALSE]
    df$.axis <- seq_along(states)
    last <- df[1, , drop = FALSE]
    last$.axis <- length(states) + 1L
    rbind(df, last)
  }))
  closed$.line <- factor(closed[[line_var]], levels = line_order)
  closed$.series <- interaction(closed[keys], drop = TRUE)
  if (is.null(limits)) {
    limits <- range(means$score, 0, na.rm = TRUE)
    if (limits[1] == limits[2]) limits <- limits + c(-0.5, 0.5)
  }
  if (length(limits) != 2L || any(!is.finite(limits)) || limits[1] >= limits[2] ||
      any(means$score < limits[1] | means$score > limits[2])) {
    stop("limits 必须覆盖所有分数且下限小于上限。", call. = FALSE)
  }
  span <- diff(limits)
  pad <- span * 0.08
  closed$.radius <- closed$score - limits[1] + pad
  ticks <- pretty(limits, n = 4)
  ticks <- ticks[ticks >= limits[1] & ticks <= limits[2]]
  p <- ggplot2::ggplot(closed, ggplot2::aes(x = .axis, y = .radius,
                                            group = .series, color = .line, fill = .line)) +
    ggplot2::geom_polygon(alpha = 0.10, linewidth = 0.7) +
    ggplot2::geom_point(data = closed[closed$.axis <= length(states), , drop = FALSE], size = 1.8) +
    ggplot2::scale_x_continuous(breaks = seq_along(states), labels = states,
                                limits = c(1, length(states) + 1), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, span + 2 * pad),
                                breaks = ticks - limits[1] + pad,
                                labels = format(ticks, trim = TRUE), expand = c(0, 0)) +
    ggplot2::scale_color_manual(values = colors, name = line_var) +
    ggplot2::scale_fill_manual(values = colors, name = line_var) +
    ggplot2::coord_polar(theta = "x", start = 0) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.background = ggplot2::element_rect(fill = "white", color = NA),
                   panel.background = ggplot2::element_rect(fill = "white", color = NA),
                   axis.text.x = ggplot2::element_text(color = "black"),
                   axis.text.y = ggplot2::element_text(color = "grey30"),
                   legend.position = "bottom")
  if (!is.null(facet_var)) p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_var)))
  p
}

.t_score_reference_radar <- function(means, states, line_var, facet_var = NULL,
                                     line_order, colors, limits = NULL, title = NULL) {
  if (length(states) < 3L || anyDuplicated(states) || !all(states %in% means$State)) {
    stop("参考样式雷达图需要至少 3 个不重复且存在的状态。", call. = FALSE)
  }
  means <- means[means$State %in% states, , drop = FALSE]
  if (is.null(limits)) limits <- c(-0.7, 0.7)
  if (length(limits) != 2L || any(!is.finite(limits)) || limits[1] >= limits[2]) {
    stop("limits 必须是长度为 2 的递增数值向量。", call. = FALSE)
  }
  if (any(means$score < limits[1] | means$score > limits[2])) {
    stop("GSVA 分数超出雷达图固定范围；请扩大 limits，避免图形截断。", call. = FALSE)
  }
  structure(list(means = means, states = states, line_var = line_var,
                 facet_var = facet_var, line_order = line_order, colors = colors,
                 limits = limits, title = title), class = "t_score_reference_radar")
}

plot.t_score_reference_radar <- function(x, ...) {
  if (!requireNamespace("fmsb", quietly = TRUE)) {
    stop("参考样式雷达图需要安装 fmsb；如暂未安装，可用 style='polar'。", call. = FALSE)
  }
  facet_values <- if (is.null(x$facet_var)) NULL else {
    values <- x$means[[x$facet_var]]
    if (is.factor(values)) levels(values)[levels(values) %in% as.character(values)] else unique(as.character(values))
  }
  if (!is.null(facet_values)) {
    ncol <- ceiling(sqrt(length(facet_values)))
    nrow <- ceiling(length(facet_values) / ncol)
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mfrow = c(nrow, ncol), mar = c(2, 2, 3, 2))
  }
  panels <- if (is.null(facet_values)) list(x$means) else
    lapply(facet_values, function(value) x$means[as.character(x$means[[x$facet_var]]) == value, , drop = FALSE])
  for (i in seq_along(panels)) {
    panel <- panels[[i]]
    lines <- x$line_order[x$line_order %in% as.character(panel[[x$line_var]])]
    if (!length(lines)) stop("雷达图面板没有可绘制的系列。", call. = FALSE)
    series <- t(vapply(lines, function(line) {
      values <- panel$score[as.character(panel[[x$line_var]]) == line]
      state_names <- panel$State[as.character(panel[[x$line_var]]) == line]
      if (anyDuplicated(state_names) || !setequal(state_names, x$states)) {
        stop("雷达图每个系列必须完整且唯一地包含全部状态。", call. = FALSE)
      }
      as.numeric(values[match(x$states, state_names)])
    }, numeric(length(x$states))))
    radar_df <- as.data.frame(rbind(
      Max = rep(x$limits[2], length(x$states)),
      Min = rep(x$limits[1], length(x$states)),
      series
    ), check.names = FALSE)
    names(radar_df) <- x$states
    rownames(radar_df) <- c("Max", "Min", lines)
    line_cols <- unname(x$colors[lines])
    fmsb::radarchart(
      radar_df,
      axistype = 0,
      seg = 4,
      pcol = line_cols,
      plwd = 2.5,
      plty = 1,
      pfcol = grDevices::adjustcolor(line_cols, alpha.f = 0.18),
      cglcol = "grey80",
      cglty = 1,
      cglwd = 0.8,
      vlcex = 0.9,
      title = if (is.null(facet_values)) x$title else facet_values[[i]]
    )
    graphics::legend("topright", legend = lines, col = line_cols,
                     lty = 1, lwd = 2.5, bty = "n", cex = 0.8)
  }
  invisible(x)
}

print.t_score_reference_radar <- function(x, ...) {
  plot(x, ...)
  invisible(x)
}

t_score_save_radar <- function(plot, filename, width, height, dpi = 300) {
  if (inherits(plot, "ggplot")) {
    ggplot2::ggsave(filename, plot = plot, width = width, height = height, units = "in", dpi = dpi)
    return(invisible(filename))
  }
  if (!inherits(plot, "t_score_reference_radar")) stop("plot 不是受支持的雷达图对象。", call. = FALSE)
  if (!requireNamespace("fmsb", quietly = TRUE)) stop("参考样式雷达图需要安装 fmsb。", call. = FALSE)
  ext <- tolower(tools::file_ext(filename))
  if (ext == "pdf") {
    grDevices::pdf(filename, width = width, height = height, onefile = FALSE)
  } else if (ext == "png") {
    grDevices::png(filename, width = width, height = height, units = "in", res = dpi, bg = "white")
  } else {
    stop("参考样式雷达图仅支持保存为 PDF 或 PNG。", call. = FALSE)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::plot(plot)
  invisible(filename)
}

t_score_plot_group_radar <- function(result, group, cycle_colors = NULL,
                                     states = rownames(result$all_samples), limits = NULL,
                                     style = c("reference", "polar")) {
  style <- match.arg(style)
  if (!group %in% names(result$by_group)) stop("未找到 Group：", group, call. = FALSE)
  scores <- result$by_group[[group]]
  means <- t_score_cycle_means(scores, result)
  cycles <- result$cycle_order[result$cycle_order %in% means$Cycle]
  colors <- .t_score_colors(cycle_colors, cycles, "Dark 3")
  renderer <- if (style == "reference") .t_score_reference_radar else .t_score_radar
  renderer(means, states, "Cycle", line_order = cycles,
           colors = colors, limits = limits, title = group)
}

t_score_plot_faceted_radar <- function(result, group_colors = NULL,
                                       states = rownames(result$all_samples), limits = NULL,
                                       style = c("reference", "polar")) {
  style <- match.arg(style)
  means <- t_score_cycle_means(result$all_samples, result)
  means$Cycle <- factor(means$Cycle, levels = result$cycle_order)
  colors <- .t_score_colors(group_colors, result$group_order, "Dark 3")
  renderer <- if (style == "reference") .t_score_reference_radar else .t_score_radar
  renderer(means, states, "Group", facet_var = "Cycle",
           line_order = result$group_order, colors = colors, limits = limits)
}

t_score_plot_trajectories <- function(result, states = rownames(result$all_samples),
                                      groups = result$group_order, ncol = 3L,
                                      group_colors = NULL, error = c("sem", "sd", "none"),
                                      show_samples = TRUE) {
  error <- match.arg(error)
  if (!all(states %in% rownames(result$all_samples)) || !all(groups %in% result$group_order)) {
    stop("states/groups 包含未知名称。", call. = FALSE)
  }
  if (ncol < 1L) stop("ncol 至少为 1。", call. = FALSE)
  df <- t_score_long(result$all_samples, result)
  df <- df[df$State %in% states & df$Group %in% groups, , drop = FALSE]
  if (!nrow(df)) stop("轨迹图无可绘制数据。", call. = FALSE)
  df$State <- factor(df$State, levels = states)
  df$Group <- factor(df$Group, levels = groups)
  df$Cycle <- factor(df$Cycle, levels = result$cycle_order)
  df$cycle_index <- as.integer(df$Cycle)
  key <- interaction(df$State, df$Group, df$Cycle, drop = TRUE)
  summary <- do.call(rbind, lapply(split(df, key), function(x) {
    data.frame(State = x$State[1], Group = x$Group[1], Cycle = x$Cycle[1],
               cycle_index = x$cycle_index[1], n = nrow(x), mean = mean(x$score),
               sd = if (nrow(x) > 1L) stats::sd(x$score) else NA_real_)
  }))
  summary$State <- factor(as.character(summary$State), levels = states)
  summary$Group <- factor(as.character(summary$Group), levels = groups)
  summary$err <- if (error == "sem") summary$sd / sqrt(summary$n) else if (error == "sd") summary$sd else 0
  summary$err[is.na(summary$err)] <- 0
  colors <- .t_score_colors(group_colors, groups, "Dark 3")
  p <- ggplot2::ggplot() + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey70")
  if (error != "none") {
    p <- p + ggplot2::geom_ribbon(data = summary,
          ggplot2::aes(x = cycle_index, ymin = mean - err, ymax = mean + err,
                       fill = Group, group = Group), alpha = 0.15, show.legend = FALSE)
  }
  p <- p + ggplot2::geom_line(data = summary,
          ggplot2::aes(x = cycle_index, y = mean, color = Group, group = Group), linewidth = 0.85) +
    ggplot2::geom_point(data = summary,
          ggplot2::aes(x = cycle_index, y = mean, color = Group, shape = Group), size = 2.4)
  if (show_samples) {
    p <- p + ggplot2::geom_point(data = df,
      ggplot2::aes(x = cycle_index, y = score, color = Group, shape = Group),
      alpha = 0.45, size = 1.3,
      position = ggplot2::position_jitter(width = 0.06, height = 0, seed = 123), show.legend = FALSE)
  }
  p + ggplot2::scale_x_continuous(breaks = seq_along(result$cycle_order),
                                 labels = result$cycle_order) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::facet_wrap(~State, ncol = ncol, scales = "free_y") +
    ggplot2::labs(x = "Cycle", y = "GSVA score", color = "Group", shape = "Group") +
    ggplot2::theme_classic(base_size = 11)
}
