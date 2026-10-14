# Configurable summary plots for DESeq2, trend and GSEA outputs.
#
# This stage is deliberately read-only with respect to analysis results: it
# consumes completed DE/trend/enrichment files and only writes figures plus a
# plotting manifest.  It can therefore be re-run freely after style changes.

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
}

result_plot_abort <- function(...) stop(paste0(...), call. = FALSE)

require_result_plot_packages <- function(config = list()) {
  load_configured_ggplot2(config)
  packages <- c("ggplot2", "patchwork")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    result_plot_abort("结果绘图缺少 R 包：", paste(missing, collapse = ", "))
  }
}

result_plot_number <- function(value, label, lower = -Inf, strict = FALSE) {
  value <- suppressWarnings(as.numeric(value))
  if (length(value) != 1L || !is.finite(value) ||
      if (strict) value <= lower else value < lower) {
    result_plot_abort(label, " 必须是", if (strict) "大于" else "大于等于", " ", lower, " 的有限数值。")
  }
  value
}

result_plot_colour <- function(colours, name, fallback) {
  value <- colours[[name]] %||% fallback
  value <- as.character(unlist(value, use.names = FALSE))
  if (length(value) != 1L || is.na(value) || !nzchar(value)) {
    result_plot_abort("颜色配置 ", name, " 必须是一个非空颜色值。")
  }
  value
}

result_plot_safe_rel_path <- function(project_root, path) {
  if (is.na(path) || !nzchar(path)) return(NA_character_)
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) return(path)
  file.path(project_root, path)
}

result_plot_size <- function(save_config, name, defaults) {
  spec <- save_config[[name]] %||% list()
  width <- result_plot_number(spec$width %||% defaults[[1]], paste0("save.", name, ".width"), 0, TRUE)
  height <- result_plot_number(spec$height %||% defaults[[2]], paste0("save.", name, ".height"), 0, TRUE)
  dpi <- result_plot_number(spec$dpi %||% save_config$dpi %||% 300, paste0("save.", name, ".dpi"), 0, TRUE)
  list(width = width, height = height, dpi = dpi)
}

save_result_plot <- function(plot, output_base, size, save_config = list(), formats_override = NULL) {
  formats <- tolower(as.character(unlist(formats_override %||% save_config$formats %||% c("pdf", "png"), use.names = FALSE)))
  formats <- unique(formats[nzchar(formats)])
  unsupported <- setdiff(formats, c("pdf", "png"))
  if (length(unsupported)) {
    result_plot_abort("save.formats 仅支持 pdf、png；不支持：", paste(unsupported, collapse = ", "))
  }
  if (!length(formats)) result_plot_abort("save.formats 至少应包含 pdf 或 png。")
  ensure_dir(dirname(output_base))
  outputs <- file.path(dirname(output_base), paste0(basename(output_base), ".", formats))
  for (index in seq_along(formats)) {
    if (formats[[index]] == "pdf") {
      if (inherits(plot, "ggbreak")) {
        # ggbreak's grid-draw method creates its own page.  ggsave() creates
        # one as well, resulting in a blank/extra PDF page on some devices.
        # Open the PDF device directly and let print.ggbreak create its sole
        # page instead.
        grDevices::pdf(outputs[[index]], width = size$width, height = size$height, useDingbats = FALSE)
        tryCatch(print(plot), finally = grDevices::dev.off())
      } else {
        ggplot2::ggsave(
          outputs[[index]], plot, width = size$width, height = size$height,
          units = "in", device = "pdf", useDingbats = FALSE, limitsize = FALSE
        )
      }
    } else {
      dimensions_px <- c(width = size$width, height = size$height) * size$dpi
      if (any(dimensions_px > 50000)) {
        result_plot_abort(
          "PNG 输出尺寸为 ",
          format(round(dimensions_px[["width"]]), big.mark = ",", trim = TRUE), " × ",
          format(round(dimensions_px[["height"]]), big.mark = ",", trim = TRUE),
          " px，超过 ragg 默认的 50,000 px 单边上限。",
          "请降低 save.width/save.height 或 dpi；对于 all 图，建议只保留 PDF（",
          "save.by_plot.all.formats: [\"pdf\"]）。"
        )
      }
      ggplot2::ggsave(
        outputs[[index]], plot, width = size$width, height = size$height,
        units = "in", dpi = size$dpi, limitsize = FALSE
      )
    }
  }
  outputs
}

read_successful_de_results_for_plot <- function(project_root) {
  manifest_path <- file.path(project_root, "data", "pipeline", "de_manifest.tsv")
  if (!file.exists(manifest_path)) {
    result_plot_abort("找不到 DE manifest：", manifest_path, "；请先完成 DESeq2 步骤。")
  }
  manifest <- read_tsv(manifest_path, c("group_id", "comparison_id", "status", "result_file"))
  manifest <- manifest[manifest$status == "success", , drop = FALSE]
  if (!nrow(manifest)) result_plot_abort("DE manifest 中没有 status=success 的比较。")
  output <- vector("list", nrow(manifest))
  for (index in seq_len(nrow(manifest))) {
    result_path <- result_plot_safe_rel_path(project_root, manifest$result_file[[index]])
    if (!file.exists(result_path)) {
      result_plot_abort("DE manifest 指向的结果文件不存在：", result_path)
    }
    table <- read_tsv(result_path, c("Gene", "log2FoldChange", "padj"))
    output[[index]] <- list(
      group_id = as.character(manifest$group_id[[index]]),
      comparison_id = as.character(manifest$comparison_id[[index]]),
      result = table
    )
  }
  output
}

classify_de_for_plot <- function(de_result, fc_cutoff, padj_cutoff, padj_floor = 1e-300) {
  log2fc_cutoff <- log2(fc_cutoff)
  if (!is.finite(log2fc_cutoff) || log2fc_cutoff < 0) {
    result_plot_abort("fc_cutoff 必须大于等于 1。")
  }
  if (!all(c("log2FoldChange", "padj") %in% names(de_result))) {
    result_plot_abort("DE 结果缺少 log2FoldChange 或 padj 列。")
  }
  out <- as.data.frame(de_result, stringsAsFactors = FALSE)
  out$log2FoldChange <- suppressWarnings(as.numeric(out$log2FoldChange))
  out$padj <- suppressWarnings(as.numeric(out$padj))
  out <- out[is.finite(out$log2FoldChange), , drop = FALSE]
  padj_plot <- out$padj
  padj_plot[!is.finite(padj_plot) | is.na(padj_plot)] <- 1
  padj_plot[padj_plot <= 0] <- padj_floor
  padj_plot[padj_plot > 1] <- 1
  out$neg_log10_padj <- -log10(padj_plot)
  out$Status <- ifelse(
    !is.na(out$padj) & is.finite(out$padj) & out$padj <= padj_cutoff & out$log2FoldChange >= log2fc_cutoff,
    "Up",
    ifelse(
      !is.na(out$padj) & is.finite(out$padj) & out$padj <= padj_cutoff & out$log2FoldChange <= -log2fc_cutoff,
      "Down", "NotSig"
    )
  )
  # Draw non-significant points first so the significant points remain visible.
  out <- out[order(match(out$Status, c("NotSig", "Down", "Up"))), , drop = FALSE]
  out$Status <- factor(out$Status, levels = c("Up", "Down", "NotSig"))
  out
}

validate_manual_limits <- function(limits, label, symmetric_x = FALSE) {
  if (is.null(limits) || !length(limits)) return(NULL)
  limits <- suppressWarnings(as.numeric(unlist(limits, use.names = FALSE)))
  if (length(limits) != 2L || any(!is.finite(limits)) || limits[[1]] >= limits[[2]]) {
    result_plot_abort(label, " 必须是递增的两个有限数值。")
  }
  if (symmetric_x) {
    maximum <- max(abs(limits))
    return(c(-maximum, maximum))
  }
  limits
}

calculate_volcano_limits <- function(plot_tables, volcano_config) {
  manual_x <- validate_manual_limits(volcano_config$x_limits %||% NULL, "volcano.x_limits", symmetric_x = TRUE)
  manual_y <- validate_manual_limits(volcano_config$y_limits %||% NULL, "volcano.y_limits")
  if (!is.null(manual_x) && !is.null(manual_y)) return(list(x = manual_x, y = manual_y))
  all_lfc <- unlist(lapply(plot_tables, function(x) x$log2FoldChange), use.names = FALSE)
  all_y <- unlist(lapply(plot_tables, function(x) x$neg_log10_padj), use.names = FALSE)
  x_max <- ceiling(max(abs(all_lfc[is.finite(all_lfc)]), na.rm = TRUE))
  y_max <- ceiling(max(all_y[is.finite(all_y)], na.rm = TRUE))
  if (!is.finite(x_max) || x_max == 0) x_max <- 1
  if (!is.finite(y_max) || y_max == 0) y_max <- 1
  list(
    x = manual_x %||% c(-x_max, x_max),
    y = manual_y %||% c(0, y_max * 1.05)
  )
}

get_named_config <- function(config, name) {
  if (is.null(config) || !length(config)) return(NULL)
  config[[as.character(name)]] %||% NULL
}

resolve_volcano_label_positions <- function(volcano_config, group_id, comparison_id, x_limits, y_limits) {
  configured <- volcano_config$label_positions %||% list()
  group_positions <- get_named_config(configured, group_id) %||% list()
  comparison_positions <- get_named_config(group_positions, comparison_id) %||% list()
  parse_position <- function(value, label) {
    if (is.null(value) || !length(value)) return(NULL)
    value <- suppressWarnings(as.numeric(unlist(value, use.names = FALSE)))
    if (length(value) != 2L || any(!is.finite(value))) {
      result_plot_abort(label, " 必须为 [x, y] 两个有限数值。")
    }
    value
  }
  x_max <- max(abs(x_limits))
  y_max <- max(y_limits)
  list(
    up = parse_position(comparison_positions$up, "volcano.label_positions") %||% c(x_max * 0.75, y_max * 0.95),
    down = parse_position(comparison_positions$down, "volcano.label_positions") %||% c(-x_max * 0.65, y_max * 0.95)
  )
}

make_volcano_plot <- function(plot_table, group_id, comparison_id, thresholds, limits, volcano_config) {
  colours <- volcano_config$colors %||% list()
  colour_map <- c(
    Up = result_plot_colour(colours, "up", "#D73027"),
    Down = result_plot_colour(colours, "down", "#4575B4"),
    NotSig = result_plot_colour(colours, "not_significant", "grey70")
  )
  positions <- resolve_volcano_label_positions(volcano_config, group_id, comparison_id, limits$x, limits$y)
  n_up <- sum(plot_table$Status == "Up", na.rm = TRUE)
  n_down <- sum(plot_table$Status == "Down", na.rm = TRUE)
  point_size <- result_plot_number(volcano_config$point_size %||% 0.8, "volcano.point_size", 0, TRUE)
  point_alpha <- result_plot_number(volcano_config$point_alpha %||% 0.7, "volcano.point_alpha", 0)
  if (point_alpha > 1) result_plot_abort("volcano.point_alpha 必须小于等于 1。")
  label_size <- result_plot_number(volcano_config$label_size %||% 5, "volcano.label_size", 0, TRUE)
  show_legend <- isTRUE(volcano_config$show_individual_legend %||% TRUE)
  title <- paste0(group_id, ": ", comparison_id)

  ggplot2::ggplot(plot_table, ggplot2::aes(x = .data$log2FoldChange, y = .data$neg_log10_padj, colour = .data$Status)) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::scale_colour_manual(values = colour_map, limits = c("Up", "Down", "NotSig"), drop = FALSE) +
    ggplot2::geom_vline(
      xintercept = c(-thresholds$log2fc_cutoff, thresholds$log2fc_cutoff),
      linetype = "dashed", color = "grey40"
    ) +
    ggplot2::geom_hline(yintercept = -log10(thresholds$padj_cutoff), linetype = "dashed", color = "grey40") +
    ggplot2::coord_cartesian(xlim = limits$x, ylim = limits$y) +
    ggplot2::labs(x = "log2 Fold Change", y = "-log10(adjusted p-value)", title = title, colour = NULL) +
    ggplot2::annotate("text", x = positions$up[[1]], y = positions$up[[2]], label = paste0("Up: ", n_up),
                      color = colour_map[["Up"]], size = label_size, fontface = "bold") +
    ggplot2::annotate("text", x = positions$down[[1]], y = positions$down[[2]], label = paste0("Down: ", n_down),
                      color = colour_map[["Down"]], size = label_size, fontface = "bold") +
    ggplot2::theme_classic(base_size = result_plot_number(volcano_config$base_size %||% 11, "volcano.base_size", 0, TRUE)) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = if (show_legend) "right" else "none"
    )
}

resolve_comparison_layout <- function(plot_config, group_id, comparison_ids) {
  layout <- get_named_config(plot_config$comparison_layout_by_group %||% list(), group_id) %||%
    plot_config$comparison_layout %||% list()
  comparison_ids <- unique(as.character(comparison_ids))
  positions <- layout$positions %||% list()
  if (!length(positions)) {
    ncol <- min(3L, max(1L, length(comparison_ids)))
    nrow <- ceiling(length(comparison_ids) / ncol)
    cells <- stats::setNames(seq_along(comparison_ids), comparison_ids)
    return(list(nrow = nrow, ncol = ncol, cells = cells, comparison_order = comparison_ids))
  }
  nrow <- suppressWarnings(as.integer(layout$nrow %||% NA_integer_))
  ncol <- suppressWarnings(as.integer(layout$ncol %||% NA_integer_))
  if (length(nrow) != 1L || is.na(nrow) || nrow < 1L || length(ncol) != 1L || is.na(ncol) || ncol < 1L) {
    result_plot_abort("comparison_layout 配置 positions 时必须同时给出正整数 nrow 和 ncol。")
  }
  missing <- setdiff(comparison_ids, names(positions))
  if (length(missing)) {
    result_plot_abort("Group '", group_id, "' 的 comparison_layout 缺少位置：", paste(missing, collapse = ", "))
  }
  cells <- vapply(comparison_ids, function(comparison_id) {
    position <- suppressWarnings(as.integer(unlist(positions[[comparison_id]], use.names = FALSE)))
    if (length(position) != 2L || any(is.na(position)) || any(position < 1L) ||
        position[[1]] > nrow || position[[2]] > ncol) {
      result_plot_abort("comparison_layout.positions.", comparison_id,
                        " 必须为 nrow/ncol 范围内的 [row, column]。")
    }
    (position[[1]] - 1L) * ncol + position[[2]]
  }, integer(1))
  if (anyDuplicated(cells)) result_plot_abort("comparison_layout 中多个 comparison 占用了同一位置。")
  cells <- stats::setNames(cells, comparison_ids)
  list(nrow = nrow, ncol = ncol, cells = cells, comparison_order = names(sort(cells)))
}

compose_comparison_grid <- function(plot_by_comparison, layout, collect_guides = FALSE) {
  panels <- replicate(layout$nrow * layout$ncol, patchwork::plot_spacer(), simplify = FALSE)
  for (comparison_id in names(plot_by_comparison)) {
    panels[[layout$cells[[comparison_id]]]] <- plot_by_comparison[[comparison_id]]
  }
  combined <- patchwork::wrap_plots(panels, ncol = layout$ncol)
  if (collect_guides) combined <- combined + patchwork::plot_layout(guides = "collect")
  combined
}

read_successful_trend_sets_for_plot <- function(project_root, selected_trend_ids = character()) {
  manifest_path <- file.path(project_root, "data", "pipeline", "trend_manifest.tsv")
  if (!file.exists(manifest_path)) {
    result_plot_abort("找不到趋势 manifest：", manifest_path, "；请先完成趋势统计步骤。")
  }
  manifest <- read_tsv(manifest_path, c("group_id", "trend_id", "status", "trend_summary_file"))
  manifest <- manifest[manifest$status == "success", , drop = FALSE]
  if (length(selected_trend_ids)) manifest <- manifest[manifest$trend_id %in% selected_trend_ids, , drop = FALSE]
  if (!nrow(manifest)) result_plot_abort("没有符合配置的成功趋势统计结果。")
  lapply(seq_len(nrow(manifest)), function(index) {
    path <- result_plot_safe_rel_path(project_root, manifest$trend_summary_file[[index]])
    if (!file.exists(path)) result_plot_abort("趋势 manifest 指向的文件不存在：", path)
    list(
      group_id = as.character(manifest$group_id[[index]]),
      trend_id = as.character(manifest$trend_id[[index]]),
      table = read_tsv(path, c("Full_Trend", "Count"))
    )
  })
}

resolve_trend_break <- function(trend_config, group_id, trend_id) {
  values <- get_named_config(get_named_config(trend_config$y_breaks %||% list(), group_id) %||% list(), trend_id)
  if (is.null(values) || !length(values)) return(NULL)
  values <- suppressWarnings(as.numeric(unlist(values, use.names = FALSE)))
  if (length(values) != 2L || any(!is.finite(values)) || values[[1]] >= values[[2]]) {
    result_plot_abort("trend_bar.y_breaks.", group_id, ".", trend_id, " 必须为递增的两个有限数值。")
  }
  values
}

make_trend_bar_plot <- function(summary_table, group_id, trend_id, trend_config, fill_colour) {
  table <- as.data.frame(summary_table, stringsAsFactors = FALSE)
  table$Count <- suppressWarnings(as.numeric(table$Count))
  if (any(!is.finite(table$Count)) || any(table$Count < 0)) {
    result_plot_abort("趋势结果 ", group_id, "/", trend_id, " 的 Count 必须为非负有限数值。")
  }
  table <- table[order(-table$Count, table$Full_Trend), , drop = FALSE]
  table$Full_Trend <- factor(table$Full_Trend, levels = table$Full_Trend)
  show_count <- isTRUE(trend_config$show_count %||% TRUE)
  top_padding <- result_plot_number(
    trend_config$label_top_padding %||% 0.10,
    "trend_bar.label_top_padding", 0, TRUE
  )
  y_max <- max(table$Count, na.rm = TRUE)
  if (!is.finite(y_max) || y_max == 0) y_max <- 1
  # Explicit scale headroom is retained after scale_y_break(), unlike a plain
  # expand() setting.  It prevents the label over the tallest bar being cut.
  y_upper <- y_max * (1 + top_padding)
  plot <- ggplot2::ggplot(table, ggplot2::aes(x = .data$Full_Trend, y = .data$Count)) +
    ggplot2::geom_col(width = 0.72, fill = fill_colour, color = "black", linewidth = 0.25) +
    ggplot2::labs(title = paste0(group_id, ": ", trend_id, " trends"), x = NULL, y = "Feature count") +
    ggplot2::scale_y_continuous(limits = c(0, y_upper), expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::theme_classic(base_size = result_plot_number(trend_config$base_size %||% 12, "trend_bar.base_size", 0, TRUE)) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 9),
      axis.text.y = ggplot2::element_text(color = "black"),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
  if (show_count) plot <- plot + ggplot2::geom_text(ggplot2::aes(label = .data$Count), vjust = -0.35, size = 3)
  y_break <- resolve_trend_break(trend_config, group_id, trend_id)
  if (!is.null(y_break)) {
    if (y_break[[2]] >= y_upper) {
      result_plot_abort(
        "trend_bar.y_breaks.", group_id, ".", trend_id,
        " 的上界必须小于最高柱加顶部标签留白后的范围（当前上界为 ",
        format(y_upper, trim = TRUE), "）。"
      )
    }
    if (!requireNamespace("ggbreak", quietly = TRUE)) {
      result_plot_abort("配置了 trend_bar.y_breaks，但缺少 R 包 'ggbreak'。")
    }
    plot <- plot + ggbreak::scale_y_break(
      breaks = y_break,
      scales = result_plot_number(trend_config$y_break_scales %||% 0.35, "trend_bar.y_break_scales", 0, TRUE),
      space = 0.08
    )
  }
  plot
}

read_gsea_table_for_plot <- function(project_root, group_id, comparison_id) {
  directory <- file.path(project_root, "data", "enrichment_results", safe_path_component(group_id), safe_path_component(comparison_id))
  rds_path <- file.path(directory, "gsea_go_bp.rds")
  tsv_path <- file.path(directory, "gsea_go_bp.tsv")
  if (file.exists(rds_path)) {
    object <- readRDS(rds_path)
    if (is.null(object)) return(data.frame())
    table <- tryCatch(as.data.frame(object), error = function(e) e)
    if (!inherits(table, "error")) return(as.data.frame(table, stringsAsFactors = FALSE))
    warning("无法读取 GSEA RDS，将尝试 TSV：", rds_path, call. = FALSE)
  }
  if (file.exists(tsv_path)) return(read_tsv(tsv_path))
  data.frame()
}

filter_significant_gsea <- function(gsea_table, p_adjust_cutoff, qvalue_cutoff = NULL) {
  if (!nrow(gsea_table)) return(data.frame(NES = numeric(), stringsAsFactors = FALSE))
  if (!"NES" %in% names(gsea_table)) result_plot_abort("GSEA 结果缺少 NES 列。")
  p_column <- intersect(c("p.adjust", "padj"), names(gsea_table))
  if (!length(p_column)) result_plot_abort("GSEA 结果缺少 p.adjust 列。")
  keep <- !is.na(gsea_table[[p_column[[1]]]]) &
    suppressWarnings(as.numeric(gsea_table[[p_column[[1]]]])) <= p_adjust_cutoff
  if (!is.null(qvalue_cutoff)) {
    q_column <- intersect(c("qvalue", "qvalues"), names(gsea_table))
    if (!length(q_column)) result_plot_abort("配置了 GSEA qvalue 阈值，但结果中没有 qvalue/qvalues 列。")
    keep <- keep & !is.na(gsea_table[[q_column[[1]]]]) &
      suppressWarnings(as.numeric(gsea_table[[q_column[[1]]]])) <= qvalue_cutoff
  }
  output <- gsea_table[keep, , drop = FALSE]
  output$NES <- suppressWarnings(as.numeric(output$NES))
  output[is.finite(output$NES), , drop = FALSE]
}

# ---- Enrichment plots -----------------------------------------------------
# These helpers always start from the untouched RDS objects when available.
# This keeps redraws reproducible while allowing a plotting-specific threshold.

read_ora_table_for_plot <- function(project_root, group_id, comparison_id, set_id) {
  set_id <- tolower(as.character(set_id))
  supported <- c("up", "down", "up_down")
  if (!set_id %in% supported) {
    result_plot_abort("ORA set_id 只能是：", paste(supported, collapse = ", "), "。")
  }
  directory <- file.path(
    project_root, "data", "enrichment_results", safe_path_component(group_id), safe_path_component(comparison_id)
  )
  rds_path <- file.path(directory, paste0("ora_", set_id, ".rds"))
  tsv_path <- file.path(directory, paste0("ora_", set_id, ".tsv"))
  if (file.exists(rds_path)) {
    object <- readRDS(rds_path)
    if (is.null(object)) return(data.frame())
    table <- tryCatch(as.data.frame(object), error = function(e) e)
    if (!inherits(table, "error")) return(as.data.frame(table, stringsAsFactors = FALSE))
    warning("无法读取 ORA RDS，将尝试 TSV：", rds_path, call. = FALSE)
  }
  if (file.exists(tsv_path)) return(read_tsv(tsv_path))
  data.frame()
}

filter_significant_ora <- function(ora_table, p_adjust_cutoff, qvalue_cutoff = NULL) {
  if (!nrow(ora_table)) return(data.frame())
  p_column <- intersect(c("p.adjust", "padj"), names(ora_table))
  if (!length(p_column)) result_plot_abort("ORA 结果缺少 p.adjust 列。")
  p_adjust <- suppressWarnings(as.numeric(ora_table[[p_column[[1]]]]))
  keep <- is.finite(p_adjust) & p_adjust <= p_adjust_cutoff
  if (!is.null(qvalue_cutoff)) {
    q_column <- intersect(c("qvalue", "qvalues"), names(ora_table))
    if (!length(q_column)) result_plot_abort("配置了 ORA qvalue 阈值，但结果中没有 qvalue/qvalues 列。")
    qvalue <- suppressWarnings(as.numeric(ora_table[[q_column[[1]]]]))
    keep <- keep & is.finite(qvalue) & qvalue <= qvalue_cutoff
  }
  output <- as.data.frame(ora_table[keep, , drop = FALSE], stringsAsFactors = FALSE)
  output$p_adjust_plot <- p_adjust[keep]
  output
}

as_gene_ratio_numeric <- function(value) {
  value <- trimws(as.character(value))
  result <- suppressWarnings(as.numeric(value))
  fraction <- grepl("/", value, fixed = TRUE)
  if (any(fraction)) {
    pieces <- strsplit(value[fraction], "/", fixed = TRUE)
    result[fraction] <- vapply(pieces, function(piece) {
      if (length(piece) != 2L) return(NA_real_)
      numerator <- suppressWarnings(as.numeric(trimws(piece[[1]])))
      denominator <- suppressWarnings(as.numeric(trimws(piece[[2]])))
      if (!is.finite(numerator) || !is.finite(denominator) || denominator <= 0) return(NA_real_)
      numerator / denominator
    }, numeric(1))
  }
  result
}

enrichment_palette <- function(name) {
  palettes <- list(
    blue = c(low = "#DEEBF7", high = "#08519C"),
    red = c(low = "#FEE5D9", high = "#A50F15"),
    green = c(low = "#E5F5E0", high = "#006D2C"),
    purple = c(low = "#EFEDF5", high = "#54278F"),
    orange = c(low = "#FEE6CE", high = "#A63603"),
    teal = c(low = "#E0F3F3", high = "#006D6F"),
    grey = c(low = "#F0F0F0", high = "#252525"),
    pink = c(low = "#FDE0DD", high = "#C51B8A")
  )
  name <- tolower(as.character(name %||% "blue"))
  if (!name %in% names(palettes)) {
    result_plot_abort("未知 ORA 气泡图渐变色：", name, "；可用值为：", paste(names(palettes), collapse = ", "), "。")
  }
  palettes[[name]]
}

resolve_enrichment_size <- function(spec, auto_width, auto_height, save_config, label, variant = NULL) {
  spec <- spec %||% list()
  # `width` and `height` are intentionally resolved independently: for
  # example, {width: 4, height: null} keeps the calculated height.  Optional
  # by_plot entries override only a named output such as top_10 or all.
  override <- list()
  if (!is.null(variant) && length(variant)) {
    configured <- get_named_config(spec$by_plot %||% list(), as.character(variant))
    if (!is.null(configured)) {
      if (is.list(configured)) {
        override <- configured
      } else {
        override <- list(width = configured)
      }
    }
  }
  width <- result_plot_number(override$width %||% spec$width %||% auto_width, paste0(label, ".width"), 0, TRUE)
  height <- result_plot_number(override$height %||% spec$height %||% auto_height, paste0(label, ".height"), 0, TRUE)
  dpi <- result_plot_number(override$dpi %||% spec$dpi %||% save_config$dpi %||% 300, paste0(label, ".dpi"), 0, TRUE)
  # Full-result (all) plots can be hundreds of inches tall. At a normal PNG
  # DPI they exceed ragg's 50,000 px limit, so PDF is the safe default. A
  # user can explicitly opt into PNG via save.by_plot.all.formats if desired.
  is_all <- !is.null(variant) && length(variant) && identical(tolower(as.character(variant[[1]])), "all")
  formats <- override$formats %||% spec$formats %||% if (is_all) "pdf" else NULL
  list(width = width, height = height, dpi = dpi, formats = formats)
}

auto_enrichment_size <- function(labels, n, width_base = 7, height_base = 1.8,
                                 width_min = 8, width_max = 24, height_min = 5, height_max = 300,
                                 width_per_character = 0.075, height_per_row = 0.22) {
  labels <- as.character(labels)
  label_width <- if (length(labels)) max(nchar(labels, type = "width"), na.rm = TRUE) else 0
  list(
    width = max(width_min, min(width_max, width_base + label_width * width_per_character)),
    height = max(height_min, min(height_max, height_base + n * height_per_row))
  )
}

normalise_top_n <- function(value, label, default = 10L) {
  value <- value %||% default
  value <- suppressWarnings(as.numeric(unlist(value, use.names = FALSE)))
  if (!length(value) || any(!is.finite(value)) || any(value < 1) || any(abs(value - round(value)) > 1e-8)) {
    result_plot_abort(label, " 必须是一个或多个正整数。")
  }
  unique(as.integer(value))
}

normalise_top_selection <- function(value, label, default = 10L) {
  value <- value %||% default
  values <- trimws(as.character(unlist(value, use.names = FALSE)))
  if (!length(values) || any(is.na(values) | !nzchar(values))) {
    result_plot_abort(label, " 必须是一个或多个正整数，也可包含 'all'。")
  }
  is_all <- tolower(values) == "all"
  numeric_values <- suppressWarnings(as.numeric(values[!is_all]))
  if (any(!is.finite(numeric_values)) || any(numeric_values < 1) ||
      any(abs(numeric_values - round(numeric_values)) > 1e-8)) {
    result_plot_abort(label, " 必须是一个或多个正整数，也可包含 'all'。")
  }
  output <- c(if (length(numeric_values)) as.character(as.integer(numeric_values)) else character(),
              if (any(is_all)) "all" else character())
  unique(output)
}

top_selection_name <- function(selection) {
  selection <- as.character(selection)
  if (identical(selection, "all")) "all" else paste0("top_", selection)
}

top_selection_label <- function(selection) {
  selection <- as.character(selection)
  if (identical(selection, "all")) "All significant pathways" else paste0("Top ", selection, " per direction")
}

make_gsea_bidirectional_bar_plot <- function(gsea_table, group_id, comparison_id, gsea_config) {
  term_column <- as.character(gsea_config$term_column %||% "Description")
  if (length(term_column) != 1L || !nzchar(term_column) || !term_column %in% names(gsea_table)) {
    result_plot_abort("GSEA 双向柱图找不到 term_column：", term_column, "。")
  }
  table <- as.data.frame(gsea_table, stringsAsFactors = FALSE)
  table$NES <- suppressWarnings(as.numeric(table$NES))
  table$Pathway <- as.character(table[[term_column]])
  table <- table[is.finite(table$NES) & !is.na(table$Pathway) & nzchar(table$Pathway), , drop = FALSE]
  if (!nrow(table)) return(NULL)
  table$GSEA_Direction <- ifelse(table$NES >= 0, "Positive", "Negative")
  table$NES_plot <- ifelse(table$GSEA_Direction == "Positive", abs(table$NES), -abs(table$NES))
  colors <- gsea_config$colors %||% list()
  color_map <- c(
    Negative = result_plot_colour(colors, "negative", "#377EB8"),
    Positive = result_plot_colour(colors, "positive", "#E41A1C")
  )
  make_one <- function(selection) {
    selected <- do.call(rbind, lapply(c("Negative", "Positive"), function(direction) {
      part <- table[table$GSEA_Direction == direction, , drop = FALSE]
      if (!nrow(part)) return(part)
      part <- part[order(-abs(part$NES_plot), part$Pathway), , drop = FALSE]
      if (!identical(selection, "all")) {
        top_n <- as.integer(selection)
        part <- part[seq_len(min(nrow(part), top_n)), , drop = FALSE]
      }
      part
    }))
    if (!nrow(selected)) return(NULL)
    selected <- selected[order(selected$NES_plot, selected$Pathway), , drop = FALSE]
    selected$Pathway_plot <- make.unique(selected$Pathway)
    selected$Pathway_plot <- factor(selected$Pathway_plot, levels = selected$Pathway_plot)
    max_abs <- max(abs(selected$NES_plot), na.rm = TRUE)
    x_limit <- max(1.2, ceiling(max_abs * 10) / 10 + 0.15)
    p <- ggplot2::ggplot(selected, ggplot2::aes(x = .data$NES_plot, y = .data$Pathway_plot, fill = .data$GSEA_Direction)) +
      ggplot2::geom_col(width = 0.76, linewidth = 0) +
      ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.55, color = "grey35") +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.5, color = "black") +
      ggplot2::scale_fill_manual(values = color_map, limits = c("Negative", "Positive"), drop = FALSE) +
      ggplot2::scale_x_continuous(limits = c(-x_limit, x_limit), expand = ggplot2::expansion(mult = c(0.015, 0.015))) +
      ggplot2::labs(
        title = paste0(group_id, ": ", comparison_id, " (", top_selection_label(selection), ")"),
        x = "Normalized Enrichment Score (NES)", y = NULL, fill = NULL
      ) +
      ggplot2::theme_classic(base_size = result_plot_number(gsea_config$base_size %||% 11, "enrichment.gsea_bar.base_size", 0, TRUE)) +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8),
        axis.line = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13),
        axis.text.y = ggplot2::element_text(size = 8, color = "black"),
        axis.text.x = ggplot2::element_text(size = 10, color = "black"),
        axis.title.x = ggplot2::element_text(size = 11, margin = ggplot2::margin(t = 8)),
        legend.position = "top", legend.justification = "center", legend.key.width = grid::unit(1.2, "cm"),
        plot.margin = ggplot2::margin(t = 8, r = 12, b = 8, l = 8)
      )
    auto_size <- auto_enrichment_size(selected$Pathway, nrow(selected))
    list(plot = p, table = selected, size = auto_size)
  }
  selections <- normalise_top_selection(gsea_config$top_n, "enrichment.gsea_bar.top_n", c(10L, 20L, 50L))
  stats::setNames(lapply(selections, make_one), vapply(selections, top_selection_name, character(1)))
}

make_ora_dot_plot <- function(ora_table, group_id, comparison_id, set_id, ora_config) {
  required <- c("Description", "GeneRatio", "Count")
  missing <- setdiff(required, names(ora_table))
  if (length(missing)) result_plot_abort("ORA 结果缺少列：", paste(missing, collapse = ", "), "。")
  table <- as.data.frame(ora_table, stringsAsFactors = FALSE)
  table$GeneRatio_numeric <- as_gene_ratio_numeric(table$GeneRatio)
  table$Count <- suppressWarnings(as.numeric(table$Count))
  table$p_adjust_plot <- suppressWarnings(as.numeric(table$p_adjust_plot %||% table[[intersect(c("p.adjust", "padj"), names(table))[[1]]]]))
  table$Description <- as.character(table$Description)
  table <- table[
    is.finite(table$GeneRatio_numeric) & is.finite(table$Count) & table$Count >= 0 &
      is.finite(table$p_adjust_plot) & !is.na(table$Description) & nzchar(table$Description), , drop = FALSE
  ]
  if (!nrow(table)) return(NULL)
  table <- table[order(table$p_adjust_plot, -table$Count, table$Description), , drop = FALSE]
  color_by <- match.arg(as.character(ora_config$color_by %||% "log10_padjust"), c("log10_padjust", "padjust"))
  palette <- enrichment_palette(ora_config$palette %||% "blue")
  size_range <- suppressWarnings(as.numeric(unlist(ora_config$point_size_range %||% c(3, 10), use.names = FALSE)))
  if (length(size_range) != 2L || any(!is.finite(size_range)) || any(size_range <= 0) || size_range[[1]] >= size_range[[2]]) {
    result_plot_abort("enrichment.ora_dot.point_size_range 必须是递增的两个正数。")
  }
  point_alpha <- result_plot_number(ora_config$point_alpha %||% 0.85, "enrichment.ora_dot.point_alpha", 0)
  if (point_alpha > 1) result_plot_abort("enrichment.ora_dot.point_alpha 必须小于等于 1。")
  if (color_by == "log10_padjust") {
    color_aes <- ggplot2::aes(size = .data$Count, color = .data$log10_padjust)
    color_scale <- ggplot2::scale_color_gradient(low = unname(palette[["low"]]), high = unname(palette[["high"]]),
                                                  name = expression(-log[10](adjusted~P)))
  } else {
    color_aes <- ggplot2::aes(size = .data$Count, color = .data$p_adjust_plot)
    color_scale <- ggplot2::scale_color_gradient(low = unname(palette[["high"]]), high = unname(palette[["low"]]), name = "Adjusted P")
  }
  set_label <- c(up = "Up", down = "Down", up_down = "Up + Down")[[set_id]]
  make_one <- function(selection) {
    selected <- table
    if (!identical(selection, "all")) {
      top_n <- as.integer(selection)
      selected <- selected[seq_len(min(top_n, nrow(selected))), , drop = FALSE]
    }
    selected <- selected[order(selected$GeneRatio_numeric, selected$Description), , drop = FALSE]
    selected$Description_plot <- make.unique(selected$Description)
    selected$Description_plot <- factor(selected$Description_plot, levels = selected$Description_plot)
    selected$log10_padjust <- -log10(pmax(selected$p_adjust_plot, .Machine$double.xmin))
    p <- ggplot2::ggplot(selected, ggplot2::aes(x = .data$GeneRatio_numeric, y = .data$Description_plot)) +
      ggplot2::geom_point(mapping = color_aes, alpha = point_alpha) + color_scale +
      ggplot2::scale_size_continuous(name = "Count", range = size_range) +
      ggplot2::scale_x_continuous(labels = scales::label_percent(accuracy = 0.1), expand = ggplot2::expansion(mult = c(0.03, 0.08))) +
      ggplot2::labs(title = paste0(group_id, ": ", comparison_id, " (", set_label, ", ", top_selection_label(selection), ")"), x = "GeneRatio", y = NULL) +
      ggplot2::guides(
        color = ggplot2::guide_colorbar(order = 1, title.position = "top", title.hjust = 0.5),
        size = ggplot2::guide_legend(order = 2, title.position = "top", title.hjust = 0.5)
      ) +
      ggplot2::theme_bw(base_size = result_plot_number(ora_config$base_size %||% 11, "enrichment.ora_dot.base_size", 0, TRUE)) +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.7),
        panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        axis.text.y = ggplot2::element_text(size = 9, color = "black"), axis.text.x = ggplot2::element_text(color = "black"),
        axis.ticks = ggplot2::element_line(color = "black"), legend.position = ora_config$legend_position %||% "right",
        legend.direction = ora_config$legend_direction %||% "vertical", legend.title = ggplot2::element_text(size = 9, face = "bold"),
        legend.text = ggplot2::element_text(size = 8), legend.box = "vertical", plot.margin = ggplot2::margin(t = 8, r = 8, b = 8, l = 8)
      )
    auto_size <- auto_enrichment_size(selected$Description, nrow(selected), width_base = 7, height_base = 1.8,
                                      width_min = 8, width_max = 24, height_min = 4, height_max = 300,
                                      width_per_character = 0.085, height_per_row = 0.42)
    list(plot = p, table = selected, size = auto_size)
  }
  selections <- normalise_top_selection(ora_config$top_n, "enrichment.ora_dot.top_n", 10L)
  stats::setNames(lapply(selections, make_one), vapply(selections, top_selection_name, character(1)))
}

read_classified_gsea_for_plot <- function(project_root, path) {
  path <- result_plot_safe_rel_path(project_root, path)
  if (!file.exists(path)) result_plot_abort("找不到分类后的 GSEA 文件：", path)
  read_tsv(path)
}

prepare_classified_gsea <- function(table, category_config) {
  required <- c("Group", "Comparison", "NES")
  category_column <- as.character(category_config$category_column %||% "Category_Merged")
  description_column <- as.character(category_config$description_column %||% "Description")
  required <- unique(c(required, category_column, description_column))
  missing <- setdiff(required, names(table))
  if (length(missing)) result_plot_abort("分类 GSEA 文件缺少列：", paste(missing, collapse = ", "), "。")
  out <- as.data.frame(table, stringsAsFactors = FALSE)
  out$Group <- trimws(as.character(out$Group))
  out$Comparison <- trimws(as.character(out$Comparison))
  out$Category_plot <- trimws(as.character(out[[category_column]]))
  out$Description_plot <- trimws(as.character(out[[description_column]]))
  out$NES <- suppressWarnings(as.numeric(out$NES))
  out <- out[
    nzchar(out$Group) & nzchar(out$Comparison) & nzchar(out$Category_plot) & nzchar(out$Description_plot) & is.finite(out$NES), , drop = FALSE
  ]
  if (!nrow(out)) return(out)
  p_column <- intersect(c("qvalue", "qvalues", "p.adjust", "padj"), names(out))
  if (!length(p_column)) result_plot_abort("分类 GSEA 文件缺少 qvalue 或 p.adjust 列。")
  out$significance_plot <- suppressWarnings(as.numeric(out[[p_column[[1]]]]))
  # qvalue is preferred, but a missing qvalue must not accidentally become
  # the most significant pathway.  Fall back to p.adjust/padj per pathway.
  fallback_column <- intersect(c("p.adjust", "padj"), names(out))
  if (length(fallback_column)) {
    fallback <- suppressWarnings(as.numeric(out[[fallback_column[[1]]]]))
    needs_fallback <- !is.finite(out$significance_plot) | out$significance_plot <= 0
    out$significance_plot[needs_fallback] <- fallback[needs_fallback]
  }
  out$significance_plot[!is.finite(out$significance_plot) | out$significance_plot <= 0] <- NA_real_
  out$GSEA_Direction <- ifelse(out$NES >= 0, "Positive", "Negative")
  out[is.finite(out$significance_plot), , drop = FALSE]
}

result_plot_ordered_categories <- function(values, configured_order = NULL) {
  present <- unique(as.character(values[!is.na(values) & nzchar(values)]))
  if (!length(present)) return(character())
  if (is.null(configured_order) || !length(configured_order)) return(sort(present))
  configured_order <- trimws(as.character(unlist(configured_order, use.names = FALSE)))
  configured_order <- configured_order[!is.na(configured_order) & nzchar(configured_order)]
  if (anyDuplicated(configured_order)) result_plot_abort("category_order 中存在重复类别。")
  c(configured_order[configured_order %in% present], sort(setdiff(present, configured_order)))
}

result_plot_category_colours <- function(categories, configured_colours = NULL) {
  categories <- unique(as.character(categories))
  categories <- categories[!is.na(categories) & nzchar(categories)]
  configured <- character()
  if (!is.null(configured_colours) && length(configured_colours)) {
    configured <- unlist(configured_colours, use.names = TRUE)
    if (is.null(names(configured)) || any(!nzchar(names(configured)))) {
      result_plot_abort("分类 GSEA 桑基图的 category_colors 必须是“类别: 颜色”的命名映射。")
    }
    configured <- stats::setNames(as.character(configured), names(configured))
  }
  output <- configured[intersect(names(configured), categories)]
  missing <- setdiff(categories, names(output))
  if (length(missing)) {
    generated <- if (length(missing) == 1L) "#4C78A8" else grDevices::hcl.colors(length(missing), palette = "Dark 3")
    output <- c(output, stats::setNames(generated, missing))
  }
  output[categories]
}

read_category_sankey_trend_sets <- function(config, project_root, sankey_config) {
  input_config <- config$input %||% list()
  path <- result_plot_safe_rel_path(project_root, input_config$trend_sets_file %||% "data/trend_sets.tsv")
  if (!file.exists(path)) result_plot_abort("找不到 trend_sets 文件：", path)
  table <- read_tsv(path, c("group_id", "trend_id", "comparison_ids"))
  for (column in c("group_id", "trend_id", "comparison_ids")) {
    table[[column]] <- trimws(as.character(table[[column]]))
  }
  table <- table[nzchar(table$group_id) & nzchar(table$trend_id) & nzchar(table$comparison_ids), , drop = FALSE]
  selected_groups <- trimws(as.character(unlist(sankey_config$groups %||% character(), use.names = FALSE)))
  selected_ids <- trimws(as.character(unlist(sankey_config$trend_ids %||% character(), use.names = FALSE)))
  if (length(selected_groups)) table <- table[table$group_id %in% selected_groups, , drop = FALSE]
  if (length(selected_ids)) table <- table[table$trend_id %in% selected_ids, , drop = FALSE]
  if (!nrow(table)) result_plot_abort("classified_gsea.category_sankey 没有匹配到 trend_sets。")
  keys <- paste(table$group_id, table$trend_id, sep = "\r")
  if (anyDuplicated(keys)) result_plot_abort("trend_sets 中同一 Group / trend_id 不能重复。")
  lapply(seq_len(nrow(table)), function(index) {
    comparison_ids <- trimws(unlist(strsplit(table$comparison_ids[[index]], ";", fixed = TRUE)))
    comparison_ids <- comparison_ids[nzchar(comparison_ids)]
    if (length(comparison_ids) < 2L || anyDuplicated(comparison_ids)) {
      result_plot_abort(table$group_id[[index]], " / ", table$trend_id[[index]],
                        " 的 comparison_ids 必须是至少两个且不重复的比较。")
    }
    list(group_id = table$group_id[[index]], trend_id = table$trend_id[[index]], comparison_ids = comparison_ids)
  })
}

build_classified_gsea_category_composition <- function(table, group_id, comparison_ids, category_order) {
  selected <- table[table$Group == group_id & table$Comparison %in% comparison_ids, , drop = FALSE]
  if (!nrow(selected)) return(NULL)
  category_order <- category_order[category_order %in% selected$Category_plot]
  if (!length(category_order)) return(NULL)
  counted <- stats::aggregate(
    rep(1L, nrow(selected)),
    by = list(Comparison = selected$Comparison, Category = selected$Category_plot), FUN = sum
  )
  names(counted)[[3]] <- "Pathway_Count"
  complete <- expand.grid(
    Comparison = comparison_ids, Category = category_order,
    stringsAsFactors = FALSE
  )
  output <- merge(complete, counted, by = c("Comparison", "Category"), all.x = TRUE, sort = FALSE)
  output$Pathway_Count[is.na(output$Pathway_Count)] <- 0L
  total <- stats::aggregate(Pathway_Count ~ Comparison, data = output, FUN = sum)
  output <- merge(output, total, by = "Comparison", suffixes = c("", "_Total"), all.x = TRUE, sort = FALSE)
  output$Proportion <- ifelse(output$Pathway_Count_Total > 0, output$Pathway_Count / output$Pathway_Count_Total, 0)
  output$Comparison <- factor(output$Comparison, levels = comparison_ids, ordered = TRUE)
  output$Category <- factor(output$Category, levels = category_order, ordered = TRUE)
  output[order(output$Comparison, output$Category), , drop = FALSE]
}

make_classified_gsea_category_sankey_plot <- function(composition, group_id, trend_id, mode, category_colours,
                                                       sankey_config) {
  mode <- match.arg(mode, c("count", "proportion"))
  value_column <- if (mode == "count") "Pathway_Count" else "Proportion"
  y_label <- if (mode == "count") "Number of significant GSEA pathways" else "Proportion of significant GSEA pathways"
  mode_label <- if (mode == "count") "Pathway count" else "Pathway proportion"
  alpha <- result_plot_number(sankey_config$flow_alpha %||% 0.78, "classified_gsea.category_sankey.flow_alpha", 0)
  if (alpha > 1) result_plot_abort("classified_gsea.category_sankey.flow_alpha 必须小于等于 1。")
  width <- result_plot_number(sankey_config$stratum_width %||% 0.48, "classified_gsea.category_sankey.stratum_width", 0, TRUE)
  if (width > 1) result_plot_abort("classified_gsea.category_sankey.stratum_width 必须小于等于 1。")
  p <- ggplot2::ggplot(
    composition,
    ggplot2::aes(
      x = .data$Comparison, stratum = .data$Category, alluvium = .data$Category,
      y = .data[[value_column]], fill = .data$Category
    )
  ) +
    # The categories themselves are the alluvia. curve_type='linear' gives
    # straight-sided transitions for composition change rather than a pathway
    # significance-trajectory Sankey.
    ggalluvial::geom_alluvium(curve_type = "linear", alpha = alpha, width = width, reverse = TRUE) +
    ggalluvial::geom_stratum(width = width, color = "grey25", linewidth = 0.25, reverse = TRUE) +
    ggplot2::scale_fill_manual(values = category_colours, limits = names(category_colours), drop = FALSE) +
    ggplot2::scale_x_discrete(drop = FALSE, expand = c(0.08, 0.08)) +
    ggplot2::labs(
      title = paste0(group_id, ": ", trend_id, " — classified GSEA category composition"),
      subtitle = mode_label,
      x = NULL, y = y_label, fill = NULL
    ) +
    ggplot2::theme_classic(base_size = result_plot_number(sankey_config$base_size %||% 11, "classified_gsea.category_sankey.base_size", 0, TRUE)) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.7),
      axis.line = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      legend.position = sankey_config$legend_position %||% "right"
    )
  if (mode == "proportion") {
    p <- p + ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1), expand = c(0, 0))
  } else {
    p <- p + ggplot2::scale_y_continuous(expand = c(0, 0))
  }
  p
}

auto_category_sankey_size <- function(composition) {
  n_comparisons <- length(unique(composition$Comparison))
  n_categories <- length(unique(composition$Category))
  longest_category <- if (n_categories) max(nchar(as.character(composition$Category)), na.rm = TRUE) else 0
  list(
    width = max(7, min(20, 3.5 + 1.25 * n_comparisons + 0.045 * longest_category)),
    height = max(5.5, min(18, 3.5 + 0.28 * n_categories))
  )
}

run_classified_gsea_category_sankey_stage <- function(config, project_root, output_dir, save_config, manifest_rows) {
  enrichment_config <- (config$result_plots %||% list())$enrichment %||% list()
  category_config <- enrichment_config$classified_gsea %||% list()
  sankey_config <- category_config$category_sankey %||% list()
  if (!isTRUE(sankey_config$enabled %||% FALSE)) return(manifest_rows)
  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    result_plot_abort("分类 GSEA 组成桑基图缺少 R 包 'ggalluvial'。")
  }
  classified <- prepare_classified_gsea(
    read_classified_gsea_for_plot(project_root, category_config$file %||% "table/enrichment_results/gsea_sig_group.merged.classified.csv"),
    category_config
  )
  if (!nrow(classified)) {
    warning("分类 GSEA 文件没有可用于组成桑基图的结果。", call. = FALSE)
    return(manifest_rows)
  }
  gsva_config <- config$gobp_gsva_plots %||% list()
  configured_order <- sankey_config$category_order %||% gsva_config$category_order %||% NULL
  configured_colours <- sankey_config$category_colors %||% gsva_config$category_colors %||% NULL
  modes <- tolower(trimws(as.character(unlist(sankey_config$modes %||% c("count", "proportion"), use.names = FALSE))))
  invalid_modes <- setdiff(modes, c("count", "proportion"))
  if (length(invalid_modes) || !length(modes)) {
    result_plot_abort("classified_gsea.category_sankey.modes 仅支持 count、proportion。")
  }
  for (entry in read_category_sankey_trend_sets(config, project_root, sankey_config)) {
    available <- classified[classified$Group == entry$group_id & classified$Comparison %in% entry$comparison_ids, , drop = FALSE]
    if (!nrow(available)) {
      warning("分类 GSEA 组成桑基图无对应结果，跳过：", entry$group_id, " / ", entry$trend_id, call. = FALSE)
      next
    }
    category_order <- result_plot_ordered_categories(available$Category_plot, configured_order)
    composition <- build_classified_gsea_category_composition(classified, entry$group_id, entry$comparison_ids, category_order)
    if (is.null(composition)) next
    category_colours <- result_plot_category_colours(category_order, configured_colours)
    auto_size <- auto_category_sankey_size(composition)
    for (mode in unique(modes)) {
      plot <- make_classified_gsea_category_sankey_plot(
        composition, entry$group_id, entry$trend_id, mode, category_colours, sankey_config
      )
      size <- resolve_enrichment_size(
        sankey_config$save %||% list(), auto_size$width, auto_size$height, save_config,
        "enrichment.classified_gsea.category_sankey.save", variant = mode
      )
      files <- save_result_plot(
        plot,
        file.path(
          output_dir, "enrichment", "classified_gsea_category_sankey",
          safe_path_component(entry$group_id), safe_path_component(entry$trend_id),
          paste0("category_composition_", mode)
        ),
        size, save_config, formats_override = size$formats
      )
      manifest_rows[[length(manifest_rows) + 1L]] <- result_plot_manifest_row(
        paste0("classified_gsea_category_sankey_", mode), entry$group_id,
        trend_id = entry$trend_id, files = files
      )
    }
  }
  manifest_rows
}

make_classified_gsea_count_plot <- function(table, group_id, comparison_id, category_config) {
  selected <- table[table$Group == group_id & table$Comparison == comparison_id, , drop = FALSE]
  if (!nrow(selected)) return(NULL)
  count_table <- stats::aggregate(
    rep(1L, nrow(selected)), by = list(Category = selected$Category_plot, GSEA_Direction = selected$GSEA_Direction), FUN = sum
  )
  names(count_table)[[3]] <- "Pathway_Count"
  total <- stats::aggregate(Pathway_Count ~ Category, data = count_table, FUN = sum)
  count_table <- merge(count_table, total, by = "Category", suffixes = c("", "_Total"), sort = FALSE)
  category_order <- total$Category[order(total$Pathway_Count, total$Category)]
  count_table$Category <- factor(count_table$Category, levels = category_order)
  count_table$Plot_Count <- ifelse(count_table$GSEA_Direction == "Negative", -count_table$Pathway_Count, count_table$Pathway_Count)
  colors <- category_config$count_bar$colors %||% list()
  color_map <- c(
    Negative = result_plot_colour(colors, "negative", "#3B75AF"),
    Positive = result_plot_colour(colors, "positive", "#D95F5F")
  )
  p <- ggplot2::ggplot(count_table, ggplot2::aes(x = .data$Plot_Count, y = .data$Category, fill = .data$GSEA_Direction)) +
    ggplot2::geom_col(width = 0.72) + ggplot2::geom_vline(xintercept = 0, linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = .data$Pathway_Count), position = ggplot2::position_stack(vjust = 0.5), size = 3) +
    ggplot2::scale_x_continuous(labels = abs, expand = ggplot2::expansion(mult = c(0.08, 0.08))) +
    ggplot2::scale_fill_manual(values = color_map, limits = c("Negative", "Positive"), drop = FALSE) +
    ggplot2::labs(title = paste0(group_id, ": ", comparison_id), x = "Number of enriched pathways", y = NULL, fill = NULL) +
    ggplot2::theme_classic(base_size = result_plot_number(category_config$count_bar$base_size %||% 12, "enrichment.classified_gsea.count_bar.base_size", 0, TRUE)) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10), plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"), legend.position = "top")
  auto_size <- auto_enrichment_size(as.character(category_order), length(category_order), width_base = 6, height_base = 1.8,
                                    width_min = 7, width_max = 20, height_min = 4, height_max = 20,
                                    width_per_character = 0.10, height_per_row = 0.42)
  list(plot = p, table = count_table, size = auto_size, category_order = category_order)
}

make_classified_gsea_bubble_plot <- function(table, group_id, comparison_id, category_config, category_order) {
  selected <- table[table$Group == group_id & table$Comparison == comparison_id, , drop = FALSE]
  if (!nrow(selected)) return(NULL)
  required <- c("setSize")
  missing <- setdiff(required, names(selected))
  if (length(missing)) result_plot_abort("分类 GSEA 气泡图缺少列：", paste(missing, collapse = ", "), "。")
  selected$setSize <- suppressWarnings(as.numeric(selected$setSize))
  selected <- selected[is.finite(selected$setSize) & selected$setSize > 0, , drop = FALSE]
  if (!nrow(selected)) return(NULL)
  top_n <- normalise_top_n(category_config$representative_bubble$top_n_per_category,
                           "enrichment.classified_gsea.representative_bubble.top_n_per_category", 1L)[[1]]
  selected <- selected[order(selected$Category_plot, selected$significance_plot, -abs(selected$NES), selected$Description_plot), , drop = FALSE]
  selected <- do.call(rbind, lapply(split(selected, selected$Category_plot), function(part) part[seq_len(min(nrow(part), top_n)), , drop = FALSE]))
  selected$Category_plot <- factor(selected$Category_plot, levels = category_order)
  selected$neg_log10_significance <- -log10(pmax(selected$significance_plot, .Machine$double.xmin))
  bubble_config <- category_config$representative_bubble %||% list()
  label_size <- result_plot_number(bubble_config$label_size %||% 3, "enrichment.classified_gsea.representative_bubble.label_size", 0, TRUE)
  point_alpha <- result_plot_number(bubble_config$point_alpha %||% 0.85, "enrichment.classified_gsea.representative_bubble.point_alpha", 0)
  if (point_alpha > 1) result_plot_abort("enrichment.classified_gsea.representative_bubble.point_alpha 必须小于等于 1。")
  palette_option <- as.character(bubble_config$viridis_option %||% "D")
  p <- ggplot2::ggplot(selected, ggplot2::aes(x = .data$NES, y = .data$Category_plot)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(size = .data$setSize, color = .data$neg_log10_significance), alpha = point_alpha) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = .data$Description_plot), size = label_size, max.overlaps = Inf,
      box.padding = 0.4, point.padding = 0.25, min.segment.length = 0
    ) +
    ggplot2::scale_color_viridis_c(option = palette_option, direction = 1, name = expression(-log[10](qvalue~or~adjusted~P))) +
    ggplot2::labs(title = paste0(group_id, ": ", comparison_id), x = "NES", y = NULL, size = "Gene set size") +
    ggplot2::theme_classic(base_size = result_plot_number(bubble_config$base_size %||% 12, "enrichment.classified_gsea.representative_bubble.base_size", 0, TRUE)) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9), plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                   legend.position = bubble_config$legend_position %||% "right", plot.margin = ggplot2::margin(t = 8, r = 20, b = 8, l = 8)) +
    ggplot2::coord_cartesian(clip = "off")
  auto_size <- auto_enrichment_size(as.character(category_order), length(category_order), width_base = 8, height_base = 1.8,
                                    width_min = 10, width_max = 28, height_min = 4, height_max = 20,
                                    width_per_character = 0.10, height_per_row = 0.50)
  list(plot = p, table = selected, size = auto_size)
}

run_enrichment_plot_stage <- function(project_root, groups, entries_by_group, plot_config, save_config, output_dir, manifest_rows,
                                      default_gsea_p_adjust_cutoff) {
  enrichment_plot_config <- plot_config$enrichment %||% list()
  if (!isTRUE(enrichment_plot_config$enabled %||% FALSE)) return(manifest_rows)
  gsea_config <- enrichment_plot_config$gsea_bar %||% list()
  ora_config <- enrichment_plot_config$ora_dot %||% list()
  category_config <- enrichment_plot_config$classified_gsea %||% list()
  gsea_p_cutoff <- result_plot_number(
    gsea_config$p_adjust_cutoff %||% default_gsea_p_adjust_cutoff,
    "enrichment.gsea_bar.p_adjust_cutoff", 0, TRUE
  )
  if (gsea_p_cutoff > 1) result_plot_abort("enrichment.gsea_bar.p_adjust_cutoff 必须小于等于 1。")
  gsea_q_cutoff <- gsea_config$qvalue_cutoff %||% NULL
  if (!is.null(gsea_q_cutoff) && length(gsea_q_cutoff)) {
    gsea_q_cutoff <- result_plot_number(gsea_q_cutoff, "enrichment.gsea_bar.qvalue_cutoff", 0, TRUE)
    if (gsea_q_cutoff > 1) result_plot_abort("enrichment.gsea_bar.qvalue_cutoff 必须小于等于 1。")
  } else gsea_q_cutoff <- NULL
  ora_p_cutoff <- result_plot_number(ora_config$p_adjust_cutoff %||% gsea_p_cutoff, "enrichment.ora_dot.p_adjust_cutoff", 0, TRUE)
  if (ora_p_cutoff > 1) result_plot_abort("enrichment.ora_dot.p_adjust_cutoff 必须小于等于 1。")
  ora_q_cutoff <- ora_config$qvalue_cutoff %||% NULL
  if (!is.null(ora_q_cutoff) && length(ora_q_cutoff)) {
    ora_q_cutoff <- result_plot_number(ora_q_cutoff, "enrichment.ora_dot.qvalue_cutoff", 0, TRUE)
    if (ora_q_cutoff > 1) result_plot_abort("enrichment.ora_dot.qvalue_cutoff 必须小于等于 1。")
  } else ora_q_cutoff <- NULL

  for (group_id in groups) {
    comparison_ids <- vapply(entries_by_group[[group_id]], `[[`, character(1), "comparison_id")
    for (comparison_id in comparison_ids) {
      if (isTRUE(gsea_config$enabled %||% TRUE)) {
        gsea <- filter_significant_gsea(read_gsea_table_for_plot(project_root, group_id, comparison_id), gsea_p_cutoff, gsea_q_cutoff)
        if (nrow(gsea)) {
          plots <- make_gsea_bidirectional_bar_plot(gsea, group_id, comparison_id, gsea_config)
          for (name in names(plots)) {
            item <- plots[[name]]
            if (is.null(item)) next
            size <- resolve_enrichment_size(gsea_config$save %||% list(), item$size$width, item$size$height, save_config,
                                            "enrichment.gsea_bar.save", variant = name)
            files <- save_result_plot(item$plot,
              file.path(output_dir, "enrichment", "gsea_bidirectional_bar", safe_path_component(group_id), safe_path_component(comparison_id), name),
              size, save_config, formats_override = size$formats
            )
            manifest_rows[[length(manifest_rows) + 1L]] <- result_plot_manifest_row("gsea_bidirectional_bar", group_id, comparison_id,
                                                                                     files = files, x_limits = NULL, y_limits = NULL)
          }
        } else warning("GSEA 双向柱图无满足阈值的通路，跳过：", group_id, "/", comparison_id, call. = FALSE)
      }
      if (isTRUE(ora_config$enabled %||% TRUE)) {
        for (set_id in as.character(unlist(ora_config$sets %||% c("up", "down", "up_down"), use.names = FALSE))) {
          set_id <- tolower(trimws(set_id))
          ora <- filter_significant_ora(read_ora_table_for_plot(project_root, group_id, comparison_id, set_id), ora_p_cutoff, ora_q_cutoff)
          plots <- make_ora_dot_plot(ora, group_id, comparison_id, set_id, ora_config)
          if (is.null(plots)) {
            warning("ORA 气泡图无满足阈值的通路，跳过：", group_id, "/", comparison_id, "/", set_id, call. = FALSE)
            next
          }
          for (name in names(plots)) {
            item <- plots[[name]]
            if (is.null(item)) next
            size <- resolve_enrichment_size(ora_config$save %||% list(), item$size$width, item$size$height, save_config,
                                            "enrichment.ora_dot.save", variant = name)
            files <- save_result_plot(item$plot,
              file.path(output_dir, "enrichment", "ora_dot", safe_path_component(group_id), safe_path_component(comparison_id),
                        paste0("ora_", safe_path_component(set_id), "_", name)),
              size, save_config, formats_override = size$formats
            )
            manifest_rows[[length(manifest_rows) + 1L]] <- result_plot_manifest_row("ora_dot", group_id, comparison_id, files = files)
          }
        }
      }
    }
  }

  if (isTRUE(category_config$enabled %||% TRUE)) {
    if (!requireNamespace("ggrepel", quietly = TRUE)) result_plot_abort("分类 GSEA 气泡图缺少 R 包 'ggrepel'。")
    classified <- prepare_classified_gsea(
      read_classified_gsea_for_plot(project_root, category_config$file %||% "table/enrichment_results/gsea_sig_group.merged.classified.csv"),
      category_config
    )
    for (group_id in groups) {
      comparison_ids <- vapply(entries_by_group[[group_id]], `[[`, character(1), "comparison_id")
      for (comparison_id in comparison_ids) {
        count_item <- make_classified_gsea_count_plot(classified, group_id, comparison_id, category_config)
        if (is.null(count_item)) {
          warning("分类 GSEA 文件无对应结果，跳过：", group_id, "/", comparison_id, call. = FALSE)
          next
        }
        count_size <- resolve_enrichment_size(category_config$count_bar$save %||% list(), count_item$size$width, count_item$size$height,
                                               save_config, "enrichment.classified_gsea.count_bar.save")
        files <- save_result_plot(count_item$plot,
          file.path(output_dir, "enrichment", "classified_gsea", safe_path_component(group_id), safe_path_component(comparison_id), "category_bidirectional_bar"),
          count_size, save_config, formats_override = count_size$formats
        )
        manifest_rows[[length(manifest_rows) + 1L]] <- result_plot_manifest_row("classified_gsea_category_count", group_id, comparison_id, files = files)
        bubble_item <- make_classified_gsea_bubble_plot(classified, group_id, comparison_id, category_config, count_item$category_order)
        if (!is.null(bubble_item)) {
          bubble_size <- resolve_enrichment_size(category_config$representative_bubble$save %||% list(), bubble_item$size$width,
                                                  bubble_item$size$height, save_config, "enrichment.classified_gsea.representative_bubble.save")
          files <- save_result_plot(bubble_item$plot,
            file.path(output_dir, "enrichment", "classified_gsea", safe_path_component(group_id), safe_path_component(comparison_id), "category_representative_bubble"),
            bubble_size, save_config, formats_override = bubble_size$formats
          )
          manifest_rows[[length(manifest_rows) + 1L]] <- result_plot_manifest_row("classified_gsea_representative_bubble", group_id, comparison_id, files = files)
        }
      }
    }
  }
  manifest_rows
}

make_de_gsea_count_data <- function(de_plot, gsea_table, p_adjust_cutoff, qvalue_cutoff) {
  gsea <- filter_significant_gsea(gsea_table, p_adjust_cutoff, qvalue_cutoff)
  values <- data.frame(
    Data_type = rep(c("DEG", "GSEA pathway"), each = 2L),
    Direction = rep(c("Positive", "Negative"), times = 2L),
    n = c(
      sum(de_plot$Status == "Up", na.rm = TRUE),
      sum(de_plot$Status == "Down", na.rm = TRUE),
      sum(gsea$NES > 0, na.rm = TRUE),
      sum(gsea$NES < 0, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  values$Data_type <- factor(values$Data_type, levels = c("DEG", "GSEA pathway"))
  values$Direction <- factor(values$Direction, levels = c("Positive", "Negative"))
  values$plot_n <- ifelse(values$Direction == "Negative", -values$n, values$n)
  values$Fill_group <- paste(values$Data_type, values$Direction, sep = "_")
  values$x <- ifelse(values$Data_type == "DEG", 0.88, 1.12)
  values
}

calculate_count_limit <- function(count_tables) {
  maximum <- max(abs(unlist(lapply(count_tables, function(x) x$plot_n), use.names = FALSE)), na.rm = TRUE)
  if (!is.finite(maximum) || maximum == 0) maximum <- 1
  ceiling(maximum * 1.15)
}

make_de_gsea_count_plot <- function(count_table, comparison_id, y_limit, count_config) {
  colours <- count_config$colors %||% list()
  colour_map <- c(
    DEG_Positive = result_plot_colour(colours, "deg_up", "#C1121F"),
    DEG_Negative = result_plot_colour(colours, "deg_down", "#669BBC"),
    `GSEA pathway_Positive` = result_plot_colour(colours, "gsea_positive", "#8E6CBE"),
    `GSEA pathway_Negative` = result_plot_colour(colours, "gsea_negative", "#56A3A6")
  )
  ggplot2::ggplot(count_table, ggplot2::aes(x = .data$x, y = .data$plot_n, fill = .data$Fill_group)) +
    ggplot2::geom_col(width = 0.16, color = "grey25", linewidth = 0.2) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$n, vjust = ifelse(.data$plot_n >= 0, -0.45, 1.45)),
      size = result_plot_number(count_config$label_size %||% 3.1, "deg_gsea_count.label_size", 0, TRUE),
      color = "grey15", fontface = "bold"
    ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.55, color = "grey20") +
    ggplot2::scale_x_continuous(
      breaks = c(0.88, 1.12), labels = c("DEG", "GSEA"), limits = c(0.68, 1.32),
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::scale_fill_manual(
      values = colour_map, limits = names(colour_map), drop = FALSE,
      labels = c(
        "Upregulated genes", "Downregulated genes",
        "Positively enriched pathways", "Negatively enriched pathways"
      ), name = NULL
    ) +
    ggplot2::scale_y_continuous(limits = c(-y_limit, y_limit), labels = abs,
                                 expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(title = comparison_id, x = NULL, y = "Count") +
    ggplot2::theme_classic(base_size = result_plot_number(count_config$base_size %||% 12, "deg_gsea_count.base_size", 0, TRUE)) +
    ggplot2::theme(
      legend.position = "top",
      legend.justification = "center",
      legend.key.width = grid::unit(0.65, "cm"),
      legend.key.height = grid::unit(0.35, "cm"),
      legend.spacing.x = grid::unit(0.15, "cm"),
      axis.text.x = ggplot2::element_text(color = "grey15", face = "bold", size = 9.5, margin = ggplot2::margin(t = 4)),
      axis.text.y = ggplot2::element_text(color = "grey20", size = 9),
      axis.title.y = ggplot2::element_text(size = 11, margin = ggplot2::margin(r = 10)),
      axis.line.x = ggplot2::element_line(color = "grey25", linewidth = 0.45),
      axis.line.y = ggplot2::element_line(color = "grey25", linewidth = 0.45),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey92", linewidth = 0.35),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 10),
      plot.margin = ggplot2::margin(t = 8, r = 10, b = 8, l = 10)
    )
}

result_plot_manifest_row <- function(plot_type, group_id, comparison_id = NA_character_, trend_id = NA_character_,
                                     files, x_limits = NULL, y_limits = NULL) {
  data.frame(
    plot_type = plot_type,
    group_id = group_id,
    comparison_id = comparison_id,
    trend_id = trend_id,
    files = paste(files, collapse = ";"),
    x_min = if (is.null(x_limits)) NA_real_ else x_limits[[1]],
    x_max = if (is.null(x_limits)) NA_real_ else x_limits[[2]],
    y_min = if (is.null(y_limits)) NA_real_ else y_limits[[1]],
    y_max = if (is.null(y_limits)) NA_real_ else y_limits[[2]],
    stringsAsFactors = FALSE
  )
}

run_result_plot_stage <- function(config, project_root) {
  plot_config <- config$result_plots %||% list()
  if (!isTRUE(plot_config$enabled %||% TRUE)) {
    message("result_plots.enabled=FALSE，跳过结果绘图。")
    return(invisible(NULL))
  }
  require_result_plot_packages(config)
  output_dir <- resolve_path(project_root, plot_config$output_dir %||% "figures/result_plots")
  ensure_dir(output_dir)
  save_config <- plot_config$save %||% list()
  volcano_config <- plot_config$volcano %||% list()
  trend_config <- plot_config$trend_bar %||% list()
  count_config <- plot_config$deg_gsea_count %||% list()
  thresholds <- list(
    fc_cutoff = result_plot_number(plot_config$fc_cutoff %||% (config$trend %||% list())$fc_cutoff %||% 2, "result_plots.fc_cutoff", 1),
    padj_cutoff = result_plot_number(plot_config$padj_cutoff %||% (config$trend %||% list())$padj_cutoff %||% 0.05, "result_plots.padj_cutoff", 0, TRUE)
  )
  if (thresholds$padj_cutoff > 1) result_plot_abort("result_plots.padj_cutoff 必须小于等于 1。")
  thresholds$log2fc_cutoff <- log2(thresholds$fc_cutoff)
  padj_floor <- result_plot_number(volcano_config$padj_floor %||% 1e-300, "volcano.padj_floor", 0, TRUE)
  if (padj_floor > 1) result_plot_abort("volcano.padj_floor 必须小于等于 1。")

  de_entries <- read_successful_de_results_for_plot(project_root)
  groups <- unique(vapply(de_entries, `[[`, character(1), "group_id"))
  entries_by_group <- split(de_entries, vapply(de_entries, `[[`, character(1), "group_id"))
  axis_scope <- volcano_config$axis_scope %||% "within_group"
  if (!axis_scope %in% c("within_group", "global")) {
    result_plot_abort("volcano.axis_scope 只能是 within_group 或 global。")
  }
  all_plot_tables <- lapply(de_entries, function(entry) {
    classify_de_for_plot(entry$result, thresholds$fc_cutoff, thresholds$padj_cutoff, padj_floor)
  })
  global_limits <- if (axis_scope == "global") calculate_volcano_limits(all_plot_tables, volcano_config) else NULL
  plot_manifest <- list()

  volcano_individual_size <- result_plot_size(save_config, "volcano_individual", c(5, 5))
  volcano_combined_size <- result_plot_size(save_config, "volcano_combined", c(12, 12))
  count_combined_size <- result_plot_size(save_config, "deg_gsea_count_combined", c(12, 10))
  trend_bar_size <- result_plot_size(save_config, "trend_bar", c(9, 6))
  gsea_p_adjust_cutoff <- result_plot_number(
    count_config$gsea_p_adjust_cutoff %||% (config$enrichment %||% list())$pvalue_cutoff %||% 0.05,
    "deg_gsea_count.gsea_p_adjust_cutoff", 0, TRUE
  )
  if (gsea_p_adjust_cutoff > 1) result_plot_abort("deg_gsea_count.gsea_p_adjust_cutoff 必须小于等于 1。")
  gsea_qvalue_cutoff <- count_config$gsea_qvalue_cutoff %||% NULL
  if (!is.null(gsea_qvalue_cutoff) && length(gsea_qvalue_cutoff)) {
    gsea_qvalue_cutoff <- result_plot_number(gsea_qvalue_cutoff, "deg_gsea_count.gsea_qvalue_cutoff", 0, TRUE)
    if (gsea_qvalue_cutoff > 1) result_plot_abort("deg_gsea_count.gsea_qvalue_cutoff 必须小于等于 1。")
  } else {
    gsea_qvalue_cutoff <- NULL
  }

  for (group_id in groups) {
    group_entries <- entries_by_group[[group_id]]
    comparison_ids <- vapply(group_entries, `[[`, character(1), "comparison_id")
    layout <- resolve_comparison_layout(plot_config, group_id, comparison_ids)
    group_plot_tables <- stats::setNames(
      lapply(group_entries, function(entry) classify_de_for_plot(entry$result, thresholds$fc_cutoff, thresholds$padj_cutoff, padj_floor)),
      comparison_ids
    )
    limits <- global_limits %||% calculate_volcano_limits(unname(group_plot_tables), volcano_config)
    volcano_plots <- lapply(layout$comparison_order, function(comparison_id) {
      make_volcano_plot(group_plot_tables[[comparison_id]], group_id, comparison_id, thresholds, limits, volcano_config)
    })
    names(volcano_plots) <- layout$comparison_order
    for (comparison_id in layout$comparison_order) {
      individual_base <- file.path(output_dir, "volcano", safe_path_component(group_id), safe_path_component(comparison_id), "volcano")
      files <- save_result_plot(volcano_plots[[comparison_id]], individual_base, volcano_individual_size, save_config)
      plot_manifest[[length(plot_manifest) + 1L]] <- result_plot_manifest_row(
        "volcano_individual", group_id, comparison_id, files = files, x_limits = limits$x, y_limits = limits$y
      )
    }
    combined_volcano <- compose_comparison_grid(volcano_plots, layout)
    combined_volcano <- combined_volcano & ggplot2::theme(legend.position = "none")
    files <- save_result_plot(
      combined_volcano, file.path(output_dir, "volcano", safe_path_component(group_id), "volcano_combined"),
      volcano_combined_size, save_config
    )
    plot_manifest[[length(plot_manifest) + 1L]] <- result_plot_manifest_row(
      "volcano_combined", group_id, files = files, x_limits = limits$x, y_limits = limits$y
    )

    count_tables <- lapply(layout$comparison_order, function(comparison_id) {
      make_de_gsea_count_data(
        group_plot_tables[[comparison_id]], read_gsea_table_for_plot(project_root, group_id, comparison_id),
        gsea_p_adjust_cutoff, gsea_qvalue_cutoff
      )
    })
    names(count_tables) <- layout$comparison_order
    count_limit <- calculate_count_limit(count_tables)
    count_plots <- lapply(layout$comparison_order, function(comparison_id) {
      make_de_gsea_count_plot(count_tables[[comparison_id]], comparison_id, count_limit, count_config)
    })
    names(count_plots) <- layout$comparison_order
    combined_counts <- compose_comparison_grid(count_plots, layout, collect_guides = TRUE)
    combined_counts <- combined_counts & ggplot2::theme(legend.position = "top")
    files <- save_result_plot(
      combined_counts, file.path(output_dir, "deg_gsea_counts", safe_path_component(group_id), "deg_gsea_count_combined"),
      count_combined_size, save_config
    )
    plot_manifest[[length(plot_manifest) + 1L]] <- result_plot_manifest_row(
      "deg_gsea_count_combined", group_id, files = files, y_limits = c(-count_limit, count_limit)
    )
  }

  plot_manifest <- run_enrichment_plot_stage(
    project_root = project_root,
    groups = groups,
    entries_by_group = entries_by_group,
    plot_config = plot_config,
    save_config = save_config,
    output_dir = output_dir,
    manifest_rows = plot_manifest,
    default_gsea_p_adjust_cutoff = gsea_p_adjust_cutoff
  )

  plot_manifest <- run_classified_gsea_category_sankey_stage(
    config = config,
    project_root = project_root,
    output_dir = output_dir,
    save_config = save_config,
    manifest_rows = plot_manifest
  )

  if (isTRUE(trend_config$enabled %||% TRUE)) {
    selected_trend_ids <- as.character(unlist(trend_config$trend_ids %||% character(), use.names = FALSE))
    trend_entries <- read_successful_trend_sets_for_plot(project_root, selected_trend_ids)
    trend_groups <- unique(vapply(trend_entries, `[[`, character(1), "group_id"))
    configured_colours <- trend_config$group_colors %||% list()
    fallback_colours <- stats::setNames(grDevices::hcl.colors(length(trend_groups), "Dark 3"), trend_groups)
    for (entry in trend_entries) {
      fill_colour <- result_plot_colour(configured_colours, entry$group_id, fallback_colours[[entry$group_id]])
      plot <- make_trend_bar_plot(entry$table, entry$group_id, entry$trend_id, trend_config, fill_colour)
      files <- save_result_plot(
        plot,
        file.path(output_dir, "trend_bars", safe_path_component(entry$group_id), safe_path_component(entry$trend_id), "trend_bar"),
        trend_bar_size, save_config
      )
      plot_manifest[[length(plot_manifest) + 1L]] <- result_plot_manifest_row(
        "trend_bar", entry$group_id, trend_id = entry$trend_id, files = files
      )
    }
  }

  manifest <- do.call(rbind, plot_manifest)
  write_tsv(manifest, file.path(output_dir, "plot_manifest.tsv"))
  message("结果绘图已保存至：", output_dir)
  invisible(manifest)
}
