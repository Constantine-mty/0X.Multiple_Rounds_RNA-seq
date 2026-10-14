# Metadata-aware plotting for categorized GO-BP GSVA pathways.
#
# This module deliberately consumes the outputs generated after the manual
# GSEA-category step.  It does not modify src/gobp_category or choose GSEA
# representative pathways itself.

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
}

require_gobp_plot_packages <- function(config = list()) {
  load_configured_ggplot2(config)
  packages <- c("ggplot2", "pheatmap", "patchwork", "readxl", "readr", "dplyr", "yaml")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("GSVA 分类绘图缺少 R 包：", paste(missing, collapse = ", "), call. = FALSE)
  }
}

gobp_plot_abort <- function(...) stop(paste0(...), call. = FALSE)

as_nonempty_character <- function(x, label) {
  out <- trimws(as.character(unlist(x, use.names = FALSE)))
  if (anyNA(out) || any(!nzchar(out))) gobp_plot_abort(label, " 不能包含空值。")
  out
}

safe_plot_component <- function(x) {
  out <- trimws(as.character(x))
  out <- gsub("[/\\\\:*?\"<>|[:cntrl:]]+", "_", out)
  out <- gsub("[[:space:]]+", "_", out)
  out <- gsub("(^_+|_+$)", "", out)
  if (!nzchar(out)) gobp_plot_abort("类别名称不能只包含文件名保留字符。")
  out
}

read_classified_gsea_table <- function(path, category_column) {
  if (!file.exists(path)) gobp_plot_abort("找不到分类后的 GSEA 表：", path)
  tab <- readr::read_csv(path, show_col_types = FALSE, name_repair = "minimal")
  unnamed <- grepl("^(Unnamed: ?[0-9]*|\\.\\.\\.[0-9]+)$", names(tab))
  if (any(unnamed)) tab <- tab[, !unnamed, drop = FALSE]
  required <- c("Group", "ID", "Description", category_column)
  absent <- setdiff(required, names(tab))
  if (length(absent)) {
    gobp_plot_abort(basename(path), " 缺少必要列：", paste(absent, collapse = ", "))
  }
  tab <- as.data.frame(tab, stringsAsFactors = FALSE)
  for (column in required) tab[[column]] <- trimws(as.character(tab[[column]]))
  tab <- tab[!is.na(tab$Group) & nzchar(tab$Group) &
             !is.na(tab$ID) & nzchar(tab$ID) &
             !is.na(tab$Description) & nzchar(tab$Description) &
             !is.na(tab[[category_column]]) & nzchar(tab[[category_column]]), , drop = FALSE]
  if (!nrow(tab)) gobp_plot_abort("分类后的 GSEA 表不含可用于绘图的通路。")
  tab
}

ordered_present_values <- function(values, configured_order = NULL) {
  present <- unique(as.character(values[!is.na(values) & nzchar(values)]))
  if (!length(present)) return(character())
  if (is.null(configured_order) || !length(configured_order)) return(sort(present))
  configured_order <- as_nonempty_character(configured_order, "category_order")
  c(configured_order[configured_order %in% present], sort(setdiff(present, configured_order)))
}

normalise_colour_map <- function(colours) {
  if (is.null(colours) || !length(colours)) return(character())
  values <- unlist(colours, use.names = TRUE)
  if (is.null(names(values)) || any(!nzchar(names(values)))) {
    gobp_plot_abort("颜色配置必须是 '名称: 颜色值' 的命名映射。")
  }
  stats::setNames(as.character(values), names(values))
}

make_colour_map <- function(values, configured_colours = NULL, palette = "Dark 3") {
  values <- unique(as.character(values))
  values <- values[!is.na(values) & nzchar(values)]
  configured <- normalise_colour_map(configured_colours)
  output <- configured[intersect(names(configured), values)]
  missing <- setdiff(values, names(output))
  if (length(missing)) {
    generated <- if (length(missing) == 1L) "#4C78A8" else grDevices::hcl.colors(length(missing), palette = palette)
    output <- c(output, stats::setNames(generated, missing))
  }
  output[values]
}

resolve_cycle_levels <- function(group_metadata, cycle_column, cycle_order = NULL) {
  observed <- unique(as.character(group_metadata[[cycle_column]]))
  if (is.null(cycle_order) || !length(cycle_order)) return(observed)
  cycle_order <- as_nonempty_character(cycle_order, "cycle_order")
  if (anyDuplicated(cycle_order)) gobp_plot_abort("cycle_order 中有重复的轮次名称。")
  absent <- setdiff(observed, cycle_order)
  if (length(absent)) {
    gobp_plot_abort("cycle_order 缺少 metadata 中实际出现的轮次：", paste(absent, collapse = ", "))
  }
  cycle_order[cycle_order %in% observed]
}

prepare_group_gsva_scores <- function(scores, metadata, group_id, input_config, cycle_order = NULL) {
  if (is.data.frame(scores)) scores <- as.matrix(scores)
  if (!is.matrix(scores) || is.null(rownames(scores)) || is.null(colnames(scores))) {
    gobp_plot_abort("GSVA RDS 必须是行名为 GO ID、列名为样本名的矩阵。")
  }
  suppressWarnings(storage.mode(scores) <- "numeric")
  if (any(!is.finite(scores))) gobp_plot_abort("GSVA score matrix 含有非数值或缺失值。")

  sample_column <- input_config$sample_column
  group_column <- input_config$group_column
  cycle_column <- input_config$cycle_column
  group_metadata <- metadata[metadata[[group_column]] == group_id, , drop = FALSE]
  if (!nrow(group_metadata)) gobp_plot_abort("metadata 中不存在 Group='", group_id, "'。")
  present <- group_metadata[[sample_column]] %in% colnames(scores)
  if (!any(present)) {
    gobp_plot_abort("Group '", group_id, "' 的 metadata 样本均不在 GSVA score matrix 中。")
  }
  missing_scores <- group_metadata[[sample_column]][!present]
  if (length(missing_scores)) {
    warning("Group '", group_id, "' 有 ", length(missing_scores),
            " 个 metadata 样本未在 GSVA score matrix 中，将忽略：",
            paste(missing_scores, collapse = ", "), call. = FALSE)
  }
  group_metadata <- group_metadata[present, , drop = FALSE]
  cycle_levels <- resolve_cycle_levels(group_metadata, cycle_column, cycle_order)
  group_metadata$.plot_cycle <- factor(group_metadata[[cycle_column]], levels = cycle_levels, ordered = TRUE)
  group_metadata <- group_metadata[order(group_metadata$.plot_cycle, group_metadata[[sample_column]]), , drop = FALSE]
  sample_order <- group_metadata[[sample_column]]
  scores <- scores[, sample_order, drop = FALSE]

  annotation <- data.frame(
    Cycle = factor(as.character(group_metadata[[cycle_column]]), levels = cycle_levels, ordered = TRUE),
    row.names = sample_order,
    check.names = FALSE
  )
  list(scores = scores, metadata = group_metadata, annotation = annotation, cycle_levels = cycle_levels)
}

make_plot_matrix <- function(scores, terms, id_column = "ID", description_column = "Description") {
  if (!nrow(terms)) return(NULL)
  terms <- terms[!duplicated(terms[[id_column]]), , drop = FALSE]
  matched <- terms[[id_column]] %in% rownames(scores)
  missing <- terms[[id_column]][!matched]
  if (length(missing)) {
    warning(length(missing), " 条分类通路不在 GSVA score matrix 中，将不绘制。", call. = FALSE)
  }
  terms <- terms[matched, , drop = FALSE]
  if (!nrow(terms)) return(NULL)
  matrix <- scores[terms[[id_column]], , drop = FALSE]
  terms$.row_label <- make.unique(as.character(terms[[description_column]]))
  rownames(matrix) <- terms$.row_label
  list(matrix = matrix, terms = terms, missing_ids = missing)
}

calculate_heatmap_breaks <- function(matrix, n_colours, symmetric = TRUE) {
  values <- as.numeric(matrix)
  values <- values[is.finite(values)]
  if (!length(values)) gobp_plot_abort("没有可用于热图的有限 GSVA 值。")
  if (symmetric) {
    maximum <- max(abs(values))
    if (maximum == 0) maximum <- 1
    return(seq(-maximum, maximum, length.out = n_colours + 1L))
  }
  limits <- range(values)
  if (diff(limits) == 0) limits <- limits + c(-0.5, 0.5)
  seq(limits[[1]], limits[[2]], length.out = n_colours + 1L)
}

write_gsva_heatmap <- function(matrix, column_annotation, cycle_colours, output_file,
                               title, heatmap_config = list(), row_annotation = NULL,
                               row_colours = NULL, gaps_row = NULL) {
  palette <- unlist(heatmap_config$palette %||% c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"), use.names = FALSE)
  n_colours <- as.integer(heatmap_config$n_colours %||% 100L)
  colours <- grDevices::colorRampPalette(palette)(n_colours)
  breaks <- calculate_heatmap_breaks(matrix, n_colours, isTRUE(heatmap_config$symmetric %||% TRUE))
  extra_width <- suppressWarnings(as.numeric(heatmap_config$extra_width %||% 0))
  if (length(extra_width) != 1L || !is.finite(extra_width) || extra_width < 0) {
    gobp_plot_abort("heatmap.extra_width 必须是大于等于 0 的有限数值。")
  }
  # `extra_width` only adds blank canvas around the pheatmap grob. It does
  # not change cellwidth/cellheight, so cell geometry remains comparable.
  width <- max(as.numeric(heatmap_config$min_width %||% 7), ncol(matrix) * as.numeric(heatmap_config$cellwidth %||% 16) / 25 + 3) + extra_width
  height <- max(as.numeric(heatmap_config$min_height %||% 4), nrow(matrix) * as.numeric(heatmap_config$cellheight %||% 10) / 25 + 2)
  annotation_colours <- list(Cycle = cycle_colours)
  if (!is.null(row_annotation)) annotation_colours$Category <- row_colours
  # pheatmap delegates row clustering to hclust(), which requires at least two
  # pathways.  A Category with one GSVA-matched pathway is still informative,
  # so draw it in its original order instead of aborting the whole batch.
  cluster_rows <- isTRUE(heatmap_config$cluster_rows %||% TRUE) && nrow(matrix) >= 2L

  ensure <- dirname(output_file)
  if (!dir.exists(ensure)) dir.create(ensure, recursive = TRUE, showWarnings = FALSE)
  pheatmap::pheatmap(
    matrix,
    color = colours,
    breaks = breaks,
    annotation_col = column_annotation,
    annotation_row = row_annotation,
    annotation_colors = annotation_colours,
    cluster_rows = cluster_rows,
    cluster_cols = FALSE,
    gaps_row = gaps_row,
    border_color = heatmap_config$border_color %||% "grey85",
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize = as.numeric(heatmap_config$fontsize %||% 9),
    fontsize_row = as.numeric(heatmap_config$fontsize_row %||% 7),
    fontsize_col = as.numeric(heatmap_config$fontsize_col %||% 8),
    cellwidth = as.numeric(heatmap_config$cellwidth %||% 16),
    cellheight = as.numeric(heatmap_config$cellheight %||% 10),
    angle_col = heatmap_config$angle_col %||% "45",
    main = title,
    filename = output_file,
    width = width,
    height = height,
    silent = TRUE
  )
  invisible(output_file)
}

matrix_to_long <- function(matrix, terms, sample_annotation, category_column) {
  long <- data.frame(
    ID = rep(terms$ID, times = ncol(matrix)),
    Description = rep(terms$Description, times = ncol(matrix)),
    Category = rep(terms[[category_column]], times = ncol(matrix)),
    Sample = rep(colnames(matrix), each = nrow(matrix)),
    Score = as.vector(matrix),
    stringsAsFactors = FALSE
  )
  long$Time <- sample_annotation[long$Sample, "Cycle"]
  long$Time <- factor(as.character(long$Time), levels = levels(sample_annotation$Cycle), ordered = TRUE)
  long
}

summarise_gsva_trajectory <- function(long) {
  dplyr::summarise(
    dplyr::group_by(long, .data$ID, .data$Description, .data$Category, .data$Time),
    Mean = mean(.data$Score, na.rm = TRUE),
    SD = stats::sd(.data$Score, na.rm = TRUE),
    N = sum(is.finite(.data$Score)),
    SE = .data$SD / sqrt(.data$N),
    .groups = "drop"
  )
}

resolve_y_limits <- function(values, requested = NULL) {
  if (!is.null(requested) && length(requested)) {
    requested <- suppressWarnings(as.numeric(unlist(requested, use.names = FALSE)))
    if (length(requested) != 2L || any(!is.finite(requested)) || requested[[1]] >= requested[[2]]) {
      gobp_plot_abort("line_plot.y_limits 必须为递增的两个有限数值。")
    }
    return(requested)
  }
  values <- values[is.finite(values)]
  if (!length(values)) gobp_plot_abort("没有可用于折线图的有限 GSVA 值。")
  limits <- range(values)
  pad <- diff(limits) * 0.05
  if (pad == 0) pad <- max(abs(limits), 1) * 0.05
  limits + c(-pad, pad)
}

resolve_symmetric_y_limits <- function(values) {
  # Representative panels are compared directly, therefore they use one
  # symmetric GSVA scale around zero and a 5% presentation margin.
  values <- values[is.finite(values)]
  if (!length(values)) gobp_plot_abort("没有可用于代表通路折线图的有限 GSVA 均值。")
  maximum <- max(abs(range(values)))
  if (!is.finite(maximum) || maximum == 0) maximum <- 1
  expansion <- maximum * 0.05
  c(-maximum - expansion, maximum + expansion)
}

save_category_trajectory_plot <- function(long, summary, y_limits, output_file, title,
                                          line_config = list()) {
  pathway_order <- unique(as.character(summary$Description))
  long$Description <- factor(long$Description, levels = pathway_order)
  summary$Description <- factor(summary$Description, levels = pathway_order)
  pathway_colours <- make_colour_map(pathway_order, palette = "Dynamic")
  show_raw <- isTRUE(line_config$show_sample_points %||% TRUE)

  plot <- ggplot2::ggplot(summary, ggplot2::aes(x = .data$Time, y = .data$Mean,
                                                group = .data$Description, colour = .data$Description))
  if (show_raw) {
    plot <- plot + ggplot2::geom_point(
      data = long,
      ggplot2::aes(x = .data$Time, y = .data$Score, colour = .data$Description),
      inherit.aes = FALSE, alpha = 0.28, size = 1.15,
      position = ggplot2::position_jitter(width = 0.06, height = 0)
    )
  }
  plot <- plot +
    ggplot2::geom_line(linewidth = 0.75) +
    ggplot2::geom_point(size = 1.9) +
    ggplot2::scale_colour_manual(values = pathway_colours, drop = FALSE, guide = "none") +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::coord_cartesian(ylim = y_limits) +
    ggplot2::labs(title = title, x = "Cycle", y = "GSVA score") +
    ggplot2::theme_classic(base_size = as.numeric(line_config$base_size %||% 10)) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    )

  ensure <- dirname(output_file)
  if (!dir.exists(ensure)) dir.create(ensure, recursive = TRUE, showWarnings = FALSE)
  # Every Category uses a fixed canvas. Without a pathway legend, neither
  # pathway count nor description length can compress the plotting panel.
  width <- as.numeric(line_config$category_width %||% line_config$min_width %||% 8)
  height <- as.numeric(line_config$category_height %||% line_config$min_height %||% 4.5)
  if (length(width) != 1L || !is.finite(width) || width <= 0 ||
      length(height) != 1L || !is.finite(height) || height <= 0) {
    gobp_plot_abort("line_plot.category_width 和 category_height 必须为正的有限数值。")
  }
  ggplot2::ggsave(output_file, plot, width = width, height = height, units = "in", useDingbats = FALSE)
  invisible(plot)
}

get_group_specific_config <- function(config, group_id) {
  if (is.null(config) || !length(config)) return(NULL)
  config[[as.character(group_id)]] %||% NULL
}

parse_summary_pathway_text <- function(pathway_text, candidate_descriptions) {
  # gsea_summary.py joins descriptions with ", ".  GO descriptions themselves
  # can contain commas, so strsplit(pathway_text, ",") would corrupt terms.
  # Recover the original list by segmenting against descriptions known for this
  # Group × Category, preferring the longest valid prefix at each position.
  text <- trimws(as.character(pathway_text))
  candidates <- unique(trimws(as.character(candidate_descriptions)))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  candidates <- candidates[order(nchar(candidates), decreasing = TRUE)]
  if (!nzchar(text) || !length(candidates)) return(NULL)

  parse_remaining <- function(remaining) {
    if (!nzchar(remaining)) return(character())
    possible <- candidates[startsWith(remaining, candidates)]
    for (candidate in possible) {
      suffix <- substr(remaining, nchar(candidate) + 1L, nchar(remaining))
      if (!nzchar(suffix)) return(candidate)
      if (!startsWith(suffix, ", ")) next
      tail <- parse_remaining(substr(suffix, 3L, nchar(suffix)))
      if (!is.null(tail)) return(c(candidate, tail))
    }
    NULL
  }
  parse_remaining(text)
}

read_representative_pathways <- function(summary_workbook, group_id, classified_gsea,
                                         summary_category_column = "English Category",
                                         summary_pathway_column = "代表通路",
                                         classified_category_column = "Category_Merged",
                                         summary_sheet = NULL) {
  if (!file.exists(summary_workbook)) gobp_plot_abort("找不到 GSEA 汇总工作簿：", summary_workbook)
  sheet <- summary_sheet %||% group_id
  sheets <- readxl::excel_sheets(summary_workbook)
  if (!sheet %in% sheets) {
    gobp_plot_abort("GSEA 汇总工作簿中找不到 Group '", group_id, "' 对应的 sheet '", sheet,
                    "'。可用 sheet：", paste(sheets, collapse = ", "))
  }
  tab <- as.data.frame(readxl::read_excel(summary_workbook, sheet = sheet), stringsAsFactors = FALSE)
  missing <- setdiff(c(summary_category_column, summary_pathway_column), names(tab))
  if (length(missing)) {
    gobp_plot_abort("汇总 sheet '", sheet, "' 缺少列：", paste(missing, collapse = ", "))
  }
  tab$.sheet_row <- seq_len(nrow(tab))
  tab <- tab[!is.na(tab[[summary_category_column]]) & nzchar(trimws(as.character(tab[[summary_category_column]]))) &
             !is.na(tab[[summary_pathway_column]]) & nzchar(trimws(as.character(tab[[summary_pathway_column]]))), , drop = FALSE]
  if (!nrow(tab)) return(data.frame())

  group_table <- classified_gsea[classified_gsea$Group == group_id, , drop = FALSE]
  output <- vector("list", nrow(tab))
  for (index in seq_len(nrow(tab))) {
    category <- trimws(as.character(tab[[summary_category_column]][[index]]))
    pathway_text <- trimws(as.character(tab[[summary_pathway_column]][[index]]))
    category_table <- group_table[group_table[[classified_category_column]] == category, , drop = FALSE]
    pathways <- parse_summary_pathway_text(pathway_text, category_table$Description)
    if (is.null(pathways)) {
      warning("无法按分类 GSEA 表还原 Group '", group_id, "' 的代表通路列表：",
              category, " / ", pathway_text, call. = FALSE)
      next
    }
    rows <- lapply(seq_along(pathways), function(pathway_index) {
      pathway <- pathways[[pathway_index]]
      candidates <- category_table[
        category_table$Description == pathway,
        c("ID", "Description", classified_category_column), drop = FALSE
      ]
      candidates <- candidates[!duplicated(candidates$ID), , drop = FALSE]
      if (!nrow(candidates)) {
        warning("Group '", group_id, "' 的汇总代表通路未能在分类 GSEA 表中回查到：",
                category, " / ", pathway, call. = FALSE)
        return(NULL)
      }
      if (nrow(candidates) > 1L) {
        gobp_plot_abort("代表通路描述映射到多个 GO ID，无法可靠绘图：Group '", group_id,
                        "' / ", category, " / ", pathway)
      }
      data.frame(
        ID = candidates$ID,
        Description = candidates$Description,
        Category = category,
        summary_row = tab$.sheet_row[[index]],
        pathway_rank = pathway_index,
        stringsAsFactors = FALSE
      )
    })
    rows <- Filter(Negate(is.null), rows)
    if (length(rows)) output[[index]] <- do.call(rbind, rows)
  }
  output <- Filter(Negate(is.null), output)
  if (!length(output)) return(data.frame())
  output <- do.call(rbind, output)
  output[!duplicated(output$ID), , drop = FALSE]
}

save_representative_trajectory_plot <- function(summary, pathway_terms, y_limits, output_file,
                                                category_colours, line_config = list()) {
  # `pathway_terms` is taken directly from the representative heatmap, so this
  # Category sequence and colour coding exactly match its top-to-bottom order.
  category_order <- unique(as.character(pathway_terms$Category))
  cycle_levels <- levels(summary$Time)
  if (!length(cycle_levels)) cycle_levels <- unique(as.character(summary$Time))
  panels <- lapply(category_order, function(category) {
    category_terms <- pathway_terms[pathway_terms$Category == category, "ID", drop = TRUE]
    category_summary <- summary[summary$Category == category & summary$ID %in% category_terms, , drop = FALSE]
    category_summary$ID <- factor(category_summary$ID, levels = unique(category_terms))
    category_summary$Time <- factor(as.character(category_summary$Time), levels = cycle_levels, ordered = TRUE)
    category_summary <- category_summary[order(category_summary$ID, category_summary$Time), , drop = FALSE]
    category_colour <- unname(category_colours[category])
    if (length(category_colour) != 1L || is.na(category_colour)) category_colour <- "#4D4D4D"

    ggplot2::ggplot(
      category_summary,
      ggplot2::aes(x = .data$Time, y = .data$Mean, group = .data$ID)
    ) +
      ggplot2::geom_hline(
        yintercept = 0, linewidth = 0.35, linetype = "dashed", color = "grey65"
      ) +
      ggplot2::geom_line(color = category_colour, linewidth = 0.8, alpha = 0.65) +
      ggplot2::geom_point(color = category_colour, size = 1.8, alpha = 0.90) +
      ggplot2::scale_x_discrete(limits = cycle_levels, drop = FALSE) +
      ggplot2::scale_y_continuous(
        limits = y_limits, expand = ggplot2::expansion(mult = c(0, 0))
      ) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(x = NULL, y = "GSVA score") +
      ggplot2::theme_classic(base_size = as.numeric(line_config$base_size %||% 10)) +
      ggplot2::theme(
        plot.title = ggplot2::element_blank(),
        legend.position = "none",
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(size = 9, color = "black", margin = ggplot2::margin(r = 5)),
        axis.text.x = ggplot2::element_text(size = 9, color = "black"),
        axis.text.y = ggplot2::element_text(size = 8, color = "black"),
        axis.ticks = ggplot2::element_line(linewidth = 0.35, color = "black"),
        axis.line = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(t = 2, r = 2, b = 2, l = 2, unit = "mm")
      )
  })
  names(panels) <- category_order
  plot <- patchwork::wrap_plots(panels, ncol = 1)

  ensure <- dirname(output_file)
  if (!dir.exists(ensure)) dir.create(ensure, recursive = TRUE, showWarnings = FALSE)
  width <- as.numeric(line_config$representative_width %||% 2)
  panel_height <- as.numeric(line_config$representative_panel_height %||% 1.35)
  if (length(width) != 1L || !is.finite(width) || width <= 0 ||
      length(panel_height) != 1L || !is.finite(panel_height) || panel_height <= 0) {
    gobp_plot_abort("line_plot.representative_width 和 representative_panel_height 必须为正的有限数值。")
  }
  # Cairo gives the cleanest text embedding on the usual Linux analysis host.
  # Some local R installations advertise Cairo but cannot dynamically load it;
  # in that case retry with the standard PDF device instead of failing a run.
  cairo_save <- try(
    suppressWarnings(ggplot2::ggsave(
      output_file, plot, width = width, height = length(panels) * panel_height,
      units = "in", device = grDevices::cairo_pdf, limitsize = FALSE
    )),
    silent = TRUE
  )
  written <- !inherits(cairo_save, "try-error") && file.exists(output_file) &&
    isTRUE(file.info(output_file)$size > 0)
  if (!written) {
    ggplot2::ggsave(
      output_file, plot, width = width, height = length(panels) * panel_height,
      units = "in", device = "pdf", useDingbats = FALSE, limitsize = FALSE
    )
  }
  invisible(plot)
}

run_gobp_gsva_plot_stage <- function(config, project_root) {
  require_gobp_plot_packages(config)
  plot_config <- config$gobp_gsva_plots %||% list()
  if (!isTRUE(plot_config$enabled %||% TRUE)) {
    message("gobp_gsva_plots.enabled=FALSE，跳过 GSVA 分类通路绘图。")
    return(invisible(NULL))
  }
  if (is.null(config$input)) gobp_plot_abort("配置缺少 input 节。")
  input_config <- config$input
  required_input <- c("metadata_file", "sample_column", "group_column", "cycle_column")
  missing_input <- setdiff(required_input, names(input_config))
  if (length(missing_input)) gobp_plot_abort("input 节缺少：", paste(missing_input, collapse = ", "))

  required_plot <- c("classified_gsea_file", "summary_workbook", "gsva_by_group_dir", "output_dir")
  missing_plot <- required_plot[vapply(required_plot, function(x) is.null(plot_config[[x]]) || !nzchar(plot_config[[x]]), logical(1))]
  if (length(missing_plot)) {
    gobp_plot_abort("gobp_gsva_plots 节缺少：", paste(missing_plot, collapse = ", "))
  }
  category_column <- plot_config$category_column %||% "Category_Merged"
  representative_category_column <- plot_config$representative_category_column %||% "Category_Merged"
  metadata_path <- resolve_path(project_root, input_config$metadata_file)
  metadata <- read_metadata(metadata_path, input_config$sample_column, input_config$group_column, input_config$cycle_column)
  classified_path <- resolve_path(project_root, plot_config$classified_gsea_file)
  summary_path <- resolve_path(project_root, plot_config$summary_workbook)
  score_root <- resolve_path(project_root, plot_config$gsva_by_group_dir)
  output_root <- resolve_path(project_root, plot_config$output_dir)
  ensure_dir(output_root)
  classified <- read_classified_gsea_table(classified_path, category_column)
  if (!representative_category_column %in% names(classified)) {
    gobp_plot_abort("分类 GSEA 表缺少代表通路对应的类别列：", representative_category_column)
  }

  configured_groups <- plot_config$groups %||% character()
  groups <- if (length(configured_groups)) as_nonempty_character(configured_groups, "groups") else unique(classified$Group)
  unavailable_groups <- setdiff(groups, unique(classified$Group))
  if (length(unavailable_groups)) {
    gobp_plot_abort("这些 groups 不在分类 GSEA 表中：", paste(unavailable_groups, collapse = ", "))
  }

  heatmap_config <- plot_config$heatmap %||% list()
  line_config <- plot_config$line_plot %||% list()
  cycle_colours_config <- plot_config$cycle_colors %||% list()
  sheet_map <- plot_config$summary_sheet_by_group %||% list()
  y_limit_map <- line_config$y_limits %||% list()
  manifest <- list()

  for (group_id in groups) {
    score_path <- file.path(score_root, safe_path_component(group_id), "gsva_scores.rds")
    if (!file.exists(score_path)) {
      gobp_plot_abort("找不到 Group '", group_id, "' 的独立 GSVA 分数：", score_path,
                      "。请先运行 gsva.run_groupwise_scoring 分支。")
    }
    prepared <- prepare_group_gsva_scores(
      readRDS(score_path), metadata, group_id, input_config, plot_config$cycle_order %||% NULL
    )
    cycle_colours <- make_colour_map(prepared$cycle_levels, cycle_colours_config, palette = "Set 2")
    group_terms <- classified[classified$Group == group_id, c("ID", "Description", category_column), drop = FALSE]
    names(group_terms)[[3]] <- "Category"
    category_order <- ordered_present_values(group_terms$Category, plot_config$category_order %||% NULL)
    group_directory <- file.path(output_root, safe_plot_component(group_id))
    category_directory <- file.path(group_directory, "categories")
    ensure_dir(category_directory)

    all_category_matrix <- make_plot_matrix(prepared$scores, group_terms)
    if (is.null(all_category_matrix)) {
      warning("Group '", group_id, "' 没有分类 GSEA 通路与 GSVA score matrix 匹配。", call. = FALSE)
      next
    }
    y_limits <- resolve_y_limits(
      as.numeric(all_category_matrix$matrix), get_group_specific_config(y_limit_map, group_id)
    )

    for (category in category_order) {
      terms <- group_terms[group_terms$Category == category, , drop = FALSE]
      matrix_info <- make_plot_matrix(prepared$scores, terms)
      if (is.null(matrix_info)) next
      category_directory_name <- file.path(category_directory, safe_plot_component(category))
      ensure_dir(category_directory_name)
      write_tsv(matrix_info$terms[, c("ID", "Description", "Category"), drop = FALSE],
                file.path(category_directory_name, "pathway_mapping.tsv"))
      write_gsva_heatmap(
        matrix_info$matrix, prepared$annotation, cycle_colours,
        file.path(category_directory_name, "gsva_heatmap.pdf"),
        paste0(group_id, " - ", category), heatmap_config
      )
      long <- matrix_to_long(matrix_info$matrix, matrix_info$terms, prepared$annotation, "Category")
      trajectory <- summarise_gsva_trajectory(long)
      write_tsv(trajectory, file.path(category_directory_name, "gsva_trajectory.tsv"))
      save_category_trajectory_plot(
        long, trajectory, y_limits, file.path(category_directory_name, "gsva_trajectory.pdf"),
        paste0(group_id, " - ", category), line_config
      )
      manifest[[length(manifest) + 1L]] <- data.frame(
        group_id = group_id, plot_type = "category", category = category,
        n_pathways_selected = nrow(terms), n_pathways_plotted = nrow(matrix_info$terms),
        y_min = y_limits[[1]], y_max = y_limits[[2]], output_directory = category_directory_name,
        stringsAsFactors = FALSE
      )
    }

    representatives <- read_representative_pathways(
      summary_workbook = summary_path, group_id = group_id, classified_gsea = classified,
      summary_category_column = plot_config$summary_category_column %||% "English Category",
      summary_pathway_column = plot_config$summary_pathway_column %||% "代表通路",
      classified_category_column = representative_category_column,
      summary_sheet = get_group_specific_config(sheet_map, group_id)
    )
    if (!nrow(representatives)) {
      warning("Group '", group_id, "' 的 GSEA 汇总表没有可绘制的代表通路。", call. = FALSE)
      next
    }
    representative_order <- ordered_present_values(representatives$Category, plot_config$category_order %||% NULL)
    representatives$Category <- factor(representatives$Category, levels = representative_order, ordered = TRUE)
    representatives <- representatives[order(representatives$Category, representatives$summary_row, representatives$pathway_rank), , drop = FALSE]
    representatives$Category <- as.character(representatives$Category)
    representative_matrix <- make_plot_matrix(prepared$scores, representatives)
    if (is.null(representative_matrix)) {
      warning("Group '", group_id, "' 的代表通路均不在 GSVA score matrix 中。", call. = FALSE)
      next
    }
    representative_directory <- file.path(group_directory, "representative_pathways")
    ensure_dir(representative_directory)
    write_tsv(representative_matrix$terms[, c("ID", "Description", "Category", "summary_row", "pathway_rank"), drop = FALSE],
              file.path(representative_directory, "representative_pathway_mapping.tsv"))
    row_annotation <- data.frame(
      Category = factor(representative_matrix$terms$Category, levels = representative_order, ordered = TRUE),
      row.names = representative_matrix$terms$.row_label,
      check.names = FALSE
    )
    sizes <- table(factor(representative_matrix$terms$Category, levels = representative_order))
    sizes <- sizes[sizes > 0]
    gaps <- cumsum(as.integer(sizes))
    if (length(gaps)) gaps <- gaps[-length(gaps)]
    representative_colours <- make_colour_map(representative_order, plot_config$category_colors %||% list(), palette = "Dark 3")
    write_gsva_heatmap(
      representative_matrix$matrix, prepared$annotation, cycle_colours,
      file.path(representative_directory, "representative_pathways_heatmap.pdf"),
      paste0(group_id, " - representative GSEA pathways"),
      utils::modifyList(
        heatmap_config,
        list(
          cluster_rows = FALSE,
          extra_width = heatmap_config$representative_extra_width %||% 4
        )
      ),
      row_annotation = row_annotation, row_colours = representative_colours, gaps_row = gaps
    )
    representative_long <- matrix_to_long(
      representative_matrix$matrix, representative_matrix$terms, prepared$annotation, "Category"
    )
    representative_trajectory <- summarise_gsva_trajectory(representative_long)
    write_tsv(representative_trajectory, file.path(representative_directory, "representative_pathways_trajectory.tsv"))
    representative_y_limits <- resolve_symmetric_y_limits(representative_trajectory$Mean)
    save_representative_trajectory_plot(
      representative_trajectory, representative_matrix$terms, representative_y_limits,
      file.path(representative_directory, "representative_pathways_trajectory.pdf"),
      representative_colours, line_config
    )
    manifest[[length(manifest) + 1L]] <- data.frame(
      group_id = group_id, plot_type = "representative", category = NA_character_,
      n_pathways_selected = nrow(representatives), n_pathways_plotted = nrow(representative_matrix$terms),
      y_min = representative_y_limits[[1]], y_max = representative_y_limits[[2]], output_directory = representative_directory,
      stringsAsFactors = FALSE
    )
  }
  manifest <- if (length(manifest)) do.call(rbind, manifest) else data.frame()
  write_tsv(manifest, file.path(output_root, "plot_manifest.tsv"))
  message("GSVA 分类通路图已保存至：", output_root)
  invisible(manifest)
}
