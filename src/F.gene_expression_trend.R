# Project-local TPM expression trends for one gene or a gene list.
# Source this file, then call plot_gene_expression_trend() or
# plot_gene_expression_trend_grid(). No sample-name parsing is used.

gene_trend_abort <- function(...) stop(paste0(...), call. = FALSE)

gene_trend_read_tsv <- function(path, required_columns) {
  if (!file.exists(path)) gene_trend_abort("找不到文件：", path)
  table <- utils::read.delim(
    path, check.names = FALSE, stringsAsFactors = FALSE,
    quote = "\"", comment.char = "", na.strings = c("", "NA")
  )
  missing <- setdiff(required_columns, names(table))
  if (length(missing)) {
    gene_trend_abort(basename(path), " 缺少列：", paste(missing, collapse = ", "))
  }
  table
}

gene_trend_load_ggplot2 <- function(project_dir, ggplot2_lib_loc = NULL) {
  if (is.null(ggplot2_lib_loc)) {
    config_path <- file.path(project_dir, "config.yml")
    if (file.exists(config_path)) {
      if (!requireNamespace("yaml", quietly = TRUE)) {
        gene_trend_abort("读取 config.yml 中的 ggplot2 库配置需要 R 包 yaml。")
      }
      config <- yaml::read_yaml(config_path)
      ggplot2_lib_loc <- config$runtime$ggplot2_lib_loc
    }
  }
  if (is.null(ggplot2_lib_loc) || !length(ggplot2_lib_loc)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) gene_trend_abort("缺少 R 包 ggplot2。")
    return(invisible(NULL))
  }
  ggplot2_lib_loc <- as.character(ggplot2_lib_loc)
  if (length(ggplot2_lib_loc) != 1L || is.na(ggplot2_lib_loc) ||
      !dir.exists(ggplot2_lib_loc)) {
    gene_trend_abort("ggplot2_lib_loc 不是有效的 R 库目录。")
  }
  wanted <- normalizePath(file.path(ggplot2_lib_loc, "ggplot2"), mustWork = FALSE)
  if (!dir.exists(wanted)) gene_trend_abort("指定库中找不到 ggplot2：", wanted)
  if ("ggplot2" %in% loadedNamespaces()) {
    loaded <- normalizePath(getNamespaceInfo(asNamespace("ggplot2"), "path"), mustWork = FALSE)
    if (!identical(loaded, wanted)) {
      gene_trend_abort("ggplot2 已从其他库加载：", loaded, "；请在新 R 会话中使用：", wanted)
    }
  } else {
    loadNamespace("ggplot2", lib.loc = ggplot2_lib_loc)
  }
  invisible(NULL)
}

gene_trend_cycle_axis <- function(cycles, cycle_order = NULL) {
  observed <- unique(as.character(cycles))
  inferred_order <- is.null(cycle_order)
  if (is.null(cycle_order)) {
    levels <- observed
    if (all(grepl("^R[0-9]+$", levels))) {
      levels <- levels[order(as.integer(sub("^R", "", levels)), levels)]
    }
  } else {
    levels <- as.character(cycle_order)
    if (!length(levels) || anyNA(levels) || any(!nzchar(levels)) || anyDuplicated(levels) ||
        !setequal(levels, observed)) {
      gene_trend_abort("cycle_order 必须无重复且恰好包含 metadata 中的全部 Cycle。")
    }
  }
  positions <- seq_along(levels) - 1
  if (all(grepl("^R[0-9]+$", levels))) {
    round_numbers <- as.integer(sub("^R", "", levels))
    if (!anyDuplicated(round_numbers) &&
        (inferred_order || identical(order(round_numbers), seq_along(round_numbers)))) {
      positions <- round_numbers
    }
  }
  names(positions) <- levels
  list(levels = levels, positions = positions)
}

gene_trend_group_colours <- function(groups, group_colors = NULL) {
  default <- stats::setNames(grDevices::hcl.colors(length(groups), "Dark 3"), groups)
  if (is.null(group_colors)) return(default)
  values <- unlist(group_colors, use.names = TRUE)
  if (is.null(names(values)) || anyDuplicated(names(values)) ||
      any(!groups %in% names(values))) {
    gene_trend_abort("group_colors 必须为包含全部所选 Group 的命名颜色向量。")
  }
  values <- as.character(values[groups])
  if (anyNA(values) || any(!nzchar(values))) gene_trend_abort("group_colors 含有空颜色值。")
  invalid <- vapply(values, function(value) {
    inherits(try(grDevices::col2rgb(value), silent = TRUE), "try-error")
  }, logical(1))
  if (any(invalid)) gene_trend_abort("group_colors 中有无效颜色：", paste(values[invalid], collapse = ", "))
  stats::setNames(values, groups)
}

gene_trend_group_shapes <- function(groups, group_shapes = NULL) {
  default <- stats::setNames(rep(c(16L, 17L, 15L, 18L, 8L), length.out = length(groups)), groups)
  if (is.null(group_shapes)) return(default)
  values <- unlist(group_shapes, use.names = TRUE)
  if (is.null(names(values)) || anyDuplicated(names(values)) ||
      any(!groups %in% names(values))) {
    gene_trend_abort("group_shapes 必须为包含全部所选 Group 的命名数字向量。")
  }
  values <- suppressWarnings(as.numeric(values[groups]))
  if (anyNA(values) || any(!is.finite(values)) || any(values < 0 | values > 25) ||
      any(values != floor(values))) {
    gene_trend_abort("group_shapes 必须使用 ggplot2 的整数点形状编号 0–25。")
  }
  stats::setNames(as.integer(values), groups)
}

gene_trend_show_legend <- function(show_legend) {
  if (!is.logical(show_legend) || length(show_legend) != 1L || is.na(show_legend)) {
    gene_trend_abort("show_legend 必须为 TRUE 或 FALSE。")
  }
  show_legend
}

# Reads only the current project's data/metadata.tsv and data/tpm.tsv.
read_gene_expression_project <- function(project_dir = ".", groups = NULL, cycle_order = NULL,
                                         group_colors = NULL, sample_column = "Sample",
                                         group_column = "Group", cycle_column = "Cycle",
                                         gene_column = "Gene", ggplot2_lib_loc = NULL) {
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  gene_trend_load_ggplot2(project_dir, ggplot2_lib_loc)
  metadata <- gene_trend_read_tsv(
    file.path(project_dir, "data", "metadata.tsv"),
    c(sample_column, group_column, cycle_column)
  )
  metadata <- data.frame(
    Sample = trimws(as.character(metadata[[sample_column]])),
    Group = trimws(as.character(metadata[[group_column]])),
    Cycle = trimws(as.character(metadata[[cycle_column]])),
    stringsAsFactors = FALSE
  )
  for (column in names(metadata)) {
    if (anyNA(metadata[[column]]) || any(!nzchar(metadata[[column]]))) {
      gene_trend_abort("metadata.tsv 的 ", column, " 列有空值。")
    }
  }
  if (anyDuplicated(metadata$Sample)) gene_trend_abort("metadata.tsv 中 Sample 不能重复。")
  available_groups <- unique(metadata$Group)
  if (is.null(groups)) groups <- available_groups
  groups <- as.character(groups)
  if (!length(groups) || anyNA(groups) || any(!nzchar(groups)) || anyDuplicated(groups) ||
      any(!groups %in% available_groups)) {
    gene_trend_abort("groups 必须是不重复的、metadata 中已有的 Group。")
  }
  axis <- gene_trend_cycle_axis(metadata$Cycle, cycle_order)
  selected <- metadata[metadata$Group %in% groups, , drop = FALSE]
  selected$Group <- factor(selected$Group, levels = groups)
  selected$Cycle <- factor(selected$Cycle, levels = axis$levels, ordered = TRUE)
  selected$Cycle_X <- unname(axis$positions[as.character(selected$Cycle)])

  tpm <- gene_trend_read_tsv(file.path(project_dir, "data", "tpm.tsv"), gene_column)
  genes <- trimws(as.character(tpm[[gene_column]]))
  if (anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes)) {
    gene_trend_abort("tpm.tsv 的 ", gene_column, " 必须是非空且唯一的基因 ID。")
  }
  if (anyDuplicated(names(tpm))) gene_trend_abort("tpm.tsv 含有重复列名。")
  missing_samples <- setdiff(metadata$Sample, names(tpm))
  if (length(missing_samples)) {
    gene_trend_abort("metadata 中的样本不在 tpm.tsv：", paste(missing_samples, collapse = ", "))
  }
  values <- as.matrix(tpm[, selected$Sample, drop = FALSE])
  suppressWarnings(storage.mode(values) <- "numeric")
  if (any(!is.finite(values)) || any(values < 0)) {
    gene_trend_abort("tpm.tsv 的所选样本含有缺失、非数值或负数 TPM。")
  }
  rownames(values) <- genes
  list(project_dir = project_dir, metadata = selected, tpm = values, groups = groups,
       cycle_levels = axis$levels, cycle_positions = axis$positions,
       group_colors = gene_trend_group_colours(groups, group_colors))
}

gene_trend_means <- function(samples, groups, cycle_levels) {
  means <- stats::aggregate(TPM ~ Group + Cycle + Cycle_X, data = samples, FUN = mean,
                            drop = TRUE)
  means <- means[order(match(as.character(means$Group), groups),
                       match(as.character(means$Cycle), cycle_levels)), , drop = FALSE]
  rownames(means) <- NULL
  means
}

gene_trend_edges <- function(means, groups) {
  pieces <- lapply(groups, function(group) {
    part <- means[as.character(means$Group) == group, , drop = FALSE]
    part <- part[order(part$Cycle_X), , drop = FALSE]
    if (nrow(part) < 2L) return(NULL)
    data.frame(
      Group = group,
      Previous_Cycle = as.character(utils::head(part$Cycle, -1L)),
      Cycle = as.character(utils::tail(part$Cycle, -1L)),
      Previous_Mean = utils::head(part$TPM, -1L),
      Mean = utils::tail(part$TPM, -1L),
      X = (utils::head(part$Cycle_X, -1L) + utils::tail(part$Cycle_X, -1L)) / 2,
      stringsAsFactors = FALSE
    )
  })
  pieces <- Filter(Negate(is.null), pieces)
  if (!length(pieces)) return(data.frame())
  do.call(rbind, pieces)
}

gene_trend_safe_component <- function(value) {
  output <- gsub("[/\\\\]+", "_", as.character(value))
  if (!nzchar(output) || output %in% c(".", "..")) gene_trend_abort("无效的结果路径组成：", value)
  output
}

gene_trend_de_lookup <- function(project_data, edges) {
  if (!nrow(edges)) return(list())
  path <- file.path(project_data$project_dir, "data", "contrasts.tsv")
  contrasts <- gene_trend_read_tsv(path, c("group_id", "comparison_id", "numerator", "denominator"))
  manifest_path <- file.path(project_data$project_dir, "data", "pipeline", "de_manifest.tsv")
  manifest <- if (file.exists(manifest_path)) {
    gene_trend_read_tsv(manifest_path, c("group_id", "comparison_id", "status"))
  } else NULL
  for (column in c("group_id", "comparison_id", "numerator", "denominator")) {
    contrasts[[column]] <- trimws(as.character(contrasts[[column]]))
  }
  if ("enabled" %in% names(contrasts)) {
    enabled <- tolower(trimws(as.character(contrasts$enabled))) %in% c("true", "t", "1", "yes", "y")
    contrasts <- contrasts[enabled, , drop = FALSE]
  }
  lookup <- list()
  unavailable <- character()
  for (index in seq_len(nrow(edges))) {
    edge <- edges[index, , drop = FALSE]
    key <- paste(edge$Group, edge$Previous_Cycle, edge$Cycle, sep = "\r")
    if (!is.null(lookup[[key]])) next
    match_rows <- contrasts[
      contrasts$group_id == edge$Group & contrasts$denominator == edge$Previous_Cycle &
        contrasts$numerator == edge$Cycle, , drop = FALSE
    ]
    label <- paste0(edge$Group, " / ", edge$Previous_Cycle, "→", edge$Cycle)
    if (nrow(match_rows) > 1L) gene_trend_abort("contrasts.tsv 中比较定义重复：", label)
    if (!nrow(match_rows)) {
      unavailable <- c(unavailable, paste0(label, "（缺少同方向 contrast）"))
      next
    }
    if (!is.null(manifest)) {
      recorded <- manifest[manifest$group_id == edge$Group &
                             manifest$comparison_id == match_rows$comparison_id[[1]], , drop = FALSE]
      if (nrow(recorded) != 1L || !identical(as.character(recorded$status[[1]]), "success")) {
        unavailable <- c(unavailable, paste0(label, "（DE manifest 未记录成功）"))
        next
      }
    }
    result_file <- file.path(
      project_data$project_dir, "data", "de_results",
      gene_trend_safe_component(edge$Group),
      gene_trend_safe_component(match_rows$comparison_id[[1]]), "deseq2_results.tsv"
    )
    if (!file.exists(result_file)) {
      unavailable <- c(unavailable, paste0(label, "（缺少 DESeq2 结果文件）"))
      next
    }
    result <- gene_trend_read_tsv(result_file, c("Gene", "log2FoldChange", "padj"))
    result$Gene <- as.character(result$Gene)
    if (anyDuplicated(result$Gene)) gene_trend_abort("DESeq2 结果 Gene 重复：", result_file)
    result$log2FoldChange <- suppressWarnings(as.numeric(result$log2FoldChange))
    result$padj <- suppressWarnings(as.numeric(result$padj))
    lookup[[key]] <- result[, c("Gene", "log2FoldChange", "padj"), drop = FALSE]
  }
  if (length(unavailable)) {
    warning("以下区间没有可用的 DESeq2 比较，图中不标注其 log2FC：",
            paste(unique(unavailable), collapse = "；"), call. = FALSE)
  }
  lookup
}

gene_trend_format_fc <- function(value, digits) {
  if (is.na(value)) return("FC=NA")
  if (is.infinite(value)) return("FC=Inf")
  paste0("FC=", formatC(value, format = "f", digits = digits))
}

gene_trend_annotations <- function(edges, gene, annotation, pseudocount, de_lookup,
                                   padj_cutoff, log2fc_cutoff, significant_rule, digits) {
  if (!nrow(edges) || annotation == "none") return(data.frame())
  output <- edges
  output$Label <- NA_character_
  for (index in seq_len(nrow(output))) {
    if (annotation == "tpm_fc") {
      fold_change <- (output$Mean[[index]] + pseudocount) /
        (output$Previous_Mean[[index]] + pseudocount)
      value <- gene_trend_format_fc(fold_change, digits)
    } else {
      key <- paste(output$Group[[index]], output$Previous_Cycle[[index]], output$Cycle[[index]], sep = "\r")
      result <- de_lookup[[key]]
      if (is.null(result)) next
      row <- match(gene, result$Gene)
      log2fc <- if (is.na(row)) NA_real_ else result$log2FoldChange[[row]]
      padj <- if (is.na(row)) NA_real_ else result$padj[[row]]
      significant <- is.finite(log2fc) && is.finite(padj) && padj <= padj_cutoff &&
        (significant_rule == "padj" || abs(log2fc) >= log2fc_cutoff)
      value <- paste0("log2FC=", if (is.finite(log2fc)) formatC(log2fc, format = "f", digits = digits) else "NA",
                      if (significant) "*" else "")
    }
    output$Label[[index]] <- paste0(output$Previous_Cycle[[index]], "→", output$Cycle[[index]], "\n", value)
  }
  output[!is.na(output$Label), , drop = FALSE]
}

gene_trend_plot_one <- function(project_data, gene, line_type, annotation, pseudocount,
                                de_lookup, padj_cutoff, log2fc_cutoff, significant_rule,
                                digits, point_size, label_size, group_shapes, show_legend) {
  if (!gene %in% rownames(project_data$tpm)) gene_trend_abort("基因不在 tpm.tsv 中：", gene)
  metadata <- project_data$metadata
  samples <- metadata
  samples$TPM <- as.numeric(project_data$tpm[gene, metadata$Sample])
  means <- gene_trend_means(samples, project_data$groups, project_data$cycle_levels)
  group_cycle_counts <- table(as.character(means$Group))
  line_means <- means[as.character(means$Group) %in%
                        names(group_cycle_counts[group_cycle_counts >= 2L]), , drop = FALSE]
  edges <- gene_trend_edges(means, project_data$groups)
  smooth_groups <- names(group_cycle_counts[group_cycle_counts >= 3L])
  smooth_samples <- samples[as.character(samples$Group) %in% smooth_groups, , drop = FALSE]
  straight_fallback <- line_means[!as.character(line_means$Group) %in% smooth_groups, , drop = FALSE]
  annotations <- gene_trend_annotations(
    edges, gene, annotation, pseudocount, de_lookup, padj_cutoff,
    log2fc_cutoff, significant_rule, digits
  )
  y_max <- max(samples$TPM, na.rm = TRUE)
  y_span <- if (y_max > 0) y_max else 1
  if (nrow(annotations)) {
    annotations$Label_Y <- y_max +
      (length(project_data$groups) - match(annotations$Group, project_data$groups) + 1L) * 0.13 * y_span
  }
  y_upper <- if (nrow(annotations)) max(annotations$Label_Y) + 0.12 * y_span else y_max + 0.08 * y_span
  if (!is.finite(y_upper) || y_upper <= 0) y_upper <- 1

  p <- ggplot2::ggplot(samples, ggplot2::aes(x = .data$Cycle_X, y = .data$TPM,
                                              colour = .data$Group, shape = .data$Group)) +
    ggplot2::geom_point(size = point_size, alpha = 0.72,
                        position = ggplot2::position_jitter(width = 0.035, height = 0, seed = 1))
  if (line_type == "straight") {
    p <- p + ggplot2::geom_line(
      data = line_means, ggplot2::aes(x = .data$Cycle_X, y = .data$TPM, colour = .data$Group,
                                 group = .data$Group), inherit.aes = FALSE, linewidth = 1.05,
      show.legend = FALSE
    )
  } else {
    if (nrow(smooth_samples)) {
      # Match demo1: fit a local LOESS curve to the replicate TPM values in
      # each Group. A global natural spline can distort individual intervals.
      p <- p + ggplot2::geom_smooth(
        data = smooth_samples,
        ggplot2::aes(x = .data$Cycle_X, y = .data$TPM, colour = .data$Group,
                     group = .data$Group),
        inherit.aes = FALSE, method = "loess", formula = y ~ x, se = TRUE,
        alpha = 0.2, linewidth = 1.05, show.legend = FALSE
      )
    }
    if (nrow(straight_fallback)) {
      p <- p + ggplot2::geom_line(
        data = straight_fallback,
        ggplot2::aes(x = .data$Cycle_X, y = .data$TPM, colour = .data$Group,
                     group = .data$Group),
        inherit.aes = FALSE, linewidth = 1.05, show.legend = FALSE
      )
    }
  }
  p <- p + ggplot2::geom_point(
    data = means, ggplot2::aes(x = .data$Cycle_X, y = .data$TPM, colour = .data$Group),
    inherit.aes = FALSE, size = point_size + 0.55, shape = 21, fill = "white", stroke = 1.0,
    show.legend = FALSE
  )
  if (nrow(annotations)) {
    annotations$Group <- factor(annotations$Group, levels = project_data$groups)
    p <- p + ggplot2::geom_text(
      data = annotations,
      ggplot2::aes(x = .data$X, y = .data$Label_Y, label = .data$Label, colour = .data$Group),
      inherit.aes = FALSE, size = label_size, lineheight = 0.9, fontface = "bold",
      show.legend = FALSE
    )
  }
  x_positions <- unname(project_data$cycle_positions)
  x_span <- max(x_positions) - min(x_positions)
  x_padding <- if (x_span > 0) max(0.15, x_span * 0.04) else 0.5
  p + ggplot2::scale_colour_manual(values = project_data$group_colors, breaks = project_data$groups,
                                    drop = FALSE) +
    ggplot2::scale_shape_manual(values = group_shapes, breaks = project_data$groups, drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = x_positions, labels = project_data$cycle_levels,
                                limits = range(x_positions) + c(-x_padding, x_padding),
                                expand = ggplot2::expansion(mult = 0)) +
    ggplot2::scale_y_continuous(limits = c(0, y_upper), expand = ggplot2::expansion(mult = 0)) +
    ggplot2::labs(title = gene, x = "Cycle", y = "TPM", colour = "Group", shape = "Group") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = if (show_legend) "top" else "none",
                   plot.title = ggplot2::element_text(face = "bold"),
                   panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.45))
}

# line_type: "straight" joins observed cycle means; "curve" locally smooths
# replicate TPM within each group using LOESS (straight with <3 cycles).
# annotation: TPM mean ratio, DESeq2 log2FC (asterisk means padj rule passed),
# or no fold-change label. A missing DE contrast is never inferred from TPM.
plot_gene_expression_trend <- function(gene, project_dir = ".", groups = NULL,
                                       group_colors = NULL, line_type = c("straight", "curve"),
                                       annotation = c("tpm_fc", "deseq2_log2fc", "none"),
                                       cycle_order = NULL, pseudocount = 0,
                                       padj_cutoff = 0.05, log2fc_cutoff = 1,
                                       significant_rule = c("padj", "padj_and_log2fc"),
                                       digits = 2, point_size = 2.7, label_size = 2.8,
                                       ggplot2_lib_loc = NULL, group_shapes = NULL,
                                       show_legend = TRUE) {
  line_type <- match.arg(line_type)
  annotation <- match.arg(annotation)
  significant_rule <- match.arg(significant_rule)
  show_legend <- gene_trend_show_legend(show_legend)
  gene <- as.character(gene)
  if (length(gene) != 1L || is.na(gene) || !nzchar(gene)) gene_trend_abort("gene 必须是一个基因 ID。")
  gene_trend_validate_options(pseudocount, padj_cutoff, log2fc_cutoff, digits, point_size, label_size)
  data <- read_gene_expression_project(
    project_dir, groups, cycle_order, group_colors, ggplot2_lib_loc = ggplot2_lib_loc
  )
  if (!gene %in% rownames(data$tpm)) gene_trend_abort("基因不在 tpm.tsv 中：", gene)
  group_shapes <- gene_trend_group_shapes(data$groups, group_shapes)
  de_lookup <- if (annotation == "deseq2_log2fc") {
    gene_trend_de_lookup(data, gene_trend_edges(gene_trend_means_for_metadata(data), data$groups))
  } else list()
  gene_trend_plot_one(data, gene, line_type, annotation, pseudocount, de_lookup,
                      padj_cutoff, log2fc_cutoff, significant_rule, digits, point_size, label_size,
                      group_shapes, show_legend)
}

gene_trend_means_for_metadata <- function(project_data) {
  dummy <- project_data$metadata
  dummy$TPM <- 0
  gene_trend_means(dummy, project_data$groups, project_data$cycle_levels)
}

gene_trend_validate_options <- function(pseudocount, padj_cutoff, log2fc_cutoff,
                                        digits, point_size, label_size) {
  if (length(pseudocount) != 1L || !is.finite(pseudocount) || pseudocount < 0) {
    gene_trend_abort("pseudocount 必须为非负有限数值。")
  }
  if (length(padj_cutoff) != 1L || !is.finite(padj_cutoff) || padj_cutoff <= 0 || padj_cutoff > 1) {
    gene_trend_abort("padj_cutoff 必须在 (0, 1] 内。")
  }
  if (length(log2fc_cutoff) != 1L || !is.finite(log2fc_cutoff) || log2fc_cutoff < 0) {
    gene_trend_abort("log2fc_cutoff 必须为非负有限数值。")
  }
  if (length(digits) != 1L || !is.finite(digits) || digits < 0 || digits != as.integer(digits)) {
    gene_trend_abort("digits 必须为非负整数。")
  }
  if (length(point_size) != 1L || !is.finite(point_size) || point_size <= 0 ||
      length(label_size) != 1L || !is.finite(label_size) || label_size <= 0) {
    gene_trend_abort("point_size 与 label_size 必须为正数。")
  }
}

gene_trend_sort_genes <- function(genes, tpm, sort) {
  if (sort == "none") return(genes)
  mean_tpm <- rowMeans(tpm[genes, , drop = FALSE])
  index <- seq_along(genes)
  if (sort == "descending") index <- order(-mean_tpm, index)
  if (sort == "ascending") index <- order(mean_tpm, index)
  genes[index]
}

# Multi-gene wrapper. TPM, metadata and DESeq2 result tables are loaded once.
plot_gene_expression_trend_grid <- function(genes, project_dir = ".", ncol = 3L,
                                            groups = NULL, group_colors = NULL,
                                            line_type = c("straight", "curve"),
                                            annotation = c("tpm_fc", "deseq2_log2fc", "none"),
                                            cycle_order = NULL, pseudocount = 0,
                                            padj_cutoff = 0.05, log2fc_cutoff = 1,
                                            significant_rule = c("padj", "padj_and_log2fc"),
                                            digits = 2, point_size = 2.7, label_size = 2.8,
                                            missing_gene = c("warn", "error"), ggplot2_lib_loc = NULL,
                                            sort = c("none", "descending", "ascending"),
                                            group_shapes = NULL, show_legend = TRUE) {
  line_type <- match.arg(line_type)
  annotation <- match.arg(annotation)
  significant_rule <- match.arg(significant_rule)
  missing_gene <- match.arg(missing_gene)
  sort <- match.arg(sort)
  show_legend <- gene_trend_show_legend(show_legend)
  if (length(ncol) != 1L || !is.finite(ncol) || ncol < 1 || ncol != as.integer(ncol)) {
    gene_trend_abort("ncol 必须为正整数。")
  }
  gene_trend_validate_options(pseudocount, padj_cutoff, log2fc_cutoff, digits, point_size, label_size)
  genes <- unique(trimws(as.character(genes)))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (!length(genes)) gene_trend_abort("genes 不能为空。")
  data <- read_gene_expression_project(
    project_dir, groups, cycle_order, group_colors, ggplot2_lib_loc = ggplot2_lib_loc
  )
  missing <- setdiff(genes, rownames(data$tpm))
  if (length(missing)) {
    if (missing_gene == "error") gene_trend_abort("基因不在 tpm.tsv 中：", paste(missing, collapse = ", "))
    warning("跳过 tpm.tsv 中不存在的基因：", paste(missing, collapse = ", "), call. = FALSE)
    genes <- setdiff(genes, missing)
  }
  if (!length(genes)) gene_trend_abort("指定基因均不在 tpm.tsv 中。")
  genes <- gene_trend_sort_genes(genes, data$tpm, sort)
  group_shapes <- gene_trend_group_shapes(data$groups, group_shapes)
  de_lookup <- if (annotation == "deseq2_log2fc") {
    gene_trend_de_lookup(data, gene_trend_edges(gene_trend_means_for_metadata(data), data$groups))
  } else list()
  plots <- lapply(genes, function(gene) {
    gene_trend_plot_one(data, gene, line_type, annotation, pseudocount, de_lookup,
                        padj_cutoff, log2fc_cutoff, significant_rule, digits, point_size, label_size,
                        group_shapes, show_legend)
  })
  names(plots) <- genes
  if (!requireNamespace("patchwork", quietly = TRUE)) gene_trend_abort("多基因拼图需要 R 包 patchwork。")
  patchwork::wrap_plots(plots, ncol = as.integer(ncol), guides = "collect") &
    ggplot2::theme(legend.position = if (show_legend) "bottom" else "none")
}
