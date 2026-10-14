# Post-classification trend tracking for genes and classified GSEA pathways.
#
# GSEA categories only exist after the user has run the external category
# workflow.  This stage therefore calculates gene and pathway trends together,
# exports every non-neutral pattern, and draws alluvial/Sankey summaries.

post_category_abort <- function(...) abort(...)

post_category_number <- function(value, label, lower = -Inf, strict = FALSE) {
  value <- suppressWarnings(as.numeric(value))
  invalid <- length(value) != 1L || !is.finite(value) || if (strict) value <= lower else value < lower
  if (invalid) {
    post_category_abort(label, " 必须是", if (strict) "大于" else "大于等于", " ", lower, " 的有限数值。")
  }
  value
}

post_category_patterns <- function(states, n_comparisons) {
  grid <- do.call(expand.grid, c(rep(list(states), n_comparisons), stringsAsFactors = FALSE))
  apply(grid, 1, paste, collapse = "-")
}

# ---- Gene trend calculation (run only after GSEA classification) ---------

post_category_prepare_gene_comparison <- function(de_result, comparison_id, fc_cutoff, padj_cutoff) {
  required <- c("Gene", "log2FoldChange", "padj")
  missing <- setdiff(required, names(de_result))
  if (length(missing)) post_category_abort(comparison_id, " 的 DE 结果缺少列：", paste(missing, collapse = ", "))
  log2fc_cutoff <- log2(fc_cutoff)
  direction <- ifelse(
    !is.na(de_result$padj) & de_result$padj <= padj_cutoff & de_result$log2FoldChange >= log2fc_cutoff,
    "Up",
    ifelse(
      !is.na(de_result$padj) & de_result$padj <= padj_cutoff & de_result$log2FoldChange <= -log2fc_cutoff,
      "Down", "Stable"
    )
  )
  output <- data.frame(
    Gene = as.character(de_result$Gene),
    log2FoldChange = suppressWarnings(as.numeric(de_result$log2FoldChange)),
    padj = suppressWarnings(as.numeric(de_result$padj)),
    Direction = direction,
    stringsAsFactors = FALSE
  )
  output <- output[!is.na(output$Gene) & nzchar(output$Gene) & !duplicated(output$Gene), , drop = FALSE]
  names(output) <- c(
    "Gene",
    paste0("log2FC_", comparison_id),
    paste0("padj_", comparison_id),
    paste0("Direction_", comparison_id)
  )
  output
}

post_category_run_one_gene_trend <- function(de_manifest, trend_row, trend_config, project_root) {
  group_id <- trimws(as.character(trend_row$group_id[[1]]))
  trend_id <- trimws(as.character(trend_row$trend_id[[1]]))
  comparison_ids <- post_category_comparison_ids(trend_row)
  rows <- de_manifest[
    de_manifest$group_id == group_id & de_manifest$comparison_id %in% comparison_ids,
    , drop = FALSE
  ]
  successful_ids <- rows$comparison_id[rows$status == "success"]
  missing_ids <- setdiff(comparison_ids, successful_ids)
  if (length(missing_ids)) {
    post_category_abort(group_id, "/", trend_id, " 所需比较没有成功完成：", paste(missing_ids, collapse = ", "))
  }
  rows <- rows[match(comparison_ids, rows$comparison_id), , drop = FALSE]
  de_tables <- lapply(rows$result_file, function(path) read_tsv(file.path(project_root, path)))
  fc_cutoff <- post_category_number(trend_config$fc_cutoff %||% 2, "trend.fc_cutoff", 1)
  padj_cutoff <- post_category_number(trend_config$padj_cutoff %||% 0.05, "trend.padj_cutoff", 0, TRUE)
  if (padj_cutoff > 1) post_category_abort("trend.padj_cutoff 必须小于等于 1。")
  prepared <- Map(
    post_category_prepare_gene_comparison, de_tables, comparison_ids,
    MoreArgs = list(fc_cutoff = fc_cutoff, padj_cutoff = padj_cutoff)
  )
  trend_table <- Reduce(function(left, right) merge(left, right, by = "Gene", all = TRUE, sort = FALSE), prepared)
  direction_columns <- paste0("Direction_", comparison_ids)
  if (isTRUE(trend_config$missing_as_stable %||% TRUE)) {
    for (column in direction_columns) trend_table[[column]][is.na(trend_table[[column]])] <- "Stable"
  }
  trend_table$Full_Trend <- apply(trend_table[, direction_columns, drop = FALSE], 1, function(value) {
    if (anyNA(value)) NA_character_ else paste(value, collapse = "-")
  })
  summary_table <- post_category_summary(trend_table, c("Up", "Stable", "Down"), length(comparison_ids))
  output_dir <- file.path(project_root, "data", "trend_results", safe_path_component(group_id), safe_path_component(trend_id))
  write_tsv(trend_table, file.path(output_dir, "trend_table.tsv"))
  write_tsv(summary_table, file.path(output_dir, "trend_summary.tsv"))
  data.frame(
    group_id = group_id,
    trend_id = trend_id,
    status = "success",
    reason = NA_character_,
    n_comparisons = length(comparison_ids),
    n_genes = nrow(trend_table),
    trend_table_file = file.path("data", "trend_results", safe_path_component(group_id), safe_path_component(trend_id), "trend_table.tsv"),
    trend_summary_file = file.path("data", "trend_results", safe_path_component(group_id), safe_path_component(trend_id), "trend_summary.tsv"),
    stringsAsFactors = FALSE
  )
}

post_category_run_gene_trend_stage <- function(config, project_root, trend_sets) {
  trend_config <- config$trend %||% list()
  de_manifest <- read_manifest(project_root)
  rows <- lapply(seq_len(nrow(trend_sets)), function(index) {
    trend_row <- trend_sets[index, , drop = FALSE]
    message("Gene trend: ", trend_row$group_id[[1]], " / ", trend_row$trend_id[[1]])
    tryCatch(
      post_category_run_one_gene_trend(de_manifest, trend_row, trend_config, project_root),
      error = function(error) data.frame(
        group_id = trend_row$group_id[[1]], trend_id = trend_row$trend_id[[1]], status = "failed",
        reason = conditionMessage(error), n_comparisons = NA_integer_, n_genes = NA_integer_,
        trend_table_file = NA_character_, trend_summary_file = NA_character_, stringsAsFactors = FALSE
      )
    )
  })
  manifest <- do.call(rbind, rows)
  write_tsv(manifest, file.path(project_root, "data", "pipeline", "trend_manifest.tsv"))
  message("基因趋势统计完成：", sum(manifest$status == "success"), "/", nrow(manifest), " 个 trend set 成功。")
  manifest
}

post_category_safe_sheet_name <- function(value, existing = character()) {
  value <- as.character(value)
  for (forbidden in c("\\", "/", ":", "*", "?", "[", "]")) {
    value <- gsub(forbidden, "_", value, fixed = TRUE)
  }
  value <- trimws(value)
  if (!nzchar(value)) value <- "pattern"
  value <- substr(value, 1L, 31L)
  candidate <- value
  index <- 1L
  while (candidate %in% existing) {
    suffix <- paste0("_", index)
    candidate <- paste0(substr(value, 1L, max(1L, 31L - nchar(suffix))), suffix)
    index <- index + 1L
  }
  candidate
}

post_category_read_trend_sets <- function(config, project_root) {
  input <- config$input %||% list()
  path <- resolve_path(project_root, input$trend_sets_file %||% "data/trend_sets.tsv")
  sets <- read_tsv(path, c("group_id", "trend_id", "comparison_ids"))
  for (column in c("group_id", "trend_id", "comparison_ids")) {
    sets[[column]] <- trimws(as.character(sets[[column]]))
  }
  if (any(!nzchar(sets$group_id)) || any(!nzchar(sets$trend_id)) || any(!nzchar(sets$comparison_ids))) {
    post_category_abort("trend_sets_file 存在空的 group_id、trend_id 或 comparison_ids。")
  }
  keys <- paste(sets$group_id, sets$trend_id, sep = "\r")
  if (anyDuplicated(keys)) post_category_abort("trend_sets_file 中同一 group_id 内的 trend_id 必须唯一。")
  sets
}

post_category_comparison_ids <- function(trend_row) {
  value <- trimws(unlist(strsplit(as.character(trend_row$comparison_ids), ";", fixed = TRUE)))
  value <- value[nzchar(value)]
  if (length(value) < 2L) {
    post_category_abort(trend_row$group_id, "/", trend_row$trend_id, " 至少需要两个 comparison_ids。")
  }
  if (anyDuplicated(value)) post_category_abort(trend_row$group_id, "/", trend_row$trend_id, " 的 comparison_ids 存在重复。")
  value
}

post_category_select_trend_sets <- function(trend_sets, post_config) {
  groups <- trimws(as.character(unlist(post_config$groups %||% character(), use.names = FALSE)))
  trend_ids <- trimws(as.character(unlist(post_config$trend_ids %||% character(), use.names = FALSE)))
  if (length(groups)) trend_sets <- trend_sets[trend_sets$group_id %in% groups, , drop = FALSE]
  if (length(trend_ids)) trend_sets <- trend_sets[trend_sets$trend_id %in% trend_ids, , drop = FALSE]
  if (!nrow(trend_sets)) post_category_abort("post_category_trends 的 groups/trend_ids 没有匹配到 trend_sets。")
  trend_sets
}

post_category_read_gene_table <- function(project_root, group_id, trend_id, comparison_ids) {
  path <- file.path(project_root, "data", "trend_results", safe_path_component(group_id), safe_path_component(trend_id), "trend_table.tsv")
  table <- read_tsv(path, c("Gene", "Full_Trend"))
  direction_columns <- paste0("Direction_", comparison_ids)
  missing <- setdiff(direction_columns, names(table))
  if (length(missing)) {
    post_category_abort("基因趋势表缺少方向列：", paste(missing, collapse = ", "), "；文件：", path)
  }
  table$Feature_ID <- as.character(table$Gene)
  table$Description <- table$Feature_ID
  table <- table[!is.na(table$Feature_ID) & nzchar(table$Feature_ID), , drop = FALSE]
  if (anyDuplicated(table$Feature_ID)) post_category_abort("基因趋势表包含重复 Gene：", path)
  for (column in direction_columns) {
    table[[column]] <- trimws(as.character(table[[column]]))
    table[[column]][is.na(table[[column]]) | !nzchar(table[[column]])] <- "Stable"
  }
  table$Full_Trend <- apply(table[, direction_columns, drop = FALSE], 1, paste, collapse = "-")
  list(table = table, source_file = path)
}

post_category_read_classified_gsea <- function(project_root, post_config) {
  path <- resolve_path(project_root, post_config$classified_gsea_file %||% "table/enrichment_results/gsea_sig_group.merged.classified.csv")
  table <- read_tsv(path)
  columns <- post_config$gsea_columns %||% list()
  id_column <- as.character(columns$id %||% "ID")
  description_column <- as.character(columns$description %||% "Description")
  direction_column <- as.character(columns$direction %||% "Direction")
  required <- c("Group", "Comparison", id_column, description_column)
  missing <- setdiff(required, names(table))
  if (length(missing)) post_category_abort("分类 GSEA 文件缺少列：", paste(missing, collapse = ", "), "；文件：", path)
  output <- as.data.frame(table, stringsAsFactors = FALSE)
  output$Group <- trimws(as.character(output$Group))
  output$Comparison <- trimws(as.character(output$Comparison))
  output$Feature_ID <- trimws(as.character(output[[id_column]]))
  output$Description <- trimws(as.character(output[[description_column]]))
  if (direction_column %in% names(output)) {
    output$State <- trimws(as.character(output[[direction_column]]))
  } else if ("NES" %in% names(output)) {
    output$State <- ifelse(suppressWarnings(as.numeric(output$NES)) >= 0, "Positive", "Negative")
  } else {
    post_category_abort("分类 GSEA 文件缺少 ", direction_column, " 和 NES，无法确定通路方向。")
  }
  output$State <- ifelse(tolower(output$State) %in% c("positive", "pos", "up"), "Positive",
                         ifelse(tolower(output$State) %in% c("negative", "neg", "down"), "Negative", NA_character_))
  output <- output[
    !is.na(output$Group) & nzchar(output$Group) & !is.na(output$Comparison) & nzchar(output$Comparison) &
      !is.na(output$Feature_ID) & nzchar(output$Feature_ID) & !is.na(output$State), , drop = FALSE
  ]
  annotation_columns <- as.character(unlist(columns$annotation_columns %||% character(), use.names = FALSE))
  annotation_columns <- annotation_columns[annotation_columns %in% names(output)]
  list(table = output, annotation_columns = annotation_columns, source_file = path)
}

# Plot-only universe: use the unfiltered @result slot, not as.data.frame(),
# which can apply clusterProfiler's significance cutoffs to a gseaResult.
post_category_read_raw_gsea_ids <- function(project_root, raw_dir, group_id, comparison_id) {
  path <- file.path(raw_dir, safe_path_component(group_id), safe_path_component(comparison_id), "gsea_go_bp.rds")
  if (!file.exists(path)) post_category_abort("完整通路桑基图缺少原始 GSEA RDS：", path)
  object <- readRDS(path)
  if (is.null(object)) return(character())
  result <- if (base::isS4(object)) {
    # slotNames() can be empty when the S4 class definition has not yet been
    # loaded in this R session, while slot(object, "result") still works.
    tryCatch(methods::slot(object, "result"), error = function(error) NULL)
  } else if (is.list(object) && is.data.frame(object$result)) {
    object$result
  } else {
    post_category_abort("原始 GSEA RDS 不是含 result 表的对象：", path)
  }
  if (!is.data.frame(result) || !"ID" %in% names(result)) {
    post_category_abort("原始 GSEA RDS 的 result 缺少 ID 列：", path)
  }
  ids <- trimws(as.character(result$ID))
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    post_category_abort("原始 GSEA RDS 的 ID 存在空值或重复：", path)
  }
  ids
}

post_category_build_full_gsea_plot_table <- function(classified, group_id, comparison_ids,
                                                     project_root, full_config = list()) {
  raw_dir <- resolve_path(project_root, full_config$raw_gsea_dir %||% "data/enrichment_results")
  selected <- classified[classified$Group == group_id & classified$Comparison %in% comparison_ids, , drop = FALSE]
  keys <- paste(selected$Comparison, selected$Feature_ID, sep = "\r")
  if (anyDuplicated(keys)) post_category_abort("完整通路桑基图的显著结果含重复 Group/Comparison/ID：", group_id)
  raw_ids <- lapply(comparison_ids, function(comparison_id) {
    post_category_read_raw_gsea_ids(project_root, raw_dir, group_id, comparison_id)
  })
  names(raw_ids) <- comparison_ids
  universe <- sort(unique(unlist(raw_ids, use.names = FALSE)))
  if (!length(universe)) post_category_abort("完整通路桑基图的原始 GSEA 结果没有受检通路：", group_id)
  table <- data.frame(Feature_ID = universe, stringsAsFactors = FALSE)
  counts <- vector("list", length(comparison_ids))
  for (index in seq_along(comparison_ids)) {
    comparison_id <- comparison_ids[[index]]
    sig <- selected[selected$Comparison == comparison_id, , drop = FALSE]
    outside <- setdiff(sig$Feature_ID, raw_ids[[index]])
    if (length(outside)) {
      post_category_abort("显著通路不在对应比较的原始 GSEA RDS 中：", group_id, "/", comparison_id,
                          "；示例 ID：", paste(utils::head(outside, 5L), collapse = ", "))
    }
    # This plot intentionally has only three states.  A pathway in the raw
    # union but absent from this comparison is displayed as NS as requested.
    states <- rep("NS", length(universe))
    states[match(sig$Feature_ID, universe)] <- sig$State
    table[[paste0("Direction_", comparison_id)]] <- states
    counts[[index]] <- data.frame(
      Comparison = comparison_id, n_union = length(universe), n_tested = length(raw_ids[[index]]),
      n_significant = nrow(sig), n_nonsig = sum(states == "NS"),
      stringsAsFactors = FALSE
    )
  }
  direction_columns <- paste0("Direction_", comparison_ids)
  table$Full_Trend <- apply(table[, direction_columns, drop = FALSE], 1, paste, collapse = "-")
  list(table = table, counts = do.call(rbind, counts))
}

post_category_first_nonempty <- function(value) {
  value <- as.character(value)
  value <- value[!is.na(value) & nzchar(trimws(value))]
  if (!length(value)) NA_character_ else value[[1]]
}

post_category_build_gsea_table <- function(classified, group_id, comparison_ids, annotation_columns) {
  selected <- classified[classified$Group == group_id & classified$Comparison %in% comparison_ids, , drop = FALSE]
  direction_columns <- paste0("Direction_", comparison_ids)
  if (!nrow(selected)) {
    output <- data.frame(Feature_ID = character(), Description = character(), stringsAsFactors = FALSE)
    for (column in direction_columns) output[[column]] <- character()
    output$Full_Trend <- character()
    return(output)
  }
  duplicate_keys <- paste(selected$Feature_ID, selected$Comparison, sep = "\r")
  if (anyDuplicated(duplicate_keys)) {
    duplicated <- unique(duplicate_keys[duplicated(duplicate_keys)])
    post_category_abort("分类 GSEA 文件中同一 Group/Comparison/ID 出现重复：", paste(utils::head(duplicated, 5L), collapse = ", "))
  }
  ids <- sort(unique(selected$Feature_ID))
  output <- data.frame(Feature_ID = ids, stringsAsFactors = FALSE)
  annotations <- c("Description", annotation_columns)
  for (column in annotations) {
    values <- tapply(selected[[column]], selected$Feature_ID, post_category_first_nonempty)
    output[[column]] <- unname(values[output$Feature_ID])
  }
  for (index in seq_along(comparison_ids)) {
    comparison_id <- comparison_ids[[index]]
    values <- selected$State[selected$Comparison == comparison_id]
    names(values) <- selected$Feature_ID[selected$Comparison == comparison_id]
    output[[direction_columns[[index]]]] <- unname(values[output$Feature_ID])
    output[[direction_columns[[index]]]][is.na(output[[direction_columns[[index]]]])] <- "NS"
  }
  output$Full_Trend <- apply(output[, direction_columns, drop = FALSE], 1, paste, collapse = "-")
  output
}

post_category_summary <- function(feature_table, states, n_comparisons) {
  patterns <- post_category_patterns(states, n_comparisons)
  counts <- table(factor(feature_table$Full_Trend, levels = patterns), useNA = "no")
  output <- data.frame(Full_Trend = names(counts), Count = as.integer(counts), stringsAsFactors = FALSE)
  total <- sum(output$Count)
  output$Percentage <- if (total) round(output$Count / total * 100, 2) else 0
  output
}

post_category_write_gsea_results <- function(feature_table, summary_table, project_root, post_config, group_id, trend_id) {
  directory <- file.path(
    project_root, post_config$gsea_output_dir %||% "data/gsea_trend_results",
    safe_path_component(group_id), safe_path_component(trend_id)
  )
  write_tsv(feature_table, file.path(directory, "trend_table.tsv"))
  write_tsv(summary_table, file.path(directory, "trend_summary.tsv"))
  list(
    trend_table_file = file.path(post_config$gsea_output_dir %||% "data/gsea_trend_results", safe_path_component(group_id), safe_path_component(trend_id), "trend_table.tsv"),
    trend_summary_file = file.path(post_config$gsea_output_dir %||% "data/gsea_trend_results", safe_path_component(group_id), safe_path_component(trend_id), "trend_summary.tsv")
  )
}

post_category_export_patterns <- function(feature_table, summary_table, project_root, post_config, group_id, trend_id,
                                          feature_type, neutral_pattern) {
  export_config <- post_config$export %||% list()
  base_directory <- file.path(
    project_root, post_config$table_output_dir %||% "table/trend_feature_subsets",
    safe_path_component(group_id), safe_path_component(trend_id), feature_type
  )
  patterns <- summary_table$Full_Trend[summary_table$Full_Trend != neutral_pattern]
  if (isTRUE(export_config$write_individual_tsv %||% TRUE)) {
    for (pattern in patterns) {
      subset <- feature_table[feature_table$Full_Trend == pattern, , drop = FALSE]
      write_tsv(subset, file.path(base_directory, "patterns", paste0(safe_path_component(pattern), ".tsv")))
    }
  }
  workbook_file <- NA_character_
  if (isTRUE(export_config$write_xlsx %||% TRUE)) {
    require_package("openxlsx")
    workbook <- openxlsx::createWorkbook()
    existing <- character()
    if (isTRUE(export_config$include_summary_sheet %||% TRUE)) {
      openxlsx::addWorksheet(workbook, "Summary")
      openxlsx::writeData(workbook, "Summary", summary_table)
      existing <- "Summary"
    }
    for (pattern in patterns) {
      sheet <- post_category_safe_sheet_name(pattern, existing)
      openxlsx::addWorksheet(workbook, sheet)
      openxlsx::writeData(workbook, sheet, feature_table[feature_table$Full_Trend == pattern, , drop = FALSE])
      existing <- c(existing, sheet)
    }
    workbook_file <- file.path(base_directory, paste0(feature_type, "_trend_patterns.xlsx"))
    ensure_dir(dirname(workbook_file))
    openxlsx::saveWorkbook(workbook, workbook_file, overwrite = TRUE)
  }
  list(directory = base_directory, workbook_file = workbook_file)
}

post_category_make_sankey_long <- function(feature_table, comparison_ids) {
  direction_columns <- paste0("Direction_", comparison_ids)
  if (!nrow(feature_table)) {
    return(data.frame(Feature_ID = character(), Full_Trend = character(), Comparison = character(),
                      State = character(), Flow_Weight = integer(), stringsAsFactors = FALSE))
  }
  # Features with the same direction pattern trace the same route. Collapse
  # them once, then draw that route with its actual count as the y weight.
  # This preserves stratum/transition totals while avoiding thousands of
  # identical polygons on every highlighted plot.
  patterns <- as.character(feature_table$Full_Trend)
  unique_patterns <- unique(patterns)
  first_rows <- match(unique_patterns, patterns)
  pattern_counts <- tabulate(match(patterns, unique_patterns), nbins = length(unique_patterns))
  do.call(rbind, lapply(seq_along(comparison_ids), function(index) {
    data.frame(
      Feature_ID = unique_patterns,
      Full_Trend = unique_patterns,
      Comparison = comparison_ids[[index]],
      State = as.character(feature_table[[direction_columns[[index]]]][first_rows]),
      Flow_Weight = pattern_counts,
      stringsAsFactors = FALSE
    )
  }))
}

post_category_sankey_labels <- function(long_table, comparison_ids, state_order, neutral_state,
                                        state_labels, transition_symbols, include_transitions = TRUE) {
  # `post_category_make_sankey_plot()` converts these columns to factors to
  # control the alluvial order.  Keep the label calculations character based:
  # base::nzchar() does not accept factors on all R installations.
  long_table <- as.data.frame(long_table, stringsAsFactors = FALSE)
  long_table$Comparison <- as.character(long_table$Comparison)
  long_table$State <- as.character(long_table$State)
  complete <- expand.grid(Comparison = comparison_ids, State = state_order, stringsAsFactors = FALSE)
  counted <- stats::aggregate(long_table$Flow_Weight, by = list(Comparison = long_table$Comparison, State = long_table$State), FUN = sum)
  names(counted)[[3]] <- "n"
  labels <- merge(complete, counted, by = c("Comparison", "State"), all.x = TRUE, sort = FALSE)
  labels$n[is.na(labels$n)] <- 0L
  labels$transition_text <- NA_character_
  if (isTRUE(include_transitions) && length(comparison_ids) > 1L) for (index in seq.int(2L, length(comparison_ids))) {
    previous <- comparison_ids[[index - 1L]]
    current <- comparison_ids[[index]]
    old <- long_table[long_table$Comparison == previous, c("Feature_ID", "State", "Flow_Weight"), drop = FALSE]
    new <- long_table[long_table$Comparison == current, c("Feature_ID", "State"), drop = FALSE]
    names(old)[[2]] <- "Previous_State"
    names(new)[[2]] <- "State"
    transitions <- merge(old, new, by = "Feature_ID", all = FALSE)
    transitions <- transitions[
      transitions$Previous_State %in% names(transition_symbols) &
        transitions$State %in% names(transition_symbols), , drop = FALSE
    ]
    if (!nrow(transitions)) next
    transition_count <- stats::aggregate(transitions$Flow_Weight,
      by = list(Previous_State = transitions$Previous_State, State = transitions$State), FUN = sum
    )
    names(transition_count)[[3]] <- "n"
    for (state in setdiff(state_order, neutral_state)) {
      pieces <- transition_count[transition_count$State == state, , drop = FALSE]
      if (!nrow(pieces)) next
      pieces <- pieces[match(setdiff(state_order, neutral_state), pieces$Previous_State, nomatch = 0L), , drop = FALSE]
      pieces <- pieces[nzchar(as.character(pieces$Previous_State)), , drop = FALSE]
      text <- paste0(
        unname(transition_symbols[pieces$Previous_State]), "\u2192", unname(transition_symbols[pieces$State]), ": ", pieces$n,
        collapse = "\n"
      )
      labels$transition_text[labels$Comparison == current & labels$State == state] <- text
    }
  }
  labels$y <- NA_real_
  for (comparison_id in comparison_ids) {
    part <- labels[labels$Comparison == comparison_id, , drop = FALSE]
    part <- part[match(state_order, part$State), , drop = FALSE]
    # The first state is drawn at the top.  This is retained as a useful
    # tabular fallback; plots themselves use stat_stratum() below so their
    # positions are always calculated by the same ggalluvial stat as strata.
    reversed_n <- rev(part$n)
    part$y <- rev(cumsum(reversed_n) - reversed_n / 2)
    labels$y[match(paste(part$Comparison, part$State, sep = "\r"), paste(labels$Comparison, labels$State, sep = "\r"))] <- part$y
  }
  labels$label <- paste0(unname(state_labels[labels$State]), "\n(n=", labels$n, ")")
  has_transition <- !is.na(labels$transition_text) & nzchar(labels$transition_text)
  labels$label[has_transition] <- paste(labels$label[has_transition], labels$transition_text[has_transition], sep = "\n")
  labels
}

post_category_make_sankey_plot <- function(feature_table, comparison_ids, feature_type, group_id, trend_id,
                                           sankey_config, highlight_pattern = NULL, with_labels = TRUE,
                                           long_table = NULL) {
  require_package("ggalluvial")
  is_gene <- identical(feature_type, "genes")
  is_full_gsea <- identical(feature_type, "gsea_pathways_full_universe")
  states <- if (is_gene) c("Up", "Stable", "Down") else c("Positive", "NS", "Negative")
  neutral_state <- if (is_gene) "Stable" else "NS"
  colors_config <- (sankey_config$colors %||% list())[[if (is_gene) "gene" else "gsea"]] %||% list()
  colors <- if (is_gene) {
    c(Up = colors_config$up %||% "#D95F59", Stable = colors_config$stable %||% "#D9D9D9", Down = colors_config$down %||% "#3B6FB6")
  } else {
    c(Positive = colors_config$positive %||% "#8E6CBE", NS = colors_config$stable %||% "#D9D9D9", Negative = colors_config$negative %||% "#56A3A6")
  }
  state_labels <- stats::setNames(states, states)
  if (is_full_gsea) {
    state_labels["NS"] <- "NotSig"
  }
  transition_symbols <- if (is_gene) c(Up = "U", Down = "D") else c(Positive = "+", Negative = "-")
  long_table <- long_table %||% post_category_make_sankey_long(feature_table, comparison_ids)
  if (!nrow(long_table)) return(NULL)
  long_table$Comparison <- factor(long_table$Comparison, levels = comparison_ids)
  long_table$State <- factor(long_table$State, levels = states)
  width <- post_category_number(sankey_config$stratum_width %||% 0.48, "post_category_trends.sankey.stratum_width", 0, TRUE)
  if (width > 1) post_category_abort("post_category_trends.sankey.stratum_width 必须小于等于 1。")
  alpha <- post_category_number(sankey_config$flow_alpha %||% 0.38, "post_category_trends.sankey.flow_alpha", 0)
  if (alpha > 1) post_category_abort("post_category_trends.sankey.flow_alpha 必须小于等于 1。")
  p <- ggplot2::ggplot(long_table, ggplot2::aes(
    x = .data$Comparison, stratum = .data$State, alluvium = .data$Feature_ID, y = .data$Flow_Weight
  ))
  if (is.null(highlight_pattern)) {
    flow_fill_values <- colors
    p <- p + ggalluvial::geom_flow(
      ggplot2::aes(fill = .data$State), alpha = alpha, width = width, reverse = TRUE
    )
    subtitle <- if (is_full_gsea) "Union of raw GSEA pathways; non-significant or absent = NotSig" else "All trend patterns"
  } else {
    highlight_pattern <- as.character(highlight_pattern)
    highlighted <- long_table[long_table$Full_Trend == highlight_pattern, , drop = FALSE]
    if (!nrow(highlighted)) {
      post_category_abort("高亮趋势不存在于当前特征表：", highlight_pattern)
    }
    # Do not split flow data into background and highlighted geom layers.
    # ggalluvial computes each layer's vertical layout independently, which
    # would move the highlighted flow away from its real branch.  One flow
    # layer preserves the full-data layout while this class controls style.
    long_table$Flow_Class <- ifelse(
      long_table$Full_Trend == highlight_pattern, "Highlighted", "Background"
    )
    long_table$Flow_Class <- factor(long_table$Flow_Class, levels = c("Background", "Highlighted"))
    flow_fill_values <- c(
      colors,
      Background = sankey_config$background_flow_color %||% "grey85",
      Highlighted = sankey_config$highlight_color %||% "#F2A541"
    )
    outline_color <- as.character(sankey_config$highlight_outline_color %||% "#8C2D04")
    if (length(outline_color) != 1L || is.na(outline_color) || !nzchar(outline_color)) {
      post_category_abort("post_category_trends.sankey.highlight_outline_color 必须是一个颜色值。")
    }
    outline_width <- post_category_number(
      sankey_config$highlight_outline_width %||% 0.55,
      "post_category_trends.sankey.highlight_outline_width", 0
    )
    background_alpha <- post_category_number(
      sankey_config$background_flow_alpha %||% min(alpha, 0.22),
      "post_category_trends.sankey.background_flow_alpha", 0
    )
    highlight_alpha <- post_category_number(
      sankey_config$highlight_flow_alpha %||% max(alpha, 0.88),
      "post_category_trends.sankey.highlight_flow_alpha", 0
    )
    if (background_alpha > 1 || highlight_alpha > 1) {
      post_category_abort("post_category_trends.sankey 的 flow alpha 必须小于等于 1。")
    }
    p <- p +
      ggalluvial::geom_flow(
        ggplot2::aes(fill = .data$Flow_Class, alpha = .data$Flow_Class, colour = .data$Flow_Class),
        width = width, reverse = TRUE, linewidth = outline_width
      ) +
      ggplot2::scale_alpha_manual(
        values = c(Background = background_alpha, Highlighted = highlight_alpha), guide = "none"
      ) +
      ggplot2::scale_colour_manual(
        values = c(Background = NA_character_, Highlighted = outline_color), guide = "none"
      )
    highlighted_n <- sum(highlighted$Flow_Weight[as.character(highlighted$Comparison) == comparison_ids[[1]]])
    subtitle <- paste0("Highlighted trend: ", highlight_pattern, " (n=", highlighted_n, ")")
  }
  # The plot was initialised before Flow_Class is added in the highlighted
  # branch.  Refresh its inherited data before the layers are evaluated.
  p$data <- long_table
  p <- p +
    ggalluvial::geom_stratum(
      ggplot2::aes(fill = .data$State), width = width, color = "grey30", linewidth = 0.3, reverse = TRUE
    ) +
    # Flow_Class and State share the fill scale only in highlighted figures.
    # Show state colours in the legend; the subtitle identifies the orange
    # target branch explicitly.
    ggplot2::scale_fill_manual(values = flow_fill_values, breaks = states,
                               labels = unname(state_labels[states]), drop = FALSE) +
    ggplot2::scale_x_discrete(labels = comparison_ids, expand = c(0.08, 0.08)) +
    ggplot2::labs(
      title = paste0(group_id, ": ", trend_id, " — ",
                     if (is_gene) "genes" else if (is_full_gsea) "GSEA pathways (full universe)" else "GSEA pathways"),
      subtitle = subtitle,
      x = NULL, y = if (is_gene) "Number of genes" else "Number of pathways", fill = NULL
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      legend.position = "top", axis.line.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 15)
    )
  if (isTRUE(with_labels)) {
    labels <- post_category_sankey_labels(
      long_table, comparison_ids, states, neutral_state, state_labels, transition_symbols,
      include_transitions = is.null(highlight_pattern)
    )
    label_key <- paste(labels$Comparison, labels$State, sep = "\r")
    long_key <- paste(as.character(long_table$Comparison), as.character(long_table$State), sep = "\r")
    long_table$Stratum_Label <- labels$label[match(long_key, label_key)]
    if (anyNA(long_table$Stratum_Label)) {
      post_category_abort("无法为桑基图 strata 匹配标签。")
    }
    # stat_stratum() yields the exact strata centroids used by geom_stratum(),
    # avoiding any hand-calculated Y positions drifting from the columns.
    p <- p + ggalluvial::stat_stratum(
      data = long_table,
      mapping = ggplot2::aes(
        x = .data$Comparison, stratum = .data$State, alluvium = .data$Feature_ID,
        y = .data$Flow_Weight, label = .data$Stratum_Label
      ),
      geom = "text", inherit.aes = FALSE, reverse = TRUE,
      size = 3.25, lineheight = 0.90
    )
  }
  p
}

post_category_save_plot <- function(plot, output_base, size, save_config = list()) {
  formats <- tolower(as.character(unlist(save_config$formats %||% c("pdf", "png"), use.names = FALSE)))
  formats <- unique(formats[formats %in% c("pdf", "png")])
  if (!length(formats)) post_category_abort("post_category_trends.sankey.save.formats 至少应包含 pdf 或 png。")
  ensure_dir(dirname(output_base))
  output <- file.path(dirname(output_base), paste0(basename(output_base), ".", formats))
  for (index in seq_along(output)) {
    if (formats[[index]] == "pdf") {
      ggplot2::ggsave(output[[index]], plot, width = size$width, height = size$height, units = "in", device = "pdf", useDingbats = FALSE, limitsize = FALSE)
    } else {
      ggplot2::ggsave(output[[index]], plot, width = size$width, height = size$height, units = "in", dpi = size$dpi, limitsize = FALSE)
    }
  }
  output
}

post_category_sankey_size <- function(sankey_config) {
  spec <- sankey_config$save %||% list()
  list(
    width = post_category_number(spec$width %||% 10, "post_category_trends.sankey.save.width", 0, TRUE),
    height = post_category_number(spec$height %||% 7, "post_category_trends.sankey.save.height", 0, TRUE),
    dpi = post_category_number(spec$dpi %||% 300, "post_category_trends.sankey.save.dpi", 0, TRUE)
  )
}

post_category_write_sankeys <- function(feature_table, summary_table, project_root, post_config, group_id, trend_id,
                                        feature_type, comparison_ids, neutral_pattern) {
  sankey_config <- post_config$sankey %||% list()
  if (!isTRUE(sankey_config$enabled %||% TRUE) || !nrow(feature_table)) return(character())
  variants <- as.character(unlist(sankey_config$variants %||% c("with_labels", "without_labels"), use.names = FALSE))
  invalid <- setdiff(variants, c("with_labels", "without_labels"))
  if (length(invalid)) post_category_abort("post_category_trends.sankey.variants 不支持：", paste(invalid, collapse = ", "))
  variants <- unique(variants)
  size <- post_category_sankey_size(sankey_config)
  directory <- file.path(project_root, post_config$figure_output_dir %||% "figures/trend_sankey",
                         safe_path_component(group_id), safe_path_component(trend_id), feature_type)
  long_table <- post_category_make_sankey_long(feature_table, comparison_ids)
  output <- character()
  for (variant in variants) {
    plot <- post_category_make_sankey_plot(feature_table, comparison_ids, feature_type, group_id, trend_id, sankey_config,
                                           with_labels = identical(variant, "with_labels"), long_table = long_table)
    output <- c(output, post_category_save_plot(plot, file.path(directory, paste0("sankey_", variant)), size, sankey_config$save %||% list()))
  }
  if (isTRUE(sankey_config$write_highlighted_patterns %||% TRUE)) {
    patterns <- summary_table$Full_Trend[summary_table$Full_Trend != neutral_pattern]
    for (pattern in patterns) {
      if (!any(feature_table$Full_Trend == pattern)) next
      for (variant in variants) {
        plot <- post_category_make_sankey_plot(feature_table, comparison_ids, feature_type, group_id, trend_id, sankey_config,
                                               highlight_pattern = pattern, with_labels = identical(variant, "with_labels"),
                                               long_table = long_table)
        output <- c(output, post_category_save_plot(
          plot, file.path(directory, "highlighted", paste0(safe_path_component(pattern), "_", variant)), size,
          sankey_config$save %||% list()
        ))
      }
    }
  }
  output
}

post_category_write_full_gsea_sankeys <- function(plot_data, project_root, post_config,
                                                   group_id, trend_id, comparison_ids) {
  sankey_config <- post_config$sankey %||% list()
  full_config <- sankey_config$full_universe %||% list()
  plot_config <- sankey_config
  plot_config$variants <- full_config$variants %||% sankey_config$variants
  plot_config$save <- utils::modifyList(sankey_config$save %||% list(), full_config$save %||% list())
  # The full universe is only a plotting denominator; the significant-pool
  # trend patterns and their highlighted figures remain unchanged.
  plot_config$write_highlighted_patterns <- isTRUE(full_config$write_highlighted_patterns %||% FALSE)
  plot_post_config <- post_config
  plot_post_config$sankey <- plot_config
  summary <- post_category_summary(plot_data$table, c("Positive", "NS", "Negative"), length(comparison_ids))
  directory <- file.path(project_root, post_config$figure_output_dir %||% "figures/trend_sankey",
                         safe_path_component(group_id), safe_path_component(trend_id), "gsea_pathways_full_universe")
  count_file <- file.path(directory, "plot_state_counts.tsv")
  write_tsv(plot_data$counts, count_file)
  files <- post_category_write_sankeys(
    plot_data$table, summary, project_root, plot_post_config, group_id, trend_id,
    "gsea_pathways_full_universe", comparison_ids, paste(rep("NS", length(comparison_ids)), collapse = "-")
  )
  list(files = files, count_file = count_file)
}

post_category_full_gsea_manifest_row <- function(group_id, trend_id, status, reason = NA_character_,
                                                  n_pathways = NA_integer_, count_file = NA_character_,
                                                  plot_files = character()) {
  data.frame(
    group_id = group_id, trend_id = trend_id, status = status, reason = reason,
    n_pathways = n_pathways, count_file = count_file,
    plot_files = paste(plot_files, collapse = ";"), stringsAsFactors = FALSE
  )
}

post_category_manifest_row <- function(group_id, trend_id, feature_type, status, comparison_ids, n_features = NA_integer_,
                                       trend_table_file = NA_character_, trend_summary_file = NA_character_,
                                       workbook_file = NA_character_, sankey_files = character(),
                                       sankey_status = "not_requested", reason = NA_character_) {
  data.frame(
    group_id = group_id,
    trend_id = trend_id,
    feature_type = feature_type,
    status = status,
    n_comparisons = length(comparison_ids),
    n_features = n_features,
    trend_table_file = trend_table_file,
    trend_summary_file = trend_summary_file,
    workbook_file = workbook_file,
    sankey_files = paste(sankey_files, collapse = ";"),
    sankey_status = sankey_status,
    reason = reason,
    stringsAsFactors = FALSE
  )
}

run_post_category_trend_stage <- function(config, project_root) {
  post_config <- config$post_category_trends %||% list()
  if (!isTRUE(post_config$enabled %||% TRUE)) {
    message("post_category_trends.enabled=FALSE，跳过分类后趋势追踪。")
    return(invisible(NULL))
  }
  load_configured_ggplot2(config)
  trend_sets <- post_category_select_trend_sets(post_category_read_trend_sets(config, project_root), post_config)
  sankey_config <- post_config$sankey %||% list()
  sankey_requested <- isTRUE(sankey_config$enabled %||% TRUE)
  sankey_dependency_error <- if (sankey_requested && !requireNamespace("ggalluvial", quietly = TRUE)) {
    "缺少 R 包 'ggalluvial'；已保留趋势表和 Excel，但跳过桑基图。"
  } else {
    NA_character_
  }
  if (!is.na(sankey_dependency_error)) warning(sankey_dependency_error, call. = FALSE)
  # The gene branch is deliberately calculated here, rather than in the core
  # pipeline, so genes and classified GSEA pathways always use the same
  # Group × trend_id definitions and are exported together.
  if (isTRUE(post_config$run_gene_trend %||% TRUE)) {
    post_category_run_gene_trend_stage(config, project_root, trend_sets)
  }
  classified_info <- post_category_read_classified_gsea(project_root, post_config)
  rows <- list()
  full_rows <- list()
  full_config <- sankey_config$full_universe %||% list()
  full_requested <- sankey_requested && isTRUE(full_config$enabled %||% TRUE)
  for (index in seq_len(nrow(trend_sets))) {
    trend_row <- trend_sets[index, , drop = FALSE]
    group_id <- trend_row$group_id[[1]]
    trend_id <- trend_row$trend_id[[1]]
    comparison_ids <- post_category_comparison_ids(trend_row)
    message("Post-category trend: ", group_id, " / ", trend_id)
    feature_specs <- list(
      genes = list(states = c("Up", "Stable", "Down"), neutral = rep("Stable", length(comparison_ids))),
      gsea_pathways = list(states = c("Positive", "NS", "Negative"), neutral = rep("NS", length(comparison_ids)))
    )
    for (feature_type in names(feature_specs)) {
      spec <- feature_specs[[feature_type]]
      result <- tryCatch({
        if (feature_type == "genes") {
          gene_info <- post_category_read_gene_table(project_root, group_id, trend_id, comparison_ids)
          feature_table <- gene_info$table
          trend_table_file <- file.path("data", "trend_results", safe_path_component(group_id), safe_path_component(trend_id), "trend_table.tsv")
          trend_summary_file <- file.path("data", "trend_results", safe_path_component(group_id), safe_path_component(trend_id), "trend_summary.tsv")
        } else {
          feature_table <- post_category_build_gsea_table(classified_info$table, group_id, comparison_ids, classified_info$annotation_columns)
          summary_temp <- post_category_summary(feature_table, spec$states, length(comparison_ids))
          gsea_files <- post_category_write_gsea_results(feature_table, summary_temp, project_root, post_config, group_id, trend_id)
          trend_table_file <- gsea_files$trend_table_file
          trend_summary_file <- gsea_files$trend_summary_file
        }
        summary_table <- post_category_summary(feature_table, spec$states, length(comparison_ids))
        neutral_pattern <- paste(spec$neutral, collapse = "-")
        export_info <- post_category_export_patterns(feature_table, summary_table, project_root, post_config, group_id, trend_id,
                                                      feature_type, neutral_pattern)
        sankey_result <- if (!is.na(sankey_dependency_error)) {
          list(files = character(), reason = sankey_dependency_error)
        } else {
          tryCatch(
            list(
              files = post_category_write_sankeys(feature_table, summary_table, project_root, post_config, group_id, trend_id,
                                                   feature_type, comparison_ids, neutral_pattern),
              reason = NA_character_
            ),
            error = function(error) list(files = character(), reason = conditionMessage(error))
          )
        }
        sankey_status <- if (!sankey_requested) {
          "not_requested"
        } else if (is.na(sankey_result$reason)) {
          "success"
        } else {
          "failed"
        }
        row_status <- if (identical(sankey_status, "failed")) "partial_success" else "success"
        if (identical(sankey_status, "failed")) {
          message(
            "Sankey skipped [", group_id, " / ", trend_id, " / ", feature_type, "]: ",
            sankey_result$reason
          )
        }
        list(
          row = post_category_manifest_row(group_id, trend_id, feature_type, row_status, comparison_ids, nrow(feature_table),
                                           trend_table_file, trend_summary_file, export_info$workbook_file,
                                           sankey_result$files, sankey_status, sankey_result$reason),
          value = NULL
        )
      }, error = function(error) {
        message("Post-category trend failed [", group_id, " / ", trend_id, " / ", feature_type, "]: ", conditionMessage(error))
        list(row = post_category_manifest_row(group_id, trend_id, feature_type, "failed", comparison_ids, reason = conditionMessage(error)), value = NULL)
      })
      rows[[length(rows) + 1L]] <- result$row
    }
    full_row <- if (!full_requested) {
      post_category_full_gsea_manifest_row(group_id, trend_id, "not_requested")
    } else if (!is.na(sankey_dependency_error)) {
      post_category_full_gsea_manifest_row(group_id, trend_id, "failed", sankey_dependency_error)
    } else {
      tryCatch({
        plot_data <- post_category_build_full_gsea_plot_table(
          classified_info$table, group_id, comparison_ids, project_root, full_config
        )
        output <- post_category_write_full_gsea_sankeys(
          plot_data, project_root, post_config, group_id, trend_id, comparison_ids
        )
        message("Full-universe GSEA Sankey: ", group_id, " / ", trend_id,
                " (", nrow(plot_data$table), " pathways)")
        post_category_full_gsea_manifest_row(
          group_id, trend_id, "success", n_pathways = nrow(plot_data$table),
          count_file = output$count_file, plot_files = output$files
        )
      }, error = function(error) {
        message("Full-universe GSEA Sankey failed [", group_id, " / ", trend_id, "]: ", conditionMessage(error))
        post_category_full_gsea_manifest_row(group_id, trend_id, "failed", conditionMessage(error))
      })
    }
    full_rows[[length(full_rows) + 1L]] <- full_row
  }
  manifest <- do.call(rbind, rows)
  write_tsv(manifest, file.path(project_root, "data", "pipeline", "post_category_trend_manifest.tsv"))
  full_manifest <- do.call(rbind, full_rows)
  write_tsv(full_manifest, file.path(project_root, "data", "pipeline", "post_category_gsea_full_universe_plot_manifest.tsv"))
  data_success <- manifest$status %in% c("success", "partial_success")
  message("分类后趋势数据与表格导出完成：", sum(data_success), "/", nrow(manifest), " 项成功。")
  if (any(manifest$sankey_status == "failed")) {
    message("桑基图未完成：", sum(manifest$sankey_status == "failed"), " 项；详见 data/pipeline/post_category_trend_manifest.tsv 的 reason 列。")
  }
  if (full_requested) {
    message("完整通路池桑基图完成：", sum(full_manifest$status == "success"), "/", nrow(full_manifest), " 项。")
    if (any(full_manifest$status == "failed")) {
      message("完整通路池桑基图失败原因见 data/pipeline/post_category_gsea_full_universe_plot_manifest.tsv。")
    }
  }
  invisible(manifest)
}
