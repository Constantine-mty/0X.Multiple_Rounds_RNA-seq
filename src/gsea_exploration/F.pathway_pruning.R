# Coverage-constrained pruning of significant classified GSEA pathways.
# Source F.common.R before this file.

gpp_category_key <- function(x) {
  value <- tolower(as.character(x))
  value[is.na(value)] <- ""
  gsub("[^[:alnum:]]", "", value)
}

gpp_target_mask <- function(df, target_background = TRUE,
                            target_categories = c("Other_Uncategorized"),
                            target_merged_categories = character(),
                            target_pairs = NULL) {
  target <- rep(FALSE, nrow(df))
  if (target_background) {
    gsea_nb_require(df, "Background", "分类 GSEA")
    target <- target | tolower(trimws(as.character(df$Background))) %in% c("true", "t", "1")
  }
  if (length(target_categories)) {
    gsea_nb_require(df, "Category", "分类 GSEA")
    target <- target | gpp_category_key(df$Category) %in% gpp_category_key(target_categories)
  }
  if (length(target_merged_categories)) {
    gsea_nb_require(df, "Category_Merged", "分类 GSEA")
    target <- target | gpp_category_key(df$Category_Merged) %in% gpp_category_key(target_merged_categories)
  }
  if (!is.null(target_pairs)) {
    gsea_nb_require(target_pairs, c("Category", "Category_Merged"), "target_pairs")
    gsea_nb_require(df, c("Category", "Category_Merged"), "分类 GSEA")
    for (i in seq_len(nrow(target_pairs))) {
      target <- target | (gpp_category_key(df$Category) == gpp_category_key(target_pairs$Category[[i]]) &
                            gpp_category_key(df$Category_Merged) == gpp_category_key(target_pairs$Category_Merged[[i]]))
    }
  }
  target
}

gpp_prune <- function(gsea, min_coverage = 0.75,
                      target_background = TRUE,
                      target_categories = c("Other_Uncategorized"),
                      target_merged_categories = character(),
                      target_pairs = NULL,
                      group_columns = c("Group", "Comparison", "Direction")) {
  if (!is.data.frame(gsea) || !is.finite(min_coverage) ||
      min_coverage <= 0 || min_coverage > 1) {
    stop("gsea 必须是 data.frame，min_coverage 必须在 (0, 1] 内。", call. = FALSE)
  }
  df <- as.data.frame(gsea, stringsAsFactors = FALSE)
  if (!"Direction" %in% names(df) && "Direction" %in% group_columns &&
      "GSEA_Direction" %in% names(df)) df$Direction <- df$GSEA_Direction
  gsea_nb_require(df, c(group_columns, "core_enrichment", "p.adjust", "NES"), "分类 GSEA")
  if (!nrow(df)) stop("分类 GSEA 没有行。", call. = FALSE)
  if (anyNA(df[, group_columns, drop = FALSE])) stop("分组列不可包含 NA。", call. = FALSE)
  df$.original_row <- seq_len(nrow(df))
  df$.target <- gpp_target_mask(df, target_background, target_categories,
                                 target_merged_categories, target_pairs)
  df$.keep <- FALSE
  df$.selection_order <- NA_integer_
  df$.new_genes <- NA_integer_
  df$.decision <- NA_character_
  group_key <- do.call(paste, c(df[group_columns], sep = "\r"))
  groups <- split(seq_len(nrow(df)), factor(group_key, levels = unique(group_key)))
  summary_rows <- vector("list", length(groups))
  for (k in seq_along(groups)) {
    idx <- groups[[k]]
    sub <- df[idx, , drop = FALSE]
    sets <- lapply(sub$core_enrichment, gsea_nb_core_genes)
    target <- which(sub$.target)
    base <- which(!sub$.target)
    all_genes <- unique(unlist(sets, use.names = FALSE))
    base_genes <- unique(unlist(sets[base], use.names = FALSE))
    keep <- !sub$.target
    order_selected <- rep(NA_integer_, nrow(sub))
    new_genes <- rep(NA_integer_, nrow(sub))
    if (!length(all_genes)) {
      keep[] <- TRUE
      base_coverage <- NA_real_
      final_coverage <- NA_real_
    } else {
      covered <- base_genes
      base_coverage <- length(covered) / length(all_genes)
      remaining <- target
      step <- 0L
      while (length(covered) / length(all_genes) < min_coverage && length(remaining)) {
        gains <- vapply(remaining, function(i) length(setdiff(sets[[i]], covered)), integer(1))
        if (!any(gains > 0L)) break
        padj <- as.numeric(sub$p.adjust[remaining]); padj[!is.finite(padj)] <- Inf
        nes <- abs(as.numeric(sub$NES[remaining])); nes[!is.finite(nes)] <- -Inf
        best_pos <- order(-gains, padj, -nes, sub$.original_row[remaining])[[1]]
        best <- remaining[[best_pos]]
        step <- step + 1L
        keep[[best]] <- TRUE
        order_selected[[best]] <- step
        new_genes[[best]] <- gains[[best_pos]]
        covered <- union(covered, sets[[best]])
        remaining <- remaining[-best_pos]
      }
      final_coverage <- length(covered) / length(all_genes)
    }
    if (!length(all_genes)) covered <- all_genes
    df$.keep[idx] <- keep
    df$.selection_order[idx] <- order_selected
    df$.new_genes[idx] <- new_genes
    df$.decision[idx] <- ifelse(!sub$.target, "base_keep",
                                ifelse(!length(all_genes), "no_core_genes_keep",
                                       ifelse(keep, "target_selected", "target_removed")))
    info <- sub[1, group_columns, drop = FALSE]
    summary_rows[[k]] <- cbind(info,
      data.frame(n_pathways_before = nrow(sub), n_pathways_after = sum(keep),
                 n_target_before = length(target), n_target_after = sum(keep[sub$.target]),
                 n_genes_before = length(all_genes), n_genes_after = length(covered),
                 n_genes_lost = length(all_genes) - length(covered),
                 base_coverage = base_coverage, final_coverage = final_coverage,
                 target_removal_rate = if (length(target)) 1 - sum(keep[sub$.target]) / length(target) else NA_real_,
                 coverage_reached = is.na(final_coverage) || final_coverage + 1e-12 >= min_coverage,
                 stringsAsFactors = FALSE))
  }
  df <- df[order(df$.original_row), , drop = FALSE]
  df$.original_row <- NULL
  rownames(df) <- NULL
  summary <- do.call(rbind, summary_rows)
  rownames(summary) <- NULL
  list(annotated = df, filtered = df[df$.keep, , drop = FALSE], summary = summary,
       settings = list(min_coverage = min_coverage, target_background = target_background,
                       target_categories = target_categories,
                       target_merged_categories = target_merged_categories,
                       target_pairs = target_pairs, group_columns = group_columns))
}

gpp_plot_composition <- function(result, group, category_column = "Category_Merged",
                                 category_colors = NULL, comparison_order = NULL) {
  df <- result$annotated
  gsea_nb_require(df, c("Group", "Comparison", "Direction", category_column, ".keep"))
  df <- df[as.character(df$Group) == group, , drop = FALSE]
  if (!nrow(df)) stop("没有 Group：", group, call. = FALSE)
  before <- df; before$Filter <- "Before"
  after <- df[df$.keep, , drop = FALSE]; after$Filter <- "After"
  plot_df <- rbind(before, after)
  plot_df$Category_plot <- as.character(plot_df[[category_column]])
  plot_df$Category_plot[is.na(plot_df$Category_plot) | !nzchar(plot_df$Category_plot)] <- "Uncategorized"
  categories <- unique(plot_df$Category_plot)
  if (is.null(category_colors)) category_colors <- setNames(grDevices::hcl.colors(length(categories), "Dynamic"), categories)
  missing <- setdiff(categories, names(category_colors))
  if (length(missing)) stop("category_colors 缺少类别：", paste(missing, collapse = ", "), call. = FALSE)
  plot_df$Category_plot <- factor(plot_df$Category_plot, levels = categories)
  if (is.null(comparison_order)) comparison_order <- unique(as.character(plot_df$Comparison))
  plot_df$Comparison <- factor(as.character(plot_df$Comparison), levels = comparison_order)
  plot_df$Direction <- factor(as.character(plot_df$Direction), levels = c("Positive", "Negative"))
  plot_df$Filter <- factor(plot_df$Filter, levels = c("Before", "After"))
  ggplot2::ggplot(plot_df, ggplot2::aes(x = Comparison, fill = Category_plot)) +
    ggplot2::geom_bar(width = 0.78) +
    ggplot2::facet_grid(Direction ~ Filter) +
    ggplot2::scale_fill_manual(values = category_colors) +
    ggplot2::labs(title = paste(group, "| GO BP pathway composition"), x = NULL,
                  y = "Number of enriched GO BP pathways", fill = "Functional category") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())
}

gpp_plot_gene_retention <- function(result, group, comparison_order = NULL,
                                    retained_color = "#2A9D8F", removed_color = "#E76F51") {
  df <- result$summary
  gsea_nb_require(df, c("Group", "Comparison", "Direction", "n_genes_before", "n_genes_after"))
  df <- df[as.character(df$Group) == group, , drop = FALSE]
  if (!nrow(df)) stop("没有 Group：", group, call. = FALSE)
  df$Retained_Percent <- ifelse(df$n_genes_before > 0,
                                100 * df$n_genes_after / df$n_genes_before, NA_real_)
  if (is.null(comparison_order)) comparison_order <- unique(as.character(df$Comparison))
  df$Comparison <- factor(as.character(df$Comparison), levels = comparison_order)
  df$Direction <- factor(as.character(df$Direction), levels = c("Positive", "Negative"))
  full <- df[df$n_genes_before > 0, , drop = FALSE]
  full$Status <- "Removed"; full$Percentage <- 100
  retained <- df[df$n_genes_before > 0, , drop = FALSE]
  retained$Status <- "Retained"; retained$Percentage <- retained$Retained_Percent
  ggplot2::ggplot() +
    ggplot2::geom_col(data = full, ggplot2::aes(x = Comparison, y = Percentage, fill = Status), width = 0.72) +
    ggplot2::geom_col(data = retained, ggplot2::aes(x = Comparison, y = Percentage, fill = Status),
                      width = 0.72, na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = 75, linetype = "dashed", color = "#333333", linewidth = 0.6) +
    ggplot2::geom_text(data = df[df$n_genes_before > 0, , drop = FALSE],
                       ggplot2::aes(x = Comparison, y = Retained_Percent / 2,
                                    label = sprintf("%.1f%%", Retained_Percent)),
                       color = "white", fontface = "bold", size = 3.5) +
    ggplot2::facet_grid(Direction ~ .) +
    ggplot2::scale_fill_manual(values = c(Retained = retained_color, Removed = removed_color)) +
    ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                                labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
    ggplot2::labs(title = paste(group, "| Core enrichment gene retention"), x = NULL,
                  y = "Unique core genes retained", fill = NULL) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = "top", panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())
}
